"""Boxes of the packed cell for the GPU simulations, read from the encoded viewer data (config.WEB_OUT).

Mirrors web/app/src/sim/extract.ts at the schematic resolution: every packed instance whose centre lies in the box
becomes a rigid body with its library sphere tree; structures built explicitly (filaments, lattices, nucleosomes,
pores; instance flag bit 1) are fixed obstacles; membrane proteins are fixed here (no lateral diffusion yet).
Hydrodynamic radii as in params.ts: Kirkwood double sum over the residue beads (subsampled to 700), x 6/5, + 0.3 nm.
"""
import gzip, json, os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pipeline"))
from config import WEB_OUT  # noqa: E402

TILE, Q_ORG, Q_RANGE = 200.0, 80.0, 360.0
Q_STEP = Q_RANGE / 65535.0
INV_SQRT2 = 1.0 / np.sqrt(2.0)
R_MEM0, R_MEM1, R_NUCLEOSOME, R_NPC = 14, 21, 22, 23
FIXED_IDS = {"F_actin", "Arp23", "CP_cap", "TPM1_coil", "MT_tubulin", "gTuRC_nuc", "gTuRC_cyto", "ATP_dimer",
             "Pil1_dimer", "Lsp1_dimer", "Vac8_Nvj1", "RVS161_RVS167_het", "NPC", "nucleosome"}
KT = 1.380649e-23 * 303.15
ETA_WATER = 0.797e-3
HYDRATION = 0.3
HYD_VOL_PER_DA = (0.73 + 0.35) / 602.214       # nm^3 per Da, hydrated


class Library:
    """types.json + residue beads + sphere trees of the encoded data"""

    def __init__(self, data=WEB_OUT):
        tj = json.load(open(os.path.join(data, "library", "types.json")))
        self.types = tj["types"]
        raw = np.fromfile(os.path.join(data, "library", "beads.bin"), np.uint8)
        nb = raw.size // 8
        i16 = raw.view(np.int16).reshape(nb, 4)
        self.beads = i16[:, :3].astype(np.float32) * tj["bead_q"]
        self.kind = raw.reshape(nb, 8)[:, 6]                       # 1 = nucleotide
        self.st = np.fromfile(os.path.join(data, "library", "spheretree.bin"), np.float32).reshape(-1, 4)
        self.by_id = {t["id"]: i for i, t in enumerate(self.types)}
        self._rh = {}

    def spheres(self, t):
        T = self.types[t]
        return self.st[T["so"]:T["so"] + T["sn"]]

    def bead_xyz(self, t):
        T = self.types[t]
        return self.beads[T["bo"]:T["bo"] + T["nb"]]

    def rh(self, t):
        """hydrodynamic radius (nm) of library type t, as hydro() in params.ts"""
        if t in self._rh:
            return self._rh[t]
        x = self.bead_xyz(t).astype(np.float64)
        n = len(x)
        rk = 0.5
        if n >= 2:
            m = min(n, 700)
            idx = np.floor(np.arange(m) * (n / m)).astype(int)
            y = x[idx]
            d = np.sqrt(((y[:, None, :] - y[None, :, :]) ** 2).sum(-1))
            iu = np.triu_indices(m, 1)
            rk = len(iu[0]) / (1.0 / np.maximum(d[iu], 0.38)).sum()
        r = max(0.5, 1.2 * rk + HYDRATION)
        self._rh[t] = r
        return r


def d_water(rh_nm):
    """Stokes-Einstein translational diffusion in water at 30 C, nm^2/ns (= 1e3 um^2/s ... 1 nm^2/ns = 1000 um^2/s)"""
    return KT / (6 * np.pi * ETA_WATER * np.asarray(rh_nm) * 1e-9) * 1e9


def dr_water(rh_nm):
    return KT / (8 * np.pi * ETA_WATER * (np.asarray(rh_nm) * 1e-9) ** 3) * 1e-9


def _quat_decode(w0, w2, w3):
    f = lambda v: ((v.astype(np.int32) << 16) >> 16).astype(np.int16).astype(np.float64) / 32767 * INV_SQRT2
    a, b, c = f(w2 >> 16), f(w3 & 0xFFFF), f(w3 >> 16)
    big = np.sqrt(np.clip(1 - a * a - b * b - c * c, 0, 1))
    qi = (w0 >> 14) & 3
    q = np.zeros((len(w0), 4))                                    # w, x, y, z
    for k in range(4):
        m = qi == k
        others = [j for j in range(4) if j != k]
        q[m, k] = big[m]
        q[m, others[0]], q[m, others[1]], q[m, others[2]] = a[m], b[m], c[m]
    return q


def read_tile(ix, iy, iz, data=WEB_OUT):
    """instances of one encoded tile: type, region, flags, position (nm), quaternion (w,x,y,z); None if absent"""
    p = os.path.join(data, "tiles", f"{ix}_{iy}_{iz}.bin.gz")
    if not os.path.exists(p):
        return None
    b = gzip.decompress(open(p, "rb").read())
    h = np.frombuffer(b, np.uint32, 16)
    assert h[0] == int.from_bytes(b"YCT1", "little"), "bad tile magic"
    n = int(h[4])
    w = np.frombuffer(b, np.uint32, 4 * n, 64).reshape(n, 4)
    lo = np.array([ix, iy, iz], np.float64) * TILE
    px, py = w[:, 1] & 0xFFFF, w[:, 1] >> 16
    pz = w[:, 2] & 0xFFFF
    pos = (lo - Q_ORG) + np.stack([px, py, pz], 1).astype(np.float64) * Q_STEP
    out = dict(type=(w[:, 0] & 0x3FFF).astype(np.int32), region=((w[:, 0] >> 16) & 31).astype(np.int32),
               flags=((w[:, 0] >> 21) & 0xFF).astype(np.int32), pos=pos, quat=_quat_decode(w[:, 0], w[:, 2], w[:, 3]))
    # polysome identity of each instance (polymer_annot.py: transcript and rank from the 5' end; -1 if none)
    tx = np.full(n, -1, np.int64); rank = np.full(n, -1, np.int64)
    a = os.path.join(data, "annot", f"{ix}_{iy}_{iz}.json")
    if os.path.exists(a):
        for loc, code, x, y in json.loads(gzip.decompress(open(a, "rb").read()))["i"]:
            if code == 1:
                tx[loc] = x; rank[loc] = y
    out["poly_tx"] = tx; out["poly_rank"] = rank
    return out


def quat_rotate(q, v):
    """rotate vectors v (n,3) by unit quaternions q (n,4) w,x,y,z"""
    w, u = q[:, :1], q[:, 1:]
    t = 2 * np.cross(u, v)
    return v + w * t + np.cross(u, t)


def cut_box(lib, center, size, data=WEB_OUT):
    """bodies with centres in the cube [center - size/2, center + size/2): dict of numpy arrays"""
    c = np.asarray(center, np.float64)
    lo, hi = c - size / 2, c + size / 2
    t0, t1 = np.floor(lo / TILE).astype(int), np.floor(hi / TILE).astype(int)
    parts = []
    for ix in range(t0[0], t1[0] + 1):
        for iy in range(t0[1], t1[1] + 1):
            for iz in range(t0[2], t1[2] + 1):
                t = read_tile(ix, iy, iz, data)
                if t is None:
                    continue
                m = np.all((t["pos"] >= lo) & (t["pos"] < hi), 1)
                parts.append({k: v[m] for k, v in t.items()})
    box = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    n = len(box["type"])
    ids = [lib.types[t]["id"] for t in box["type"]]
    membrane = (box["region"] >= R_MEM0) & (box["region"] <= R_MEM1)
    scaffold = ((box["flags"] & 2) != 0) | np.array([i in FIXED_IDS for i in ids], bool) | \
        (box["region"] == R_NUCLEOSOME) | (box["region"] == R_NPC)
    box["mobile"] = ~membrane & ~scaffold
    box["rh"] = np.array([lib.rh(t) for t in box["type"]])
    box["mass"] = np.array([lib.types[t]["mass"] for t in box["type"]], np.float64)
    box["lo"], box["hi"], box["size"] = lo, hi, float(size)
    box["phi"] = float(min(0.5, box["mass"].sum() * HYD_VOL_PER_DA / size ** 3))
    box["n"] = n
    return box


def subset(box, keep):
    """the box with only the bodies of the boolean mask keep (e.g. a share of the free ribosomes removed)"""
    out = dict(box)
    for k, v in box.items():
        if isinstance(v, np.ndarray) and v.shape[:1] == (box["n"],):
            out[k] = v[keep]
    out["n"] = int(keep.sum())
    out["phi"] = float(min(0.5, out["mass"].sum() * HYD_VOL_PER_DA / box["size"] ** 3))
    return out


def polysome_bonds(tx, rank, pos, quat, channel):
    """links between neighbouring ribosomes of one mRNA that are both in the box: (i, j, local point on i, local
    point on j, spring constant kT/nm^2). The mRNA leaves ribosome i (rank r + 1 from the 5' end) at its channel
    exit and enters ribosome j (rank r) at its entry (polysome.py: the walk goes 3' to 5'). The linker is a Gaussian
    entropic spring whose rms extension is the placed exit-to-entry distance d0 (at least 3 nm): k = 3 / d0^2."""
    ex, en = channel[0], channel[-1]
    bi, bj, k = [], [], []
    idx = np.where(tx >= 0)[0]
    by = {}
    for i in idx:
        by[(int(tx[i]), int(rank[i]))] = i
    for (t, r), j in by.items():
        i = by.get((t, r + 1))
        if i is None:
            continue
        d0 = np.linalg.norm(pos[j] + quat_rotate(quat[j][None], en[None])[0] - pos[i] - quat_rotate(quat[i][None], ex[None])[0])
        bi.append(i); bj.append(j); k.append(3.0 / max(d0, 3.0) ** 2)
    return (np.array(bi, np.int64), np.array(bj, np.int64), np.tile(ex, (len(bi), 1)), np.tile(en, (len(bi), 1)),
            np.array(k))
