"""Pairwise proximity statistics of the packed cell (IMPROVEMENTS I3): how often is a copy of protein A within a
surface gap g of a copy of protein B?

Surface gaps are measured between the molecules' sphere trees (the simulation's schematic bodies, library
spheretree.bin), not between centres, so large and elongated complexes are treated by their shape. For every
instance i (one placed copy of a library type: a monomer, oligomer or complex) and every protein p:
  C_g[type(i), p] += 1  if any other instance containing p has a sphere within gap g of a sphere of i
  M_g[type(i), p] += number of other instances containing p within gap g
Gaps g = 2, 5, 10 nm (contact / cross-linking and PCA range / BioID labelling radius). Protein-level statistics
(proximity(): P(A copy has B within g), mean number of B neighbours) sum over the types that contain the protein.
Subunits of one complex instance are not counted as neighbours of each other (they are one body).
Not included: mRNA, DNA linkers, lipids (not instances).

  python contacts.py            # decode tiles (64 processes), then GPU pass on all visible GPUs -> work/contacts/
"""
import os, sys, json, time, argparse
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "simgpu"))
from config import WEB_OUT, WORK, NPROC

OUT = f"{WORK}/contacts"
GAPS = (2.0, 5.0, 10.0)
TILE = 200.0
H = 25.0                   # cell size >= 2 * max sphere radius (7.3 nm) + max gap


def proteins_of(types):
    """protein (gene) list and CSR membership type -> proteins"""
    genes = [sorted({g.strip() for g in t["genes"].split(";") if g.strip()}) for t in types]
    names = sorted({g for gs in genes for g in gs})
    idx = {g: k for k, g in enumerate(names)}
    off = np.zeros(len(types) + 1, np.int64)
    off[1:] = np.cumsum([len(g) for g in genes])
    mem = np.array([idx[g] for gs in genes for g in gs], np.int32)
    return names, off, mem


def _decode(name):
    import boxes
    ix, iy, iz = (int(v) for v in name.split(".")[0].split("_"))
    t = boxes.read_tile(ix, iy, iz, WEB_OUT)
    return (ix, iy, iz), t["type"].astype(np.int16), t["region"].astype(np.int8), t["pos"].astype(np.float32), \
        t["quat"].astype(np.float32)


def decode_all(log=print):
    """all primary instances of the encoded tiles -> OUT/instances.npz (sorted by tile)"""
    import multiprocessing as mp
    names = sorted(f for f in os.listdir(f"{WEB_OUT}/tiles") if f.endswith(".bin.gz"))
    t0 = time.time()
    with mp.get_context("fork").Pool(NPROC) as pool:
        parts = pool.map(_decode, names, chunksize=16)
    tiles = np.array([p[0] for p in parts], np.int32)
    n = np.array([len(p[1]) for p in parts])
    off = np.zeros(len(n) + 1, np.int64); off[1:] = np.cumsum(n)
    np.savez(f"{OUT}/instances.npz", tiles=tiles, off=off, type=np.concatenate([p[1] for p in parts]),
             region=np.concatenate([p[2] for p in parts]), pos=np.concatenate([p[3] for p in parts]),
             quat=np.concatenate([p[4] for p in parts]))
    log(f"decoded {len(names)} tiles, {off[-1] / 1e6:.2f} M instances, {time.time() - t0:.0f} s")


_SRC = r"""
extern "C" {
__global__ void cellid(const float* __restrict__ w, int* __restrict__ cid, int ns, int nc, float h) {
  int s = blockIdx.x * blockDim.x + threadIdx.x; if (s >= ns) return;
  int cx = min(nc - 1, max(0, (int)(w[3*s] / h))), cy = min(nc - 1, max(0, (int)(w[3*s+1] / h))),
      cz = min(nc - 1, max(0, (int)(w[3*s+2] / h)));
  cid[s] = (cx * nc + cy) * nc + cz;
}
// one thread per query sphere: records (body, other body, gap) for spheres of other bodies within gmax
__global__ void near(const float* __restrict__ w, const float* __restrict__ rad, const int* __restrict__ sb,
                     const int* __restrict__ order, const int* __restrict__ cstart, const int* __restrict__ cend,
                     const int* __restrict__ cid, int nq, int nc, float gmax,
                     int* __restrict__ cnt, int cap, int* __restrict__ ra, int* __restrict__ rb, float* __restrict__ rg) {
  int s = blockIdx.x * blockDim.x + threadIdx.x; if (s >= nq) return;
  int b = sb[s];
  float x = w[3*s], y = w[3*s+1], z = w[3*s+2], r = rad[s];
  int c = cid[s]; int cz = c % nc, cy = (c / nc) % nc, cx = c / (nc * nc);
  for (int dx = -1; dx <= 1; dx++) for (int dy = -1; dy <= 1; dy++) for (int dz = -1; dz <= 1; dz++) {
    int ox = cx + dx, oy = cy + dy, oz = cz + dz;
    if (ox < 0 || oy < 0 || oz < 0 || ox >= nc || oy >= nc || oz >= nc) continue;
    int cc = (ox * nc + oy) * nc + oz;
    for (int k = cstart[cc]; k < cend[cc]; k++) {
      int t = order[k]; int bt = sb[t]; if (bt == b) continue;
      float ex = x - w[3*t], ey = y - w[3*t+1], ez = z - w[3*t+2];
      float g = sqrtf(ex*ex + ey*ey + ez*ez) - r - rad[t];
      if (g > gmax) continue;
      int i = atomicAdd(cnt, 1);
      if (i < cap) { ra[i] = b; rb[i] = bt; rg[i] = g; }
    }
  }
}
}
"""


def _gpu_worker(dev, jobs, q):
    os.environ["CUDA_VISIBLE_DEVICES"] = str(dev)
    import cupy as cp
    import cupyx as cpx
    import boxes
    mod = cp.RawModule(code=_SRC)
    k_cell, k_near = mod.get_function("cellid"), mod.get_function("near")
    lib = boxes.Library(WEB_OUT)
    types = lib.types
    NT = len(types)
    names, poff, pmem = proteins_of(types)
    NP = len(names)
    so = np.array([int(t["so"]) for t in types]); sn = np.array([int(t["sn"]) for t in types])
    br = np.array([float(t["br"]) for t in types])
    st = lib.st.astype(np.float64)
    Z = np.load(f"{OUT}/instances.npz")
    tiles, off, TY, POS, QU = Z["tiles"], Z["off"], Z["type"].astype(np.int64), Z["pos"], Z["quat"]
    tindex = {tuple(t): k for k, t in enumerate(tiles)}
    poff_g, pmem_g = cp.asarray(poff), cp.asarray(pmem)
    C = [cp.zeros(NT * NP, cp.int32) for _ in GAPS]
    M = [cp.zeros(NT * NP, cp.int32) for _ in GAPS]
    gmax = max(GAPS); margin = H + float(br.max())
    cap = 40_000_000
    cnt = cp.zeros(1, cp.int32); ra = cp.empty(cap, cp.int32); rb = cp.empty(cap, cp.int32); rg = cp.empty(cap, cp.float32)

    def spheres(idx):
        """world spheres of instances idx: centres (n,3), radii, owning instance (position in idx)"""
        t = TY[idx]; k = sn[t]
        rep = np.repeat(np.arange(len(idx)), k)
        first = np.repeat(so[t] - np.concatenate([[0], np.cumsum(k)[:-1]]), k) + np.arange(k.sum())
        loc = st[first, :3]
        qq = QU[idx][rep].astype(np.float64)
        w_, u = qq[:, :1], qq[:, 1:]
        tt = 2 * np.cross(u, loc)
        xyz = POS[idx][rep] + loc + w_ * tt + np.cross(u, tt)
        return xyz, st[first, 3], rep

    done = 0
    for ti in jobs:
        lo = tiles[ti].astype(np.float64) * TILE
        own = np.arange(off[ti], off[ti + 1])
        if len(own) == 0:
            continue
        ctx = []
        for d in np.ndindex(3, 3, 3):
            nb = tindex.get(tuple(tiles[ti] + np.array(d) - 1))
            if nb is None or nb == ti:
                continue
            j = np.arange(off[nb], off[nb + 1])
            p = POS[j]
            m = np.all((p > lo - margin) & (p < lo + TILE + margin), 1)
            ctx.append(j[m])
        ctx = np.concatenate([own] + ctx)                    # owner instances first
        xyz, rad, rep = spheres(ctx)
        qmask = rep < len(own)
        o = xyz[qmask].min(0) - H; hi = xyz[qmask].max(0) + H      # window: query spheres +- one cell
        keep = np.all((xyz > o) & (xyz < hi), 1) | qmask
        xyz, rad, rep = xyz[keep], rad[keep], rep[keep]
        qmask = rep < len(own)
        perm = np.concatenate([np.where(qmask)[0], np.where(~qmask)[0]])  # query spheres first
        xyz, rad, rep = xyz[perm], rad[perm], rep[perm]
        nq, ns = int(qmask.sum()), len(rep)
        nc = int(np.ceil((hi - o).max() / H)) + 1
        w = cp.asarray((xyz - o).astype(np.float32))
        rad_g = cp.asarray(rad.astype(np.float32)); sb = cp.asarray(rep.astype(np.int32))
        cid = cp.empty(ns, cp.int32)
        bs = 128
        k_cell(((ns + bs - 1) // bs,), (bs,), (w, cid, np.int32(ns), np.int32(nc), np.float32(H)))
        order = cp.argsort(cid).astype(cp.int32); sc = cid[order]
        cells = cp.arange(nc ** 3, dtype=cp.int32)
        cstart = cp.searchsorted(sc, cells, side="left").astype(cp.int32)
        cend = cp.searchsorted(sc, cells, side="right").astype(cp.int32)
        while True:
            cnt.fill(0)
            k_near(((nq + bs - 1) // bs,), (bs,), (w, rad_g, sb, order, cstart, cend, cid, np.int32(nq), np.int32(nc),
                                                  np.float32(gmax), cnt, np.int32(cap), ra, rb, rg))
            nrec = int(cnt.get()[0])
            if nrec <= cap:
                break
            cap = int(nrec * 1.3); ra = cp.empty(cap, cp.int32); rb = cp.empty(cap, cp.int32); rg = cp.empty(cap, cp.float32)
        if nrec == 0:
            continue
        a, b, g = ra[:nrec].astype(cp.int64), rb[:nrec].astype(cp.int64), rg[:nrec]
        # min gap per instance pair
        key = a * len(ctx) + b
        o_ = cp.lexsort(cp.stack([g, key])); key, g = key[o_], g[o_]
        first = cp.concatenate([cp.ones(1, bool), key[1:] != key[:-1]])
        key, g = key[first], g[first]
        a, b = key // len(ctx), key % len(ctx)
        ctx_t = cp.asarray(TY[ctx])
        ta, tb = ctx_t[a], ctx_t[b]
        # expand the neighbour instance to its proteins
        nper = (poff_g[tb + 1] - poff_g[tb]).astype(cp.int64)
        csum = cp.cumsum(nper); tot = int(csum[-1])
        rowid = cp.searchsorted(csum, cp.arange(tot), side="right")
        within = cp.arange(tot) - (csum[rowid] - nper[rowid])
        p = pmem_g[poff_g[tb[rowid]] + within].astype(cp.int64)
        ea, eg = a[rowid], g[rowid]
        cell_key = ta[rowid] * NP + p
        for k, gk in enumerate(GAPS):
            m = eg <= gk
            cpx.scatter_add(M[k], cell_key[m], 1)
            uk = cp.unique(ea[m] * NP + p[m])                    # (instance, protein) presence
            cpx.scatter_add(C[k], ctx_t[uk // NP] * NP + uk % NP, 1)
        done += 1
    q.put((dev, [c.get() for c in C], [m_.get() for m_ in M], done))


def run(log=print):
    import multiprocessing as mp
    import cupy as cp
    os.makedirs(OUT, exist_ok=True)
    if not os.path.exists(f"{OUT}/instances.npz"):
        decode_all(log)
    Z = np.load(f"{OUT}/instances.npz")
    nt = len(Z["tiles"]); sizes = np.diff(Z["off"])
    ngpu = cp.cuda.runtime.getDeviceCount()
    order = np.argsort(-sizes)
    jobs = [order[k::ngpu] for k in range(ngpu)]
    t0 = time.time()
    ctx = mp.get_context("spawn"); q = ctx.Queue()
    ps = [ctx.Process(target=_gpu_worker, args=(k, jobs[k], q)) for k in range(ngpu)]
    for p in ps:
        p.start()
    res = [q.get() for _ in ps]
    for p in ps:
        p.join()
    import boxes
    types = boxes.Library(WEB_OUT).types
    names, poff, pmem = proteins_of(types)
    NT, NP = len(types), len(names)
    n_inst = np.bincount(Z["type"].astype(np.int64), minlength=NT)
    out = dict(proteins=np.array(names), type_ids=np.array([t["id"] for t in types]), n_inst=n_inst,
               poff=poff, pmem=pmem, gaps=np.array(GAPS))
    for k, gk in enumerate(GAPS):
        out[f"C{int(gk)}"] = sum(r[1][k].astype(np.int64) for r in res).reshape(NT, NP)
        out[f"M{int(gk)}"] = sum(r[2][k].astype(np.int64) for r in res).reshape(NT, NP)
    np.savez_compressed(f"{OUT}/contacts.npz", **out)
    log(f"contacts: {sum(r[3] for r in res)} tiles on {ngpu} GPUs, {time.time() - t0:.0f} s -> {OUT}/contacts.npz")


class Proximity:
    """protein-level proximity from contacts.npz"""

    def __init__(self, path=f"{OUT}/contacts.npz"):
        Z = np.load(path)
        self.proteins = list(Z["proteins"]); self.pidx = {p: k for k, p in enumerate(self.proteins)}
        self.n_inst, self.poff, self.pmem = Z["n_inst"], Z["poff"], Z["pmem"]
        self.gaps = list(Z["gaps"])
        self.C = {int(g): Z[f"C{int(g)}"] for g in self.gaps}
        self.M = {int(g): Z[f"M{int(g)}"] for g in self.gaps}
        NT = len(self.n_inst)
        rows = np.repeat(np.arange(NT), np.diff(self.poff))
        # protein x type incidence (instances of type t contain protein p)
        from scipy.sparse import csr_matrix
        self.PT = csr_matrix((np.ones(len(self.pmem)), (self.pmem, rows)), shape=(len(self.proteins), NT))
        self.copies = self.PT @ self.n_inst               # instances containing each protein

    def matrix(self, g=10, kind="C"):
        """protein x protein: kind C -> P(an instance with A has B within g); M -> mean number of B neighbours"""
        X = (self.C if kind == "C" else self.M)[int(g)]
        S = np.asarray(self.PT @ X)
        with np.errstate(invalid="ignore", divide="ignore"):
            return S / self.copies[:, None]

    def pair(self, a, b, g=10):
        ia, ib = self.pidx[a], self.pidx[b]
        row = self.PT[ia]
        ts = row.indices
        n = self.n_inst[ts].sum()
        return dict(A=a, B=b, gap_nm=g, copies_A=int(n),
                    p_near=float(self.C[int(g)][ts, ib].sum() / max(n, 1)),
                    mean_B_neighbours=float(self.M[int(g)][ts, ib].sum() / max(n, 1)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--decode-only", action="store_true")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    if a.decode_only:
        decode_all()
    else:
        run()
