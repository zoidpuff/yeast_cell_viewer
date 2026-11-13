"""Per-protein feature table for the data-driven layers (localization benchmark, diffusion model, condensate sampler).

Output: work/features/protein_features.parquet, one row per UniProt entry of the inventory, with
  abundance and compartment weights (inventory.build), sequence composition, disordered regions (IDRs: AlphaFold
  residues with pLDDT < 50 in runs of >= 10, the residues the structural model trims), IDR composition (charge,
  aromatic and arginine "stickers", Q/N prion-like content, RG/RGG motifs), hydrodynamic radius of the folded part and
  of the full chain with its IDRs as random coils (Kirkwood, as the simulation's hydro()), and AlphaFold 3 pool
  interaction degrees at several PAE thresholds.
usage: python pipeline/features.py
"""
import os, io, re, gzip, json, tarfile
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

from config import WORK, DATA, AFDB_TAR, NPROC

OUT = f"{WORK}/features"
os.makedirs(OUT, exist_ok=True)
PLDDT_IDR = 50.0          # the structural model deletes AlphaFold residues below this (structlib)
IDR_MIN_RUN = 10          # shorter low-confidence runs are loops, not disordered regions
HYDRATION = 0.3           # nm, as params.ts
BOND = 0.38               # nm, C-alpha virtual bond


def read_fasta(path=f"{DATA}/uniprot/UP000002311.fasta"):
    seqs, acc, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if acc:
                seqs[acc] = "".join(buf)
            acc = line.split("|")[1]; buf = []
        else:
            buf.append(line.strip())
    if acc:
        seqs[acc] = "".join(buf)
    return seqs


def _parse_pdb(data):
    """CA coordinates (nm) and pLDDT per residue from an AlphaFold DB PDB file"""
    xyz, pl = [], []
    for line in gzip.decompress(data).decode().splitlines():
        if line.startswith("ATOM") and line[12:16] == " CA ":
            xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            pl.append(float(line[60:66]))
    return np.array(xyz) / 10.0, np.array(pl)


def afdb_models():
    """uniprot -> (CA xyz nm, pLDDT); cached in work/features/afdb_ca.npz"""
    cache = f"{OUT}/afdb_ca.npz"
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        return z["models"].item()
    out = {}
    with tarfile.open(AFDB_TAR) as tf:
        for m in tf:
            if m.name.endswith(".pdb.gz") and "-F1-" in m.name:
                acc = m.name.split("-")[1]
                out[acc] = _parse_pdb(tf.extractfile(m).read())
    np.savez(cache, models=np.array(out, dtype=object))
    return out


def idr_segments(pl):
    """[(start, end)) runs of pLDDT < PLDDT_IDR of at least IDR_MIN_RUN residues"""
    low = pl < PLDDT_IDR
    segs, i, n = [], 0, len(pl)
    while i < n:
        if low[i]:
            j = i
            while j < n and low[j]:
                j += 1
            if j - i >= IDR_MIN_RUN:
                segs.append((i, j))
            i = j
        else:
            i += 1
    return segs


def kirkwood_rh(x, nsub=700):
    """R_h (nm) of a bead model: 1/R_K = <1/r_ij>, x 6/5, + hydration (params.ts hydro())"""
    n = len(x)
    if n < 2:
        return 0.5
    m = min(n, nsub)
    y = x[np.floor(np.arange(m) * (n / m)).astype(int)]
    d = np.sqrt(((y[:, None] - y[None]) ** 2).sum(-1))
    iu = np.triu_indices(m, 1)
    rk = len(iu[0]) / (1.0 / np.maximum(d[iu], 0.38)).sum()
    return max(0.5, 1.2 * rk + HYDRATION)


def grow_coil(start, direction, n, obstacles, rng, theta=2.15, excl=0.42):
    """self-avoiding C-alpha chain of n beads from `start`, first step along `direction`: fixed virtual bond angle
    theta (rad), random dihedral, beads kept `excl` nm from each other and from the obstacle beads. Returns (n,3).
    theta = 2.15 gives R_g and R_h of free coils within ~10% of the IDP scaling laws R_g = 0.254 N^0.522 nm and
    R_h(Kirkwood, without hydration) = 0.249 N^0.509 nm for N = 50-400 (unchecked against the sources)."""
    pts = []
    prev = start - direction * BOND
    cur = start
    obs = obstacles
    for k in range(n):
        for _ in range(60):
            b = cur - prev; b /= max(np.linalg.norm(b), 1e-9)
            # random direction at angle pi - theta to the previous bond
            a = rng.normal(size=3); a -= a.dot(b) * b; a /= np.linalg.norm(a)
            ct = np.cos(np.pi - theta)
            v = ct * b + np.sqrt(1 - ct * ct) * a
            nxt = cur + BOND * v
            ok = True
            if len(pts) > 2 and np.min(np.sum((np.array(pts[:-2]) - nxt) ** 2, 1)) < excl ** 2:
                ok = False
            if ok and len(obs) and np.min(np.sum((obs - nxt) ** 2, 1)) < excl ** 2:
                ok = False
            if ok:
                break
        pts.append(nxt); prev, cur = cur, nxt
    return np.array(pts)


def full_chain_rh(xyz, pl, segs, n_conf=4, seed=0):
    """R_h of the folded part and mean R_h of the full chain with IDRs regrown as self-avoiding coils"""
    keep = np.ones(len(pl), bool)
    for s, e in segs:
        keep[s:e] = False
    fold = xyz[keep]
    rh_f = kirkwood_rh(fold) if keep.sum() >= 2 else 0.5
    if not segs:
        return rh_f, rh_f
    rng = np.random.default_rng(seed)
    rhs = []
    cen = fold.mean(0) if len(fold) else xyz.mean(0)
    for _ in range(n_conf):
        chain = xyz.copy()
        obs = fold.copy() if len(fold) else np.zeros((0, 3))
        for s, e in segs:
            n = e - s
            if s > 0:                       # grow from the preceding ordered residue toward the C-terminus
                anchor = chain[s - 1]; direction = anchor - cen
                seg = grow_coil(anchor, direction / (np.linalg.norm(direction) + 1e-9), n, obs, rng)
                chain[s:e] = seg
            else:                           # N-terminal IDR: grow backwards from the first ordered residue
                anchor = chain[e] if e < len(chain) else cen
                direction = anchor - cen
                seg = grow_coil(anchor, direction / (np.linalg.norm(direction) + 1e-9), n, obs, rng)
                chain[s:e] = seg[::-1]
            obs = np.vstack([obs, seg]) if len(obs) else seg
        rhs.append(kirkwood_rh(chain))
    return rh_f, float(np.mean(rhs))


AROM = set("FYW")


def seq_features(seq, segs):
    L = len(seq)
    pos, neg = sum(seq.count(c) for c in "KR"), sum(seq.count(c) for c in "DE")
    idr = "".join(seq[s:e] for s, e in segs)
    n_idr = len(idr)
    f = dict(length=L, net_charge=pos - neg, fcr=(pos + neg) / max(1, L), n_idr=n_idr, f_idr=n_idr / max(1, L),
             n_idr_segments=len(segs), longest_idr=max([e - s for s, e in segs], default=0))
    if n_idr:
        ip, ineg = sum(idr.count(c) for c in "KR"), sum(idr.count(c) for c in "DE")
        f.update(idr_net_charge=ip - ineg, idr_fcr=(ip + ineg) / n_idr,
                 idr_arom=sum(idr.count(c) for c in AROM), idr_arg=idr.count("R"),
                 idr_qn=(idr.count("Q") + idr.count("N")) / n_idr,
                 idr_rg=len(re.findall("RG", idr)), idr_rgg=len(re.findall("RGG", idr)),
                 idr_sticker=sum(idr.count(c) for c in AROM) + idr.count("R"),
                 idr_g=idr.count("G") / n_idr, idr_s=idr.count("S") / n_idr, idr_p=idr.count("P") / n_idr)
    else:
        f.update(idr_net_charge=0, idr_fcr=0.0, idr_arom=0, idr_arg=0, idr_qn=0.0, idr_rg=0, idr_rgg=0,
                 idr_sticker=0, idr_g=0.0, idr_s=0.0, idr_p=0.0)
    # prion-like composition: maximum Q+N+S+G+Y fraction in an 80-residue window (Alberti et al. 2009-style PrLD proxy)
    w = 80
    if L >= w:
        a = np.array([c in "QNSGY" for c in seq], float)
        cs = np.concatenate([[0], np.cumsum(a)])
        f["prld_max80"] = float((cs[w:] - cs[:-w]).max() / w)
        qn = np.array([c in "QN" for c in seq], float); cq = np.concatenate([[0], np.cumsum(qn)])
        f["qn_max80"] = float((cq[w:] - cq[:-w]).max() / w)
    else:
        f["prld_max80"] = f["qn_max80"] = 0.0
    return f


def _one(args):
    acc, seq, model = args
    if model is None:
        return acc, seq_features(seq, []), np.nan, np.nan
    xyz, pl = model
    segs = idr_segments(pl)
    f = seq_features(seq, segs) if len(seq) == len(pl) else seq_features(seq, [])
    try:
        rh_f, rh_full = full_chain_rh(xyz, pl, segs, seed=hash(acc) % 2 ** 31)
    except Exception:
        rh_f = rh_full = np.nan
    f["plddt_mean"] = float(pl.mean())
    return acc, f, rh_f, rh_full


def pool_degrees(accs, abundance):
    """AlphaFold 3 pool pair scores (work/interactions/pair_summary.parquet): partners per protein at PAE thresholds
    (best chain-pair PAE min-min), and partner copies weighted"""
    ps = pd.read_parquet(f"{WORK}/interactions/pair_summary.parquet", columns=["p1", "p2", "pae_min", "iptm_max"])
    ps["p1"] = ps.p1.str.upper(); ps["p2"] = ps.p2.str.upper()
    ps = ps[ps.p1 != ps.p2]
    out = {a: {} for a in accs}
    for lo, hi, name in ((0, 8, "pool_deg_le8"), (8, 12, "pool_deg_8_12"), (12, 16, "pool_deg_12_16")):
        m = (ps.pae_min > lo) & (ps.pae_min <= hi) if lo else ps.pae_min <= hi
        sub = ps[m]
        for col in ("p1", "p2"):
            c = sub[col].value_counts()
            for a, v in c.items():
                if a in out:
                    out[a][name] = out[a].get(name, 0) + int(v)
        # partner copies (abundance-weighted degree)
        ab = abundance
        w = sub.p2.map(ab).fillna(0.0).groupby(sub.p1).sum().add(sub.p1.map(ab).fillna(0.0).groupby(sub.p2).sum(), fill_value=0)
        for a, v in w.items():
            if a in out:
                out[a][name + "_copies"] = float(v)
    return pd.DataFrame.from_dict(out, orient="index").fillna(0.0)


def build():
    import inventory as I
    inv, _ = I.build()
    seqs = read_fasta()
    models = afdb_models()
    accs = [a for a in inv.uniprot if isinstance(a, str)]
    jobs = [(a, seqs.get(a, ""), models.get(a)) for a in accs]
    rows = {}
    with ProcessPoolExecutor(NPROC) as ex:
        for acc, f, rh_f, rh_full in ex.map(_one, jobs, chunksize=8):
            f["rh_folded"] = rh_f; f["rh_full"] = rh_full
            rows[acc] = f
    F = pd.DataFrame.from_dict(rows, orient="index")
    inv = inv.drop_duplicates("uniprot").set_index("uniprot")
    F = inv[["orf", "gene", "copies_per_cell", "mass_da", "n_tm", "compartments", "comp_weights", "loc_source"]].join(F, how="right")
    deg = pool_degrees(list(F.index), F.copies_per_cell.to_dict())
    F = F.join(deg).fillna({c: 0.0 for c in deg.columns})
    F.index.name = "uniprot"
    F.to_parquet(f"{OUT}/protein_features.parquet")
    return F


if __name__ == "__main__":
    F = build()
    print(F.shape)
    print(F[["length", "f_idr", "rh_folded", "rh_full", "net_charge", "pool_deg_le8", "pool_deg_8_12"]].describe().T)
