"""Per-bead data for the viewer's simulation mode (web/app/src/sim): residue charges and native contacts of predicted
complexes. Both are recovered from the source structures, since the structure library keeps only bead coordinates.

Outputs (in <data>/library, next to beads.bin, in its bead order: Morton order within each type):
  charges.bin     int8 per bead: +1 Lys/Arg, -1 Asp/Glu (cytosolic pH 7.2: His neutral), -1 per nucleotide (the
                  simulation scales nucleotides for counterion condensation), 0 otherwise or when unmatched
  contacts.json   {"cutoff_nm", "min_iptm", "pairs": [[typeA, typeB, offset, count, model], ...], stats}
  contacts.bin    uint16 pairs (bead of typeA, bead of typeB): inter-chain residue contacts (C-alpha within cutoff_nm)
                  of the confident AlphaFold DB yeast dimers (YC_DIMERS heterodimers/ and homodimers/, ipTM >= MIN_IPTM)
                  whose two proteins are both monomer types of the library

Matching library beads to source residues: a library type's protein beads are its source's C-alpha atoms after a
rigid motion (centring; membrane orientation) and dropping residues (low pLDDT, unused chains). Four-residue windows of
consecutive C-alphas are matched by their six internal distances; a window pose (Kabsch) that brings at least 90% of
the library beads within 0.03 nm of a source C-alpha is accepted, and every bead takes the residue of its nearest one.
Composite structures (a ribosome with a translocon) take the part found in the structure of the type their id extends
(40% of the beads suffice); unmatched beads stay neutral.

  python simdata.py [--data <viewer data dir>]      (default: config.WEB_OUT)
"""
import argparse, gzip, json, os, tarfile, time, glob, re
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from config import LIB, PDB, AFDB_TAR, DIMERS, WEB_OUT, NPROC, WORK

BEAD_Q = 0.005
CHARGE = {"ASP": -1, "GLU": -1, "LYS": 1, "ARG": 1}
CONTACT_NM = 0.8
MIN_IPTM = 0.5


def morton_order(b):
    """as webenc.morton_order (bead order within a type in beads.bin)"""
    lo, hi = b.min(0), b.max(0)
    g = ((b - lo) / np.maximum(hi - lo, 1e-6) * 1023).astype(np.uint64)

    def spread(x):
        x = (x | (x << 16)) & 0x030000FF; x = (x | (x << 8)) & 0x0300F00F
        x = (x | (x << 4)) & 0x030C30C3; x = (x | (x << 2)) & 0x09249249
        return x
    code = spread(g[:, 0]) | (spread(g[:, 1]) << np.uint64(1)) | (spread(g[:, 2]) << np.uint64(2))
    return np.argsort(code, kind="stable")


def ca_table(st):
    """C-alpha atoms of model 0: xyz (nm), residue name, chain, sequence number"""
    xyz, name, chain, num = [], [], [], []
    for ch in st[0]:
        for res in ch:
            a = res.find_atom("CA", "*")
            if a is None:
                continue
            xyz.append([a.pos.x / 10, a.pos.y / 10, a.pos.z / 10]); name.append(res.name); chain.append(ch.name); num.append(res.seqid.num)
    return dict(xyz=np.array(xyz, float).reshape(-1, 3), name=np.array(name), chain=np.array(chain), num=np.array(num, int))


def _kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, S, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T)); R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def _windows(X, chain=None):
    """6 internal distances of every window of 4 consecutive C-alphas (consecutive residues only)"""
    n = len(X)
    if n < 4:
        return np.zeros((0, 6)), np.zeros(0, int)
    step = np.linalg.norm(np.diff(X, axis=0), axis=1)
    ok = (step < 0.42)
    if chain is not None:
        ok &= chain[1:] == chain[:-1]
    good = np.where(ok[:-2] & ok[1:-1] & ok[2:])[0]
    F = np.stack([np.linalg.norm(X[good + a] - X[good + b], axis=1) for a, b in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))], 1)
    return F, good


def match(lib_xyz, src, tol=0.03, frac=0.9):
    """index of the source C-alpha of every library protein bead (-1: unmatched), or None if no pose fits"""
    S = src["xyz"]
    if len(lib_xyz) == 0 or len(S) == 0:
        return None
    tree = cKDTree(S)
    if len(lib_xyz) == len(S):
        R, t = _kabsch(lib_xyz, S)
        if np.sqrt(((lib_xyz @ R.T + t - S) ** 2).sum(1).mean()) < 0.01:
            return np.arange(len(S))
    Fs, gs = _windows(S, src["chain"])
    Fl, gl = _windows(lib_xyz)
    if not len(Fs) or not len(Fl):
        return None
    ftree = cKDTree(Fs)
    cand = ftree.query_ball_point(Fl, r=0.01)
    order = sorted([i for i in range(len(gl)) if cand[i]], key=lambda i: len(cand[i]))
    rng = np.random.default_rng(0)
    sample = rng.choice(len(lib_xyz), min(300, len(lib_xyz)), replace=False)
    for i in order[:60]:
        a = gl[i]
        for c in cand[i][:12]:
            b = gs[c]
            R, t = _kabsch(lib_xyz[a:a + 4], S[b:b + 4])
            d, _ = tree.query(lib_xyz[sample] @ R.T + t)
            if (d < tol).mean() >= frac:
                d, j = tree.query(lib_xyz @ R.T + t)
                return np.where(d < tol, j, -1)
    return None


def afdb_tables(accs):
    """C-alpha tables of the AlphaFold DB monomer models (fragment F1) of the given accessions, from the proteome tar"""
    import gemmi
    from joblib import Parallel, delayed

    def parse(acc, data):
        st = gemmi.make_structure_from_block(gemmi.cif.read_string(data)[0])
        return acc, ca_table(st)
    jobs = []
    with tarfile.open(AFDB_TAR) as tf:
        for m in tf:
            name = os.path.basename(m.name)
            if name.endswith(".cif.gz") and "-F1-" in name and name.split("-")[1] in accs:
                jobs.append((name.split("-")[1], gzip.decompress(tf.extractfile(m).read()).decode()))
    return dict(Parallel(n_jobs=NPROC, batch_size=16)(delayed(parse)(a, d) for a, d in jobs))


def pdb_files(pid):
    """cached assembly files of a PDB entry; fetched from RCSB when none is cached (assembly 1, then the full entry)"""
    import requests
    fs = sorted(glob.glob(f"{PDB}/{pid}-assembly*.cif.gz")) + sorted(glob.glob(f"{PDB}/{pid}.cif*"))
    if fs:
        return fs
    for name in (f"{pid}-assembly1.cif.gz", f"{pid}.cif.gz"):
        r = requests.get(f"https://files.rcsb.org/download/{name}", timeout=120)
        if r.status_code == 200:
            open(f"{PDB}/{name}", "wb").write(r.content); fs.append(f"{PDB}/{name}")
    return fs


def structure_table(path):
    import gemmi
    st = gemmi.read_structure(path); st.remove_hydrogens()
    return ca_table(st)


def build(data_dir=WEB_OUT, log=print):
    """write charges.bin and contacts.{json,bin} into data_dir/library; returns a summary"""
    a = argparse.Namespace(data=data_dir)
    print = log
    t0 = time.time()
    tj = json.load(open(f"{a.data}/library/types.json"))
    types, nbt = tj["types"], tj["n_beads"]
    raw = np.fromfile(f"{a.data}/library/beads.bin", dtype=[("x", "<i2"), ("y", "<i2"), ("z", "<i2"), ("k", "u1"), ("p", "u1")])
    assert len(raw) == nbt
    charges = np.zeros(nbt, np.int8)
    web_of = {}        # type index -> (npz index -> web bead index)
    res_of = {}        # monomer type index -> {sequence number: web bead index}
    # sources
    acc_needed = {t["uniprot"] for t in types if t["src"] == "AF" and t["uniprot"]}
    print(f"{len(types)} types; reading {len(acc_needed)} AlphaFold DB monomer models")
    af = afdb_tables(acc_needed)
    print(f"  {len(af)} read ({time.time() - t0:.0f} s)")
    stats = dict(types=len(types), protein_beads=0, protein_beads_matched=0, types_unmatched=[])
    for ti, t in enumerate(types):
        z = np.load(f"{LIB}/{t['id']}.npz")
        b = z["beads"].astype(np.float64); k = z["kind"].astype(np.int8)
        o = morton_order(b)
        w = raw[t["bo"]:t["bo"] + t["nb"]]
        assert len(o) == t["nb"] and np.array_equal(np.round(b[o] / BEAD_Q).astype(np.int16), np.c_[w["x"], w["y"], w["z"]]), t["id"]
        inv = np.empty(len(o), np.int64); inv[o] = np.arange(len(o))
        web = t["bo"] + inv                                       # npz bead -> global web bead
        web_of[ti] = web
        charges[web[k == 1]] = -1                                 # phosphates
        prot = np.where(k == 0)[0]
        stats["protein_beads"] += len(prot)
        if not len(prot):
            continue
        srcs = []
        if t["src"] == "AF" and t["uniprot"] in af:
            srcs = [af[t["uniprot"]]]
        elif t["src"] == "AF-multimer" and t.get("model"):
            srcs = [structure_table(p) for p in glob.glob(f"{DIMERS}/*/{t['model']}.pdb")[:1]]
        elif t["src"] == "PDB" and t["pdb"]:
            srcs = (structure_table(p) for p in pdb_files(t["pdb"]))
        idx = None
        for src in srcs:
            idx = match(b[prot], src)
            if idx is not None:
                break
        if idx is None:
            # composite structures (e.g. ribosome_80S_Sec61: a ribosome plus a translocon): the part taken from the
            # structure of a type whose id prefixes this one
            for tj, u in enumerate(types[:ti]):
                if t["id"].startswith(u["id"] + "_") and u["src"] == "PDB" and u["pdb"]:
                    for src in (structure_table(p) for p in pdb_files(u["pdb"])):
                        idx = match(b[prot], src, frac=0.4)
                        if idx is not None:
                            break
                if idx is not None:
                    break
        if idx is None:
            stats["types_unmatched"].append(t["id"]); continue
        ok = idx >= 0
        stats["protein_beads_matched"] += int(ok.sum())
        names = src["name"][idx[ok]]
        charges[web[prot[ok]]] = np.array([CHARGE.get(r, 0) for r in names], np.int8)
        if t["src"] == "AF" and not t["id"].endswith(("_oligo", "_het")):
            res_of[ti] = dict(zip(src["num"][idx[ok]].tolist(), web[prot[ok]].tolist()))
        if ti % 500 == 0:
            print(f"  {ti} types ({time.time() - t0:.0f} s)")
    charges.tofile(f"{a.data}/library/charges.bin")
    print(f"charges: {stats['protein_beads_matched']}/{stats['protein_beads']} protein beads matched; "
          f"{len(stats['types_unmatched'])} types unmatched: {stats['types_unmatched'][:30]}")
    # native contacts of confident predicted dimers whose partners are both monomer types
    mono = {}
    for ti, t in enumerate(types):
        if ti in res_of:
            mono.setdefault(t["uniprot"], ti)
    pairs, binv, n_rows, n_skip = [], [], 0, 0
    for sub in ("heterodimers", "homodimers"):
        cf = f"{DIMERS}/{sub}/candidates.csv"
        if not os.path.exists(cf):
            continue
        for _, r in pd.read_csv(cf).iterrows():
            if r["ipTM"] < MIN_IPTM:
                continue
            n_rows += 1
            ta, tb = mono.get(r["uniprotAccession_1"]), mono.get(r["uniprotAccession_2"])
            f = f"{DIMERS}/{sub}/{r['modelEntityId']}.pdb"
            if ta is None or tb is None or not os.path.exists(f):
                n_skip += 1; continue
            s = structure_table(f)
            ch = list(dict.fromkeys(s["chain"]))
            if len(ch) != 2:
                n_skip += 1; continue
            A, B = s["chain"] == ch[0], s["chain"] == ch[1]
            ia, ib = np.where(A)[0], np.where(B)[0]
            hits = cKDTree(s["xyz"][ib]).query_ball_point(s["xyz"][ia], r=CONTACT_NM)
            ca = []
            for p, lst in enumerate(hits):
                na = s["num"][ia[p]]
                if na not in res_of[ta]:
                    continue
                for q in lst:
                    nb_ = s["num"][ib[q]]
                    if nb_ in res_of[tb]:
                        ca.append((res_of[ta][na] - types[ta]["bo"], res_of[tb][nb_] - types[tb]["bo"]))
            if len(ca) < 3:
                n_skip += 1; continue
            pairs.append([ta, tb, len(binv), len(ca), r["modelEntityId"]]); binv.extend(ca)
    np.array(binv, np.uint16).reshape(-1, 2).tofile(f"{a.data}/library/contacts.bin")
    cn = [p[3] for p in pairs]
    json.dump(dict(cutoff_nm=CONTACT_NM, min_iptm=MIN_IPTM, pairs=pairs,
                   note="native inter-chain C-alpha contacts of AlphaFold DB yeast dimer models between library monomer types",
                   stats=dict(dimer_models=n_rows, used=len(pairs), skipped=n_skip, contacts=len(binv),
                              median_contacts=float(np.median(cn)) if cn else 0)),
              open(f"{a.data}/library/contacts.json", "w"), separators=(",", ":"))
    print(f"contacts: {len(pairs)} of {n_rows} confident dimer models used ({n_skip} without both monomer types or contacts), "
          f"{len(binv)} contacts, median {np.median(cn) if cn else 0:.0f} per interface; {time.time() - t0:.0f} s")
    rh_scale(types, a.data)
    import diffusion_model
    diffusion_model.write(a.data)
    import membrane_diffusion
    membrane_diffusion.write(a.data)
    return dict(protein_beads=stats["protein_beads"], protein_beads_matched=stats["protein_beads_matched"],
                types_unmatched=len(stats["types_unmatched"]), contact_pairs=len(pairs), contacts=len(binv))


def rh_scale(types, data_dir):
    """library/rh_scale.json: per AlphaFold monomer type, the factor by which its trimmed IDRs (pLDDT < 50) would
    enlarge its hydrodynamic radius: R_h(full chain, IDRs regrown as self-avoiding coils) / R_h(folded part), from
    work/features/protein_features.parquet (features.py). Types without features or factor < 1.02 are left out."""
    fp = f"{WORK}/features/protein_features.parquet"
    if not os.path.exists(fp):
        print("rh_scale: no feature table (run pipeline/features.py); skipped"); return
    F = pd.read_parquet(fp)
    f = (F.rh_full / F.rh_folded).dropna()
    out = {}
    for ti, t in enumerate(types):
        if t["src"] == "AF" and not t["id"].endswith(("_oligo", "_het")) and t["uniprot"] in f.index:
            v = float(f[t["uniprot"]])
            if v >= 1.02:
                out[str(ti)] = round(min(v, 3.0), 3)
    json.dump(dict(note="R_h(full chain with IDRs as coils) / R_h(folded part), AlphaFold monomer types (features.py)",
                   scale=out), open(f"{data_dir}/library/rh_scale.json", "w"), separators=(",", ":"))
    print(f"rh_scale: {len(out)} types with trimmed IDRs (median factor {np.median(list(out.values())) if out else 1:.2f})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=WEB_OUT)
    build(ap.parse_args().data)


if __name__ == "__main__":
    main()
