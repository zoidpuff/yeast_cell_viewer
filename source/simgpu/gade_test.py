"""Polysomes and the diffusion of 20 and 40 nm particles (IMPROVEMENTS P10/Q4; test against Gade et al. 2025 Cell Rep
44:116204, preprint 10.1101/2024.11.14.623620, Fig. 2; nothing here was fitted to it).

Measured in budding yeast: disassembling polysomes with hippuristanol (ribosomes run off, mRNAs stay) raises the
diffusion coefficient of 40 nm GEMs about 1.5x (0.135 -> 0.20 and 0.21 -> 0.29 um^2/s) and of 20 nm GEMs about 1.3x
(0.19 -> 0.25); cycloheximide (polysomes frozen) changes nothing; stm1-delta (20-25% fewer ribosomes, by the paper's
account idle ones) changes nothing.
Simulated: boxes of 200 nm cut from the packed cytosol (about 7,500 bodies, 85 ribosomes each), four GEM tracers each at the same
positions in every condition, three seeds, mean-field hydrodynamics, 10 us each:
  wt     neighbouring ribosomes of one mRNA tied by entropic springs (engine links; the mRNA itself is not a body)
  hip    no links: every ribosome diffuses on its own (how all earlier runs were done)
  stm1   links, and 22% of the box's ribosomes removed, taken from those not on an mRNA (all of them where fewer)
A single tracer's 10 us give its D only to about 25%, so the time-averaged MSDs of all tracers of a condition
(72: 6 boxes x 4 x 3 seeds) are pooled: D = MSD(lag) / (6 lag) at a lag of 300 ns (1 us also reported), with a
bootstrap over tracers. The lags are far shorter than the 10-100 ms of particle tracking, so only ratios between
conditions are compared.

  python simgpu/gade_test.py centres                  # choose and save the box centres
  python simgpu/gade_test.py run GPU JOB [JOB ...]    # JOB = <centre index>:<condition>:<gem nm>:<seed>
  python simgpu/gade_test.py summary
Output: work/sweeps/gade/*.json, work/validation/gade_polysomes.json
"""
import os, sys, json, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline")); sys.path.insert(0, HERE)
from config import WORK, LIB  # noqa: E402

OUT = f"{WORK}/sweeps/gade"
SIZE, TOTAL_NS, N_CENTRES, STM1_REMOVE = 200.0, 10000.0, 6, 0.22
LAGS_NS, SEEDS = (100, 300, 1000), (1, 2, 3)
GEM = {40: (18.5, 3.2, 20.0), 20: (8.5, 2.2, 10.0)}          # shell radius, bead radius, hydrodynamic radius (nm)
GEM_AT = ([0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.25, 0.75, 0.75], [0.75, 0.25, 0.75])


def centres():
    """box centres with cytosol only (< 1% of bodies elsewhere, < 2% fixed) and >= 30 polysome ribosomes, at random (seed 5)"""
    import boxes as bx
    lib = bx.Library(); rng = np.random.default_rng(5); out = []
    while len(out) < N_CENTRES:
        c = rng.uniform(-1900, 1900, 3)
        try:
            b = bx.cut_box(lib, c, SIZE)
        except Exception:
            continue
        if b["n"] < 5000 or (b["region"] != 0).mean() > 0.01 or (~b["mobile"]).mean() > 0.02 or (b["poly_tx"] >= 0).sum() < 30:
            continue
        if all(np.linalg.norm(c - o) > 400 for o in out):
            out.append(c.round(1)); print(len(out), c.round(1), b["n"], int((b["poly_tx"] >= 0).sum()), round(b["phi"], 3), flush=True)
    os.makedirs(OUT, exist_ok=True)
    json.dump([list(map(float, c)) for c in out], open(f"{OUT}/centres.json", "w"))


def run(job):
    import boxes as bx
    from engine import BD
    from sweep import gem_shell
    ci, cond, nm, seed = job.split(":"); ci, nm, seed = int(ci), int(nm), int(seed)
    lib = bx.Library()
    c = json.load(open(f"{OUT}/centres.json"))[ci]
    box = bx.cut_box(lib, c, SIZE)
    ids = np.array([lib.types[t]["id"] for t in box["type"]])
    ribo = np.char.startswith(ids, "ribosome_80S")
    if cond == "stm1":
        idle = np.where(ribo & (box["poly_tx"] < 0))[0]
        drop = np.random.default_rng(11).choice(idle, min(len(idle), int(round(STM1_REMOVE * ribo.sum()))), replace=False)
        keep = np.ones(box["n"], bool); keep[drop] = False
        box = bx.subset(box, keep)
    R, r, rh = GEM[nm]
    sh = gem_shell(R=R, r=r)
    extra = [(sh, rh, np.array(p) * SIZE) for p in GEM_AT]
    links = None if cond == "hip" else np.load(f"{LIB}/ribosome_80S_mrna_channel.npy")
    sim = BD(box, lib, hydro="mf", seed=seed, extra_bodies=extra, links=links)
    t0 = time.time()
    dtf = 10.0
    times, tr = sim.run(TOTAL_NS, dtf)
    tracer = np.where(sim.tracer)[0]
    rib = np.where(np.array([t >= 0 and lib.types[t]["id"].startswith("ribosome_80S") for t in sim.types]))[0]
    msd = lambda idx: [((tr[int(l / dtf):, idx] - tr[:-int(l / dtf), idx]) ** 2).sum(-1).mean(0).tolist() for l in LAGS_NS]
    res = dict(job=job, n=int(sim.n), ribosomes=int(len(rib)), ribosomes_removed=int(ribo.sum() - len(rib)) if cond == "stm1" else 0,
               polysome_ribosomes=int((box["poly_tx"] >= 0).sum()), bonds=int(getattr(sim, "n_bonds", 0)), phi=box["phi"],
               lags_ns=list(LAGS_NS), gem_msd_nm2=msd(tracer), ribosome_msd_nm2=np.mean(msd(rib), axis=1).tolist(), wall_s=time.time() - t0)
    json.dump(res, open(f"{OUT}/{job.replace(':', '_')}.json", "w"))
    print(job, "bonds", res["bonds"], "GEM D(300 ns) um2/s", (np.array(res["gem_msd_nm2"][1]) / (6 * LAGS_NS[1]) * 1e3).round(2).tolist(),
          f"{res['wall_s']:.0f}s", flush=True)


def summary():
    import glob
    R = {}
    for f in glob.glob(f"{OUT}/[0-9]*_*_*_*.json"):
        r = json.load(open(f)); ci, cond, nm, seed = r["job"].split(":"); R.setdefault((cond, int(nm)), []).append(r)
    out = dict(box_nm=SIZE, total_ns=TOTAL_NS, lags_ns=list(LAGS_NS),
               measured=dict(hip_over_wt_40nm=[1.5, 1.4], hip_over_wt_20nm=1.3, stm1_over_wt_40nm=1.0))
    rng = np.random.default_rng(0)
    M = {k: np.concatenate([np.array(r["gem_msd_nm2"]) for r in v], axis=1) for k, v in R.items()}     # (lags, tracers)
    D = lambda m, li: m[li].mean() / (6 * LAGS_NS[li]) * 1e3
    for nm in (40, 20):
        if ("wt", nm) not in M:
            continue
        w = M[("wt", nm)]
        out[f"wt_{nm}nm"] = dict(tracers=int(w.shape[1]), D_um2_s={str(l): float(D(w, i)) for i, l in enumerate(LAGS_NS)})
        for cond in ("hip", "stm1"):
            if (cond, nm) not in M:
                continue
            a = M[(cond, nm)]; res = dict(tracers=int(a.shape[1]))
            for li, l in enumerate(LAGS_NS):
                bs = [D(a[:, rng.integers(0, a.shape[1], a.shape[1])], li) / D(w[:, rng.integers(0, w.shape[1], w.shape[1])], li) for _ in range(2000)]
                res[f"ratio_lag_{l}ns"] = dict(ratio=float(D(a, li) / D(w, li)), ci95=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))])
            out[f"{cond}_over_wt_{nm}nm"] = res
    rb = {c: np.mean([r["ribosome_msd_nm2"][1] for r in R[(c, 40)]]) / (6 * LAGS_NS[1]) * 1e3 for c in ("wt", "hip", "stm1") if (c, 40) in R}
    out["ribosome_D_um2_s_lag_300ns"] = {c: float(v) for c, v in rb.items()}
    out["bonds_per_box"] = sorted({r["job"].split(":")[0]: r["bonds"] for r in R.get(("wt", 40), [])}.items())
    out["stm1_ribosomes_removed_fraction"] = float(np.mean([r["ribosomes_removed"] / (r["ribosomes"] + r["ribosomes_removed"]) for r in R.get(("stm1", 40), [])])) if ("stm1", 40) in R else None
    json.dump(out, open(f"{WORK}/validation/gade_polysomes.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "centres":
        centres()
    elif sys.argv[1] == "run":
        os.environ["CUDA_VISIBLE_DEVICES"] = sys.argv[2]
        for j in sys.argv[3:]:
            run(j)
    else:
        summary()
