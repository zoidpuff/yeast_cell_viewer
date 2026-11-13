"""Placement artefacts that the stage-4 numbers cannot see (stack review of 2026-11-13): checks, not comparisons with data.

Read from the packed tiles of the current build (tiles_raw, geometry, tile measures, structure library):

  species    placed copies / copies of the inventory (expected_totals) per (species, region). The packer drops what does
             not fit, and over-full compartments are thinned before packing; this shows who is missing: per region the placed share by copies and by volume, and how many species fall below half and
             below 0.8 of their region's share (species with >= MIN_REQ requested copies)
  seams      tiles are packed in 8 passes (tile parity). In fully cytosolic tiles: summed sphere volume of the packed
             molecules per unit volume, against the distance of the molecule centre from the nearest tile face, per pass.
             A flat profile is the aim; first-pass tiles over-fill the 5 nm at their faces and last-pass tiles under-fill it
  shells     membrane-face proteins (regions *_peripheral) whose centre is further from their membrane than the
             shell width + SHELL_TOL
  halos      number density of molecule centres against distance from each membrane on its cytosolic side, for
             cytosol-region and membrane-face molecules: an empty or depleted layer shows as a ratio below 1 between
             the layer just outside the shell and the cytosol further out

  python validation/placement_checks.py          # about 3 min
Output: work/validation/placement_checks.json, placement_checks.png
"""
import os, sys, re, json
for _v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
import pandas as pd
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WORK, W3D, RAW  # noqa: E402
import cell3d as c3, tiles3d as t3, wc3d  # noqa: E402

MIN_REQ = 100          # species with fewer requested copies in a region are not counted in the species statistics
SHELL_TOL = 4.0        # nm beyond the shell width that still counts as inside
N_SAMPLE = 400         # tiles sampled for shells and halos
EDGE = np.arange(0, 101, 5.0)       # nm from the tile face
HALO = np.arange(0, 64, 4.0)        # nm from the membrane
MEMS = dict(PM=c3.M_PM, ER=c3.M_ER, NE=c3.M_NE, VAC=c3.M_VAC, MITO_OM=c3.M_OM, GOLGI=c3.M_GOLGI, PEROX=c3.M_PEROX, LD=c3.M_LD)


@nb.njit(cache=True)
def _dist_lab(P, G, pores, pidx, PT, PQ, PC):
    out = np.zeros(c3.NMEM); d = np.empty((P.shape[0], c3.NMEM)); lab = np.empty(P.shape[0], np.int8)
    for i in range(P.shape[0]):
        lab[i] = c3.label_of(P[i, 0], P[i, 1], P[i, 2], G, pores, pidx, PT, PQ, PC, out)
        c3.mem_sdfs(P[i, 0], P[i, 1], P[i, 2], G, pores, pidx, PT, PQ, PC, out)
        for m in range(c3.NMEM):
            d[i, m] = abs(out[m])
    return d, lab


def load():
    cell = c3.Cell.load(f"{W3D}/geometry.json")
    meas = pd.read_parquet(f"{W3D}/tile_measure.parquet"); names = list(meas.index)
    tiles = np.array([[int(x) for x in re.findall(r"[+-]\d+", n)] for n in names])
    lib_ids = [str(x) for x in np.load(f"{W3D}/lib_ids.npy")]
    lp = wc3d.build_libpack(lib_ids)
    vol = np.array([(4 / 3 * np.pi * lp["STR"][a:b].astype(float) ** 3).sum() for a, b in zip(lp["st_off"][:-1], lp["st_off"][1:])])
    return cell, meas, names, tiles, lib_ids, vol


def species(names, lib_ids, vol):
    rows = []
    for n in names:
        f = f"{RAW}/{n}.npz"
        if not os.path.exists(f):
            continue
        z = np.load(f)
        rows.append(np.c_[z["job_lib"], z["job_region"], z["job_n"], z["job_placed"]])
    J = pd.DataFrame(np.concatenate(rows), columns=["lib", "region", "n", "placed"]).groupby(["lib", "region"], as_index=False).sum()
    # against the inventory, not the packing request: the request of an over-full compartment may have been thinned
    # (wc3d.pilot_capacity), and a species the packer was never asked to place in full must still show here
    fe = f"{W3D}/expected_totals.parquet"
    if os.path.exists(fe):
        e = pd.read_parquet(fe); ix = {k: i for i, k in enumerate(lib_ids)}; rx = {k: i for i, k in enumerate(t3.REGION_NAMES)}
        e = e[e.lib_id.isin(ix) & e.region.isin(rx)]
        e = pd.DataFrame(dict(lib=e.lib_id.map(ix).values, region=e.region.map(rx).values, expected=e.expected.values)).groupby(["lib", "region"], as_index=False).sum()
        J = J.merge(e, on=["lib", "region"], how="outer").fillna(0)
        J["n"] = np.maximum(J.expected.round(), J.placed); J = J[J.n > 0].drop(columns="expected")
        J["lib"] = J.lib.astype(int); J["region"] = J.region.astype(int)
    J["vol"] = vol[J.lib.values]
    out, worst = {}, []
    for r, g in J.groupby("region"):
        if g.n.sum() == 0:
            continue
        name = t3.REGION_NAMES[r]; share = g.placed.sum() / g.n.sum()
        b = g[g.n >= MIN_REQ]; fr = b.placed / b.n
        out[name] = dict(requested=int(g.n.sum()), placed_share=float(share),
                         placed_share_by_volume=float((g.placed * g.vol).sum() / (g.n * g.vol).sum()),
                         species=int(len(b)), species_below_half_of_region_share=int((fr < 0.5 * share).sum()),
                         species_below_0p8_of_region_share=int((fr < 0.8 * share).sum()),
                         species_share_p05_p50=[float(fr.quantile(0.05)), float(fr.median())] if len(b) else None)
        for _, x in b[fr < 0.5 * share].iterrows():
            worst.append((lib_ids[int(x.lib)], name, int(x.n), float(x.placed / x.n)))
    worst.sort(key=lambda w: -w[2])
    tot = J.groupby("lib")[["n", "placed"]].sum(); tot = tot[tot.n >= MIN_REQ]; f = tot.placed / tot.n
    summary = dict(placed_share=float(J.placed.sum() / J.n.sum()), placed_share_by_volume=float((J.placed * J.vol).sum() / (J.n * J.vol).sum()),
                   species=int(len(tot)), species_below_50pct=int((f < 0.5).sum()), species_below_80pct=int((f < 0.8).sum()),
                   copies_requested_of_species_below_50pct=int(tot.n[f < 0.5].sum()))
    return dict(summary=summary, regions=out, largest_species_below_half_of_region_share=[dict(species=w[0], region=w[1], requested=w[2], placed_share=round(w[3], 3)) for w in worst[:25]]), J


def seams(meas, names, tiles, vol):
    full = np.where((meas["cytosol"].values / t3.TILE ** 3) > 0.98)[0]
    phase = (tiles[:, 0] & 1) + 2 * (tiles[:, 1] & 1) + 4 * (tiles[:, 2] & 1)
    shell = np.array([(t3.TILE - 2 * a) ** 3 - (t3.TILE - 2 * b) ** 3 for a, b in zip(EDGE[:-1], EDGE[1:])])
    HV = np.zeros((8, len(EDGE) - 1)); NT = np.zeros(8)
    for ti in full:
        f = f"{RAW}/{names[ti]}.npz"
        if not os.path.exists(f):
            continue
        z = np.load(f); pos = z["pos"].astype(float) - tiles[ti] * t3.TILE
        d = np.minimum(pos, t3.TILE - pos).min(1)
        HV[phase[ti]] += np.histogram(d, EDGE, weights=vol[z["types"]])[0]; NT[phase[ti]] += 1
    prof = HV / np.maximum(NT[:, None], 1) / shell
    inner = prof[:, 3:7].mean(1)                       # 15-35 nm from the face
    return dict(tiles=int(NT.sum()), tiles_per_pass=NT.astype(int).tolist(), bin_nm=EDGE.tolist(),
                volume_fraction_by_pass=np.round(prof, 4).tolist(),
                face_layer_0_5nm=dict(first_pass=float(prof[0, 0]), last_pass=float(prof[7, 0]), all=float(HV[:, 0].sum() / NT.sum() / shell[0])),
                interior_15_35nm=dict(first_pass=float(inner[0]), last_pass=float(inner[7]), all=float(HV[:, 3:7].sum() / NT.sum() / shell[3:7].sum())),
                face_over_interior_range=[float((prof[:, 0] / inner).min()), float((prof[:, 0] / inner).max())]), prof


def shells_halos(cell, meas, names, tiles):
    rng = np.random.default_rng(1)
    peri = [c for c in meas.columns if c.endswith("_peripheral") or c.startswith("mem_")]
    cand = np.where((meas[peri].sum(1) > 0).values)[0]
    sel = rng.permutation(cand)[:N_SAMPLE]
    sh = {r: [] for r in t3.BANDS}
    cnt = {m: np.zeros((2, len(HALO) - 1)) for m in MEMS}; volh = {m: np.zeros(len(HALO) - 1) for m in MEMS}
    NQ = 40000
    for ti in sel:
        f = f"{RAW}/{names[ti]}.npz"
        if not os.path.exists(f):
            continue
        z = np.load(f); reg = np.repeat(z["job_region"], z["job_placed"]); pos = z["pos"].astype(np.float64)
        lo = tiles[ti].astype(float) * t3.TILE; a = cell.args(lo, lo + t3.TILE, margin=80.0)
        d, lab = _dist_lab(pos, *a)
        Q = lo + rng.uniform(0, t3.TILE, (NQ, 3)); dq, labq = _dist_lab(Q, *a)
        for r, (m, w) in t3.BANDS.items():
            k = reg == r
            if k.any():
                sh[r].append(d[k, m])
        soluble = reg == t3.R_CYTO; face = np.isin(reg, list(t3.BANDS))
        for name, mi in MEMS.items():
            other = np.delete(np.arange(c3.NMEM), mi)
            near = (lab == c3.CYTO) & (d[:, mi] <= d[:, other].min(1)); nearq = (labq == c3.CYTO) & (dq[:, mi] <= dq[:, other].min(1))
            cnt[name][0] += np.histogram(d[near & soluble, mi], HALO)[0]; cnt[name][1] += np.histogram(d[near & face, mi], HALO)[0]
            volh[name] += np.histogram(dq[nearq, mi], HALO)[0] * (t3.TILE ** 3 / NQ)
    shells = {}
    for r, (m, w) in t3.BANDS.items():
        if sh[r]:
            d = np.concatenate(sh[r])
            shells[t3.REGION_NAMES[r]] = dict(shell_nm=w, copies=int(len(d)), outside=float((d > w + SHELL_TOL).mean()),
                                              beyond_60nm=float((d > 60).mean()), median_nm=float(np.median(d)))
    halos, prof = {}, {}
    for name in MEMS:
        with np.errstate(divide="ignore", invalid="ignore"):
            dens = cnt[name] / volh[name] * 1e4                      # centres per 1e4 nm3
        if volh[name][9:].sum() < 5e6:                                # too little cytosol sampled next to this membrane
            continue
        tot = dens.sum(0); far = np.nanmean(tot[9:])                  # 36-60 nm
        halos[name] = dict(cytosol_region=np.round(np.nan_to_num(dens[0]), 2).tolist(), membrane_face=np.round(np.nan_to_num(dens[1]), 2).tolist(),
                           all_4_20nm_over_far=float(np.nanmean(tot[1:5]) / far), all_20_28nm_over_far=float(np.nanmean(tot[5:7]) / far),
                           cytosol_region_4_20nm_over_far=float(np.nanmean(dens[0, 1:5]) / np.nanmean(dens[0, 9:])))
        prof[name] = dens
    return shells, dict(bin_nm=HALO.tolist(), unit="centres per 1e4 nm3, cytosolic side", membranes=halos), prof


def figure(seam_prof, halo_prof, J, lib_ids, fn):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, muted, grid = "#1f2328", "#6b7280", "#e5e7eb"
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    x = (EDGE[:-1] + EDGE[1:]) / 2
    blues = plt.cm.Blues(np.linspace(0.35, 0.95, 8))
    for p in range(8):
        ax[0].plot(x, seam_prof[p], color=blues[p], lw=1.6)
    ax[0].annotate("pass 1", (x[0], seam_prof[0, 0]), xytext=(6, 0), textcoords="offset points", color=ink, fontsize=9, va="center")
    ax[0].annotate("pass 8", (x[0], seam_prof[7, 0]), xytext=(6, 0), textcoords="offset points", color=ink, fontsize=9, va="center")
    ax[0].set_xlabel("distance of the molecule centre from the nearest tile face (nm)"); ax[0].set_ylabel("sphere volume fraction")
    ax[0].set_title("Tile seams: cytosolic tiles, by packing pass (light to dark)", loc="left", fontsize=10, color=ink)
    ax[0].set_ylim(0, None)
    xh = (HALO[:-1] + HALO[1:]) / 2
    cols = dict(GOLGI="#b45309", VAC="#0f766e", ER="#6d28d9", MITO_OM="#be123c", PM="#1d4ed8", NE="#4d7c0f", PEROX="#64748b", LD="#a16207")
    for name, dens in halo_prof.items():
        ax[1].plot(xh, dens.sum(0), color=cols.get(name, muted), lw=1.6)
        ax[1].annotate(name, (xh[-1], dens.sum(0)[-1]), xytext=(4, 0), textcoords="offset points", color=ink, fontsize=8, va="center")
    ax[1].set_xlabel("distance from the membrane, cytosolic side (nm)"); ax[1].set_ylabel("molecule centres per 10$^4$ nm$^3$")
    ax[1].set_title("Layer next to each membrane (all molecules)", loc="left", fontsize=10, color=ink); ax[1].set_ylim(0, None)
    ax[1].set_xlim(0, HALO[-1] + 9)
    b = J[J.n >= MIN_REQ]
    for k, (name, col) in enumerate((("mito_matrix", "#be123c"), ("nucleolus", "#6d28d9"), ("GOLGI_peripheral", "#b45309"), ("cytosol", "#1d4ed8"))):
        g = b[b.region == t3.REGION_NAMES.index(name)]
        if len(g):
            ax[2].scatter(g.vol ** (1 / 3), g.placed / g.n, s=np.clip(np.sqrt(g.n) / 6, 4, 60), color=col, alpha=0.55, lw=0, label=name)
    ax[2].set_xlabel("molecule size: cube root of sphere volume (nm)"); ax[2].set_ylabel("placed / requested")
    ax[2].set_title("Who is dropped: one dot per species (area ~ copies)", loc="left", fontsize=10, color=ink)
    ax[2].legend(frameon=False, fontsize=8, loc="lower left"); ax[2].set_ylim(0, 1.03)
    for a in ax:
        a.grid(color=grid, lw=0.6); a.set_axisbelow(True)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        a.tick_params(colors=muted, labelsize=8); a.xaxis.label.set_color(muted); a.yaxis.label.set_color(muted)
        a.xaxis.label.set_size(9); a.yaxis.label.set_size(9)
    fig.tight_layout(); fig.savefig(fn, dpi=130); plt.close(fig)


def main():
    cell, meas, names, tiles, lib_ids, vol = load()
    sp, J = species(names, lib_ids, vol)
    se, seam_prof = seams(meas, names, tiles, vol)
    sh, ha, halo_prof = shells_halos(cell, meas, names, tiles)
    res = dict(species=sp, seams=se, shells=sh, halos=ha)
    if os.path.exists(f"{WORK}/build.json"):
        res["build_id"] = json.load(open(f"{WORK}/build.json")).get("build_id")
    os.makedirs(f"{WORK}/validation", exist_ok=True)
    json.dump(res, open(f"{WORK}/validation/placement_checks.json", "w"), indent=1)
    figure(seam_prof, halo_prof, J, lib_ids, f"{WORK}/validation/placement_checks.png")
    s = sp["summary"]
    print(f"placed {s['placed_share']:.3f} by copies, {s['placed_share_by_volume']:.3f} by volume; species below 50% placed: "
          f"{s['species_below_50pct']} of {s['species']} ({s['copies_requested_of_species_below_50pct']} requested copies); below 80%: {s['species_below_80pct']}")
    print(f"seams ({se['tiles']} cytosolic tiles): volume fraction within 5 nm of a tile face {se['face_layer_0_5nm']['first_pass']:.2f} (first pass) "
          f"to {se['face_layer_0_5nm']['last_pass']:.2f} (last), interior {se['interior_15_35nm']['all']:.2f}")
    for k, v in sh.items():
        print(f"shell {k:17s} {v['outside']:.3f} outside ({v['copies']} copies sampled), {v['beyond_60nm']:.3f} beyond 60 nm")
    for k, v in ha["membranes"].items():
        print(f"halo  {k:8s} all molecules 4-20 nm / far {v['all_4_20nm_over_far']:.2f}; 20-28 nm / far {v['all_20_28nm_over_far']:.2f}; "
              f"cytosol-region only 4-20 nm / far {v['cytosol_region_4_20nm_over_far']:.2f}")
    return res


if __name__ == "__main__":
    main()
