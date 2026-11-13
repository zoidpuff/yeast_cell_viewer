"""Figures for model_validation_and_tuning.md -> docs/figures/*.png (static, light surface).
Reads the outputs of the validation scripts in work/; run them first.
Palette: the validated categorical order blue, orange, aqua, yellow; baselines in neutral gray; text in neutral inks.
"""
import os, sys, json, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WORK, W3D

FIG = os.path.join(HERE, "..", "docs", "figures")
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
GRAY, INK, INK2, GRID, SURF = "#a3a29c", "#0b0b0b", "#52514e", "#e6e5e0", "#fcfcfb"
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF, "font.size": 10,
                     "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
                     "axes.axisbelow": True, "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left"})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(f"{FIG}/{name}.png", dpi=150)
    plt.close(fig)
    print("wrote", name)


def localization():
    # the "before" bars come from the round-H evaluation, a file of the bundle's own work folder; isolated builds lack it
    fb = f"{WORK}/validation/h7_test_eval.json"
    base = json.load(open(fb)) if os.path.exists(fb) else None
    reb = json.load(open(f"{WORK}/validation/latest/localization_test.json"))
    pick = lambda name, pre: next(r for r in reb if r["name"] == name and r["prefix"] == pre and r["mode"] == "share")
    refs = [("SWAT_native", "SWAT, native promoter"), ("SWAT_NOP1", "SWAT, NOP1 promoter"), ("hyperLOPIT", "hyperLOPIT")]
    fig, ax = plt.subplots(figsize=(7, 3.4))
    y = np.arange(len(refs)); h = 0.26
    series = ([("before", GRAY, [base[f"baseline|test|{k}"]["top1"] for k, _ in refs])] if base else []) + [
              ("nucleus rule: inventory", BLUE, [pick(k, "assigned")["top1"] for k, _ in refs]),
              ("nucleus rule: placed copies", ORANGE, [pick(k, "placed")["top1"] for k, _ in refs])]
    for i, (lab, col, v) in enumerate(series):
        i += 3 - len(series)
        ax.barh(y + (i - 1) * h, np.array(v) * 100, h * 0.9, color=col, label=lab)
        for yy, vv in zip(y + (i - 1) * h, v):
            ax.text(vv * 100 + 0.5, yy, f"{vv * 100:.1f}", va="center", fontsize=8, color=INK2)
    ax.axvline(74, color=INK2, lw=1, ls="--")
    ax.set_yticks(y, [r for _, r in refs]); ax.invert_yaxis(); ax.set_xlim(50, 92); ax.set_ylim(2.5, -0.85)
    ax.text(73.6, -0.62, "reference ceiling: SWAT vs CYCLoPs agree 74%", fontsize=8, color=INK2, va="center", ha="right")
    ax.set_xlabel("top-1 agreement on held-out proteins (%)")
    ax.set_title("Localization vs independent references (test half)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.2), ncol=3, fontsize=8)
    save(fig, "localization_top1")
    # per-class F1, SWAT NOP1 and hyperLOPIT, placed copies
    fig, ax = plt.subplots(figsize=(7, 3.4))
    cls = ["cytosol", "nucleus", "nucleolus", "ER", "mitochondrion", "golgi_endosome", "vacuole", "periphery",
           "peroxisome", "lipid_droplet"]
    x = np.arange(len(cls)); w = 0.38
    for i, (k, lab, col) in enumerate([("SWAT_NOP1", "SWAT NOP1", BLUE), ("hyperLOPIT", "hyperLOPIT", ORANGE)]):
        pc = pick(k, "placed")["per_class"]
        v = [pc[c]["f1"] if pc[c]["n_ref"] > 0 else np.nan for c in cls]
        ax.bar(x + (i - 0.5) * w, v, w * 0.9, color=col, label=lab)
    ax.set_xticks(x, [c.replace("_", "/").replace("golgi/endosome", "Golgi/endo").replace("lipid/droplet", "lipid dropl.")
                      for c in cls], rotation=35, ha="right")
    ax.set_ylabel("F1"); ax.set_ylim(0, 1); ax.legend(fontsize=8)
    ax.set_title("Per-compartment agreement (placed copies, test half)")
    save(fig, "localization_per_class")


def diffusion():
    d = pd.read_csv(f"{WORK}/validation/diffusion_benchmark.tsv", sep="\t")
    pcol = "pred::+ crowding(R_h) from BD"
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.8))
    ax = axs[0]
    for split, col, mk in (("train", BLUE, "o"), ("test", ORANGE, "s")):
        m = d.split == split
        ax.scatter(d.loc[m, pcol], d.loc[m, "D"], s=36, color=col, marker=mk, edgecolor=SURF, lw=1.5,
                   label=f"{split} papers (n={m.sum()})", zorder=3)
    lim = [0.05, 40]
    ax.plot(lim, lim, color=INK2, lw=1); ax.fill_between(lim, [l / 2 for l in lim], [l * 2 for l in lim], color=GRID, alpha=0.6, lw=0)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("model D (µm²/s)"); ax.set_ylabel("measured D (µm²/s)")
    for _, r in d[d.molecule.str.contains("GEM|Ede1|Srv2", regex=True)].drop_duplicates("molecule").iterrows():
        ax.annotate(r.molecule.split(" ")[0], (r[pcol], r.D), xytext=(5, -3), textcoords="offset points", fontsize=7, color=INK2)
    ax.legend(fontsize=8, loc="upper left"); ax.set_title("Diffusion: model vs literature (band = 2×)")
    ax = axs[1]
    js = json.load(open(f"{WORK}/validation/diffusion_benchmark.json"))
    names = list(js)
    labels = ["Stokes–Einstein\n× constant", "+ crowding(R_h)\nfrom BD", "+ binding\n(rejected)"]
    def get(v, key):
        t = v.get("test", v)
        for k in (key, key.replace("rmse_", "")):
            if k in t:
                return t[k]
        return np.nan
    ratio = [js[n]["test"]["rmse_log10_ratio"] for n in names]; absr = [js[n]["test"]["rmse_log10_abs"] for n in names]
    x = np.arange(len(names)); w = 0.38
    ax.bar(x - w / 2, ratio, w * 0.9, color=BLUE, label="within-paper ratio")
    ax.bar(x + w / 2, absr, w * 0.9, color=ORANGE, label="absolute")
    for xx, v in zip(x - w / 2, ratio):
        ax.text(xx, v + 0.005, f"{v:.2f}", ha="center", fontsize=8, color=INK2)
    for xx, v in zip(x + w / 2, absr):
        ax.text(xx, v + 0.005, f"{v:.2f}", ha="center", fontsize=8, color=INK2)
    ax.set_xticks(x, labels, fontsize=8); ax.set_ylabel("RMSE on held-out papers (log10)")
    ax.legend(fontsize=8); ax.set_title("Held-out error")
    save(fig, "diffusion_benchmark")


def crowding():
    fig, ax = plt.subplots(figsize=(6, 3.6))
    from diffusion_model import CROWD
    for job, comp, col, lab in (("cyto120_mf", "cytosol", BLUE, "cytoplasm box (120 nm)"),
                                ("nuc100_mf", "nucleus", ORANGE, "nucleoplasm box (100 nm, fixed chromatin)")):
        s = json.load(open(f"{WORK}/sweeps/{job}.json"))
        rh, D, D0 = map(np.array, (s["rh"], s["D"], s["D0"]))
        ok = D > 0
        ax.scatter(rh[ok], np.log(D[ok] / D0[ok]), s=6, color=col, alpha=0.25, lw=0)
        a, k = CROWD[comp]; xx = np.linspace(1, 16, 50)
        ax.plot(xx, a + k * xx, color=col, lw=2, label=f"{lab}: {a:.2f} {k:+.3f}·R_h")
    g = json.load(open(f"{WORK}/sweeps/gem200_mf.json"))
    tr = np.array(g["tracer"]); D, D0 = np.array(g["D"])[tr], np.array(g["D0"])[tr]
    ax.scatter(np.full(tr.sum(), 20.0), np.log(D / D0), s=40, marker="D", color=AQUA, edgecolor=SURF, lw=1.5,
               label="40 nm GEM tracers", zorder=3)
    ax.set_xlabel("hydrodynamic radius R_h (nm)"); ax.set_ylabel("ln(D / D_water)")
    ax.set_ylim(-5, 0.5); ax.legend(fontsize=8, loc="lower left"); ax.set_title("Crowding from Brownian dynamics of packed boxes")
    save(fig, "crowding_bd")


def condensates():
    c = json.load(open(f"{WORK}/condensates/test_cytosol.json"))
    fc = f"{WORK}/condensates/calibration_nucleus.json"     # written by cond_calibrate.py, which the build does not run
    cal = json.load(open(fc)) if os.path.exists(fc) else None
    lp = c["log phase"]
    nuc = None if cal is None else max(r["auroc_nucleolus"] for r in cal) if isinstance(cal, list) else cal.get("auroc_nucleolus")
    items = ([("nucleolus members\n(calibration, nucleus)", nuc, GRAY)] if nuc is not None else []) + [("P-body members\n(held-out test)", lp["auroc_pb"], BLUE),
             ("P-body, RNA-binding\nlabels permuted (null)", float(np.mean(lp["null_pb_rbp_permuted"])), GRAY),
             ("stress-granule-only\n(negative control)", lp["auroc_sg_only"], ORANGE),
             ("glycolytic G-bodies\n(negative control)", lp["auroc_gbody"], ORANGE)]
    fig, ax = plt.subplots(figsize=(7, 3.4))
    x = np.arange(len(items))
    ax.bar(x, [v for _, v, _ in items], 0.6, color=[c_ for *_, c_ in items])
    for xx, (_, v, _) in zip(x, items):
        ax.text(xx, v + 0.01, f"{v:.2f}", ha="center", fontsize=8, color=INK2)
    ax.axhline(0.5, color=INK2, lw=1, ls="--"); ax.text(len(items) - 0.55, 0.51, "chance", fontsize=8, color=INK2, ha="right")
    ax.set_xticks(x, [i for i, *_ in items], fontsize=8); ax.set_ylim(0.4, 1); ax.set_ylabel("AUROC of enrichment")
    ax.set_title("Emergent condensates: controls fail (negative result)")
    save(fig, "condensates")


def proximity():
    js = json.load(open(f"{WORK}/validation/proximity_test_g10.json"))["sets"]
    sets = [("prox_label", "proximity labelling\n(BioID/TurboID/APEX)"), ("pca", "DHFR PCA\n(Tarassov 2008)"),
            ("xlms", "cross-linking MS"), ("kinase_invivo", "kinase→substrate\n(in vivo, curated)"),
            ("kinase_invitro", "kinase→substrate\n(in vitro arrays)")]
    scores = [("abundance", "abundance of B", GRAY), ("compartment", "shared compartments", BLUE),
              ("snapshot", "packed-cell snapshot", ORANGE), ("snapshot, no shared complex", "snapshot, complex partners removed", AQUA)]
    fig, ax = plt.subplots(figsize=(9, 3.8))
    x = np.arange(len(sets)); w = 0.2
    for i, (k, lab, col) in enumerate(scores):
        v = [js[s]["median_auroc"][k] for s, _ in sets]
        ci = [js[s]["median_auroc_ci95"].get(k, [np.nan, np.nan]) for s, _ in sets]
        err = np.array([[vv - c[0], c[1] - vv] for vv, c in zip(v, ci)]).T
        ax.bar(x + (i - 1.5) * w, v, w * 0.9, color=col, label=lab, yerr=err, error_kw=dict(lw=1, ecolor=INK2, capsize=0))
    ax.axhline(0.5, color=INK2, lw=1, ls="--")
    ax.set_xticks(x, [s for _, s in sets], fontsize=8); ax.set_ylim(0.4, 1.02); ax.set_ylabel("median per-bait AUROC")
    ax.legend(fontsize=8, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.2)); ax.set_title("Does model proximity predict measured proximity? (10 nm, 95% CI)")
    save(fig, "proximity_auroc")


def kinases():
    t = pd.read_parquet(f"{WORK}/contacts/encounters_kinases.parquet", columns=["A", "B", "p_near_10", "c_comp_uM", "k_enc_per_s", "c_local_uM", "R_A", "R_B"])
    s = pd.read_csv(f"{WORK}/validation/proximity_sets.tsv", sep="\t", dtype=str)
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.6))
    ax = axs[0]
    bg = t.c_comp_uM.values
    groups = [("all kinase–protein pairs", GRAY, t)]
    for name, lab, col in (("kinase_invivo", "in vivo (curated)", BLUE), ("kinase_invitro", "in vitro arrays", ORANGE)):
        p = set(zip(s.loc[s.set == name, "A"], s.loc[s.set == name, "B"]))
        groups.append((lab, col, t[[(a, b) in p for a, b in zip(t.A, t.B)]]))
    for lab, col, g in groups:
        v = np.sort(np.log10(np.clip(g.c_comp_uM.values, 1e-4, None)))
        ax.plot(v, np.linspace(0, 1, len(v)), color=col, lw=2, label=f"{lab} (n={len(v):,})")
    ax.set_xlabel("log10 concentration of B in A's compartments (µM)"); ax.set_ylabel("cumulative fraction")
    ax.legend(fontsize=8, loc="upper left"); ax.set_title("Kinase–substrate pairs share compartments")
    ax = axs[1]
    m = (t.c_comp_uM > 0) & np.isfinite(t.R_A) & np.isfinite(t.R_B)
    tt = t[m].sample(min(200000, m.sum()), random_state=0)
    Rs = tt.R_A + tt.R_B
    shell = 4 / 3 * np.pi * ((Rs + 10) ** 3 - Rs ** 3)
    pred = 1 - np.exp(-tt.c_comp_uM * 6.02214e-7 * shell)
    hb = ax.hexbin(pred, tt.p_near_10, gridsize=45, bins="log", cmap="Blues", mincnt=1, extent=(0, 1, 0, 1))
    ax.plot([0, 1], [0, 1], color=INK2, lw=1)
    r = np.corrcoef(pred, tt.p_near_10)[0, 1]
    ax.text(0.03, 0.93, f"Pearson r = {r:.2f}", fontsize=8, color=INK2, transform=ax.transAxes)
    ax.set_xlabel("well-mixed estimate 1 − exp(−c·V_shell)"); ax.set_ylabel("snapshot P(B within 10 nm)")
    ax.set_title("Mean-field vs packed cell")
    save(fig, "kinase_proximity")


def packing():
    v = json.load(open(f"{WORK}/results/wholecell3d_validation.json"))
    reg = pd.DataFrame(v["regions"]).set_index("region")
    reg = reg[reg.requested > 0].sort_values("frac_placed")
    fig, ax = plt.subplots(figsize=(7, 3.8))
    col = [ORANGE if f < 0.9 else BLUE for f in reg.frac_placed]
    ax.barh(range(len(reg)), reg.frac_placed * 100, 0.7, color=col)
    for i, f in enumerate(reg.frac_placed):
        ax.text(f * 100 + 0.5, i, f"{f * 100:.0f}%", va="center", fontsize=7, color=INK2)
    ax.set_yticks(range(len(reg)), [r.replace("_", " ") for r in reg.index], fontsize=8)
    ax.set_xlim(0, 108); ax.set_xlabel("placed / requested copies (%)")
    ax.set_title("Packing: copies placed per region (orange < 90%)")
    save(fig, "packing_placed")


def copy_images():
    for f in ("gallery_midplane.png", "mip_punctate.png", "spot_counts.png"):
        p = f"{WORK}/microscopy/{f}"
        if os.path.exists(p):
            shutil.copy(p, f"{FIG}/microscopy_{f}")


def idr0078():
    res = json.load(open(f"{WORK}/validation/microscopy_idr0078.json"))
    names = {"SAC6": ("Sac6", "actin patches"), "SLA1": ("Sla1", "endocytic sites"), "SNF7": ("Snf7", "late endosomes"),
             "VPH1": ("Vph1", "vacuoles")}
    fig, ax = plt.subplots(figsize=(7, 3.8))
    for i, (g, v) in enumerate(res.items()):
        ax.plot([i - 0.18, i - 0.18], v["real_p10_p90"], color=GRAY, lw=6, solid_capstyle="butt")
        ax.plot(i - 0.18, v["real_median"], "o", color=INK, ms=6)
        ax.plot(i + 0.05 + np.linspace(-0.06, 0.06, 9), v["model_objects"], "s", color=BLUE, ms=5)
        ax.plot(i + 0.25 + np.linspace(-0.06, 0.06, 9), v["model_spots"], "D", color=ORANGE, ms=5)
    ax.plot([], [], color=GRAY, lw=6, label="screen, wild-type cells: 10th–90th percentile (dot: median)")
    ax.plot([], [], "s", color=BLUE, label="model: threshold objects (9 views)")
    ax.plot([], [], "D", color=ORANGE, label="model: contrast spots (9 views)")
    ax.set_xticks(range(len(res)), [f"{names[g][0]}\n({names[g][1]})" for g in res])
    ax.set_ylabel("objects per cell (max projection)")
    ax.legend(fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=2)
    ax.set_title("Foci per cell: simulated vs an image screen (Mattiazzi Usaj 2020)")
    save(fig, "microscopy_idr0078")


def _saved(name, key, default):
    """a number from a validation JSON of the current build, else the value of the last build that had it"""
    try:
        return json.load(open(f"{WORK}/validation/{name}"))[key]
    except (OSError, KeyError):
        return default


# model value, literature low, high, label, kind (input = the model was built from or set to this value)
LIT = [
    (3.5, 4.0, 6.0, "protein per cell (pg; cultures 4–6, model: one G1 cell)", "input"),
    (211500, 1.87e5, 2.0e5, "80S ribosomes per cell (from rRNA)", "check"),
    (17.0, 23.0, 23.0, "cytosolic ribosomes (µM, cryo-ET)", "test"),
    (28.6, 29.4, 30.7, "ribosome nearest-neighbour distance (nm, cryo-ET)", "test"),
    (round(_saved("er_ribosome_density.json", "per_um2", 782)), 550, 900, "ribosomes per µm² on cortical ER (cryo-ET)", "partly tuned"),
    (7.0, 5.9, 10.5, "nucleus (% of cell volume; incl. soft X-ray 5.9)", "input"),
    (17.1, 17.2, 17.2, "vacuole (% of cell volume; soft X-ray G1 cells)", "input"),
    (0.34, 0.3, 1.0, "lipid droplets (% of cell volume)", "input"),
    (2.9, 1.0, 4.77, "mitochondria (% of cell volume; 4.8 in glucose)", "check"),
    (248, 200, 225, "cytoplasmic dry mass incl. estimates (mg/mL)", "check"),
    (68.2, 92, 92, "peroxisome lumen protein copies (thousands)", "check"),
    (22, 30, 30, "eisosomes per cell", "input"),
    (333, 2500, 2500, "Pil1/Lsp1 per eisosome coat (implied by copies)", "check"),
    (1500, 3000, 3000, "actin per endocytic patch at peak", "assumption"),
    (5, 5, 10, "peroxisomes per cell (simulated foci)", "test"),
    (5507, 7715, 7715, "translocon-engaged ER ribosomes per cell (cryo-ET)", "test"),
    (51, 200, 200, "TORC1 per cell (GFP calibration)", "test"),
    (12.7, 13.6, 13.6, "free GFP diffusion, held-out paper (µm²/s)", "test"),
    (0.66, 0.90, 0.90, "nucleoplasm / cytosol density (soft X-ray)", "test"),
    (0.90, 1.0, 1.0, "core histones per nucleosome need", "check"),
    (1.36, 1.0, 1.0, "nucleoporins per pore need (median)", "check"),
]


def literature():
    fig, ax = plt.subplots(figsize=(8.5, 6.6))
    col = {"test": ORANGE, "check": BLUE, "input": GRAY, "assumption": GRAY, "partly tuned": AQUA}
    for i, (m, lo, hi, lab, kind) in enumerate(LIT):
        ax.plot([np.log2(lo / m), np.log2(hi / m)], [i, i], color=GRID, lw=7, solid_capstyle="round")
        ax.plot(0, i, "o", color=col[kind], ms=7, zorder=3)
        mid = np.sqrt(lo * hi)
        ax.plot(np.log2(mid / m), i, "|", color=INK2, ms=10, mew=1.5)
    ax.axvline(0, color=INK2, lw=1)
    ax.set_yticks(range(len(LIT)), [f"{l}  [{k}]" for _, _, _, l, k in LIT], fontsize=8)
    ax.invert_yaxis()
    t = [-3, -2, -1, 0, 1, 2, 3]
    ax.set_xticks(t, [f"{2.0 ** v:g}×" for v in t]); ax.set_xlim(-3.6, 3.6)
    ax.set_xlabel("literature ÷ model (band: reported range, tick: centre)")
    ax.set_title("The model's quantities against literature values")
    for k, lab in (("test", "independent test"), ("check", "consistency check"), ("partly tuned", "partly tuned"), ("input", "input or assumption")):
        ax.plot([], [], "o", color=col[k], label=lab)
    ax.legend(fontsize=8, loc="lower left")
    save(fig, "literature_numbers")


def chromatin():
    """R(s) against the Bystricky 2004 worm-like-chain range, and Hi-C P(s) of the model against G1 and log-phase data"""
    cs = json.load(open(f"{WORK}/validation/chromatin_stats.json"))
    hc = json.load(open(f"{WORK}/validation/hic_test.json"))
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.7))
    ax = axs[0]
    s = np.array([float(k) for k in cs["rms_distance_nm"]]); r = np.array(list(cs["rms_distance_nm"].values()))
    band = cs["bystricky2004_wlc_rms_range_nm"]
    bs = np.array([float(k) for k in band]); lo = np.array([v[0] for v in band.values()]); hi = np.array([v[1] for v in band.values()])
    ax.fill_between(bs, lo, hi, color=GRAY, alpha=0.35, lw=0, label="worm-like chain, Bystricky 2004 range")
    if "arbona2017_wlc_rms_range_nm" in cs:
        ab = cs["arbona2017_wlc_rms_range_nm"]; xs = np.array([float(k) for k in ab])
        ax.fill_between(xs, [v[0] for v in ab.values()], [v[1] for v in ab.values()], color=YELLOW, alpha=0.3, lw=0,
                        label="worm-like chain, Arbona 2017 range")
    ax.plot(s, r, "o-", color=BLUE, lw=2, ms=6, label="model")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("genomic distance (kb)"); ax.set_ylabel("rms 3D distance (nm)")
    ax.set_title("Chromatin compaction"); ax.legend(fontsize=8, loc="upper left")
    ax = axs[1]
    for k, lab, col in (("model", "model (nucleosome pairs < 30 nm)", BLUE), ("lazarstefanita2017_G1", "Hi-C, G1 (Lazar-Stefanita 2017)", ORANGE),
                        ("swygert2019_log", "Micro-C XL, log (Swygert 2019)", AQUA)):
        ps = hc[k]["ps"]; x = np.array(ps["s_bp"]) / 1000; y = np.array(ps["p"], float)
        ok = (x >= 10) & (x <= 200) & np.isfinite(y) & (y > 0)
        y = y / y[ok][0]
        ax.plot(x[ok], y[ok], "o-", color=col, lw=2, ms=5, label=f"{lab}, slope {hc[k]['ps_slope']:.2f}")
        if k + "_ice" in hc:          # the same after ICE balancing (hic_test.ice), dashed
            ps = hc[k + "_ice"]["ps"]; y2 = np.array(ps["p"], float); ok2 = (x >= 10) & (x <= 200) & np.isfinite(y2) & (y2 > 0)
            ax.plot(x[ok2], y2[ok2] / y2[ok2][0], "--", color=col, lw=1.5, label=f"  balanced, slope {hc[k + '_ice']['ps_slope']:.2f}")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("genomic distance (kb)"); ax.set_ylabel("contact frequency (relative)")
    ax.set_title("Contact probability P(s), 10-200 kb"); ax.legend(fontsize=7.5, loc="lower left")
    save(fig, "chromatin")


def mrna_localization():
    f = f"{WORK}/validation/mrna_localization.json"
    if not os.path.exists(f):
        return
    js = json.load(open(f))
    er = f"{WORK}/validation/mrna_localization_er.tsv"
    if os.path.exists(er) and "er_spearman" in js:
        fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 3.4), gridspec_kw=dict(width_ratios=[1.3, 1]))
        d = pd.read_csv(er, sep="\t")
        d["bin"] = pd.qcut(d.chartron_log2_membrane, 10, labels=False, duplicates="drop")
        g = d.groupby("bin").agg(x=("chartron_log2_membrane", "median"), y=("share", "mean"))
        ax2.plot(g.x, 100 * g.y, "o-", color=ORANGE, lw=2, ms=8)
        ax2.set_xlabel("membrane enrichment, log2 (Chartron 2016)")
        ax2.set_ylabel("model: mRNAs at the ER (%)")
        ax2.set_title("mRNAs at the ER")
        ax2.text(0.03, 0.95, f"Spearman {js['er_spearman']:.2f}\nAUROC top decile {js['er_auroc_top_decile']:.2f}\n{js['er_genes_compared']} genes",
                 transform=ax2.transAxes, ha="left", va="top", fontsize=9, color=INK2)
    else:
        fig, ax = plt.subplots(figsize=(7, 3.4))
    labs = ["model: tethered genes", "model: other genes", "MLR class I", "MLR class III-IV"]
    v = [100 * js["model_share_median_tethered_genes"], 100 * js["model_share_median_other"], js["mlr_median_class1"], js["mlr_median_class3_4"]]
    ax.barh(labs[::-1], v[::-1], color=[BLUE, BLUE, ORANGE, ORANGE], height=0.6)       # data blue, model orange
    for i, x in enumerate(v[::-1]):
        ax.text(x + 1, i, f"{x:.0f}%", va="center", fontsize=9, color=INK2)
    ax.set_xlabel("median share of mRNAs at mitochondria (%)")
    ax.set_title("mRNAs at the mitochondrial surface")
    ax.text(0.98, 0.5, f"model vs MLR %:\nSpearman {js['spearman_mlr']:.2f}\nAUROC class I {js['auroc_class1_vs_3_4']:.2f}",
            transform=ax.transAxes, ha="right", va="center", fontsize=9, color=INK2)
    save(fig, "mrna_localization")


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    for f in (localization, diffusion, crowding, condensates, proximity, kinases, packing, idr0078, literature, copy_images,
              chromatin, mrna_localization):
        try:
            f()
        except Exception as e:
            print("FAILED", f.__name__, repr(e))
