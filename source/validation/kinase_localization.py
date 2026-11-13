"""Curated localization test of signalling proteins in log phase (IMPROVEMENTS I12).

Each entry states where a protein (all placed instances containing it) should be in an unstressed log-phase cell, from
dedicated studies that are not model inputs (literature review 2026-11-11, Part D; citations there), with a minimum
share of placed copies. Negative controls state where it must NOT be (stress-only locations): at most the given share.
Locations are groups of packing regions:
  vac_face   vacuolar membrane or its cytosolic face        nucleus  nucleoplasm, nucleolus
  cytosol    cytosol and cytosolic shells                    pm       plasma membrane and its peripheral shell
  wall       cell wall
Output: work/validation/kinase_localization.json (table printed)
"""
import os, sys, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WORK, WEB_OUT
from contacts import proteins_of

GROUPS = {"vac_face": {"mem_VAC", "VAC_peripheral"}, "nucleus": {"nucleoplasm", "nucleolus"},
          "cytosol": {"cytosol", "PM_peripheral", "NE_peripheral", "ER_peripheral", "LD_peripheral", "VAC_peripheral",
                      "GOLGI_peripheral", "PEROX_peripheral"},
          "pm": {"mem_PM", "PM_peripheral"}, "wall": {"cell_wall"}}

# (protein, group, ">=" or "<=", share, kind, evidence)
TESTS = [
    ("KOG1", "vac_face", ">=", 0.8, "positive", "TORC1 on the vacuolar membrane (Berchtold & Walther 2009; Prouteau 2017)"),
    ("TCO89", "vac_face", ">=", 0.8, "positive", "TORC1 subunit, Rag tether (Nicastro 2025)"),
    ("SCH9", "vac_face", ">=", 0.3, "positive", "partly on the vacuolar surface via PI(3,5)P2 (Novarina 2021; Takeda 2018)"),
    ("BCY1", "nucleus", ">=", 0.8, "positive", "almost exclusively nuclear in glucose (Griffioen 2001; Tudisca 2010)"),
    ("TPK2", "nucleus", ">=", 0.6, "positive", "nuclear in glucose (Tudisca 2010)"),
    ("PHO85", "nucleus", ">=", 0.5, "positive", "Pho80-Pho85 nuclear (O'Neill 1996; Kaffman 1998)"),
    ("RPN1", "nucleus", ">=", 0.6, "positive", "26S proteasome mostly nuclear in growing cells (Laporte 2008)"),
    ("PHO4", "cytosol", ">=", 0.7, "positive", "cytoplasmic in phosphate-rich medium (O'Neill 1996)"),
    ("RIM15", "cytosol", ">=", 0.7, "positive", "cytoplasmic in log phase (Wanke 2005; Pedruzzi 2003)"),
    ("HOG1", "cytosol", ">=", 0.9, "positive", "cytoplasmic without osmotic stress (Ferrigno 1998)"),
    ("PBS2", "cytosol", ">=", 0.9, "positive", "cytoplasmic (Ferrigno 1998)"),
    ("GAL83", "cytosol", ">=", 0.8, "positive", "Snf1-Gal83 cytoplasmic in glucose (Vincent 2001)"),
    ("YCK2", "pm", ">=", 0.8, "positive", "palmitoylated, plasma membrane (Babu 2004)"),
    ("PKH1", "pm", ">=", 0.5, "positive", "at eisosomes (Walther 2007; Froehlich 2009)"),
    ("SLM1", "pm", ">=", 0.8, "positive", "eisosomes / plasma membrane (Kamble 2011; Berchtold 2012)"),
    ("HOG1", "nucleus", "<=", 0.1, "negative control", "nuclear only after osmotic stress"),
    ("PHO4", "nucleus", "<=", 0.2, "negative control", "nuclear only on phosphate starvation"),
    ("GAL83", "nucleus", "<=", 0.1, "negative control", "nuclear only on glucose limitation"),
    ("STE5", "pm", "<=", 0.05, "negative control", "at the plasma membrane only with pheromone (Mahanty 1999)"),
    ("FUS3", "wall", "<=", 0.0, "negative control", "a cytoplasmic kinase cannot be in the wall"),
    ("SSA1", "wall", "<=", 0.0, "negative control", "a cytosolic chaperone cannot be in the wall"),
]


def main():
    sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
    import tiles3d as t3
    Z = np.load(f"{WORK}/contacts/instances.npz")
    types = json.load(open(f"{WEB_OUT}/library/types.json"))["types"]
    names, off, mem = proteins_of(types)
    pidx = {n: k for k, n in enumerate(names)}
    rn = t3.REGION_NAMES
    cnt = np.zeros((len(types), len(rn)), np.int64)
    np.add.at(cnt, (Z["type"].astype(np.int64), Z["region"].astype(np.int64)), 1)
    out = []
    for prot, grp, op, thr, kind, ev in TESTS:
        ts = [t for t in range(len(types)) if pidx.get(prot) in set(mem[off[t]:off[t + 1]])]
        c = cnt[ts].sum(0)
        n = int(c.sum())
        share = float(sum(c[i] for i, r in enumerate(rn) if r in GROUPS[grp]) / n) if n else float("nan")
        ok = (share >= thr) if op == ">=" else (share <= thr)
        top = {rn[i]: round(float(c[i] / n), 3) for i in np.argsort(-c)[:3] if c[i] > 0} if n else {}
        out.append(dict(protein=prot, location=grp, rule=f"{op} {thr}", kind=kind, placed=n, share=round(share, 3),
                        passed=bool(ok), top_regions=top, evidence=ev))
        print(f"{prot:6s} {kind:16s} {grp:8s} {op}{thr:<4} share {share:5.2f}  {'PASS' if ok else 'FAIL'}  {top}")
    res = dict(tests=out, positive_passed=sum(o["passed"] for o in out if o["kind"] == "positive"),
               positive_total=sum(o["kind"] == "positive" for o in out),
               negative_passed=sum(o["passed"] for o in out if o["kind"] != "positive"),
               negative_total=sum(o["kind"] != "positive" for o in out))
    print({k: v for k, v in res.items() if k != "tests"})
    json.dump(res, open(f"{WORK}/validation/kinase_localization.json", "w"), indent=1)
    return res


if __name__ == "__main__":
    main()
