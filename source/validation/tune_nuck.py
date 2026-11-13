"""Fit inventory.NUC_K on the dev half of the SWAT (native promoter) comparison (assigned shares)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pipeline"))
import inventory as I
import localization as L
ref = {o: v for o, v in L.load_swat().items() if L.protein_split(o) == "dev"}
nop = {o: v for o, v in L.load_swat("NOP1pr-GFP localization ").items() if L.protein_split(o) == "dev"}
for k in [float(x) for x in sys.argv[1:]]:
    I.NUC_K = k
    M = L.model_shares()
    r, _ = L.score(M, ref, "assigned"); r2, _ = L.score(M, nop, "assigned")
    print(f"NUC_K={k}: SWAT_native top1 {r['top1']:.3f} macroF1 {r['macro_f1']:.3f} | NOP1 top1 {r2['top1']:.3f} macroF1 {r2['macro_f1']:.3f}", flush=True)
