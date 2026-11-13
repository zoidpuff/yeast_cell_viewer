"""Fit inventory.BOUND_F on the dev half of the SWAT (native promoter) comparison (assigned shares)."""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pipeline"))
import inventory as I
import membrane_anchor as MA
import localization as L
ref = {o: v for o, v in L.load_swat().items() if L.protein_split(o) == "dev"}
anc = MA.anchors()
sub = {o: v for o, v in ref.items() if o in anc}
print("anchored proteins in dev SWAT set:", len(sub), "of", len(anc), "anchored overall")
for fl, fd in itertools.product([0.0, 0.3, 0.5, 0.7, 0.9], [0.0, 0.2, 0.4, 0.6]):
    I.BOUND_F = {"lipid": fl, "domain": fd}
    M = L.model_shares()
    r, _ = L.score(M, ref, "assigned"); rs, _ = L.score(M, sub, "assigned")
    pc = r["per_class"]
    print(f"lipid {fl} domain {fd}: all top1 {r['top1']:.4f} macroF1 {r['macro_f1']:.4f} periphery F1 {pc['periphery']['f1']:.3f} golgi F1 {pc['golgi_endosome']['f1']:.3f} | anchored top1 {rs['top1']:.3f}", flush=True)
