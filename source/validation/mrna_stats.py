"""mRNA and polysomes of the packed cell (IMPROVEMENTS M1-M3): numbers to compare with measurements.

From prepass.npz: per polysome chain its ribosomes (poly_chain), its gene (poly_chain_gene, gene-specific builds),
and the mRNA threads (mrna_pts / mrna_off / mrna_kind: 1 = 3' tail leaving the first ribosome's entry, 2 = 5' tail
leaving the last ribosome's exit, 0 = free mRNA; one 3' and one 5' tail per chain, in chain order).
  ribosomes per mRNA (distribution), fraction of cytosolic ribosomes on mRNAs, mRNAs per cell
  5'-3' end-to-end distance of each mRNA (5' tail end to 3' tail end; free mRNAs first to last bead) against its
  length, the quantity smFISH measures for mRNP compaction
Output: work/validation/mrna_stats.json
"""
import os, sys, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WORK, W3D, WEB_OUT


def main():
    z = np.load(f"{W3D}/prepass.npz")
    chn = z["poly_chain"]; P, off, kind = z["mrna_pts"], z["mrna_off"], z["mrna_kind"]
    target = z["mrna_target"]
    n_chain = len(np.unique(chn))
    per = np.bincount(np.unique(chn, return_inverse=True)[1])
    tails = np.where(kind > 0)[0]; free = np.where(kind == 0)[0]
    e2e, length = [], []
    for c in range(n_chain):
        t3, t5 = tails[2 * c], tails[2 * c + 1]
        if off[t3 + 1] > off[t3] and off[t5 + 1] > off[t5]:
            e2e.append(np.linalg.norm(P[off[t5 + 1] - 1] - P[off[t3 + 1] - 1]))
            length.append(int(target[t3] + target[t5]))
    for f in free:
        if off[f + 1] - off[f] > 1:
            e2e.append(np.linalg.norm(P[off[f + 1] - 1] - P[off[f]])); length.append(int(target[f]))
    e2e, length = np.array(e2e), np.array(length)
    types = json.load(open(f"{WEB_OUT}/library/types.json"))["types"]
    Z = np.load(f"{WORK}/contacts/instances.npz")
    ribo = np.array([t["id"].startswith("ribosome_80S") for t in types])
    cyto80 = int((ribo[Z["type"].astype(int)] & (Z["region"] == 0)).sum())
    bins = [0, 500, 1000, 2000, 4000, 20000]
    res = dict(mrnas=int(n_chain + len(free)), polysome_mrnas=int(n_chain), free_mrnas=int(len(free)),
               ribosomes_on_mrna=int(len(chn)), fraction_cytosolic_80S_on_mrna=float(len(chn) / max(cyto80, 1)),
               ribosomes_per_polysome=dict(mean=float(per.mean()), median=float(np.median(per)),
                                           p10_p90=[float(np.percentile(per, 10)), float(np.percentile(per, 90))]),
               end_to_end_nm_by_length={f"{a}-{b} nt": dict(n=int(((length >= a) & (length < b)).sum()),
                                                            median=float(np.median(e2e[(length >= a) & (length < b)])) if ((length >= a) & (length < b)).any() else None)
                                        for a, b in zip(bins[:-1], bins[1:])},
               gene_specific="poly_chain_gene" in z.files)
    json.dump(res, open(f"{WORK}/validation/mrna_stats.json", "w"), indent=1)
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    main()
