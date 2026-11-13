"""Per-gene ledger of where each protein's copies are placed, build by build (IMPROVEMENTS N7; after the provenance
ledger of the Jurkat cell model).

  python validation/ledger.py snapshot BUILD [TYPES_JSON]   # placed coarse shares of every gene -> work/ledger/BUILD.parquet
  python validation/ledger.py compare A B                   # what moved between two snapshots -> work/ledger/A_to_B.*

A snapshot holds, per ORF: inventory copies and mass, the assigned shares (inventory weights) and the placed shares
(placed copies per packing region, from the library's types.json, reduced to the coarse classes of
validation/localization.py). A comparison lists genes whose largest placed share changed by >= MOVE (0.1) and the
mass that moved between classes: per gene the change in placed shares times its placed mass, netted per gene, split
from losing to gaining classes in proportion to the gains (Jurkat ledger rule). Shifts are not attributed to rules
automatically: the Progress entries in IMPROVEMENTS.md name the changes of each build.
"""
import os, sys, json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline")); sys.path.insert(0, HERE)
from config import WORK, WEB_OUT  # noqa: E402

OUT = f"{WORK}/ledger"
MOVE = 0.1


def snapshot(build, types_json=None):
    import localization as L
    if types_json:                      # model_shares reads WEB_OUT; point it at an archived types.json
        L.WEB_OUT = os.path.dirname(os.path.dirname(os.path.abspath(types_json)))
        assert os.path.exists(f"{L.WEB_OUT}/library/types.json"), "TYPES_JSON must sit in <dir>/library/types.json"
    df = L.model_shares()
    import inventory as I
    inv, _ = I.build()
    m = inv.drop_duplicates("orf").set_index("orf")
    df["gene"] = m.gene.reindex(df.index).values
    df["mass_da"] = m.mass_da.reindex(df.index).values
    df["inv_copies"] = m.copies_per_cell.reindex(df.index).values
    os.makedirs(OUT, exist_ok=True)
    df.to_parquet(f"{OUT}/{build}.parquet")
    print(f"{build}: {len(df)} genes, {int((df.get('placed_copies', 0) > 0).sum())} with placed copies -> {OUT}/{build}.parquet")
    return df


def compare(a, b):
    import localization as L
    A = pd.read_parquet(f"{OUT}/{a}.parquet"); B = pd.read_parquet(f"{OUT}/{b}.parquet")
    cls = L.CLASSES
    pa = A.reindex(columns=[f"placed_{c}" for c in cls]).fillna(0.0)
    pb = B.reindex(columns=[f"placed_{c}" for c in cls]).fillna(0.0)
    common = pa.index.intersection(pb.index)
    pa, pb = pa.loc[common], pb.loc[common]
    pa.columns = pb.columns = cls
    has = (pa.sum(1) > 0) & (pb.sum(1) > 0)
    pa, pb = pa[has], pb[has]
    d = pb - pa
    mass_b = (B.loc[pb.index, "placed_copies"].fillna(0) * B.loc[pb.index, "mass_da"].fillna(0)) / 6.022e23 * 1e12   # pg
    moved = d.abs().max(1) >= MOVE
    rows = pd.DataFrame(dict(gene=B.loc[pb.index, "gene"], top_before=pa.idxmax(1), top_after=pb.idxmax(1),
                             max_shift=d.abs().max(1), placed_pg_after=mass_b))
    for c in cls:
        rows[f"d_{c}"] = d[c].round(3)
    rows = rows[moved].sort_values("placed_pg_after", ascending=False)
    # mass flow: losses to gains in proportion to the gains, per gene
    flow = pd.DataFrame(0.0, index=cls, columns=cls)
    for orf in d.index:
        dd = d.loc[orf].values * mass_b.loc[orf]
        gain = np.clip(dd, 0, None); loss = np.clip(-dd, 0, None)
        if gain.sum() <= 0:
            continue
        flow += np.outer(loss, gain / gain.sum())
    tag = f"{a}_to_{b}"
    rows.to_csv(f"{OUT}/{tag}_genes.tsv", sep="\t")
    flow.round(4).to_csv(f"{OUT}/{tag}_mass_flow_pg.tsv", sep="\t")
    summary = dict(genes_compared=int(has.sum()), genes_moved=int(moved.sum()), move_threshold=MOVE,
                   placed_pg_before=float((A.loc[pa.index, "placed_copies"].fillna(0) * A.loc[pa.index, "mass_da"].fillna(0)).sum() / 6.022e11),
                   placed_pg_after=float(mass_b.sum()),
                   largest_flows_pg={f"{i}->{j}": round(float(flow.loc[i, j]), 3) for i, j in
                                     sorted(((i, j) for i in cls for j in cls if i != j), key=lambda x: -flow.loc[x])[:8]})
    json.dump(summary, open(f"{OUT}/{tag}.json", "w"), indent=1)
    print(json.dumps(summary, indent=1))
    print(rows.head(15)[["gene", "top_before", "top_after", "max_shift", "placed_pg_after"]].to_string())
    return summary


if __name__ == "__main__":
    if sys.argv[1] == "snapshot":
        snapshot(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
    else:
        compare(sys.argv[2], sys.argv[3])
