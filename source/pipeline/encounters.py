"""Pair statistics for kinetic models (IMPROVEMENTS I5): how much of protein B does a molecule of protein A see, and how
often does it meet one?

Per ordered pair (A, B):
  p_near_2, p_near_10   fraction of A copies with a B copy within a 2 / 10 nm surface gap in the packed cell
                        (contacts.py; separate instances; p_same_complex adds B in the same complex instance)
  c_comp_uM             compartment-level concentration of B particles (placed instances containing B: an oligomer or
                        complex counts once) around A: sum_c w_A(c) N_B(c) / V_c (inventory weights,
                        built volumes; membranes as area x 10 nm plus peripheral shell). This is the concentration a
                        well-mixed kinetic model of each compartment would use.
  c_local_uM            concentration of B in the 10 nm shell around A in the packed cell: mean number of B neighbours
                        within 10 nm / shell volume (spheres of the two hydrodynamic radii). c_local / c_comp > 1 means B
                        is enriched next to A beyond compartment mixing (shared membrane, scaffold, polysome, shell).
  k_enc_per_s           diffusion-limited encounter rate of one A molecule with any B (Smoluchowski):
                        4 pi (D_A + D_B) (R_A + R_B) c_comp, with in-vivo D of the calibrated diffusion model (MODEL.md
                        14.5; soluble species in the cytosol and nucleus only; NaN otherwise) and R = hydrodynamic radii.
                        An upper bound on productive encounters (no orientation or reaction-probability factor).
Units: uM = umol/L; 1 uM = 6.022e-7 molecules/nm^3.

  python encounters.py [--a GENE ...] [--kinases]   -> work/contacts/encounters_*.tsv
"""
import os, sys, json, argparse
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "simgpu")); sys.path.insert(0, os.path.join(HERE, "..", "validation"))
from config import WORK, WEB_OUT
from contacts import Proximity, OUT

UM = 6.02214e-7          # molecules / nm^3 per uM


class Encounters:
    def __init__(self):
        import boxes
        from proximity_test import ZONES, zone_volumes
        self.P = P = Proximity()
        lib = boxes.Library(WEB_OUT)
        F = pd.read_parquet(f"{WORK}/features/protein_features.parquet")
        key = [g if isinstance(g, str) and g else o for o, g in zip(F.orf, F.gene)]
        F = F.assign(key=key).drop_duplicates("key").set_index("key")
        self.names = [n for n in P.proteins if n in F.index]
        self.loc = {n: k for k, n in enumerate(self.names)}
        self.idx = np.array([P.pidx[n] for n in self.names])
        zones = list(ZONES); V = zone_volumes()
        W = np.array([[json.loads(w).get(z, 0.0) for z in zones] for w in F.loc[self.names, "comp_weights"]])
        self.W = W / np.maximum(W.sum(1, keepdims=True), 1e-12)
        self.copies = F.loc[self.names, "copies_per_cell"].values.astype(float)        # chains per cell (inventory)
        self.inst = P.copies[self.idx].astype(float)                                   # placed particles containing it
        self.conc = self.W * self.inst[:, None] / np.array([V[z] for z in zones])[None, :]
        # diffusing unit of each protein: the library type holding most of its instances
        sc = json.load(open(f"{WEB_OUT}/library/rh_scale.json"))["scale"]
        dj = json.load(open(f"{WEB_OUT}/library/diffusion.json"))["types"]
        R = np.full(len(self.names), np.nan); D = np.full(len(self.names), np.nan); unit = [""] * len(self.names)
        for k, n in enumerate(self.names):
            ts = P.PT[P.pidx[n]].indices
            if len(ts) == 0:
                continue
            t = int(ts[np.argmax(P.n_inst[ts])])
            R[k] = lib.rh(t) * sc.get(str(t), 1.0); unit[k] = lib.types[t]["id"]
            if str(t) in dj:
                D[k] = dj[str(t)][0]
        self.R, self.D, self.unit = R, D, unit
        self.M10 = P.matrix(10, "M"); self.C10 = P.matrix(10, "C"); self.C2 = P.matrix(2, "C")
        from scipy.sparse import diags
        self.same = np.asarray((P.PT @ diags(P.n_inst) @ P.PT.T).todense()) / np.maximum(P.copies[:, None], 1)

    def table(self, A, B=None):
        A = [a for a in A if a in self.loc]
        B = self.names if B is None else [b for b in B if b in self.loc]
        ib = np.array([self.loc[b] for b in B]); pb = self.idx[ib]
        rows = []
        for a in A:
            ia = self.loc[a]; pa = self.idx[ia]
            c_comp = (self.W[ia][None, :] * self.conc[ib]).sum(1)
            Rs = self.R[ia] + self.R[ib]
            shell = 4 / 3 * np.pi * ((Rs + 10) ** 3 - Rs ** 3)
            c_loc = np.nan_to_num(self.M10[pa, pb]) / shell
            k = 4 * np.pi * (self.D[ia] + self.D[ib]) * 1e6 * Rs * c_comp
            rows.append(pd.DataFrame(dict(
                A=a, B=B, copies_A=self.copies[ia], copies_B=self.copies[ib], particles_B=self.inst[ib], unit_A=self.unit[ia],
                unit_B=[self.unit[i] for i in ib], p_near_2=np.nan_to_num(self.C2[pa, pb]),
                p_near_10=np.nan_to_num(self.C10[pa, pb]), p_same_complex=np.clip(self.same[pa, pb], 0, 1),
                c_comp_uM=c_comp / UM, c_local_uM=c_loc / UM,
                enrichment=np.where(c_comp > 0, c_loc / np.maximum(c_comp, 1e-30), np.nan),
                D_A=self.D[ia], D_B=self.D[ib], R_A=self.R[ia], R_B=self.R[ib], k_enc_per_s=k)))
        out = pd.concat(rows, ignore_index=True)
        return out[out.A != out.B]


def kinases():
    s = pd.read_csv(f"{WORK}/validation/proximity_sets.tsv", sep="\t", dtype=str)
    return sorted(set(s.loc[s.set.str.startswith("kinase"), "A"]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--a", nargs="*"); ap.add_argument("--kinases", action="store_true")
    a = ap.parse_args()
    E = Encounters()
    if a.kinases:
        t = E.table(kinases())
        t.to_parquet(f"{OUT}/encounters_kinases.parquet", index=False)
        print(f"{len(t)} rows, {t.A.nunique()} kinases -> {OUT}/encounters_kinases.parquet")
    if a.a:
        t = E.table(a.a)
        p = f"{OUT}/encounters_{'_'.join(a.a[:3])}.tsv"
        t.sort_values("p_near_10", ascending=False).to_csv(p, sep="\t", index=False, float_format="%.4g")
        print(t.sort_values("k_enc_per_s", ascending=False).head(15).to_string())
