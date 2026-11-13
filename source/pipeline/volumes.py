"""Compartment volumes and membrane areas of the built geometry (stage 2: tile_labvol.parquet, tile_measure.parquet),
so that no script carries a copy of them (until 2026-11-12 several validation scripts had 26.62 um3 of cytosol
written in, which went stale when the mitochondria grew, and would be wrong for any other cell size)."""
import pandas as pd
from config import W3D

_V = None


def measured():
    """{name: um^3} for cytosol, nucleoplasm, nucleolus, ER (cortical + perinuclear lumen), mitochondrion (matrix +
    intermembrane space), mito_matrix, mito_IMS, vacuole, golgi_endosome, peroxisome, lipid_droplet, cell_wall, cell;
    and {name: um^2} under 'area' for the membranes (PM, ER, NE, VAC, MITO_OM, MITO_IM, GOLGI, PEROX)"""
    global _V
    if _V is None:
        lv = pd.read_parquet(f"{W3D}/tile_labvol.parquet").sum() / 1e9
        me = pd.read_parquet(f"{W3D}/tile_measure.parquet").sum()
        v = dict(cytosol=lv.cytosol, nucleoplasm=lv.nucleoplasm, nucleolus=lv.nucleolus, ER=lv.cER_lumen + lv.perinuclear,
                 mito_matrix=lv.mito_matrix, mito_IMS=lv.mito_IMS, mitochondrion=lv.mito_matrix + lv.mito_IMS,
                 vacuole=lv.vacuole, golgi_endosome=lv.golgi_endosome, peroxisome=lv.peroxisome,
                 lipid_droplet=lv.lipid_droplet, cell_wall=lv.wall_mannan + lv.wall_glucan)
        v["cell"] = float(lv.drop("exterior").sum())
        v = {k: float(x) for k, x in v.items()}
        v["area"] = {k[4:]: float(me[k]) / 1e6 for k in me.index if k.startswith("mem_")}
        _V = v
    return _V
