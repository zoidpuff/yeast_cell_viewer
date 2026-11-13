"""Ribosomes on the cortical ER: centres of all 80S ribosome instances (docked on translocons, packed in the ER band and
on ER-tethered polysomes) within NEAR_NM of the ER membrane, per µm² of the ER's cytoplasmic face, against cryo-ET
of the cortical ER (West et al. 2011 J Cell Biol 193:333: 550-900 per µm²). Partly tuned: the ER tether of round N
was chosen with this number in view.

  python validation/er_ribosome_density.py [data_dir]     (default: config.WEB_OUT; any encoded build's data folder)

The area is one face of the cortical ER sheet: half the closed M_ER membrane surface of the tiles' meshes, which
wraps the thin lumen on both sides (the rims are negligible). Round N quoted about 790 per µm² from a one-off area
estimate; this script gives 782 for the same build. Output: work/validation/er_ribosome_density.json (only for the
default data folder).
"""
import os, sys, json, gzip, glob
for _v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):    # one BLAS thread per worker
    os.environ.setdefault(_v, "1")
import multiprocessing as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WORK, WEB_OUT, W3D, NPROC  # noqa: E402

NEAR_NM = 30.0
_G = {}


def _tile(fn):
    import webenc as we, cell3d as c3
    D = we.decode_tile(gzip.decompress(open(fn, "rb").read()))
    types, pos = D["inst"][0], D["inst"][1]
    k = np.isin(types, _G["ribo"])
    out = np.zeros(c3.NMEM); n = 0
    if k.any():
        cell = _G["cell"]; lo = D["lo"]
        a = cell.args(lo - 50, lo + 250, margin=80.0)
        for p in pos[k]:
            c3.mem_sdfs(p[0], p[1], p[2], *a, out)
            n += abs(out[c3.M_ER]) <= NEAR_NM
    # ER membrane area from this tile's mesh triangles of membrane M_ER
    V, F, FM = D["verts"], D["faces"], D["face_mem"]
    area = 0.0
    if len(F):
        m = FM == c3.M_ER
        if m.any():
            T = V[F[m]]
            area = 0.5 * np.linalg.norm(np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]), axis=1).sum()
    return int(n), float(area)


def main(data=WEB_OUT):
    import cell3d as c3
    types = json.loads(gzip.decompress(b) if (b := open(f"{data}/library/types.json", "rb").read())[:2] == b"\x1f\x8b" else b)["types"]
    _G.update(ribo=np.array([i for i, t in enumerate(types) if t["id"].startswith("ribosome_80S")]),
              cell=c3.Cell.load(f"{W3D}/geometry.json"))
    fs = sorted(glob.glob(f"{data}/tiles/*.bin.gz"))
    with mp.get_context("fork").Pool(NPROC) as pool:
        r = pool.map(_tile, fs, chunksize=16)
    n = sum(x[0] for x in r); area_um2 = 0.5 * sum(x[1] for x in r) / 1e6        # one face of the sheet
    res = dict(near_nm=NEAR_NM, ribosomes_near_er=n, er_area_um2=round(area_um2, 2), per_um2=round(n / area_um2, 1),
               ribosome_types=[types[i]["id"] for i in _G["ribo"]], west2011_per_um2=[550, 900], data=data)
    print(json.dumps(res, indent=1))
    if data == WEB_OUT:
        json.dump(res, open(f"{WORK}/validation/er_ribosome_density.json", "w"), indent=1)
    return res


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else WEB_OUT)
