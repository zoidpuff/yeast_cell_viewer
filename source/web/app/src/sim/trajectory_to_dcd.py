#!/usr/bin/env python
"""Expand the rigid-body trajectory of a browser simulation into a bead-level DCD (VMD, ChimeraX, MDAnalysis, mdtraj).

Files (from the viewer's simulation mode, "Export MD files"): yeast_box.cif, yeast_box.json, yeast_box_trajectory.json.
Each molecule's beads are moved rigidly with its recorded pose (position, quaternion x,y,z,w); membrane beads stay put.

  python trajectory_to_dcd.py [--out yeast_box_browser.dcd]
  then e.g.  vmd yeast_box.cif yeast_box_browser.dcd   or   ChimeraX: open yeast_box.cif; open yeast_box_browser.dcd
Requires: openmm (for the mmCIF reader and DCD writer), numpy.
"""
import argparse, json
import numpy as np
import openmm.app as app
import openmm.unit as u


def rotmat(q):
    """rotation matrices from quaternions (..., 4) as x, y, z, w"""
    x, y, z, w = q[..., 0], q[..., 1], q[..., 2], q[..., 3]
    return np.stack([np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)], -1),
                     np.stack([2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)], -1),
                     np.stack([2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)], -1)], -2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", default="yeast_box.cif"); ap.add_argument("--json", default="yeast_box.json")
    ap.add_argument("--traj", default="yeast_box_trajectory.json"); ap.add_argument("--out", default="yeast_box_browser.dcd")
    a = ap.parse_args()
    cif = app.PDBxFile(a.cif)
    x = np.array(cif.getPositions(asNumpy=True).value_in_unit(u.nanometer))
    meta = json.load(open(a.json)); tr = json.load(open(a.traj))["frames"]
    bodies = meta["bodies"]
    chain_of = np.array([at.residue.chain.id for at in cif.topology.atoms()])
    # bead -> body (manifest order = chain order); local coordinates from the exported pose
    idx = {b["chain"]: i for i, b in enumerate(bodies)}
    body = np.array([idx.get(c, -1) for c in chain_of])
    P = np.array([b["pos_nm"] for b in bodies]); Q = np.array([b["quat_xyzw"] for b in bodies])
    R0 = rotmat(Q)
    m = body >= 0
    local = np.zeros_like(x)
    local[m] = np.einsum("nji,nj->ni", R0[body[m]], x[m] - P[body[m]])      # R^T (x - p)
    lo, hi = np.array(meta["box"]["lo"]), np.array(meta["box"]["hi"])
    dt_ps = (tr[1]["t_ns"] - tr[0]["t_ns"]) * 1e3 if len(tr) > 1 else 1.0
    with open(a.out, "wb") as fh:
        dcd = app.DCDFile(fh, cif.topology, dt_ps * u.picoseconds)
        for f in tr:
            pf = np.array(f["pos"]).reshape(-1, 3); qf = np.array(f["quat"]).reshape(-1, 4)
            Rf = rotmat(qf)
            y = x.copy()
            y[m] = pf[body[m]] + np.einsum("nij,nj->ni", Rf[body[m]], local[m])
            dcd.writeModel(y * u.nanometer, periodicBoxVectors=None)
    print(f"{len(tr)} frames, {len(x)} beads ({m.sum()} in {len(bodies)} molecules), "
          f"{tr[0]['t_ns']:.1f}-{tr[-1]['t_ns']:.1f} ns -> {a.out}; box {lo} to {hi} nm")


if __name__ == "__main__":
    main()
