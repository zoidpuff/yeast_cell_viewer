#!/usr/bin/env python
"""Run a box exported from the yeast whole-cell viewer in OpenMM, as a coarse-grained model.

Files (from the viewer's simulation mode, "Export MD files"):
  yeast_box.cif    one pseudo-atom per residue (CA) / nucleotide (P) bead, one chain per molecule; chain ZZZZ: membrane
                   midplane beads; occupancy 1 = mobile, 0 = fixed (membrane proteins, filaments, chromatin)
  yeast_box.json   manifest: molecule sources (PDB / AlphaFold DB), poses, hydrodynamic radii, settings of the browser run

Model (deliberately simple; edit to taste):
  - each molecule: elastic network over its beads (pairs closer than ENM_CUT, force constant ENM_K), so it stays close
    to its structure; inter-molecular: Weeks-Chandler-Andersen excluded volume (EPS_REP, sigma = bead diameters) plus
    the optional short-range attraction of the browser run (depth "eps" kT over 0.6 nm, cos^2 well)
  - fixed beads (occupancy 0, membrane beads) have zero mass, which OpenMM keeps in place
  - box faces as in the browser run: periodic boundary conditions, or soft walls keeping mobile beads inside the box;
    "embedded" boxes also tether the shell molecules (harmonic restraint of each bead to its start, the body's tether
    shared over its beads)
  - interactions of the browser run, when exported: screened electrostatics between charged beads (Debye-Hueckel,
    manifest "interactions.charges"), polysome mRNA links (worm-like chain between the beads at the mRNA exit of one
    ribosome and the entry of the next), native contacts of predicted complexes (cos^2 wells between residue pairs)
  - Langevin dynamics at 303 K after an energy minimisation

  python run_openmm.py [--steps 50000] [--dt 0.02] [--out traj.dcd] [--platform OpenCL]
Requires: openmm (pip install openmm), numpy.
"""
import argparse, json, sys
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as u

ENM_CUT, ENM_K = 1.2, 500.0          # nm, kJ/mol/nm^2
ENM_MIN = 0.1                         # nm; shorter pairs get no bond (bonds of near-zero length blew up on CUDA)
WCA_CORE = 0.8                        # below this fraction of sigma the excluded volume continues linearly (capped force)
EPS_REP = 2.5                         # kJ/mol (~1 kT)
SIGMA = {"CA": 0.50, "P": 0.80, "M": 3.0}   # nm, bead diameters (membrane midplane bead: most of the bilayer)
MASS = {"CA": 110.0, "P": 330.0, "M": 0.0}  # amu


def pairs_within(x, cut):
    """index pairs (i < j) closer than cut, with their distances (scipy k-d tree if available, else a cell list)"""
    try:
        from scipy.spatial import cKDTree
        p = cKDTree(x).query_pairs(cut, output_type="ndarray")
        return np.c_[p, np.linalg.norm(x[p[:, 0]] - x[p[:, 1]], axis=1)] if len(p) else np.zeros((0, 3))
    except ImportError:
        pass
    keys = np.floor(x / cut).astype(np.int64)
    cells = {}
    for i, k in enumerate(map(tuple, keys)):
        cells.setdefault(k, []).append(i)
    out = []
    for k, idx in cells.items():
        idx = np.array(idx)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    nb = cells.get((k[0] + dx, k[1] + dy, k[2] + dz))
                    if nb is None:
                        continue
                    nb = np.array(nb)
                    d = np.linalg.norm(x[idx][:, None] - x[nb][None], axis=2)
                    a, b = np.nonzero((d < cut) & (idx[:, None] < nb[None]))
                    out.append(np.c_[idx[a], nb[b], d[a, b]])
    return np.vstack(out) if out else np.zeros((0, 3))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", default="yeast_box.cif"); ap.add_argument("--json", default="yeast_box.json")
    ap.add_argument("--steps", type=int, default=50000); ap.add_argument("--dt", type=float, default=0.02, help="ps")
    ap.add_argument("--out", default="yeast_box.dcd"); ap.add_argument("--every", type=int, default=500)
    ap.add_argument("--platform", default=None)
    a = ap.parse_args()
    meta = json.load(open(a.json))
    eps_att = float(meta["settings"].get("eps", 0.0)) * 2.519            # kT -> kJ/mol at 303 K
    cif = app.PDBxFile(a.cif)
    top = cif.topology
    x = np.array(cif.getPositions(asNumpy=True).value_in_unit(u.nanometer))
    atoms = list(top.atoms())
    kind = np.array(["M" if at.residue.chain.id == "ZZZZ" else at.name for at in atoms])
    chain = np.array([at.residue.chain.index for at in atoms])
    mobile_chain = {b["chain"]: b["mobile"] for b in meta["bodies"]}
    mobile = np.array([kind[i] != "M" and mobile_chain.get(at.residue.chain.id, False) for i, at in enumerate(atoms)])
    print(f"{len(atoms)} beads in {top.getNumChains()} chains; {mobile.sum()} mobile", flush=True)

    system = mm.System()
    for i in range(len(atoms)):
        system.addParticle(MASS[kind[i]] if mobile[i] else 0.0)
    # elastic network within each mobile molecule
    bonds = mm.HarmonicBondForce(); nb_ = 0
    for c in np.unique(chain[mobile]):
        idx = np.where((chain == c) & mobile)[0]
        if len(idx) < 2:
            continue
        p = pairs_within(x[idx], ENM_CUT)
        p = p[p[:, 2] >= ENM_MIN] if len(p) else p      # near-duplicate beads (overlapping RNA/protein beads of
        for i, j, d in p:                               # ribosome models) would get zero-length bonds
            bonds.addBond(int(idx[int(i)]), int(idx[int(j)]), float(d), ENM_K); nb_ += 1
    system.addForce(bonds)
    # excluded volume (+ attraction) between molecules
    # pairs inside one molecule are left to the elastic network; pairs of two fixed beads are switched off (m = mobile)
    inter = meta.get("interactions") or {}
    el = inter.get("electrostatics")
    q = np.zeros(len(atoms))
    for i, z in (inter.get("charges") or []) if el else []:
        q[int(i)] = z
    kappa = float(el["kappa_per_nm"]) if el else 1.0
    kt = 2.519                                                           # kJ/mol at 303 K
    expr = ("(1 - delta(c1 - c2)) * step(m1 + m2 - 0.5) * (wca + att + elec);"
            "wca = step(rc - r) * (4 * er * ((s / rw)^12 - (s / rw)^6) + er + fw * (rw - r));"
            "fw = 4 * er * (12 * (s / rw)^12 - 6 * (s / rw)^6) / rw; rw = max(r, core * s);"
            "att = -ea * step(r - s) * step(s + 0.6 - r) * cos(1.5707963 * (r - s) / 0.6)^2;"
            "elec = lb * q1 * q2 * (exp(-kappa * rr) / rr - exp(-kappa * rcut) / rcut) * step(rcut - r); rr = max(r, s);"
            "rc = s * 1.122462; s = 0.5 * (s1 + s2)")
    nbf = mm.CustomNonbondedForce(expr)
    nbf.addGlobalParameter("er", EPS_REP); nbf.addGlobalParameter("ea", eps_att); nbf.addGlobalParameter("core", WCA_CORE)
    nbf.addGlobalParameter("lb", kt * float(el["bjerrum_nm"]) if el else 0.0); nbf.addGlobalParameter("kappa", kappa)
    nbf.addPerParticleParameter("c"); nbf.addPerParticleParameter("s"); nbf.addPerParticleParameter("m"); nbf.addPerParticleParameter("q")
    for i in range(len(atoms)):
        nbf.addParticle([float(chain[i]), SIGMA[kind[i]], 1.0 if mobile[i] else 0.0, float(q[i])])
    lo, hi = np.array(meta["box"]["lo"]), np.array(meta["box"]["hi"])
    boundary = meta["box"].get("boundary", "walls")
    if boundary == "periodic":
        L = hi - lo
        system.setDefaultPeriodicBoxVectors(mm.Vec3(L[0], 0, 0), mm.Vec3(0, L[1], 0), mm.Vec3(0, 0, L[2]))
        nbf.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    else:
        nbf.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
    s_max = 0.5 * (max(SIGMA["CA"], SIGMA["P"]) + max(SIGMA.values()))    # mobile bead against anything
    rcut = max(s_max * 1.123, s_max + 0.6) + 0.05
    if el:
        rcut = max(rcut, 2 * 0.4 + 3.0 / kappa)                              # three Debye lengths beyond contact
    nbf.addGlobalParameter("rcut", rcut)
    nbf.setCutoffDistance(rcut)
    nbf.setUseLongRangeCorrection(False)
    mob_ix = [int(i) for i in np.where(mobile)[0]]
    system.addForce(nbf)
    # box walls (closed and embedded boxes)
    if boundary != "periodic":
        wall = mm.CustomExternalForce("kw * (max(0, xl - x)^2 + max(0, x - xh)^2 + max(0, yl - y)^2 + max(0, y - yh)^2 + "
                                      "max(0, zl - z)^2 + max(0, z - zh)^2)")
        wall.addGlobalParameter("kw", 500.0)
        for nm_, v in zip(("xl", "yl", "zl", "xh", "yh", "zh"), list(lo) + list(hi)):
            wall.addGlobalParameter(nm_, float(v))
        for i in mob_ix:
            wall.addParticle(i, [])
        system.addForce(wall)
    # embedded boundary: shell molecules tethered to where they started
    tk = float(meta["box"].get("tether_kT_per_nm2", 0.0)) * 2.519        # kT/nm^2 -> kJ/mol/nm^2
    teth = {b["chain"] for b in meta["bodies"] if b.get("tethered")}
    if tk > 0 and teth:
        tf = mm.CustomExternalForce("0.5 * kt * ((x - x0)^2 + (y - y0)^2 + (z - z0)^2)")
        tf.addPerParticleParameter("kt"); tf.addPerParticleParameter("x0"); tf.addPerParticleParameter("y0"); tf.addPerParticleParameter("z0")
        nb_body = {}
        for at in atoms:
            nb_body[at.residue.chain.id] = nb_body.get(at.residue.chain.id, 0) + 1
        for i, at in enumerate(atoms):
            c = at.residue.chain.id
            if c in teth and mobile[i]:
                tf.addParticle(i, [tk / nb_body[c], *map(float, x[i])])
        system.addForce(tf)
    # polysome mRNA links: worm-like chain (Marko-Siggia), extension capped at 0.95 of the contour
    links = inter.get("mrna_links") or []
    if links:
        wlc = mm.CustomBondForce("kt * lc / lp * (0.5 * u^2 - 0.25 * u + 0.25 / (1 - u) - 0.25); u = min(r / lc, 0.95)")
        wlc.addGlobalParameter("kt", kt); wlc.addPerBondParameter("lc"); wlc.addPerBondParameter("lp")
        wlc.setUsesPeriodicBoundaryConditions(boundary == "periodic")
        for l in links:
            wlc.addBond(int(l["atoms"][0]), int(l["atoms"][1]), [float(l["contour_nm"]), float(l["persistence_nm"])])
        system.addForce(wlc)
    # native contacts of predicted complexes: each residue pair a cos^2 well of its share of the interface energy
    nat = inter.get("native_contacts") or []
    if nat:
        e_if = float(inter["native_contact_model"]["interface_energy_kT"]) * kt
        go = mm.CustomBondForce("-e * select(step(0.5 - r), 1, step(1.1 - r) * cos(1.5707963 * (r - 0.5) / 0.6)^2)")
        go.addPerBondParameter("e"); go.setUsesPeriodicBoundaryConditions(boundary == "periodic")
        for i, j, share in nat:
            if mobile[int(i)] or mobile[int(j)]:
                go.addBond(int(i), int(j), [e_if * float(share)])
        system.addForce(go)
    print(f"elastic-network bonds {nb_}; attraction {eps_att:.2f} kJ/mol; boundary {boundary}; charged beads {int((q != 0).sum())}; "
          f"mRNA links {len(links)}; native contact pairs {len(nat)}", flush=True)

    integ = mm.LangevinMiddleIntegrator(meta.get("temperature_K", 303.15) * u.kelvin, 1.0 / u.picosecond, a.dt * u.picoseconds)
    plat = mm.Platform.getPlatformByName(a.platform) if a.platform else None
    sim = app.Simulation(top, system, integ, plat) if plat else app.Simulation(top, system, integ)
    sim.context.setPositions(x * u.nanometer)
    print("platform", sim.context.getPlatform().getName(), flush=True)
    sim.minimizeEnergy(maxIterations=500)
    sim.reporters.append(app.DCDReporter(a.out, a.every))
    sim.reporters.append(app.StateDataReporter(sys.stdout, a.every * 10, step=True, time=True, potentialEnergy=True,
                                               temperature=True, speed=True))
    sim.step(a.steps)
    print("done:", a.out)


if __name__ == "__main__":
    main()
