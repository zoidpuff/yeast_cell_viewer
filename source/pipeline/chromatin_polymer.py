"""Chromosome arms as a coarse polymer relaxed in the nucleus (IMPROVEMENTS N2, after the Jurkat model v7).

The fibre walk (pack3d.fiber_walks) places each arm in one pass, steered from its kinetochore to a telomere target;
chromosomes never relax against each other, and built cells came out too territorial with too flat a contact decay
(Hi-C test). Here every arm is first a chain of BEAD_BP beads, simulated with Langevin dynamics in OpenMM:
  chain        worm-like: bond BEAD_BP / BP_PER_NM nm, bending stiffness for persistence LP_NM (Arbona et al. 2017
               Genome Biol 18:81: 53-65 bp/nm, 52-85 nm; middle values)
  excluded     soft repulsion between beads closer than SIGMA_NM (E_REP kT at full overlap; strand passage is rare)
  confinement  the nucleoplasm: inside the inner nuclear membrane, outside the nucleolus sphere, clear of the pre-placed
               obstacles (nuclear pore baskets, microtubules, kinetochores); a tabulated violation depth
  tethers      both arms of a chromosome start at its kinetochore (flat-bottom, CEN_R nm); telomere beads stay within
               TEL_R nm of the inner nuclear membrane (no foci are imposed); chromosome XII-R ends at the nucleolus and
               restarts from it (the rDNA is built in the nucleolus, wc3d.rdna_chromatin)
The relaxed chain is the guide for the nucleosomes (pack3d.guided_walks): each measured nucleosome is placed near
the guide point at its bp. Nothing here is fitted to Hi-C; the Hi-C test compares contacts (Tjong et al. 2012 Genome
Res 22:1295: tethering and volume exclusion alone predict much of the yeast inter-chromosome contact pattern).
Several replicas with different seeds give an ensemble for population statistics; the built cell uses replica 0.
"""
import os, time
import numpy as np

BEAD_BP = 2000
BP_PER_NM = 60.0
LP_NM = 70.0
SIGMA_NM, E_REP = 30.0, 5.0
CEN_R, TEL_R = 20.0, 60.0
STEPS, DT_PS, GAMMA = 400_000, 5.0, 0.01
GRID_NM, CLEAR_NM = 15.0, 15.0
KT = 2.494                     # kJ/mol at 300 K


def kappa_for(lp, b):
    """bending constant k of E = k (1 + cos theta) for persistence lp with bond b: <cos> = coth k - 1/k = exp(-b/lp)"""
    target = np.exp(-b / lp)
    lo, hi = 1e-3, 100.0
    for _ in range(80):
        k = 0.5 * (lo + hi)
        if 1 / np.tanh(k) - 1 / k < target:
            lo = k
        else:
            hi = k
    return 0.5 * (lo + hi)


def chain_layout(segs):
    """beads per arm segment (genome_map.segments order): bp of each bead in walk order, segment index"""
    bp, seg = [], []
    for s, (c, kind, s0, s1) in enumerate(segs):
        L = abs(s1 - s0); n = max(2, int(np.ceil(L / BEAD_BP)))
        x = s0 + np.sign(s1 - s0) * (np.arange(n) + 0.5) * L / n
        bp.append(x); seg.append(np.full(n, s))
    return np.concatenate(bp), np.concatenate(seg)


def violation_grid(cell, obstacles, clear=CLEAR_NM, h=GRID_NM):
    """tabulated depth (nm, >= 0) by which a bead centre violates the nucleoplasm: outside the inner nuclear membrane
    less clear, or within clear of an obstacle sphere (the nucleolus is analytic, in the force). Returns origin,
    spacing, shape, values (x fastest) and the signed distance below the inner nuclear membrane (for telomeres)."""
    import cell3d as c3
    from scipy.spatial import cKDTree
    G = cell.G; cn = G[3:6]; R = float(max(G[11:14].max(), G[6])) + 150.0
    o = cn - R; n = int(np.ceil(2 * R / h)) + 1
    ax = [o[i] + h * np.arange(n) for i in range(3)]
    Z, Y, X = np.meshgrid(ax[2], ax[1], ax[0], indexing="ij")          # x fastest when flattened
    P = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)
    sd = _sdf_many(P, G)                                             # < 0 inside the ONM midsurface
    below_inm = -(sd + c3.PERI_W)                                    # depth below the inner nuclear membrane (nm)
    v = np.maximum(0.0, clear - below_inm)
    OC, OR = obstacles
    if len(OC):
        k = np.all(np.abs(OC - cn) < R, axis=1); OC, OR = OC[k], OR[k]
        inside = below_inm > -h
        t = cKDTree(OC); rmax = float(OR.max())
        d, i = t.query(P[inside], distance_upper_bound=rmax + clear + h)
        ok = np.isfinite(d)
        gap = np.full(inside.sum(), np.inf); gap[ok] = d[ok] - OR[i[ok]]
        v_ob = np.maximum(0.0, clear - gap)
        v[inside] = np.maximum(v[inside], v_ob)
    return o, h, (n, n, n), v, below_inm


def _sdf_many(P, G):
    import numba as nb
    import cell3d as c3

    @nb.njit(parallel=True, cache=False)
    def f(P, G, out):
        for i in nb.prange(P.shape[0]):
            out[i] = c3.nuc_sdf(P[i, 0], P[i, 1], P[i, 2], G)
    out = np.empty(len(P)); f(P, G, out)
    return out


def start_config(cell, bp, seg, segs, starts, t_in, t_out, rng):
    """a random worm-like start inside the nucleus: each arm from its anchor, bond b, steps kept inside 85% of the
    nuclear radius and out of the nucleolus (overlaps are resolved by the dynamics)"""
    G = cell.G; cn = G[3:6]; rn = float(G[6]); ncl, rncl = G[7:10], G[10]
    b = BEAD_BP / BP_PER_NM
    X = np.zeros((len(bp), 3))
    for s, (c, kind, s0, s1) in enumerate(segs):
        I = np.where(seg == s)[0]
        p = (t_out if kind == "XIIR_out" else starts[c]).copy()
        u = rng.normal(size=3); u /= np.linalg.norm(u)
        for i in I:
            for _ in range(50):
                v = u + 0.6 * rng.normal(size=3); v /= np.linalg.norm(v)
                q = p + b * v
                if np.linalg.norm(q - cn) < 0.85 * rn and np.linalg.norm(q - ncl) > rncl + 20:
                    break
                v = (cn - p) / np.linalg.norm(cn - p) + 0.5 * rng.normal(size=3); v /= np.linalg.norm(v); q = p + b * v
            X[i] = q; p = q; u = v
    return X


def simulate(cell, segs, starts, t_in, t_out, obstacles, seed=0, steps=STEPS, platform="CUDA", log=print, grid=None):
    """relaxed bead positions (nm), bead bp, bead segment; one replica"""
    import openmm as mm
    import openmm.unit as u
    rng = np.random.default_rng(seed)
    bp, seg = chain_layout(segs)
    N = len(bp); b = BEAD_BP / BP_PER_NM
    o, h, shape, v, below = grid if grid is not None else violation_grid(cell, obstacles)
    X0 = start_config(cell, bp, seg, segs, starts, t_in, t_out, rng)
    sysm = mm.System()
    for _ in range(N):
        sysm.addParticle(100.0)
    bonds = mm.HarmonicBondForce(); kb = KT / (0.1 * b) ** 2
    ang = mm.CustomAngleForce("ka*(1+cos(theta))"); ang.addGlobalParameter("ka", kappa_for(LP_NM, b) * KT)
    pairs = []
    first = {}
    for s, (c, kind, s0, s1) in enumerate(segs):
        I = np.where(seg == s)[0]
        for a_, b_ in zip(I[:-1], I[1:]):
            bonds.addBond(int(a_), int(b_), b, kb); pairs.append((a_, b_))
        for a_, b_, c_ in zip(I[:-2], I[1:-1], I[2:]):
            ang.addAngle(int(a_), int(b_), int(c_))
        if kind in ("L", "R", "XIIR_in"):
            first.setdefault(c, []).append(I[0])
    for c, (i, j) in [(c, v_[:2]) for c, v_ in first.items() if len(v_) >= 2]:      # the two arms meet at the centromere
        bonds.addBond(int(i), int(j), 0.5 * b, kb); pairs.append((i, j))
    sysm.addForce(bonds); sysm.addForce(ang)
    rep = mm.CustomNonbondedForce("erep*(1-(r/sig)^2)^2")
    rep.addGlobalParameter("erep", E_REP * KT); rep.addGlobalParameter("sig", SIGMA_NM)
    rep.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic); rep.setCutoffDistance(SIGMA_NM)
    for _ in range(N):
        rep.addParticle([])
    for i, j in pairs:
        rep.addExclusion(int(i), int(j))
    sysm.addForce(rep)
    # confinement: tabulated violation + analytic nucleolus sphere
    ncl, rncl = cell.G[7:10], float(cell.G[10])
    xs = [o[k] + h * np.arange(shape[k]) for k in range(3)]
    # single-particle compound bonds: CustomExternalForce takes no tabulated functions
    conf = mm.CustomCompoundBondForce(1,
        "kc*(vgrid(x1,y1,z1)^2 + max(0, rn + cl - sqrt((x1-nx)^2+(y1-ny)^2+(z1-nz)^2))^2)")
    conf.addGlobalParameter("kc", KT); conf.addGlobalParameter("rn", rncl); conf.addGlobalParameter("cl", CLEAR_NM)
    for k_, val in zip(("nx", "ny", "nz"), ncl):
        conf.addGlobalParameter(k_, float(val))
    conf.addTabulatedFunction("vgrid", mm.Continuous3DFunction(shape[0], shape[1], shape[2], v.astype(float).tolist(),
                                                               xs[0][0], xs[0][-1], xs[1][0], xs[1][-1], xs[2][0], xs[2][-1]))
    for i in range(N):
        conf.addBond([i], [])
    sysm.addForce(conf)
    # tethers: centromere beads, XII-R ends at the nucleolus (flat-bottom wells); telomeres near the envelope
    teth = mm.CustomExternalForce("kt*max(0, sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2) - r0)^2")
    teth.addGlobalParameter("kt", KT)
    for p_ in ("x0", "y0", "z0", "r0"):
        teth.addPerParticleParameter(p_)
    tel = mm.CustomCompoundBondForce(1, "kl*max(0, dgrid(x1,y1,z1) - tr)^2")
    tel.addGlobalParameter("kl", 0.1 * KT); tel.addGlobalParameter("tr", TEL_R)
    tel.addTabulatedFunction("dgrid", mm.Continuous3DFunction(shape[0], shape[1], shape[2], below.astype(float).tolist(),
                                                              xs[0][0], xs[0][-1], xs[1][0], xs[1][-1], xs[2][0], xs[2][-1]))
    n_tel = 0
    for s, (c, kind, s0, s1) in enumerate(segs):
        I = np.where(seg == s)[0]
        if kind in ("L", "R", "XIIR_in"):
            teth.addParticle(int(I[0]), [*map(float, starts[c]), CEN_R])
        if kind == "XIIR_in":
            teth.addParticle(int(I[-1]), [*map(float, t_in), 30.0])
        if kind == "XIIR_out":
            teth.addParticle(int(I[0]), [*map(float, t_out), 30.0])
        if kind in ("L", "R", "XIIR_out"):
            tel.addBond([int(I[-1])], []); n_tel += 1
    sysm.addForce(teth); sysm.addForce(tel)
    integ = mm.LangevinMiddleIntegrator(300 * u.kelvin, GAMMA / u.picosecond, DT_PS * u.picosecond)
    integ.setRandomNumberSeed(int(seed) + 1)
    try:
        plat = mm.Platform.getPlatformByName(platform)
        ctx = mm.Context(sysm, integ, plat, {"Precision": "mixed"} if platform == "CUDA" else {})
    except Exception:
        ctx = mm.Context(sysm, integ, mm.Platform.getPlatformByName("CPU"))
    ctx.setPositions(X0 * u.nanometer)
    mm.LocalEnergyMinimizer.minimize(ctx, 1.0, 2000)
    t0 = time.time()
    chunk = 20_000
    for k in range(0, steps, chunk):
        integ.step(min(chunk, steps - k))
    st = ctx.getState(getPositions=True, getEnergy=True)
    X = np.asarray(st.getPositions(asNumpy=True).value_in_unit(u.nanometer))
    log(f"polymer replica {seed}: {N} beads ({BEAD_BP} bp), {steps} steps of {DT_PS} ps in {time.time() - t0:.0f}s, "
        f"{n_tel} telomeres tethered near the envelope; E = {st.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole) / KT:.0f} kT")
    return X, bp, seg


def guide_points(X, bp, seg, nuc_bp_per_seg):
    """per nucleosome (segments in order, nucleosomes in walk order): the polymer position at its bp (linear
    interpolation between bead centres)"""
    out = []
    for s, p in enumerate(nuc_bp_per_seg):
        I = np.where(seg == s)[0]
        b = bp[I]
        if b[0] > b[-1]:                     # walk runs towards lower bp
            b = -b; q = -np.asarray(p, np.float64)
        else:
            q = np.asarray(p, np.float64)
        out.append(np.stack([np.interp(q, b, X[I, k]) for k in range(3)], 1))
    return np.concatenate(out) if out else np.zeros((0, 3))
