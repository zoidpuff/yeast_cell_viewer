"""Spindle pole body plaques (IMPROVEMENTS N4; the Jurkat model's pericentriolar region, made specific).

The SPB of a G1 cell is a layered cylinder in the nuclear envelope (cryo-EM of isolated SPBs: Bullitt et al. 1997
Cell 89:1077; integrative model of the core: Viswanath et al. 2017 Mol Biol Cell 28:3298). Before this module its
proteins were spread through the nucleoplasm or along the envelope, so a simulated Spc42-GFP image showed no spot.
Here the copies of each SPB protein are placed in its layer, a disc of radius PLAQUE_R (the microtubule plaques,
microtubules.PLAQUE_R) centred on the SPB axis, at a height range above the outer nuclear membrane (outward normal;
negative: below the inner nuclear membrane):
  outer plaque         Spc72, Nud1            cytoplasm, +45..+75 nm
  IL2 / spacer          Cnm67                  cytoplasm, +20..+50 nm
  central plaque        Spc42, Spc29, Cmd1     cytoplasm, +3..+25 nm   (in cells: inside the envelope fenestra)
  inner plaque          Spc110                 nucleoplasm, 3..60 nm below the inner membrane
  half-bridge           Sfi1, Kar1             cytoplasm, +3..+25 nm, 90..200 nm from the axis on one side
Simplifications: the envelope has no fenestra, so the central plaque sits on its cytoplasmic face; the layer heights
are approximate (unchecked). The rest of each pool stays where it was.
Copies (2026-11-12, IMPROVEMENTS P3): the core proteins come in the stoichiometry of the core's unit cell,
Spc42 : Spc29 : Spc110 : Cmd1 : Cnm67 = 2 : 1 : 1 : 1 : 1 (Muller et al. 2005 Mol Biol Cell 16:3341, as used for the integrative model of
Viswanath et al. 2017), scaled to the Spc110 at the SPB (CORE_RATIO x n: the 14 per nuclear gamma-TuRC, ring_spc110,
plus those placed in the inner-plaque layer, which is filled first); calmodulin binds the Spc110 C-terminus in the
central plaque. The absolute
number per SPB is not measured here. Spc42, Spc29 and Cnm67 are added together in units of 2 : 1 : 1 until units stop
fitting (3 in a row), so the plaque disc caps them in their ratio; calmodulin takes one copy per Spc110. In the
90 nm disc the units run out before Spc42 reaches twice the Spc110 (build P: about 0.8 Spc42 per Spc110), so the
core is short of Spc42, Spc29 and Cnm67 relative to Spc110. Before, every protein was filled from its whole pool up
to what fit, which gave Spc42 : Spc29 : Cnm67 about 1 : 1 : 1 and one calmodulin per inner-plaque Spc110 only. Nud1, Spc72, Sfi1 and Kar1 keep their
pools (no stoichiometry in hand).
"""
import numpy as np
from config import seed as _seed
import pack3d as p3
import cell3d as c3

PLAQUE_R = 90.0
LAYERS = {"SPC72": ("cyto", 45.0, 75.0), "NUD1": ("cyto", 45.0, 75.0), "CNM67": ("cyto", 20.0, 50.0),
          "SPC42": ("cyto", 3.0, 25.0), "SPC29": ("cyto", 3.0, 25.0), "SPC110": ("nuc", 3.0, 60.0),
          "SFI1": ("bridge", 3.0, 25.0), "KAR1": ("bridge", 3.0, 25.0), "CMD1": ("cyto", 3.0, 25.0)}
ORDER = ["SPC110", "SPC42", "SPC29", "CNM67", "NUD1", "SPC72", "SFI1", "KAR1", "CMD1"]
CORE_RATIO = {"SPC42": 2.0, "SPC29": 1.0, "CNM67": 1.0, "CMD1": 1.0}     # per Spc110 at the SPB
CORE_UNIT = ["SPC42", "SPC42", "SPC29", "CNM67"]                            # added together, in this ratio
TRIES = 200


def place(cell, lp, centre_dir, pools, obstacles, seed=71, log=print, ring_spc110=None):
    """centre_dir: direction from the nucleus centre to the SPB; pools: {gene: copies}; ring_spc110: Spc110 copies
    in the nuclear gamma-TuRCs (None: every protein from its whole pool, as before 2026-11-12). Returns pre-pass
    group 'spb' (spb_type, spb_pos, spb_q), per-gene counts and the placed spheres."""
    from scipy.spatial.transform import Rotation
    ids = list(lp["ids"])
    seed = _seed(seed); rng = np.random.default_rng(seed)
    G = cell.G
    onm, n_out = c3.nuc_surface(G, centre_dir, 0.5 * c3.PERI_W)          # outer membrane surface (approx.)
    inm, _ = c3.nuc_surface(G, centre_dir, -1.5 * c3.PERI_W)             # below the inner membrane
    e1 = np.cross(n_out, [0.3, 0.5, 0.8]); e1 /= np.linalg.norm(e1); e2 = np.cross(n_out, e1)
    OC, OR = obstacles
    k = np.linalg.norm(OC - onm, axis=1) < 600
    OC, OR = OC[k].astype(np.float64), OR[k].astype(np.float64)
    glo, ghi = onm - 700, onm + 700
    a = cell.args(glo, ghi, margin=200.0)
    out = np.zeros(c3.NMEM)
    placed_c, placed_r, types, pos, quat, counts = [], [], [], [], [], {}

    grid = p3.LocalGrid(glo, ghi, OC, OR)

    # fill all layers at once, interleaving the proteins (a copy of each in turn, in proportion to its pool), so that
    # no protein crowds the others out of a shared layer; a protein stops after 3 failed copies in a row
    genes = [g for g in ORDER if g in ids and g in LAYERS and g != "CMD1"]
    want = {g: int(pools.get(g, 0)) for g in genes}
    fails = {g: 0 for g in genes}
    core_spc110 = None

    def put(g):
        side, h0, h1 = LAYERS[g]
        li = ids.index(g)
        C = lp["STC"][lp["st_off"][li]:lp["st_off"][li + 1]].astype(np.float64)
        R = lp["STR"][lp["st_off"][li]:lp["st_off"][li + 1]].astype(np.float64)
        mask = 1 << (c3.NUCP if side == "nuc" else c3.CYTO)
        for _t in range(TRIES):
            ph = rng.uniform(0, 2 * np.pi)
            if side == "bridge":
                r = rng.uniform(PLAQUE_R, 200.0); ph = rng.uniform(-0.4, 0.4)
            else:
                r = PLAQUE_R * np.sqrt(rng.random())
            h = rng.uniform(h0, h1)
            base = (inm - h * n_out) if side == "nuc" else (onm + h * n_out)
            ctr = base + r * (np.cos(ph) * e1 + np.sin(ph) * e2)
            Rm = Rotation.random(random_state=rng).as_matrix()
            t = ctr - Rm @ C.mean(0); S = C @ Rm.T + t
            if grid.free(S, R) and p3.check_volume(S, R, t[0], t[1], t[2], float(lp["BR"][li]), mask, *a, out):
                grid.add(S, R); placed_c.append(S); placed_r.append(R)
                types.append(li); pos.append(t); quat.append(np.roll(Rotation.from_matrix(Rm).as_quat(), 1))
                counts[g] = counts.get(g, 0) + 1
                return True
        return False

    if ring_spc110 is not None and "SPC110" in genes:       # the inner plaque first (nuclear side, no competition)
        for _ in range(want["SPC110"]):
            fails["SPC110"] = 0 if put("SPC110") else fails["SPC110"] + 1
            if fails["SPC110"] >= 3:
                break
        core_spc110 = ring_spc110 + counts.get("SPC110", 0)
        units = min(int(pools.get(g, 0)) // CORE_UNIT.count(g) for g in set(CORE_UNIT))
        miss = 0
        for _ in range(min(units, int(round(core_spc110)))):          # at most one unit (2 Spc42) per Spc110
            ok = all([put(g) for g in CORE_UNIT])
            miss = 0 if ok else miss + 1
            if miss >= 3:
                break
        genes = [g for g in genes if g != "SPC110" and g not in CORE_UNIT]
    queue = np.concatenate([np.full(want[g], i) for i, g in enumerate(genes)]) if genes else np.zeros(0, int)
    queue = queue[np.argsort(rng.random(len(queue)))]
    for i in queue:
        g = genes[i]
        if fails[g] >= 3:
            continue
        fails[g] = 0 if put(g) else fails[g] + 1
    if "CMD1" in ids:                               # calmodulin on the Spc110 C-termini (one each)
        f = 0
        n_cmd = counts.get("SPC110", 0) if core_spc110 is None else int(round(CORE_RATIO["CMD1"] * core_spc110))
        for _ in range(min(int(pools.get("CMD1", 0)), n_cmd)):
            f = 0 if put("CMD1") else f + 1
            if f >= 3:
                break
    log("spindle pole body: " + ", ".join(f"{g} {counts.get(g, 0)}/{int(pools.get(g, 0))}" for g in ORDER if g in ids)
        + ("" if core_spc110 is None else f"; core scaled to {core_spc110:.0f} Spc110 at the SPB (2:1:1:1:1)"))
    if not types:
        return {}, counts, (np.zeros((0, 3)), np.zeros(0))
    return (dict(spb_type=np.array(types, np.int32), spb_pos=np.array(pos, np.float32), spb_q=np.array(quat, np.float32)),
            counts, (np.concatenate(placed_c), np.concatenate(placed_r)))
