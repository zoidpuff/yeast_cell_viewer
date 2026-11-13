"""Stage 5: encode the packed cell for the WebGPU viewer -> config.WEB_OUT. Notebook cells 38-45.

  python s5_encode.py               # full encode (library, tiles, index, geometry, meshes, bookmarks)
  python s5_encode.py --tiles-only  # tiles only (baseline reproduction check; needs no inventory)
  python s5_encode.py --annot-only  # only library/genome.json and annot/ (polymer identities) for the encoded tiles
"""
import argparse, gzip, shutil, multiprocessing as mp
from _common import *
import cell3d as c3, tiles3d as t3, webenc as we, wc3d, organisation as og, cellshape

MRNA_PIECE = 64          # nt per encoded mRNA piece (pieces are assigned to tiles by their midpoint)
IDX_COLS = ["ti", "n_inst", "n_ghost", "n_mrna", "n_mrna_beads", "n_link", "n_link_beads", "n_vert", "n_tri", "raw_bytes", "gz_bytes"]


def pre_regions(pre):
    """region code per pre-placed instance, in wc3d.pre_instances order"""
    nucl = t3.REGION_NAMES.index("nucleolus")
    fixed = {"npc": we.R_NPC, "nuc": we.R_NUCLEOSOME, "poly": t3.R_CYTO, "act": t3.R_CYTO, "tpm": t3.R_CYTO, "pat": t3.R_CYTO,
             "atp": t3.REGION_NAMES.index("mem_MITO_IM"), "eis": t3.REGION_NAMES.index("mem_PM"), "cvt": t3.R_CYTO, "pol2": t3.R_NUCP, "nvj": t3.R_CYTO, "rvs": t3.R_CYTO, "coat": t3.R_PM_PERI, "mcc": t3.REGION_NAMES.index("mem_PM"), "kt": t3.R_NUCP, "mcm": t3.R_NUCP, "orc": t3.R_NUCP, "tf": t3.R_NUCP, "rnp": t3.R_CYTO, "kto": t3.R_NUCP, "spb": t3.R_CYTO, "lm": t3.R_NUCP, "kap": t3.REGION_NAMES.index("NE_pore_NPC"), "rpol": nucl, "r90": nucl, "rnuc": we.R_NUCLEOSOME}
    return np.concatenate([pre[f"{k}_reg"].astype(np.int64) if f"{k}_reg" in pre else np.full(len(pre[f"{k}_pos"]), fixed[k])
                           for k, _ in wc3d.PRE_GROUPS if f"{k}_pos" in pre])


def pre_anchored(pre):
    """per pre-placed instance: anchored structure (filament, lattice, membrane or DNA assembly) held fixed in the
    viewer's simulations; polysome ribosomes are not (they move, without their mRNA tether)"""
    return np.concatenate([np.full(len(pre[f"{k}_pos"]), k != "poly") for k, _ in wc3d.PRE_GROUPS if f"{k}_pos" in pre])


MARK_ON = 0.5       # histone-mark bit set where the nucleosome's Weiner 2015 value (log2 enrichment) is >= this


def pre_marks(pre, gmap, log=print):
    """per pre-placed instance: bits of the histone marks (genome_map.MARKS order) enriched on that nucleosome; 0 for
    everything else. Shipped in the instance flags from bit 2 (the viewer's 'Histone mark' colour modes)."""
    import genome_map
    out = []
    for k, _ in wc3d.PRE_GROUPS:
        if f"{k}_pos" not in pre:
            continue
        b = np.zeros(len(pre[f"{k}_pos"]), np.uint16)
        if k == "nuc" and gmap is not None and all(m in gmap for m in genome_map.MARKS):
            g = gmap.sort_values("nuc")
            assert len(g) == len(b)
            for i, m in enumerate(genome_map.MARKS):
                on = (g[m].values >= MARK_ON)
                b |= (on.astype(np.uint16) << i)
                log(f"  {m}: {on.mean():.0%} of nucleosomes enriched (>= {MARK_ON} log2)")
        out.append(b)
    return np.concatenate(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tiles-only", action="store_true")
    ap.add_argument("--annot-only", action="store_true", help="only the polymer identities (polymer_annot.py) for the tiles already encoded")
    ap.add_argument("--bookmarks-only", action="store_true", help="only the guided views (bookmarks.json) for a build already encoded")
    args = ap.parse_args()
    meta = dict(np.load(f"{C.W3D}/libpack_meta.npz"))
    ids = [str(x) for x in meta["ids"]]; NT = len(ids)
    cell = c3.Cell.load(f"{C.W3D}/geometry.json")
    tiles = we.tile_index_list(); tindex = {t: i for i, t in enumerate(tiles)}
    if args.bookmarks_only:
        bookmarks(tiles, ids, cell, pd.read_parquet(f"{C.W3D}/tile_labvol.parquet")); log("bookmarks ->", we.OUT)
        return
    pre = dict(np.load(f"{C.W3D}/prepass.npz"))
    gt, pp, gq = wc3d.pre_instances(pre)        # NPCs, nucleosomes, polysome ribosomes, actin, microtubules
    n_act = len(pre["act_pos"]) if "act_pos" in pre else 0
    pre_tile = np.array([tindex[tuple(int(v) for v in np.floor(p / we.TILE))] for p in pp])
    pl = polymer_lists(pre, tindex, log)
    import genome_map           # genome coordinates of the chromatin (marks for the instance flags; gene loci below)
    import polymer_annot        # identities of the polymer chains and chromatin instances (viewer hover)
    if args.annot_only:
        polymer_annot.write(we.OUT, pre, tiles, pl, pre_tile, gt, pre_regions(pre), log=log)
        return
    gmap = genome_map.build(log=lambda *a: log(*a)) if "nuc_pos" in pre else None
    shutil.rmtree(C.WEB_TMP, ignore_errors=True); os.makedirs(C.WEB_TMP)
    shutil.rmtree(we.OUT, ignore_errors=True); os.makedirs(f"{we.OUT}/tiles")
    we._S.update(tiles=tiles, tset=set(tiles), cell=cell, BR=meta["BR"].astype(np.float64), pre_tile=pre_tile,
                 pre_type=gt.astype(np.int64), pre_pos=pp.astype(np.float64), pre_q=gq.astype(np.float64),
                 pre_reg=pre_regions(pre), pre_fixed=pre_anchored(pre), pre_marks=pre_marks(pre, gmap, log=log),
                 poly_mrna=pl["mrna"], poly_mrna_tile=pl["mrna_tile"], extra=pl["extra"], extra_tile=pl["extra_tile"],
                 extra_kind=pl["extra_kind"], links=pl["links"], link_tile=pl["link_tile"], tmp=C.WEB_TMP, out=we.OUT)
    t0 = time.time()
    with mp.get_context("fork").Pool(C.NPROC) as pool:
        st1 = pool.map(we._stage1, range(len(tiles)), chunksize=8)
    counts_by_type = np.zeros(NT, np.int64); region_counts = {}
    for x in st1:
        for k, c in zip(x[3], x[4]):
            t, rg = divmod(int(k), 32); counts_by_type[t] += c; region_counts.setdefault(t, {})[rg] = region_counts.get(t, {}).get(rg, 0) + int(c)
    log(f"stage1 {time.time() - t0:.0f}s; instances {counts_by_type.sum() / 1e6:.3f}M")

    t0 = time.time()
    with mp.get_context("fork").Pool(C.NPROC) as pool:
        st2 = pool.map(we._stage2, range(len(tiles)), chunksize=4)
    st2 = pd.DataFrame(st2, columns=IDX_COLS).sort_values("ti").reset_index(drop=True)
    log(f"stage2 {time.time() - t0:.0f}s; gz {st2.gz_bytes.sum() / 1e6:.0f} MB")
    roundtrip(tiles, pl["links"], pl["link_tile"])
    if args.tiles_only:
        return

    ing, _ = load_inventory()
    log("library", we.encode_library(meta, ing, counts_by_type, region_counts))
    import trna_species         # predicted tRNA models, opened by the viewer's "Load atomic model" (src 'pred')
    pred = [l for l, s in zip(ids, meta["src"]) if str(s) == "pred"]
    if pred:
        os.makedirs(f"{we.OUT}/models", exist_ok=True)
        for l in pred:
            shutil.copy(f"{trna_species.MODELS}/{l}.cif", f"{we.OUT}/models/{l}.cif")
        log(f"models: {len(pred)} predicted structures copied for the viewer")
    import simdata      # residue charges and native contacts for the viewer's simulation mode
    log("simulation data", simdata.build(we.OUT, log=lambda *a: log(*a)))
    if gmap is not None:        # gene loci for the viewer's search
        genome_map.write_loci(we.OUT, gmap, log=lambda *a: log(*a))
    polymer_annot.write(we.OUT, pre, tiles, pl, pre_tile, gt, pre_regions(pre), log=log)
    if os.path.exists(f"{og.ODIR}/network.npz"):          # organisation layer of s3 (network colouring in the viewer)
        og.encode_viewer(we.OUT, log=log)
    labvol = pd.read_parquet(f"{C.W3D}/tile_labvol.parquet")
    index = tile_index(tiles, st2, labvol)
    types_tiles(tiles, NT)
    web_geometry(cell)
    organelle_meshes(tiles)
    bookmarks(tiles, ids, cell, labvol)
    index["wall"] = wall_budget(tiles, cell, ing, labvol)
    for i, row in enumerate(index["tiles"]):
        row += [round(float(labvol.iloc[i]["wall_glucan"]), 1), round(float(labvol.iloc[i]["wall_mannan"]), 1)]
    index["tile_columns"] += ["vol_glucan_nm3", "vol_mannan_nm3"]
    json.dump(index, open(f"{we.OUT}/tiles/index.json", "w"), separators=(",", ":"))
    log("done ->", we.OUT)


def polymer_lists(pre, tindex, log=print):
    """The polylines written into the tiles, in tile order, each with its tile (midpoint; -1 outside) and its source
    (polymer_annot.py turns sources into identities for the viewer):
      mrna   polysome mRNA, one piece per ribosome (channel + linker to the next): src (0, ribosome, 0, 0); full-length
             mRNA threads in MRNA_PIECE-nt pieces: src (1, thread, first bead, end bead)
      extra  FG repeats (kind 5), nascent pre-rRNA (4) and nascent Pol II transcripts (0) in MRNA_PIECE-bead pieces:
             src (kind, thread, first bead, end bead)
      links  DNA: linkers between nucleosomes (0, pair i) and between rDNA nucleosomes (1, pair i); active rRNA genes
             (2, gene), DNA through Pol II (3, piece) and mtDNA (4, piece) in pieces of MRNA_PIECE + 1 beads (+ a, b)"""
    def tile_of(p):
        return tindex.get(tuple(int(v) for v in np.floor(p / we.TILE)), -1)
    mrna, src = [], []
    if "poly_pos" in pre and len(pre["poly_pos"]):
        off = pre["poly_mrna_off"]
        mrna = [pre["poly_mrna"][off[i]:off[i + 1]].astype(np.float64) for i in range(len(off) - 1)]
        src = [(0, i, 0, 0) for i in range(len(mrna))]
    if "mrna_pts" in pre:   # full-length mRNA threads: cut into MRNA_PIECE-nt pieces, each assigned to its midpoint's tile
        mpts, mo = pre["mrna_pts"].astype(np.float64), pre["mrna_off"]
        n0 = len(mrna)
        for i in range(len(mo) - 1):
            for a in range(mo[i], mo[i + 1], MRNA_PIECE):
                b = min(mo[i + 1], a + MRNA_PIECE)
                if b - a >= 2:
                    mrna.append(mpts[a:b]); src.append((1, i, a - mo[i], b - mo[i]))
        log(f"mRNA threads: {len(mo) - 1} threads -> {len(mrna) - n0} pieces")
    mrna_tile = np.array([tile_of(p[len(p) // 2]) for p in mrna], np.int64)
    extra, extra_kind, xsrc = [], [], []
    for key, kind in (("fg", 5), ("prerrna", 4), ("pol2rna", 0)):
        if f"{key}_pts" in pre:
            P_, O_ = pre[f"{key}_pts"].astype(np.float64), pre[f"{key}_off"]
            for i in range(len(O_) - 1):
                for a in range(O_[i], O_[i + 1], MRNA_PIECE):
                    b = min(O_[i + 1], a + MRNA_PIECE)
                    if b - a >= 2:
                        extra.append(P_[a:b]); extra_kind.append(kind); xsrc.append((kind, i, a - O_[i], b - O_[i]))
    extra_tile = np.array([tile_of(p[len(p) // 2]) for p in extra], np.int64)
    log(f"mRNA pieces {len(mrna)} ({(mrna_tile < 0).sum()} outside tiles); extra polylines {len(extra)} ({(extra_tile < 0).sum()} outside tiles)")
    zn = np.load(f"{we.LIB}/nucleosome.npz")
    links, mids, pair = we.linker_chains(pre, zn["beads"].astype(np.float64), zn["kind"])
    lsrc = [(0, i, 0, 0) for i in pair]
    if "rnuc_pos" in pre:      # nucleosomal rDNA in the nucleolus
        l2, m2, p2 = we.linker_chains({k[1:]: v for k, v in pre.items() if k.startswith("rnuc_")}, zn["beads"].astype(np.float64), zn["kind"])
        links += l2; mids = np.concatenate([mids, m2]) if len(m2) else mids; lsrc += [(1, i, 0, 0) for i in p2]
    for key, code in (("rdna", 2), ("pol2dna", 3), ("mtdna", 4)):
        if f"{key}_pts" not in pre:
            continue
        # nucleosome-free DNA of the active rRNA genes, DNA through the transcribing Pol II, mtDNA nucleoids
        P_, O_ = pre[f"{key}_pts"].astype(np.float64), pre[f"{key}_off"]
        n0 = len(links)
        for i in range(len(O_) - 1):
            for a in range(O_[i], O_[i + 1] - 1, MRNA_PIECE):
                b = min(O_[i + 1], a + MRNA_PIECE + 1)
                links.append(P_[a:b]); lsrc.append((code, i, a - O_[i], b - O_[i]))
        mids = np.concatenate([mids, np.array([p[len(p) // 2] for p in links[n0:]])])
        log(f"{key}: {len(O_) - 1} pieces -> {len(links) - n0} DNA pieces")
    link_tile = np.array([tile_of(p) for p in mids], np.int64)
    as_src = lambda s: np.array(s, np.int64).reshape(-1, 4)
    return dict(mrna=mrna, mrna_tile=mrna_tile, mrna_src=as_src(src), extra=extra, extra_kind=np.array(extra_kind, np.int64),
                extra_tile=extra_tile, extra_src=as_src(xsrc), links=links, link_tile=link_tile, link_src=as_src(lsrc))


def roundtrip(tiles, links, link_tile):
    """decode 40 tiles with the reference decoder and compare against the raw packing output"""
    rng = np.random.default_rng(3); perr, aerr = [0.0], [0.0]
    for ti in rng.choice(len(tiles), 40, replace=False):
        t = tiles[ti]; name = t3.tile_name(t)
        d = we.decode_tile(gzip.decompress(open(f"{we.OUT}/tiles/{t[0]}_{t[1]}_{t[2]}.bin.gz", "rb").read()))
        z = np.load(f"{C.WEB_TMP}/{name}.npz"); o = np.lexsort((z["region"], z["types"]))
        ty, po, qu, rg, fl = d["inst"]
        assert (ty == z["types"][o]).all() and (rg == z["region"][o]).all() and (fl == z["flags"][o]).all()
        perr.append(np.abs(po - z["pos"][o]).max(initial=0))
        if len(qu):
            q0 = z["quat"][o]; dots = np.abs((qu * q0 / np.linalg.norm(q0, axis=1, keepdims=True)).sum(1))
            aerr.append(np.degrees(2 * np.arccos(np.clip(dots, 0, 1))).max())
    log(f"round trip: max pos err {max(perr):.4f} nm, max rot err {max(aerr):.4f} deg")
    assert max(perr) < 0.01 and max(aerr) < 0.1


def tile_index(tiles, st2, labvol):
    labmask = ((labvol.values > 0).astype(np.int64) << np.arange(labvol.shape[1])).sum(1)
    return dict(tile_nm=we.TILE, q_org=we.Q_ORG, q_range=we.Q_RANGE, mrna_q=we.MR_Q, link_q=we.LK_Q,
                label_names=c3.LABEL_NAMES, region_names=t3.REGION_NAMES,
                mem_names=c3.MEM_NAMES, label_columns=list(labvol.columns),
                tiles=[[*tiles[i], int(r.n_inst), int(r.n_ghost), int(r.n_mrna_beads), int(r.n_link_beads), int(r.n_vert),
                        int(r.n_tri), int(r.gz_bytes), int(labmask[i])] for i, r in st2.iterrows()],
                tile_columns=["ix", "iy", "iz", "n_inst", "n_ghost", "n_mrna_beads", "n_link_beads", "n_vert", "n_tri", "gz_bytes", "label_mask"])


def types_tiles(tiles, NT):
    tt = [[] for _ in range(NT)]
    for i, t in enumerate(tiles):
        for ty in np.unique(np.load(f"{C.WEB_TMP}/{t3.tile_name(t)}.npz")["types"]): tt[ty].append(i)
    offs = np.zeros(NT + 1, np.uint32); offs[1:] = np.cumsum([len(x) for x in tt])
    with open(f"{we.OUT}/types_tiles.bin", "wb") as f:
        f.write(offs.astype("<u4").tobytes()); f.write(np.concatenate([np.array(x, "<u2") for x in tt]).tobytes())


WEB_BANDS = [(180.0, 40.0), (40.0, 40.0), (100.0, 100.0)]   # viewer bands (in, out) of the cell, vacuole, nucleus fields
FIELD_Q = 0.01                                               # nm per int16 step


def web_fields(cell):
    """the real-shape distance fields for the viewer (geometry_fields.bin; layout in web/app/src/fields.ts):
    split into bricks of half the edge (less padding around the thin bands), re-classified for narrower bands, node
    values as int16 in FIELD_Q steps. None for an analytic cell."""
    if not cell.fields:
        return None
    heads, n_index, n_data, parts = [], [], [], []
    off = 1 + 12 * len(cell.fields)
    for (h, idx, data), (bin_, bout) in zip(cell.fields, WEB_BANDS):
        B = int(h[7]); n1 = B + 1; nbx, nby, nbz = int(h[4]), int(h[5]), int(h[6]); b2 = B // 2
        blocks = data.reshape(-1, n1, n1, n1)
        old = idx.astype(np.int64).reshape(nbx, nby, nbz)
        new = np.empty((2 * nbx, 2 * nby, 2 * nbz), np.int64); keep = []
        for bx in range(2 * nbx):
            for by in range(2 * nby):
                for bz in range(2 * nbz):
                    o = old[bx // 2, by // 2, bz // 2]
                    if o < 0:
                        new[bx, by, bz] = o; continue
                    i0, j0, k0 = (bx % 2) * b2, (by % 2) * b2, (bz % 2) * b2
                    blk = blocks[o, i0:i0 + b2 + 1, j0:j0 + b2 + 1, k0:k0 + b2 + 1]
                    if blk.max() <= -bin_:
                        new[bx, by, bz] = -1
                    elif blk.min() >= bout:
                        new[bx, by, bz] = -2
                    else:
                        new[bx, by, bz] = len(keep); keep.append(np.clip(blk, -bin_, bout).ravel())
        idx = new.ravel()
        dat = np.round(np.concatenate(keep) / FIELD_Q).astype(np.int16) if keep else np.zeros(0, np.int16)
        hd = [float(h[0]), float(h[1]), float(h[2]), float(h[3]), 2.0 * nbx, 2.0 * nby, 2.0 * nbz, float(b2),
              bin_, bout, float(off), float(off + len(idx))]
        off += len(idx) + len(dat)
        heads.append(hd); n_index.append(int(len(idx))); n_data.append(int(len(dat)))
        pad = (-2 * len(dat)) % 4
        parts += [idx.astype(np.int32).tobytes(), dat.tobytes() + b"\0" * pad]
    open(f"{we.OUT}/geometry_fields.bin", "wb").write(b"".join(parts))
    log(f"geometry fields for the viewer: {sum(n_data) / 1e6:.1f} M nodes, "
        f"{os.path.getsize(f'{we.OUT}/geometry_fields.bin') / 1e6:.1f} MB")
    return dict(file="geometry_fields.bin", head=heads, n_index=n_index, n_data=n_data, q=FIELD_Q)


def web_geometry(cell):
    """primitives + uniform acceleration grid for the WGSL SDF"""
    geo = json.load(open(f"{C.W3D}/geometry.json"))
    GN, EXT = 36, 2600.0; h = 2 * EXT / GN       # encloses the real cell outline (viewer tile grid: +-2.6 um)
    bc, br = cell.bc, cell.br
    pc = cell.pores[:, :3]; pr = c3.PORE_R + c3.PERI_W + 12
    cells_p, cells_q = [], []
    for ix in range(GN):
        for iy in range(GN):
            for iz in range(GN):
                lo = -EXT + h * np.array([ix, iy, iz]); hi = lo + h
                q = np.clip(bc, lo, hi); cells_p.append(np.where(((q - bc) ** 2).sum(1) <= (br + 12) ** 2)[0])
                q = np.clip(pc, lo, hi); cells_q.append(np.where(((q - pc) ** 2).sum(1) <= pr ** 2)[0])
    goff = np.zeros(GN ** 3 + 1, np.int64); goff[1:] = np.cumsum([len(a) for a in cells_p])
    qoff = np.zeros(GN ** 3 + 1, np.int64); qoff[1:] = np.cumsum([len(a) for a in cells_q])
    fields = web_fields(cell)
    gweb = dict(G=geo["G"], fields=fields, pores=cell.pores.round(6).ravel().tolist(), PT=geo["PT"], PQ=np.round(cell.PQ, 4).ravel().tolist(), PC=geo["PC"],
                consts=geo["consts"], label_names=c3.LABEL_NAMES, mem_names=c3.MEM_NAMES, npc_pos=geo["npc_pos"], mito_segments=geo["mito_segments"],
                grid=dict(n=GN, ext=EXT, prim_off=goff.tolist(), prim_idx=np.concatenate(cells_p).tolist(),
                          pore_off=qoff.tolist(), pore_idx=np.concatenate(cells_q).tolist()))
    json.dump(gweb, open(f"{we.OUT}/geometry.json", "w"), separators=(",", ":"))


def organelle_meshes(tiles):
    import fast_simplification as fsimp
    Vs = {m: [] for m in range(10)}; Fs = {m: [] for m in range(10)}; nv = {m: 0 for m in range(10)}
    for t in tiles:
        g = np.load(f"{we.GDIR}/{t3.tile_name(t)}.npz")
        if not len(g["faces"]): continue
        V, F, FM = g["verts"], g["faces"], g["face_mem"]
        for m in np.unique(FM):
            f = F[FM == m]; u, inv = np.unique(f.ravel(), return_inverse=True)
            Vs[m].append(V[u]); Fs[m].append(inv.reshape(-1, 3) + nv[m]); nv[m] += len(u)
    mesh_out = []
    for m in [m for m in range(10) if Fs[m]]:
        V = np.concatenate(Vs[m]).astype(np.float64); F = np.concatenate(Fs[m])
        _, first, inv = np.unique(np.round(V / 0.01).astype(np.int64), axis=0, return_index=True, return_inverse=True)
        V2 = V[first]; F2 = inv.ravel()[F]
        F2 = F2[(F2[:, 0] != F2[:, 1]) & (F2[:, 1] != F2[:, 2]) & (F2[:, 0] != F2[:, 2])]
        Vd, Fd = fsimp.simplify(V2.astype(np.float32), F2.astype(np.int32), target_reduction=0.95 if m in (0, 3) else 0.9)
        mesh_out.append((m, Vd.astype("<f4"), Fd.astype("<u4")))
    with gzip.open(f"{we.OUT}/organelle_meshes.bin.gz", "wb", compresslevel=9) as f:
        hdr = np.zeros(1 + 3 * len(mesh_out), np.uint32); hdr[0] = len(mesh_out)
        for i, (m, Vd, Fd) in enumerate(mesh_out): hdr[1 + 3 * i: 4 + 3 * i] = [m, len(Vd), len(Fd)]
        f.write(hdr.tobytes())
        for m, Vd, Fd in mesh_out: f.write(Vd.tobytes()); f.write(Fd.tobytes())


def bookmarks(tiles, ids, cell, labvol):
    """hand-placed guided views (tied to geometry seed 11)"""
    unit = lambda v: np.asarray(v, float) / np.linalg.norm(v)
    perp = lambda a, hint=(0.3, 0.5, 0.81): unit(np.cross(a, hint))
    geo = json.load(open(f"{C.W3D}/geometry.json"))
    li = {k: ids.index(k) for k in ["Pma1_hexamer", "ribosome_80S", "nucleosome", "NPC"]}
    tcount = {k: np.zeros(len(tiles), int) for k in li}
    for i, t in enumerate(tiles):
        ty = np.load(f"{C.WEB_TMP}/{t3.tile_name(t)}.npz")["types"]
        for k, v in li.items(): tcount[k][i] = (ty == v).sum()
    lv = labvol.values; lcol = list(labvol.columns)
    tc = lambda i: np.array(tiles[i], float) * 200 + 100
    G = cell.G; cn = G[3:6]
    vac = [i for i, c in enumerate(cell.PC) if c == c3.C_VAC]
    if vac:                                                   # analytic cell: the largest vacuole lobe
        big = max(vac, key=lambda i: cell.PQ[i, 3]); vc = cell.PQ[big, :3]; rv = cell.PQ[big, 3]
        vac_target, vac_u = vc + unit(-vc) * rv, unit(-vc)
    else:                                                     # real shapes: vacuole membrane point nearest the centre
        Pr = np.random.default_rng(3).uniform(-G[0], G[0], (400000, 3))
        dv = np.array([c3.field_sdf(p[0], p[1], p[2], G, c3.F_VAC) for p in Pr])
        vc = Pr[dv < 0].mean(0)
        near = Pr[(np.abs(dv) < 4) & (np.array([c3.outer_depth(*p, G) for p in Pr]) > 250)]     # (cell field clamped at 260)
        clear = np.array([c3.field_sdf(p[0], p[1], p[2], G, c3.F_NUC) for p in near]) > 125.0   # not at the junction with the nucleus (field clamped at 130)
        near = near[clear] if clear.any() else near
        vac_target = near[np.argmin(np.linalg.norm(near, axis=1))]
        vac_mem = Pr[np.abs(dv) < 4]                                                            # for the Cvt view
        vac_u = -c3._grad(lambda x, y, z: c3.field_sdf(x, y, z, G, c3.F_VAC), vac_target)
    n_inst = sum(t["copies"] for t in json.load(open(f"{we.OUT}/library/types.json"))["types"])
    bm = [dict(name="Whole cell", caption=f"Equatorial section through nucleus and vacuole. {n_inst / 1e6:.1f} million macromolecules, every one a real PDB or AlphaFold structure.",
               target=[0, 0, 0], normal=unit(np.cross(cn, vc)), up=unit(cn), width=5200, slab=60)]
    scars = [np.asarray(cell.PQ[i][:3], float) for i in range(len(cell.PT)) if cell.PC[i] == c3.C_SCAR]
    off_scar = lambda p: all(np.linalg.norm(np.asarray(p) - q) > 900.0 for q in scars)      # clear of the chitin ring
    u0 = unit(np.cross(cn, vc)); k0 = unit(cn); u = u0
    for ang in np.radians([0, 40, -40, 80, -80, 120, -120, 160]):       # turn about the nucleus direction until clear
        u = u0 * np.cos(ang) + np.cross(k0, u0) * np.sin(ang) + k0 * (k0 @ u0) * (1 - np.cos(ang))
        if off_scar(u * c3.ray_depth(G, u, 55.0)):
            break
    bm.append(dict(name="Cell wall", caption="Mannoprotein brush (outer 35 nm) over the beta-glucan layer, plasma membrane below.", target=u * c3.ray_depth(G, u, 55.0), normal=perp(u), up=u, width=260, slab=40))
    cer_col = lcol.index("cER_lumen")
    cand = [i for i in np.argsort(-tcount["Pma1_hexamer"]) if lv[i, cer_col] == 0][:1]
    u = unit(tc(cand[0]))
    bm.append(dict(name="Plasma membrane (face-on)", caption="Cytoplasmic face of the plasma membrane: Pma1 H+-ATPase hexamers among ~1,000 other membrane proteins; lipids generated around them.",
                   target=u * c3.ray_depth(G, u, G[0] - G[2] + 18.0), normal=u, up=perp(u), width=220, slab=30))
    cer = [i for i, c in enumerate(cell.PC) if c == c3.C_CER]
    cer = [i for i in cer if off_scar(cell.PQ[i][:3])] or cer
    q = cell.PQ[cer[0]]; u = unit(q[:3])
    tgt = u * c3.ray_depth(G, u, G[0] - G[2] + c3.CER_GAP + 0.5 * c3.CER_LUMEN)    # on the sheet, under the real outline
    bm.append(dict(name="Cortical ER", caption="Cortical ER sheet 27 nm beneath the plasma membrane, with Sec61 translocons and ER-bound ribosomes.", target=tgt, normal=perp(u), up=u, width=320, slab=40))
    cyto = lcol.index("cytosol")
    i = [i for i in np.argsort(-tcount["ribosome_80S"]) if lv[i, cyto] > 0.99 * lv[i].sum()][0]
    bm.append(dict(name="Cytoplasm", caption="Crowded cytoplasm (~93 mg/mL): 80S ribosomes (pink), most of them threaded on mRNAs as polysomes, with glycolytic enzymes and tRNAs.", target=tc(i), normal=unit([0.2, -0.4, 0.89]), up=unit([0, 1, 0.45]), width=160, slab=30))
    pre = dict(np.load(f"{C.W3D}/prepass.npz"))
    if "act_pos" in pre:   # longest actin filament: view the cortex from inside, cut plane through the cable
        f = np.bincount(pre["act_fil"]).argmax(); P = pre["act_pos"][pre["act_fil"] == f]
        mid = P[len(P) // 2]; u = unit(mid); tdir = unit(P[min(len(P) - 1, len(P) // 2 + 20)] - P[max(0, len(P) // 2 - 20)])
        bm.append(dict(name="Actin cable", caption="Actin cable beneath the cortex: a bundle of F-actin filaments (orange; 9GO5 subunit and helix: 167 deg, 2.76 nm per subunit) coated by tropomyosin TPM1 (teal) along both strands, running past polysomes and cytosolic proteins.",
                       target=mid, normal=u, up=tdir, width=400, slab=40))
    if "mt_pos" in pre:    # SPB: section through the plaque axis, nuclear MT bouquet below, cytoplasmic MTs above
        spb, sd = c3.nuc_surface(G, cn - G[7:10])
        nuc_mt = pre["mt_reg"] == t3.R_NUCP
        side = unit(np.cross(sd, unit(pre["mt_pos"][nuc_mt].mean(0) - spb) + [0.1, 0.2, 0.3]))
        bm.append(dict(name="Spindle pole body", caption="Microtubules (green; 8QV3 lattice, 13 protofilaments) grow from gamma-tubulin rings on the spindle pole body: "
                       "a bouquet of short nuclear microtubules (16 reach the centromeres, Rabl configuration) and a few long cytoplasmic ones.",
                       target=spb, normal=side, up=sd, width=900, slab=60))
        cyt = pre["mt_kind"] == "cMT"
        if cyt.any():
            j = int(np.argmax(np.where(cyt, pre["mt_len"], 0))); sub = pre["mt_id"] == j
            P = pre["mt_pos"][sub]; mid = P[len(P) // 2]; tdir = unit(P[-1] - P[0])
            bm.append(dict(name="Cytoplasmic microtubule", caption="A cytoplasmic microtubule, 25 nm wide and about a micrometre long, curving through the crowded cytoplasm.",
                           target=mid, normal=unit(np.cross(tdir, unit(mid))), up=tdir, width=700, slab=40))
    if "pat_pos" in pre and len(pre["pat_pos"]):   # largest actin patch, seen from the side (membrane at the top)
        j = int(np.bincount(pre["pat_id"]).argmax())
        es = cell.endocytic_sites()
        if es:              # the largest patch around one of the deeper invaginations
            S_ = np.array([p for p, _, _ in es]); D_ = np.array([d for _, _, d in es]); size = np.bincount(pre["pat_id"])
            cen = np.array([pre["pat_pos"][pre["pat_id"] == k].mean(0) for k in range(len(size))])
            deep = np.array([D_[np.linalg.norm(S_ - c_, axis=1).argmin()] for c_ in cen])
            j = int(np.argmax(np.where(deep > 60.0, size, 0))) if (deep > 60.0).any() else j
        P = pre["pat_pos"][pre["pat_id"] == j]; c = P.mean(0); u = -c3._grad(lambda x, y, z: c3.outer_depth(x, y, z, G), c)
        if es:              # cut plane through the invagination axis
            s0 = S_[np.linalg.norm(S_ - c, axis=1).argmin()]; c = s0 + ((c - s0) @ u) * u
        bm.append(dict(name="Actin patch", caption="Endocytic site: a plasma-membrane invagination (a tube about 40 nm wide that grows to ~100 nm before "
                       "the vesicle pinches off; coat proteins Sla2, Sla1, End3 in gold at its tip, amphiphysin Rvs161/167 in violet on its tubular part) inside an actin patch: short actin filaments (orange) branching at 73 deg from "
                       "Arp2/3 complexes (red; 8UXW junction) and capped by capping protein (purple; 8F8Q), barbed ends toward the membrane.",
                       target=c, normal=perp(u), up=u, width=360, slab=50))
    if "eis_pos" in pre and len(pre["eis_pos"]):   # eisosome furrow with the most coat, in cross-section
        j = int(np.bincount(pre["eis_id"]).argmax()); i = cell.furrows()[j]; q = cell.PQ[i]
        a, b = q[:3], q[3:6]; u = unit(a + b); e = unit(b - a)
        bm.append(dict(name="Eisosome", caption="Eisosome: a plasma-membrane furrow about 50 nm deep and 300 nm long (MCC), its half-pipe bottom lined "
                       "by a lattice of Pil1 (dark blue) and Lsp1 (light blue) BAR-domain dimers in their native arrangement (8QB7, 8QB8; "
                       "Kefauver et al. 2024), amphipathic helices inserted into the inner leaflet. About 30 per cell.",
                       target=0.5 * (a + b) + u * 15.0, normal=e, up=u, width=160, slab=60))
    if "cvt_pos" in pre and len(pre["cvt_pos"]):    # Cvt complex next to the vacuole
        c = pre["cvt_centre"].astype(np.float64); u = unit(c)
        if c3.n_fields(G) > c3.F_VAC:      # towards the vacuole: up in the picture, so that the membrane is in view
            gv = np.asarray(c3._grad(lambda x, y, z: c3.field_sdf(x, y, z, G, c3.F_VAC), c), float)
            if np.linalg.norm(gv) > 0.3:
                u = unit(-gv)
        width = 320
        if not vac:         # real shapes: frame the complex and the nearest stretch of vacuolar membrane together
            pm = vac_mem[np.argmin(np.linalg.norm(vac_mem - c, axis=1))]; dvm = float(np.linalg.norm(pm - c))
            if dvm < 600.0:
                u = unit(pm - c); c = c + 0.5 * (pm - c); width = max(320, int(dvm + 220))
        bm.append(dict(name="Cvt complex", caption="Cvt complex: precursor aminopeptidase I dodecamers (Ape1, amber) packed with alpha-mannosidase Ams1 (brown) "
                       "and coated by the receptor Atg19 (red), waiting next to the vacuole to be wrapped into a Cvt vesicle "
                       "(the cytoplasm-to-vacuole targeting pathway of growing cells).",
                       target=c, normal=perp(u), up=u, width=width, slab=60))
    if "nvj_pos" in pre and len(pre["nvj_pos"]):   # nucleus-vacuole junction, membranes edge-on
        c = pre["nvj_centre"].astype(np.float64); u = unit(c - cn)
        bm.append(dict(name="Nucleus-vacuole junction", caption="Nucleus-vacuole junction: where the vacuole lies against the nucleus, Vac8 on "
                       "the vacuolar membrane binds the outer-nuclear-membrane protein Nvj1 (teal; Vac8-Nvj1 structure 5XJG). "
                       "Nvj1's membrane anchor and linker are not resolved and not drawn.",
                       target=c, normal=perp(u), up=u, width=300, slab=40))
    sc = [i for i in range(len(cell.PT)) if cell.PC[i] == c3.C_SCAR]
    if sc:                  # bud scar, sectioned through its axis: the chitin rim on both sides, the scar's face beyond
        q = cell.PQ[sc[0]]; ax = unit(q[3:6])
        bm.append(dict(name="Bud scar", caption="Bud scar: the chitin ring a mother cell keeps from each division, about 1 um across, its rim (pale "
                       "microfibrils, half the cell's chitin) raised above the mannoprotein surface. Shown in section through the ring; "
                       "a daughter's birth scar has no chitin rim and is not modelled.",
                       target=q[:3] - ax * 60.0, normal=perp(ax), up=ax, width=1300, slab=60))
    pcen = cell.pores[:, :3]
    k = int(np.argmax(np.linalg.norm(pcen, axis=1))); ax = unit(cell.pores[k, 3:6])
    bm.append(dict(name="Nuclear pore complex", caption="NPC (PDB 7N9F) seated in a toroidal pore membrane joining the outer and inner nuclear membranes.", target=pcen[k], normal=perp(ax), up=ax, width=260, slab=50))
    nucp = lcol.index("nucleoplasm")
    i = [i for i in np.argsort(-tcount["nucleosome"]) if lv[i, nucp] > 0.95 * lv[i].sum()][0]
    bm.append(dict(name="Chromatin", caption="Nucleoplasm: nucleosome chains (73,150 nucleosomes, linker DNA) with RNA polymerase II and nuclear proteins.", target=tc(i), normal=unit([0.5, 0.3, -0.81]), up=unit([0, 1, 0.37]), width=200, slab=40))
    if "seg_start" in pre:   # Rabl layout: a thick slab through the SPB axis, chromatin and microtubules only
        spb = pre["seg_start"][[k for k, x in enumerate(pre["seg_kind"]) if str(x) == "L"]].mean(0); ax = unit(spb - cn)
        bm.append(dict(name="Chromosome arms", caption="Rabl configuration: the 16 centromeres held at the spindle pole body by kinetochore "
                       "microtubules (green), each chromosome arm running from there to its telomere at the nuclear envelope; telomeres "
                       "gather in 7 foci, longer arms reaching farther from the pole. Chromosome XII passes through the nucleolus "
                       "(its rDNA). Shown as a thick slab, nucleosomes, microtubules and pores only.",
                       target=cn, normal=perp(ax), up=ax, width=2600, slab=500, show=["Nucleosomes", "Microtubules", "Nuclear pores", "DNA"]))
    nucl = lcol.index("nucleolus")
    # the crescent on its axis: nuclear envelope (and the vacuole beyond it) at the top, nucleoplasm at the bottom
    uu = unit(G[7:10] - cn); env = c3.nuc_surface(G, uu)[0]; inner = G[7:10] - uu * G[10]
    span = float(np.linalg.norm(env - inner))
    bm.append(dict(name="Nucleolus", caption="Nucleolar crescent against the nuclear envelope on the vacuole side: ribosomes being assembled (90S and pre-60S particles with their pre-rRNA), snoRNPs and Pol I, bordering the nucleoplasm.",
                   target=0.5 * (env + inner), normal=perp(uu), up=uu, width=int(max(450, 1.5 * span)), slab=40))
    if "rpol_pos" in pre and len(pre["rpol_pos"]):   # the densest stretch of active rRNA genes, seen side-on in a thick slab
        from scipy.spatial import cKDTree
        P = pre["rpol_pos"].astype(np.float64); nn = cKDTree(P).query_ball_point(P, 120.0, return_length=True)
        near = P[np.linalg.norm(P - P[int(np.argmax(nn))], axis=1) < 120.0]; c = near.mean(0)
        ax = np.linalg.eigh(np.cov((near - c).T))[1][:, 2]; nrm = perp(ax)
        bm.append(dict(name="rRNA genes", caption="Active rRNA genes in the nucleolus, other molecules hidden: RNA polymerase I (purple) every 45 nm along "
                       "the nucleosome-free DNA (orange), nascent pre-rRNA (magenta) growing longer along each gene, and 90S pre-ribosomes (coral) "
                       "on the 5' ends of uncleaved transcripts: the 'Christmas trees' of Miller spreads. Inactive repeats are nucleosomal.",
                       target=c - nrm * 40.0, normal=nrm, up=ax, width=360, slab=80,
                       show=["DNA", "pre-rRNA", "RNA polymerases", "Ribosomes", "Nucleosomes"]))
    seg = np.array(geo["mito_segments"][len(geo["mito_segments"]) // 2]); a = unit(seg[1] - seg[0])
    m_target, m_normal = seg.mean(0), perp(a)
    cr = [i for i in range(len(cell.PT)) if cell.PC[i] == c3.C_CRISTA]
    if cr:
        cc = np.array([cell.PQ[i][:3] for i in cr], float); sg = np.array(geo["mito_segments"], np.float64)
        nb = (np.linalg.norm(cc[:, None] - cc[None], axis=2) < 150.0).sum(1)
        k = int(np.argmax(nb)); a = unit(np.asarray(cell.PQ[cr[k]][3:6], float))
        j = int(np.argmin(np.linalg.norm(sg.mean(1) - cc[k], axis=1)))
        pt = sg[j, 0] + a * ((cc[k] - sg[j, 0]) @ a)                 # on the tubule axis, level with the crista
        o = cc[k] - pt; o -= a * (o @ a)
        m_target = pt
        m_normal = unit(np.cross(a, o)) if np.linalg.norm(o) > 5.0 else perp(a)
    bm.append(dict(name="Mitochondrion", caption="Mitochondrial tubule: outer and inner membranes, short lamellar cristae reaching in from two sides, lined at their rims by rows of ATP synthase dimers (8F39 monomers in the 4B2Q dimer arrangement), mitoribosomes tethered to the inner membrane with their exit tunnels facing it, and the protein-dense matrix (~100 mg/mL).", target=m_target, normal=m_normal, up=a, width=420, slab=40))
    if "mtdna_centres" in pre:   # the nucleoid with the most genomes, in a thick slab across its tubule
        k = int(np.argmax(pre["mtdna_genomes"])); c = pre["mtdna_centres"][k].astype(np.float64)
        sg = np.array(geo["mito_segments"], np.float64); m = sg.mean(1)
        j = int(np.argmin(np.linalg.norm(m - c, axis=1))); a = unit(sg[j, 1] - sg[j, 0])
        bm.append(dict(name="Mitochondrial nucleoid", caption=f"mtDNA nucleoid in the matrix: {int(pre['mtdna_genomes'][k])} copies of the 85.8 kb "
                       "mitochondrial genome (orange) folded between the cristae; a haploid cell has 10-20 nucleoids holding 25-50 genome copies. "
                       "The packaging protein Abf2 is in the surrounding matrix.",
                       target=c, normal=perp(a), up=a, width=360, slab=80, show=["DNA", "Membranes", "Ribosomes"]))
    u = vac_u
    bm.append(dict(name="Vacuole", caption="Vacuolar membrane with V-ATPases. The lumen holds little protein (~4 mg/mL) but is filled with polyphosphate (purple; short chains, about 300 mM phosphate units, drawn procedurally), the cell's phosphate store made by the VTC complex. " + ("A G1 cell has one near-spherical vacuole (soft X-ray tomography)." if not vac and cellshape.VAC_MODE == "single" else "The vacuole is a cluster of lobes, as in cryo-FIB/SEM of log-phase cells."), target=vac_target, normal=perp(u), up=u, width=300, slab=40))
    q = cell.PQ[[i for i, c in enumerate(cell.PC) if c == c3.C_LD][0]]
    bm.append(dict(name="Lipid droplet", caption="Lipid droplet: neutral-lipid core with a phospholipid monolayer and peripheral proteins, touching ER.", target=q[:3], normal=unit([0, 0, 1]), up=unit([0, 1, 0]), width=560, slab=40))
    q = cell.PQ[[i for i, c in enumerate(cell.PC) if c == c3.C_GOLGI][0]]
    bm.append(dict(name="Golgi / endosomes", caption="Golgi cisternae and endosomes with their membrane and lumenal proteins.", target=q[:3], normal=unit([0, 0, 1]), up=unit([0, 1, 0]), width=500, slab=40))
    for b in bm:
        for k in ("target", "normal", "up"): b[k] = [round(float(x), 3) for x in b[k]]
    json.dump(bm, open(f"{we.OUT}/bookmarks.json", "w"), indent=1)


def wall_budget(tiles, cell, ing, labvol):
    """wall polysaccharide budgets (protein 45% of dry mass, wall 20% of dry mass, Klis 2002 composition)"""
    prot_da = float((ing.copies_per_cell * ing.mass_da).sum())
    dry = prot_da / 0.45; wall = dry * 0.20
    glc_res = wall * (0.525 + 0.075) / 162.14; man_res = wall * 0.375 * 0.85 / 162.14
    vol_glucan = labvol["wall_glucan"].sum(); vol_mannan = labvol["wall_mannan"].sum()
    lv = labvol.values; lcol = list(labvol.columns); n_wp = 0
    for i, t in enumerate(tiles):
        if not (lv[i, lcol.index("wall_mannan")] > 0): continue
        z = np.load(f"{C.WEB_TMP}/{t3.tile_name(t)}.npz")
        k = z["region"] == t3.REGION_NAMES.index("cell_wall")
        n_wp += int(sum(c3.outer_depth(p[0], p[1], p[2], cell.G) < cell.G[0] - cell.G[1] + 5 for p in z["pos"][k]))
    # chitin (1.5% of the wall): the bud scars' share, in their rims (torus volume), 36 GlcNAc per microfibril bead
    chi_res = wall * 0.015 / 203.2
    tor = [cell.PQ[i] for i in range(len(cell.PT)) if cell.PC[i] == c3.C_SCAR]
    vol_scar = sum(2 * np.pi ** 2 * q[6] * (q[7] - 1.2) ** 2 for q in tor)
    chi = dict(chitin_glcnac=chi_res, chitin_in_scars=c3.SCAR_CHITIN, bud_scars=len(tor),
               chitin_beads_per_nm3_scar=chi_res * c3.SCAR_CHITIN / 36 / vol_scar if vol_scar else 0.0)
    return dict(**chi, dry_mass_pg=dry * 1.66054e-12, wall_mass_pg=wall * 1.66054e-12, glucan_glucose=glc_res, mannose=man_res,
                glucan_beads_per_nm3=glc_res / 36 / vol_glucan, glucan_layer_um3=vol_glucan / 1e9, mannan_layer_um3=vol_mannan / 1e9,
                wall_proteins_in_mannan_layer=n_wp, mannan_beads_per_protein=man_res / 3 / max(n_wp, 1))


if __name__ == "__main__":
    main()
