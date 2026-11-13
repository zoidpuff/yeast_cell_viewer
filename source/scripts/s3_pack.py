"""Stage 3: expected totals per (library type, region), NPC + microtubule + chromatin + actin cable + polysome pre-pass, distribution over tiles,
checkerboard-phased parallel tile packing. Notebook cells 27-29.

  python s3_pack.py                      # from s1 inventory + s2 geometry, with the organisation layer (organisation.py)
  python s3_pack.py --random             # the same with uniform random placement
  python s3_pack.py --from-intermediates DIR
        # re-pack the shipped model: reuse DIR/{tile_counts.parquet, prepass.npz, mrna_counts.npy, libpack_meta.npz,
        # partners.npz (organisation layer; absent: uniform placement)}
"""
import argparse, shutil
from _common import *
import organisation as og
import cell3d as c3, mrna_genes as mg, mrna_localization as mloc, mrnp, tiles3d as t3, pack3d as p3, wc3d, polysome as ps, cytoskeleton as cs, microtubules as mt, actin_patches as patches, mrna, atp_dimers as atp, fg_mesh as fgm, inventory as inv, rdna, eisosomes as eis, cvt, rac, mtdna, pol2, origins, kinetochores as kin, chromatin_sites as chs, nvj, amphiphysin as amph, mcc, spb, trna_telomeres as ttl, pore_receptors as prc


def consume(exp3, ing_ids, n):
    """remove n instances from the expected totals of these ingredients (all regions, proportionally); returns n taken"""
    r = exp3.ing_id.isin(ing_ids); tot = float(exp3.loc[r, "expected"].sum()); take = min(max(n, 0), tot)
    if tot > 0:
        exp3.loc[r, "expected"] *= 1 - take / tot
    log(f"  {'/'.join(ing_ids)}: {take:.0f} of {tot:.0f} used ({n:.0f} needed)")
    return take


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from-intermediates")
    ap.add_argument("--random", action="store_true", help="uniform random placement (no organisation layer)")
    args = ap.parse_args()
    org = org_prm = None
    cell = c3.Cell.load(f"{C.W3D}/geometry.json")
    tiles = canonical_tiles()
    measure = pd.read_parquet(f"{C.W3D}/tile_measure.parquet")
    assert list(measure.index) == [t3.tile_name(t) for t in tiles]

    if args.from_intermediates:
        d = args.from_intermediates
        lib_ids = [str(x) for x in np.load(f"{d}/libpack_meta.npz")["ids"]]
        lp = wc3d.build_libpack(lib_ids)
        counts = pd.read_parquet(f"{d}/tile_counts.parquet")
        mrna_counts = np.load(f"{d}/mrna_counts.npy")
        pre = dict(np.load(f"{d}/prepass.npz"))
        for k in ("li_npc", "li_nuc"): pre[k] = int(pre[k])
        exp3 = None
        if os.path.exists(f"{d}/partners.npz") and not args.random:      # organisation layer of the shipped build
            z = np.load(f"{d}/partners.npz")
            assert [str(x) for x in z["ids"]] == [str(x) for x in lp["ids"]]
            org = {k: z[k] for k in ("off", "idx", "strength", "kind")}; org_prm = og.params()
    else:
        ing, oligo_map = load_inventory()
        exp3 = wc3d.expected_totals(ing, oligo_map)
        log("expected per region:", exp3.groupby("region").expected.sum().round(0).to_dict())
        dropped = exp3[~exp3.region.isin(t3.REGION_NAMES)].groupby("region").expected.sum()
        log("not packed:", dropped.round(0).to_dict())
        cs.filament_model(); cs.tropomyosin_model(log=log); mt.model(); atp.model(log=log); eis.model(log=log); rac.build(log=log); nvj.model(log=log); amph.model(log=log)
        lib_ids = sorted(set(exp3[exp3.region.isin(t3.REGION_NAMES)].lib_id) | {"NPC", "nucleosome", cs.LIB_ID, cs.TPM_ID, atp.LIB_ID, *mt.LIB_IDS, *patches.LIB_IDS, *eis.LIB_IDS, rac.LIB_ID, nvj.LIB_ID, amph.LIB_ID}
                         # actin that the filaments give back is placed as monomers, also when the heterodimers
                         # (cofilin, profilin, Srv2) have taken the whole monomer pool
                         | ({"ACT1"} if os.path.exists(f"{C.LIB}/ACT1.npz") else set()))
        lp = wc3d.build_libpack(lib_ids)
        n_nuc = int(round(ing.set_index("ing_id").loc["nucleosome", "copies_per_cell"]))
        n_rnuc = rdna.nucleosomes_in_rdna() if "preribo_90S" in set(lp["ids"]) else 0
        n_nuc -= n_rnuc                     # the rest are on the rDNA in the nucleolus (below)
        # microtubules from the SPB (alpha/beta-tubulin heterodimer pool), placed before chromatin: kinetochore-MT
        # plus ends seed the chromosome walks
        tub = (exp3.ing_id == mt.TUBULIN_ING) & (exp3.region == "cytosol")
        mtd = {}

        def microtubules(npc_spheres):
            mtd.update(mt.place(cell, lp, float(exp3.loc[tub, "expected"].sum()), obstacles=npc_spheres, log=log))
            nucl = mtd["mt_reg"] == t3.R_NUCP
            C, R = p3.emit_spheres(mtd["mt_type"][nucl], mtd["mt_pos"][nucl], mtd["mt_q"][nucl], lp["st_off"], lp["STC"],
                                   lp["STR"], lp["BR"], np.full(3, -1e4), np.full(3, 1e4))
            # kinetochores on the plus ends; the chromosome arms start beyond them
            kt = kin.place(cell, lp, mtd["kmt_tips"], mtd["kmt_dirs"], (np.concatenate([C, npc_spheres[0]]), np.concatenate([R, npc_spheres[1]])),
                           float(exp3.loc[pcc, "expected"].sum()), float(exp3.loc[pmi, "expected"].sum()), log=log)
            starts = mtd["kmt_tips"]
            if kt:
                mtd.update({k: v for k, v in kt.items() if k.startswith("kt_")}); mtd.update(n_ccan=kt["n_ccan"], n_mind=kt["n_mind"])
                C, R = np.concatenate([C, kt["spheres"][0]]), np.concatenate([R, kt["spheres"][1]])
                starts = kt["starts"]
            # outer kinetochore (DASH ring, Ndc80 complexes) and the SPB plaques (IMPROVEMENTS N4)
            pool = lambda genes: {g: float(exp3.loc[exp3.lib_id == g, "expected"].sum()) for g in genes}
            ko, ko_n, ko_sph = kin.outer(cell, lp, mtd["kmt_tips"], mtd["kmt_dirs"], (np.concatenate([C, npc_spheres[0]]),
                                         np.concatenate([R, npc_spheres[1]])), pool(kin.DASH + kin.NDC80C), log=log)
            C, R = np.concatenate([C, ko_sph[0]]), np.concatenate([R, ko_sph[1]])
            CA, RA = p3.emit_spheres(mtd["mt_type"], mtd["mt_pos"], mtd["mt_q"], lp["st_off"], lp["STC"], lp["STR"], lp["BR"],
                                     np.full(3, -1e4), np.full(3, 1e4))
            nrm = cell.G[7:10] - cell.G[3:6]
            ps_ = pool(spb.LAYERS); ps_["SPC110"] = max(0.0, ps_["SPC110"] - 14 * mtd["n_rings_nuc"])     # minus the ring share (below)
            sp, sp_n, sp_sph = spb.place(cell, lp, -nrm / np.linalg.norm(nrm), ps_,
                                         (np.concatenate([CA, C, npc_spheres[0]]), np.concatenate([RA, R, npc_spheres[1]])), log=log,
                                         ring_spc110=14 * mtd["n_rings_nuc"])
            C, R = np.concatenate([C, sp_sph[0]]), np.concatenate([R, sp_sph[1]])
            mtd.update(ko); mtd.update(sp); mtd["struct_counts"] = {**ko_n, **sp_n}
            return dict(obstacles=(C, R), starts=starts)
        pcc = (exp3.lib_id == kin.CCAN_ID) & (exp3.region == "nucleoplasm"); pmi = (exp3.lib_id == kin.MIND_ID) & (exp3.region == "nucleoplasm")
        pre = wc3d.prepass(cell, lp, n_nuc, log=log, before_chromatin=microtubules)
        pre.update({k: v for k, v in mtd.items() if k.split("_")[0] in ("mt", "kt", "kto", "spb")})
        for g_, n_ in mtd.get("struct_counts", {}).items():
            consume(exp3, list(exp3.loc[exp3.lib_id == g_, "ing_id"].unique()), n_)
        for m_, k_ in ((pcc, "n_ccan"), (pmi, "n_mind")):
            if k_ in mtd:
                exp3.loc[m_, "expected"] = (exp3.loc[m_, "expected"] - mtd[k_]).clip(lower=0)
        # transport receptors in the pore channels, then the FG-repeat domains threaded around them
        gt0, gp0, gq0 = wc3d.pre_instances(pre)
        KC, KR = p3.emit_spheres(gt0, gp0, gq0, lp["st_off"], lp["STC"], lp["STR"], lp["BR"], np.full(3, -1e4), np.full(3, 1e4))
        rk, k_counts, _ = prc.place(cell, lp, exp3.groupby("lib_id").expected.sum().to_dict(), (KC, KR), log=log)
        pre.update(rk)
        for g_, n_ in k_counts.items():
            consume(exp3, list(exp3.loc[exp3.lib_id == g_, "ing_id"].unique()), n_)
        # FG-repeat domains of the nucleoporins, threaded into the pores (avoiding the NPC scaffold, chromatin, MTs)
        gt0, gp0, gq0 = wc3d.pre_instances(pre)
        FC, FR = p3.emit_spheres(gt0, gp0, gq0, lp["st_off"], lp["STC"], lp["STR"], lp["BR"], np.full(3, -1e4), np.full(3, 1e4))
        pre.update(fgm.place(cell, pre, (FC, FR), log=log))
        FGC = pre["fg_pts"].astype(np.float64); FGR = np.full(len(FGC), fgm.BEAD_R)
        if n_rnuc:
            # rDNA in the nucleolus: active genes with Pol I and nascent pre-rRNA (90S on uncleaved transcripts),
            # then nucleosomal rDNA around them; engaged Pol I and tethered 90S come out of the nucleolar pools
            cn, rn = cell.G[3:6], float(cell.G[11:14].max())
            RC1, RR1 = p3.emit_spheres(gt0, gp0, gq0, lp["st_off"], lp["STC"], lp["STR"], lp["BR"], cn - rn - 100, cn + rn + 100)
            pol = (exp3.lib_id == "PolI") & (exp3.region == "nucleolus")
            r90 = (exp3.lib_id == "preribo_90S") & (exp3.region == "nucleolus")
            rd = rdna.place(cell, lp, float(exp3.loc[pol, "expected"].sum()), float(exp3.loc[r90, "expected"].sum()),
                            (np.concatenate([RC1, FGC]), np.concatenate([RR1, FGR])), log=log)
            pre.update(rd)
            exp3.loc[pol, "expected"] = (exp3.loc[pol, "expected"] - len(rd["rpol_pos"])).clip(lower=0)
            exp3.loc[r90, "expected"] = (exp3.loc[r90, "expected"] - len(rd["r90_pos"])).clip(lower=0)
            gt1, gp1, gq1 = wc3d.pre_instances(pre)
            RC2, RR2 = p3.emit_spheres(gt1, gp1, gq1, lp["st_off"], lp["STC"], lp["STR"], lp["BR"], cn - rn - 100, cn + rn + 100)
            pre.update(wc3d.rdna_chromatin(cell, lp, n_rnuc, (np.concatenate([RC2, FGC, rd["rdna_pts"], rd["prerrna_pts"]]),
                                                          np.concatenate([RR2, FGR, np.full(len(rd["rdna_pts"]), 1.0),
                                                                          np.full(len(rd["prerrna_pts"]), rdna.RNA_R)])), log=log))
        # transcribing Pol II on its genes (nucleosomes evicted around each), with DNA through the cleft and nascent RNA
        gtx, gpx, gqx = wc3d.pre_instances({k: v for k, v in pre.items() if not k.startswith("nuc_")})
        XC, XR = p3.emit_spheres(gtx, gpx, gqx, lp["st_off"], lp["STC"], lp["STR"], lp["BR"], np.full(3, -1e4), np.full(3, 1e4))
        pr2 = (exp3.ing_id == "PolII") & (exp3.region == "nucleoplasm")
        r2 = pol2.place(cell, lp, pre, float(exp3.loc[pr2, "expected"].sum()), (np.concatenate([XC, FGC]), np.concatenate([XR, FGR])), log=log)
        k2 = r2.pop("keep"); nb2 = pre["nuc_brk"].copy(); nb2[np.r_[False, ~k2[:-1]] & k2] = 1
        for key in ("nuc_pos", "nuc_q", "nuc_chain", "nuc_bp"):
            if key not in pre:
                continue
            pre[key] = pre[key][k2]
        pre["nuc_brk"] = nb2[k2]; pre.update(r2)
        exp3.loc[pr2, "expected"] = (exp3.loc[pr2, "expected"] - len(r2["pol2_pos"])).clip(lower=0)
        # licensed origins: an MCM double hexamer and ORC on each confirmed origin (needs the measured nucleosome map)
        pmc = (exp3.lib_id == origins.MCM_ID) & (exp3.region == "nucleoplasm")
        por = (exp3.lib_id == origins.ORC_ID) & (exp3.region == "nucleoplasm")
        gtx, gpx, gqx = wc3d.pre_instances({k: v for k, v in pre.items() if not k.startswith("nuc_")})
        XC, XR = p3.emit_spheres(gtx, gpx, gqx, lp["st_off"], lp["STC"], lp["STR"], lp["BR"], np.full(3, -1e4), np.full(3, 1e4))
        ro = origins.place(cell, lp, pre, float(exp3.loc[pmc, "expected"].sum()), float(exp3.loc[por, "expected"].sum()),
                           (np.concatenate([XC, FGC, pre["pol2dna_pts"], pre["pol2rna_pts"]]),
                            np.concatenate([XR, FGR, np.full(len(pre["pol2dna_pts"]) + len(pre["pol2rna_pts"]), 1.0)])), log=log)
        k3 = ro.pop("keep"); nb3 = pre["nuc_brk"].copy(); nb3[np.r_[False, ~k3[:-1]] & k3] = 1
        for key in ("nuc_pos", "nuc_q", "nuc_chain", "nuc_bp"):
            if key in pre:
                pre[key] = pre[key][k3]
        pre["nuc_brk"] = nb3[k3]; pre.update(ro)
        for m_, k_ in ((pmc, "mcm_pos"), (por, "orc_pos")):
            if k_ in ro:
                exp3.loc[m_, "expected"] = (exp3.loc[m_, "expected"] - len(ro[k_])).clip(lower=0)
        # chromatin proteins on their ChIP-exo sites (needs the measured nucleosome map)
        nucp = exp3.region == "nucleoplasm"
        gtx, gpx, gqx = wc3d.pre_instances(pre)
        XC, XR = p3.emit_spheres(gtx, gpx, gqx, lp["st_off"], lp["STC"], lp["STR"], lp["BR"], np.full(3, -1e4), np.full(3, 1e4))
        beads = np.concatenate([pre[k] for k in ("pol2dna_pts", "pol2rna_pts", "rdna_pts", "prerrna_pts") if k in pre])
        rt, tf_counts = chs.place(cell, lp, pre, exp3[nucp].groupby("lib_id").expected.sum().to_dict(),
                                  (np.concatenate([XC, FGC, beads]), np.concatenate([XR, FGR, np.full(len(beads), 1.0)])), log=log)
        pre.update(rt)
        for g_, n_ in tf_counts.items():
            m_ = nucp & (exp3.lib_id == g_)
            exp3.loc[m_, "expected"] = (exp3.loc[m_, "expected"] - n_ * exp3.loc[m_, "expected"] / max(exp3.loc[m_, "expected"].sum(), 1e-9)).clip(lower=0)
        # tRNA genes (Pol III, TFIIIC), telomere ends (Rap1), silent chromatin (Sir2/3/4): sites the ChIP-exo peaks leave out
        if "tf_pos" in pre and "nuc_bp" in pre:
            nuc_any = exp3.region.isin(["nucleoplasm", "nucleolus"])
            TC, TR = p3.emit_spheres(pre["tf_type"], pre["tf_pos"], pre["tf_q"], lp["st_off"], lp["STC"], lp["STR"], lp["BR"],
                                     np.full(3, -1e4), np.full(3, 1e4))
            rl, lm_counts = chs.place(cell, lp, pre, exp3[nuc_any].groupby("lib_id").expected.sum().to_dict(),
                                      (np.concatenate([XC, TC, FGC, beads]), np.concatenate([XR, TR, FGR, np.full(len(beads), 1.0)])),
                                      sites=ttl.sites(pre), prefix="lm", what="tRNA genes, telomeres, silent loci", seed=59, log=log)
            pre.update(rl)
            for g_, n_ in lm_counts.items():
                m_ = nuc_any & (exp3.lib_id == g_)
                exp3.loc[m_, "expected"] = (exp3.loc[m_, "expected"] - n_ * exp3.loc[m_, "expected"] / max(exp3.loc[m_, "expected"].sum(), 1e-9)).clip(lower=0)
        exp3.loc[tub, "expected"] = (exp3.loc[tub, "expected"] - mtd["n_tub"]).clip(lower=0)
        # ring subunits: gamma-TuSC complexes first, then the free subunit pools (SPC110 from its CMD1 heterodimer)
        n_all = mtd["n_rings_nuc"] + mtd["n_rings_cyto"]
        m = consume(exp3, [mt.TUSC_ING], 6 * n_all)
        for g, k in mt.RING_GENES.items():
            consume(exp3, [g], k * n_all - (2 if g == "TUB4" else 1) * m)
        consume(exp3, ["SPC110", "CMD1_SPC110_het"], 14 * mtd["n_rings_nuc"])
        cyto_mt = mtd["mt_reg"] == t3.R_CYTO
        MC, MR = p3.emit_spheres(mtd["mt_type"][cyto_mt], mtd["mt_pos"][cyto_mt], mtd["mt_q"][cyto_mt], lp["st_off"], lp["STC"],
                                 lp["STR"], lp["BR"], np.full(3, -1e4), np.full(3, 1e4))
        # eisosome coats in the PM furrows, from the Pil1/Lsp1 pools (monomers, homodimers, heterodimer); the coat takes
        # the pools' Pil1 : Lsp1 ratio, so every form is reduced by the same fraction
        genes = ing.set_index("ing_id").genes.fillna("").str.split(";")
        units = {}
        for g in ("PIL1", "LSP1"):
            r = exp3.ing_id.map(lambda i: g in genes.get(i, [])).astype(bool)
            mult = exp3.ing_id.map(lambda i: float(oligo_map[i][1]) if i in oligo_map else 1.0)   # homo-oligomers: k per instance
            units[g] = float((exp3.loc[r, "expected"] * mult[r]).sum())
        eisd, used = eis.place(cell, lp, units["PIL1"] / 2, units["LSP1"] / 2, obstacles=(MC, MR), log=log)
        pre.update(eisd)
        f_eis = 2.0 * (used["Pil1_dimer"] + used["Lsp1_dimer"]) / max(units["PIL1"] + units["LSP1"], 1.0)
        r = exp3.ing_id.map(lambda i: bool({"PIL1", "LSP1"} & set(genes.get(i, []))))
        exp3.loc[r, "expected"] *= 1.0 - f_eis
        log(f"  Pil1/Lsp1: {units['PIL1']:.0f} / {units['LSP1']:.0f} monomers, {f_eis:.1%} in eisosome coats")
        EC, ER_ = p3.emit_spheres(eisd["eis_type"], eisd["eis_pos"], eisd["eis_q"], lp["st_off"], lp["STC"], lp["STR"], lp["BR"],
                                  np.full(3, -1e4), np.full(3, 1e4))
        MC, MR = np.concatenate([MC, EC]), np.concatenate([MR, ER_])          # obstacles for everything placed after
        # MCC marker proteins (Sur7, Nce102) on the furrow walls
        mr_ = exp3.ing_id.isin(mcc.GENES) & (exp3.region == "mem_PM")
        want = {lib: mcc.MCC_FRACTION * n for lib, n in exp3.loc[mr_].groupby("lib_id").expected.sum().items()}
        md_, used_m = mcc.place(cell, lp, want, obstacles=(MC, MR), log=log)
        if len(md_["mcc_pos"]):
            pre.update(md_)
            for lib, k in used_m.items():
                r_ = mr_ & (exp3.lib_id == lib); exp3.loc[r_, "expected"] = (exp3.loc[r_, "expected"] - k).clip(lower=0)
            QC, QR = p3.emit_spheres(md_["mcc_type"], md_["mcc_pos"], md_["mcc_q"], lp["st_off"], lp["STC"], lp["STR"], lp["BR"],
                                     np.full(3, -1e4), np.full(3, 1e4))
            MC, MR = np.concatenate([MC, QC]), np.concatenate([MR, QR])
        # F-actin reserved in s1: cables, then Arp2/3-branched patches; whatever does not fit returns to the monomer pool
        fa = json.load(open(F_ACTIN_JSON))
        cab = cs.place(cell, lp, int(fa["cables"]), obstacles=(MC, MR), log=log)
        pre.update(cab)
        # tropomyosin (TPM1 coiled coils) along both long-pitch strands of the cable filaments
        tpm = cs.decorate(cell, lp, cab, obstacles=(MC, MR), log=log)
        pre.update(tpm)
        consume(exp3, [cs.TPM_ING], 2 * len(tpm["tpm_pos"]))
        CC, CR = p3.emit_spheres(np.r_[np.full(len(cab["act_pos"]), cab["li_actin"]), np.full(len(tpm["tpm_pos"]), tpm["li_tpm"])].astype(np.int32),
                                 np.r_[cab["act_pos"], tpm["tpm_pos"]].astype(np.float32), np.r_[cab["act_q"], tpm["tpm_q"]].astype(np.float32),
                                 lp["st_off"], lp["STC"], lp["STR"], lp["BR"], np.full(3, -1e4), np.full(3, 1e4))
        # endocytic coat at the invagination tips (Sla2, Sla1, End3), then amphiphysin (Rvs161-Rvs167) on the deeper tubes,
        # before the actin grows around them
        pools = {g: float(exp3.loc[exp3.ing_id == g, "expected"].sum()) for g in amph.COAT}
        co, used_c = amph.coat(cell, lp, pools, obstacles=(np.concatenate([MC, CC]), np.concatenate([MR, CR])), log=log)
        if len(co["coat_pos"]):
            pre.update(co)
            for g, k in used_c.items():
                consume(exp3, [g], k)
            KC2, KR2 = p3.emit_spheres(co["coat_type"], co["coat_pos"], co["coat_q"], lp["st_off"], lp["STC"], lp["STR"], lp["BR"],
                                       np.full(3, -1e4), np.full(3, 1e4))
            CC, CR = np.concatenate([CC, KC2]), np.concatenate([CR, KR2])
        rv_pool = min(float(exp3.loc[exp3.ing_id == "RVS161", "expected"].sum()), float(exp3.loc[exp3.ing_id == "RVS167", "expected"].sum()))
        rv, n_rv = amph.place(cell, lp, rv_pool, obstacles=(np.concatenate([MC, CC]), np.concatenate([MR, CR])), log=log)
        if n_rv:
            pre.update(rv); consume(exp3, ["RVS161"], n_rv); consume(exp3, ["RVS167"], n_rv)
            VC, VR = p3.emit_spheres(np.full(n_rv, rv["li_rvs"], np.int32), rv["rvs_pos"], rv["rvs_q"], lp["st_off"], lp["STC"], lp["STR"],
                                     lp["BR"], np.full(3, -1e4), np.full(3, 1e4))
            CC, CR = np.concatenate([CC, VC]), np.concatenate([CR, VR])
        patches.model(dict(zip(*[inv.build()[0][c] for c in ("gene", "uniprot")])))
        avail = lambda i: float(exp3.loc[(exp3.ing_id == i) & (exp3.region == "cytosol"), "expected"].sum())
        pat, used = patches.place(cell, lp, fa["patches"], avail(patches.ARP23), avail(patches.CP_ING),
                             obstacles=(np.concatenate([MC, CC]), np.concatenate([MR, CR])), log=log,
                             sites=[(p, u) for p, u, _ in cell.endocytic_sites()] or None)
        pre.update(pat)
        consume(exp3, [patches.ARP23], used[patches.ARP23]); consume(exp3, [patches.CP_ING], used[patches.CP])
        back = fa["cables"] + fa["patches"] - len(cab["act_pos"]) - used[patches.ACTIN]
        act = (exp3.ing_id == "ACT1") & (exp3.region == "cytosol")
        if act.any():
            exp3.loc[act, "expected"] += back
        else:
            exp3.loc[len(exp3)] = ["ACT1", "ACT1", "cytosol", back]
        log(f"  F-actin: {len(cab['act_pos'])} in cables, {used[patches.ACTIN]} in patches, {back:.0f} returned to the monomer pool")
        PC_, PR_ = p3.emit_spheres(pat["pat_type"], pat["pat_pos"], pat["pat_q"], lp["st_off"], lp["STC"], lp["STR"], lp["BR"],
                                   np.full(3, -1e4), np.full(3, 1e4))
        OC, OR = np.concatenate([CC, PC_]), np.concatenate([CR, PR_])
        # ATP synthase dimer rows on cristae rims (from the monomer pool)
        im = (exp3.ing_id == atp.MONOMER) & (exp3.region == "mem_MITO_IM")
        dim = atp.place(cell, lp, int(exp3.loc[im, "expected"].sum() // 2), log=log)
        pre.update(dim)
        consume(exp3, [atp.MONOMER], 2 * len(dim["atp_pos"]))
        # mtDNA nucleoids in the matrix (before the matrix proteins are packed)
        AC, AR = p3.emit_spheres(np.full(len(dim["atp_pos"]), dim["li_atp"], np.int32), dim["atp_pos"], dim["atp_q"], lp["st_off"],
                                 lp["STC"], lp["STR"], lp["BR"], np.full(3, -1e4), np.full(3, 1e4))
        pre.update(mtdna.place(cell, obstacles=(AC, AR), log=log))
        # Cvt complex next to the vacuole: precursor Ape1 dodecamers and Ams1 (from the vacuolar pools), Atg19 coat
        tot = lambda i: float(exp3.loc[exp3.ing_id == i, "expected"].sum())
        cv, used = cvt.place(cell, lp, tot("APE1"), tot("AMS1"), int(tot("ATG19")),
                             obstacles=(np.concatenate([OC, MC, FGC]), np.concatenate([OR, MR, FGR])), log=log)
        if cv:
            pre.update(cv)
            consume(exp3, ["APE1"], used[cvt.APE1]); consume(exp3, ["AMS1"], used[cvt.AMS1]); consume(exp3, ["ATG19"], used[cvt.ATG19])
            KC, KR = p3.emit_spheres(cv["cvt_type"], cv["cvt_pos"], cv["cvt_q"], lp["st_off"], lp["STC"], lp["STR"], lp["BR"],
                                     np.full(3, -1e4), np.full(3, 1e4))
            OC, OR = np.concatenate([OC, KC]), np.concatenate([OR, KR])
        # nucleus-vacuole junction: Vac8-Nvj1 on the vacuolar membrane facing the nucleus (all Nvj1, one Vac8 each)
        nv, n_nvj = nvj.place(cell, lp, tot("NVJ1"), obstacles=(np.concatenate([OC, MC, FGC]), np.concatenate([OR, MR, FGR])), log=log)
        if n_nvj:
            pre.update(nv); consume(exp3, ["NVJ1"], n_nvj); consume(exp3, ["VAC8"], n_nvj)
            JC, JR = p3.emit_spheres(np.full(n_nvj, nv["li_nvj"], np.int32), nv["nvj_pos"], nv["nvj_q"], lp["st_off"], lp["STC"], lp["STR"],
                                     lp["BR"], np.full(3, -1e4), np.full(3, 1e4))
            OC, OR = np.concatenate([OC, JC]), np.concatenate([OR, JR])
        # polysomes: a share of the cytosolic ribosomes, threaded on mRNAs, placed before packing
        rib = (exp3.lib_id == ps.RIBO) & (exp3.region == "cytosol")
        tx = None
        if mg.GENE_MRNA:      # gene-specific mRNAs: each translated mRNA is one polysome chain of its own size
            tx = mg.transcripts(np.random.default_rng(C.seed(23)), ps.POLY_FRACTION * exp3.loc[rib, "expected"].sum(), log=log)
            tx_chains, tx_free = tx[tx.ribosomes > 0].reset_index(drop=True), tx[tx.ribosomes == 0].reset_index(drop=True)
        poly = ps.place(cell, lp, ps.POLY_FRACTION * exp3.loc[rib, "expected"].sum(),
                        obstacles=(np.concatenate([OC, MC, FGC]), np.concatenate([OR, MR, FGR])), log=log, n_mrna=int(wc3d.N_MRNA),
                        sizes=tx_chains.ribosomes.values if tx is not None else None,
                        linker_nt=(tx_chains.cds / tx_chains.ribosomes - ps.FOOTPRINT_NT).values if tx is not None else None,
                        tether=mloc.tether(cell, tx_chains, np.random.default_rng(C.seed(29)), log=log) if tx is not None else None)
        pre.update(poly)
        # RAC rides on ribosomes (all RAC complexes, spread over polysome and free ribosomes in proportion)
        n_rib = float(exp3.loc[rib, "expected"].sum())
        racr = exp3.ing_id == rac.RAC_ING; n_rac = float(exp3.loc[racr, "expected"].sum()); f_rac = min(1.0, n_rac / max(n_rib, 1.0))
        pre["poly_type"] = rac.decorate(cell, lp, poly, f_rac * len(poly["poly_pos"]),
                                        obstacles=(np.concatenate([OC, MC, FGC]), np.concatenate([OR, MR, FGR])), log=log)
        on_poly = int((pre["poly_type"] != poly["li_ribo"]).sum())
        exp3.loc[rib, "expected"] = exp3.loc[rib, "expected"] - len(poly["poly_pos"])
        n_free = float(exp3.loc[rib, "expected"].sum())
        free_rac = min(n_rac - on_poly, f_rac * n_free)
        new = exp3.loc[rib].assign(lib_id=rac.LIB_ID, ing_id=rac.LIB_ID, expected=lambda d: d.expected * free_rac / max(n_free, 1.0))
        exp3.loc[rib, "expected"] *= 1.0 - free_rac / max(n_free, 1.0)
        consume(exp3, [rac.RAC_ING], on_poly + free_rac)
        exp3 = pd.concat([exp3, new], ignore_index=True)        # free ribosomes with RAC, in the same regions
        log(f"  RAC: {n_rac:.0f} complexes for {n_rib:.0f} cytosolic ribosomes; {on_poly} on polysome ribosomes, {free_rac:.0f} on free ones")
        # full-length mRNAs, threaded before packing: polysome tails + free transcripts (one mRNA per polysome)
        from scipy.spatial.transform import Rotation
        import structlib as sl
        zr = np.load(f"{C.LIB}/{ps.RIBO}.npz")["beads"].astype(np.float64)
        rc, rr = sl.sphere_tree(zr, per=100)
        Rm = Rotation.from_quat(np.c_[poly["poly_q"][:, 1:], poly["poly_q"][:, :1]].astype(np.float64)).as_matrix()
        RC = (np.einsum("nij,sj->nsi", Rm, rc) + poly["poly_pos"][:, None, :]).reshape(-1, 3)
        RAC_C, RAC_R = rac.rac_spheres(poly["poly_pos"], poly["poly_q"], pre["poly_type"], lp)     # RAC at the tunnel exits
        RC = np.concatenate([RC, RAC_C]); rr_all = np.concatenate([np.tile(rr, len(Rm)), RAC_R])
        npc_t = np.full(len(pre["npc_pos"]), pre["li_npc"], np.int32)
        NC, NR = p3.emit_spheres(npc_t, pre["npc_pos"], pre["npc_q"], lp["st_off"], lp["STC"], lp["STR"], lp["BR"], np.full(3, -1e4), np.full(3, 1e4))
        mr = mrna.place(cell, poly, ps.channel_path().astype(np.float64), int(wc3d.N_MRNA) - poly["poly_chains"],
                        (np.concatenate([OC, MC, NC, RC, poly["poly_link"], FGC]),
                         np.concatenate([OR, MR, NR, rr_all, np.full(len(poly["poly_link"]), 0.5), FGR])), log=log,
                        tx_chains=tx_chains if tx is not None else None, tx_free=tx_free if tx is not None else None)
        pre.update(mr)
        # proteins on the mRNA ends: Pab1 on the poly(A) tails, eIF4E on the caps (out of the free pools)
        pab = (exp3.lib_id == "PAB1") & (exp3.region == "cytosol"); cap = (exp3.lib_id == "CDC33") & (exp3.region == "cytosol")
        pa = mg.gene_table().drop_duplicates("gene").set_index("gene").polya.to_dict() if tx is not None else None
        rnp = mrnp.place(cell, lp, pre, float(exp3.loc[pab, "expected"].sum()), float(exp3.loc[cap, "expected"].sum()),
                         (np.concatenate([OC, MC, NC, RC, poly["poly_link"], FGC]),
                          np.concatenate([OR, MR, NR, rr_all, np.full(len(poly["poly_link"]), 0.5), FGR])),
                         polya=pa, log=log)
        if rnp:
            pre.update({k: v for k, v in rnp.items() if k.startswith("rnp_")})
            exp3.loc[pab, "expected"] = (exp3.loc[pab, "expected"] - rnp["n_pab"]).clip(lower=0)
            exp3.loc[cap, "expected"] = (exp3.loc[cap, "expected"] - rnp["n_cap"]).clip(lower=0)
        if wc3d.SHELL_CAP:
            exp3 = wc3d.shell_capacity(exp3, measure, lp, log=log)
        logmod = None
        if not args.random:
            centres = (np.array(tiles, float) + 0.5) * t3.TILE
            logmod, org = og.build(exp3, ing, [str(x) for x in lp["ids"]], centres, log=log)
            org_prm = og.params()
        vol = {str(k): float((4 / 3 * np.pi * lp["STR"][lp["st_off"][i]:lp["st_off"][i + 1]] ** 3).sum()) for i, k in enumerate(lp["ids"])}
        counts = wc3d.distribute(exp3, measure, {k: i for i, k in enumerate(lp["ids"])}, logmod=logmod, vol=vol)
        mrna_counts = np.zeros(len(tiles), np.int64)          # no in-tile threads any more

    halo = float(max(lp["BR"][i] for i, k in enumerate(lp["ids"]) if k != "NPC")) + 1.0
    log(f"{len(lp['ids'])} library types, {len(lp['STR'])} sphere-tree spheres, halo {halo:.1f} nm; "
        f"{counts.n.sum() / 1e6:.2f}M instances to place")
    counts.to_parquet(f"{C.W3D}/tile_counts.parquet")
    np.save(f"{C.W3D}/mrna_counts.npy", mrna_counts)
    np.savez(f"{C.W3D}/prepass.npz", **pre)
    np.save(f"{C.W3D}/lib_ids.npy", np.array([str(i) for i in lp["ids"]]))
    if exp3 is not None:
        exp3.to_parquet(f"{C.W3D}/expected_totals.parquet")

    wc3d.setup_shared(cell, lp, counts, tiles, pre, halo, partners=org, prm=org_prm)
    if wc3d.THIN and exp3 is not None:
        thin_f = wc3d.pilot_capacity(counts, measure, lp, log=log)
        if thin_f:
            counts = wc3d.distribute(wc3d.thin(exp3, thin_f), measure, {k: i for i, k in enumerate(lp["ids"])}, logmod=logmod, vol=vol)
            counts.to_parquet(f"{C.W3D}/tile_counts.parquet")
            json.dump(thin_f, open(f"{C.W3D}/thinning.json", "w"), indent=1)
            log(f"thinned request: {counts.n.sum() / 1e6:.2f}M instances to place")
            wc3d.setup_shared(cell, lp, counts, tiles, pre, halo, partners=org, prm=org_prm)
    wc3d._S["round_edges"] = wc3d.round_edges(counts, lp, wc3d.N_ROUNDS)
    log(f"{wc3d.N_ROUNDS} packing rounds by size class; bounding-radius thresholds (nm) {np.round(wc3d._S['round_edges'], 1).tolist()}")
    shutil.rmtree(wc3d.RAW, ignore_errors=True); os.makedirs(wc3d.RAW)
    phase = np.array([wc3d.phase_of(t) for t in tiles])
    allres = []; T0 = time.time()
    for rnd in range(wc3d.N_ROUNDS):
        last = rnd == wc3d.N_ROUNDS - 1; t0 = time.time(); r = []
        for ph in (range(8) if rnd % 2 == 0 else range(7, -1, -1)):      # reversed order in every second round
            r += wc3d.run_phase(list(np.where(phase == ph)[0]), mrna_counts, rnd=rnd, last=last)
        allres += r
        log(f"round {rnd}: {time.time() - t0:.0f}s, placed {sum(x[2] for x in r)}/{sum(x[1] for x in r)}"
            f"{' (with the squeeze)' if last else ''}, overflow {sum(x[3] for x in r)}")
    res = pd.DataFrame(allres, columns=["ti", "requested", "placed", "overflow", "sec", "n_mrna", "attract_tried", "attract_placed"])
    res = res.groupby("ti", as_index=False).agg(requested=("requested", "sum"), placed=("placed", "sum"), overflow=("overflow", "max"),
                                                sec=("sec", "sum"), n_mrna=("n_mrna", "max"), attract_tried=("attract_tried", "sum"),
                                                attract_placed=("attract_placed", "sum"))
    res.to_parquet(f"{C.W3D}/pack_summary.parquet")
    log(f"total {time.time() - T0:.0f}s; placed {res.placed.sum() / 1e6:.2f}M / {res.requested.sum() / 1e6:.2f}M "
        f"= {res.placed.sum() / res.requested.sum():.4f}; next to a partner {res.attract_placed.sum() / 1e6:.2f}M "
        f"of {res.attract_tried.sum() / 1e6:.2f}M that tried")


if __name__ == "__main__":
    main()
