"""Identity of the drawn polymers and of the chromatin-bound instances, for the viewer's hover and inspect panel.

The tiles hold mRNA, DNA and the other polymers as anonymous polylines (a picked chain is known only by its index in
its tile) and nucleosomes, polysome ribosomes and Pol II as plain instances of their type. This writes what each of
them is, per tile, in the order the tiles hold them (chains: s5_encode.polymer_lists; instances: webenc._stage1 and
_stage2, raw instances then pre-placed ones, sorted stably by type and region):
  mRNA chains    transcript (gene), the part of it (5' tail, ribosome of the polysome, 3' tail, free mRNA) and the
                 nucleotides of the piece counted from the end its part is measured from
  DNA chains     chromosome and bp range of linkers between measured nucleosome dyads and of the DNA through Pol II
                 (between its two flanking nucleosomes); nucleosomal rDNA, active rRNA genes and mtDNA by kind only
  extra chains   nascent Pol II transcripts: chromosome and bp of their polymerase
  instances      nucleosomes (chromosome, dyad bp), rDNA nucleosomes, polysome ribosomes (transcript, rank from the 5'
                 end), engaged Pol II (chromosome, bp)
Transcripts are the placed polysome mRNAs (in mrna.place chain order) followed by the free mRNAs (threads of kind 0),
each of a named gene with that gene's layout (5' UTR, CDS, 3' UTR, poly(A); mrna_genes.gene_table).
Threads are drawn at 1 nt per bead. A 5' tail runs from the 5'-most ribosome to the cap, a 3' tail from the 3'-most
ribosome to the 3' end, and a free mRNA from 5' to 3'; positions are counted from the drawn end where eIF4E or Pab1
sits (mrnp.py), as the model does: a thread that did not reach its target length (20% of 3' tails, 5% of 5' tails and
free mRNAs in build N) ends early, and its drawn end is its end. The mRNA between the ribosomes of a polysome is drawn
shorter than the nucleotides it stands for (linkers are end-to-end paths, polysome.py) and mrna.place gives the rest to
the tails, so pieces there are identified by their ribosome only. A Pol II position is the midpoint of its flanking nucleosome dyads
(pol2.py placed it at the nucleosome nearest its drawn position, evicting one on each side), so near gene ends the ORF
the viewer finds there can be the neighbour of the transcribed gene.
Outputs in <data>:
  library/genome.json     SGD R64 features (ORFs, RNA genes, centromeres, ARS, telomeres, Ty elements), the gene
                          layouts and the transcripts, and the indices of the annotated tiles
  annot/<x>_<y>_<z>.json  per annotated tile (gzip): m (mRNA chains), l (DNA chains), x (extra chains), i (instances)
"""
import os, json, gzip, shutil
import multiprocessing as mp
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import webenc as we, wc3d, genome_map, mrna_genes as mg
from config import DATA, RAW, NPROC

NUC_HALF = 73           # bp of DNA wrapped on each side of a dyad (147 bp core)
FEATURES = {"ORF": "ORF", "tRNA gene": "tRNA gene", "snoRNA gene": "snoRNA gene", "snRNA gene": "snRNA gene",
            "ncRNA gene": "ncRNA gene", "rRNA gene": "rRNA gene", "pseudogene": "pseudogene", "centromere": "centromere",
            "ARS": "replication origin (ARS)", "telomere": "telomere", "transposable element gene": "Ty element gene",
            "long terminal repeat": "LTR"}
L_CODES = ["linker DNA", "rDNA (nucleosomal)", "rDNA (active rRNA gene)", "DNA through Pol II", "mtDNA"]
I_CODES = ["nucleosome", "polysome ribosome", "Pol II", "rDNA nucleosome"]


def features():
    """per chromosome 1-16: [start, end, strand (+1/-1), feature type index, name, systematic name], sorted by start"""
    f = pd.read_csv(f"{DATA}/sgd/SGD_features.tab", sep="\t", header=None, dtype=str)
    f = f[f[1].isin(FEATURES) & f[8].isin([str(i) for i in range(1, 17)])]
    f = f[(f[1] != "ORF") | f[2].isin(["Verified", "Uncharacterized"])]
    types = list(dict.fromkeys(FEATURES.values()))
    out = [[] for _ in range(16)]
    for r in f.itertuples(index=False):
        a, b = int(r[9]), int(r[10])
        name = r[4] if isinstance(r[4], str) else r[3]
        if r[1] == "tRNA gene" and isinstance(r[5], str):
            parts = r[5].split("|"); name = f"{parts[0]} {parts[-1]}"
        out[int(r[8]) - 1].append([min(a, b), max(a, b), 1 if r[11] == "W" else -1, types.index(FEATURES[r[1]]), name, r[3]])
    return types, [sorted(x) for x in out]


def transcripts(pre, log=print):
    """gene layouts and transcripts; returns (genes [name, orf, utr5, cds, utr3, polya], tx gene index, tx ribosomes,
    tx length, per thread its transcript, per polysome chain index its transcript, 5'/3' tail targets per transcript)"""
    chn, kind, tg, tgt = pre["poly_chain"], pre["mrna_kind"], pre["mrna_gene"], pre["mrna_target"]
    chains = np.unique(chn); nch = len(chains)
    assert (kind[:2 * nch:2] == 1).all() and (kind[1:2 * nch:2] == 2).all() and (kind[2 * nch:] == 0).all()
    assert (tg[:2 * nch:2] == pre["poly_chain_gene"]).all() and (tg[1:2 * nch:2] == pre["poly_chain_gene"]).all()
    names = list(pre["poly_chain_gene"]) + list(tg[2 * nch:])
    gt = mg.gene_table().drop_duplicates("gene").set_index("gene")
    gl = sorted(set(names))
    genes = [[g, gt.orf[g], *[int(round(float(gt[c][g]))) for c in ("utr5", "cds", "utr3", "polya")]] for g in gl]
    gi = {g: i for i, g in enumerate(gl)}
    tx_g = np.array([gi[g] for g in names], np.int64)
    tx_len = np.array([int(round(float(gt.utr5[g] + gt.cds[g] + gt.utr3[g] + gt.polya[g]))) for g in names], np.int64)
    tx_rib = np.r_[np.bincount(np.searchsorted(chains, chn), minlength=nch), np.zeros(len(names) - nch, np.int64)]
    thread_tx = np.r_[np.repeat(np.arange(nch), 2), nch + np.arange(len(kind) - 2 * nch)]
    t5 = np.r_[tgt[1:2 * nch:2], np.zeros(len(names) - nch, np.int64)]
    t3 = np.r_[tgt[:2 * nch:2], np.zeros(len(names) - nch, np.int64)]
    bad = int((tgt[2 * nch:] != tx_len[nch:]).sum())
    assert bad == 0, (f"{bad} free mRNAs have a target length different from their gene's: the gene layouts "
                      "(mrna_genes.gene_table) differ from the ones this build was made with; re-run s3 first")
    log(f"transcripts: {nch} polysome mRNAs and {len(names) - nch} free mRNAs of {len(gl)} genes")
    return genes, tx_g, tx_rib, tx_len, thread_tx, chains, t5, t3


def mrna_records(pre, src, tx, log=print):
    """per mRNA piece [transcript, part, first, last, ribosome rank from the 5' end]: part 0 a polysome piece (one per
    ribosome: its channel and the linker to the next ribosome; first = last = 0), 1 a 3' tail piece (nt counted from the
    drawn 3' end, 0-based, first <= last), 2 a 5' tail piece (from the drawn 5' end), 3 a free mRNA piece (from the 5'
    end). Also per ribosome [transcript, rank] and per transcript its drawn nt (free mRNAs: the thread)."""
    genes, tx_g, tx_rib, tx_len, thread_tx, chains, t5, t3 = tx
    chn, kind, mo = pre["poly_chain"], pre["mrna_kind"], pre["mrna_off"]
    rib = np.zeros((len(chn), 2), np.int64)
    ci_of = np.searchsorted(chains, chn)
    for ci in range(len(chains)):
        idx = np.where(ci_of == ci)[0]                                # chain order: the first ribosome is 3'-most
        rib[idx] = np.c_[np.full(len(idx), ci), len(idx) - np.arange(len(idx))]
    n_thread = np.diff(mo)
    drawn = np.zeros(len(tx_len), np.int64)
    np.add.at(drawn, thread_tx, n_thread)
    out = np.zeros((len(src), 5), np.int64)
    for j, (s, i, a, b) in enumerate(src):
        if s == 0:
            out[j] = (rib[i, 0], 0, 0, 0, rib[i, 1])
        elif kind[i] == 0:
            out[j] = (thread_tx[i], 3, a, b - 1, 0)
        else:                                                          # tails grow from their ribosome to their end
            n = n_thread[i]
            out[j] = (thread_tx[i], 1 if kind[i] == 1 else 2, n - b, n - 1 - a, 0)
    return out, rib, drawn


def pol2_loci(pre, segs, log=print):
    """(chromosome, bp) per engaged Pol II: midpoint of the dyads of its flanking nucleosomes (the first and last
    points of its three DNA pieces), else of the one still placed, else of the nucleosome nearest the polymerase"""
    if "pol2_pos" not in pre or not len(pre["pol2_pos"]):
        return np.zeros((0, 2), np.int64)
    N = pre["nuc_pos"].astype(np.float64); tree = cKDTree(N)
    chrom = np.array([segs[s][0] + 1 for s in pre["nuc_chain"]]); bp = pre["nuc_bp"]
    D, O = pre["pol2dna_pts"].astype(np.float64), pre["pol2dna_off"]
    out = np.zeros((len(pre["pol2_pos"]), 2), np.int64); how = [0, 0, 0]
    for k in range(len(out)):
        ends = [D[O[3 * k]], D[O[3 * k + 3] - 1]]
        hit = [j for d, j in zip(*tree.query(ends)) if d < 0.01]
        if hit:
            out[k] = (chrom[hit[0]], int(np.mean(bp[hit]))); how[2 - len(hit)] += 1
        else:
            j = tree.query(pre["pol2_pos"][k])[1]; out[k] = (chrom[j], bp[j]); how[2] += 1
    log(f"Pol II loci: {how[0]} between both flanking nucleosomes, {how[1]} next to one, {how[2]} from the nearest nucleosome")
    return out


def pol2rna_pol(pre):
    """per nascent-transcript thread: its polymerase (threads start at the RNA exit of their polymerase, in order)"""
    if "pol2rna_off" not in pre or len(pre["pol2rna_off"]) < 2:
        return np.zeros(0, np.int64)
    P, O = pre["pol2rna_pts"].astype(np.float64), pre["pol2rna_off"]
    tree = cKDTree(pre["pol2_pos"].astype(np.float64))
    # Threads are stored in polymerase order, but not every polymerase has one. Nearest polymerase to the thread's first
    # bead, taken among those after the previous thread's: with two polymerases close together the plain nearest one
    # can be the neighbour (this stopped the encoder of a build with tighter packing, 2026-11-13).
    d, cand = tree.query(P[O[:-1]], k=min(4, len(pre["pol2_pos"])))
    cand = np.atleast_2d(cand); k = np.empty(len(cand), np.int64); last = -1; fixed = 0
    for i in range(len(cand)):
        ok = [c for c in cand[i] if c > last]
        k[i] = ok[0] if ok else cand[i][0]
        fixed += int(k[i] != cand[i][0]); last = max(last, int(k[i]))
    if fixed:
        print(f"nascent transcripts: {fixed} of {len(k)} assigned to the next polymerase in order instead of the nearest")
    return k


_G = {}


def _locals(ti):
    """local index (in the encoded tile) of each pre-placed instance of tile ti that has a record"""
    G = _G
    g = G["by_tile"].get(ti)
    if g is None:
        return ti, None
    t = G["tiles"][ti]
    r = np.load(f"{RAW}/{G['name'](t)}.npz")
    types = np.r_[r["types"].astype(np.int64), G["pre_type"][G["in_tile"][ti]]]
    region = np.r_[np.repeat(r["job_region"], r["job_placed"]).astype(np.int64), G["pre_reg"][G["in_tile"][ti]]]
    order = np.lexsort((region, types)); inv = np.empty_like(order); inv[order] = np.arange(len(order))
    n_raw = len(r["types"])
    pos_in = np.searchsorted(G["in_tile"][ti], g)
    return ti, np.c_[inv[n_raw + pos_in], g]


def write(out_dir, pre, tiles, pl, pre_tile, pre_type, pre_reg, log=print, check_tiles=300):
    import tiles3d as t3
    segs = genome_map.segments()
    tx = transcripts(pre, log=log)
    genes, tx_g, tx_rib, tx_len = tx[:4]
    m_rec, rib, drawn = mrna_records(pre, pl["mrna_src"], tx, log=log)
    ploc = pol2_loci(pre, segs, log=log)
    rna_pol = pol2rna_pol(pre)
    # DNA chains
    nchrom = np.array([segs[s][0] + 1 for s in pre["nuc_chain"]]); nbp = pre["nuc_bp"]
    l_rec = np.zeros((len(pl["link_src"]), 4), np.int64)
    for j, (s, i, a, b) in enumerate(pl["link_src"]):
        if s == 0:
            lo, hi = sorted((int(nbp[i]), int(nbp[i + 1])))
            lo, hi = (lo + NUC_HALF, hi - NUC_HALF) if hi - lo > 2 * NUC_HALF else ((lo + hi) // 2,) * 2
            l_rec[j] = (0, nchrom[i], lo, hi)
        elif s == 3:
            k = i // 3; l_rec[j] = (3, ploc[k, 0], ploc[k, 1], ploc[k, 1])
        else:
            l_rec[j] = (s, 12 if s in (1, 2) else 0, 0, 0)
    x_rec = np.zeros((len(pl["extra_src"]), 3), np.int64)
    for j, (s, i, a, b) in enumerate(pl["extra_src"]):
        x_rec[j] = (0, *ploc[rna_pol[i]]) if s == 0 else (s, 0, 0)
    # instances with a record: global pre index -> (code, a, b)
    off, o = {}, 0
    for k, _ in wc3d.PRE_GROUPS:
        if f"{k}_pos" in pre:
            off[k] = o; o += len(pre[f"{k}_pos"])
    recs = {}
    if "nuc" in off:
        g = off["nuc"] + np.arange(len(nbp)); recs.update(zip(g, zip([0] * len(g), nchrom, nbp)))
    if "poly" in off:
        g = off["poly"] + np.arange(len(rib)); recs.update(zip(g, zip([1] * len(g), rib[:, 0], rib[:, 1])))
    if "pol2" in off:
        g = off["pol2"] + np.arange(len(ploc)); recs.update(zip(g, zip([2] * len(g), ploc[:, 0], ploc[:, 1])))
    if "rnuc" in off:
        g = off["rnuc"] + np.arange(len(pre["rnuc_pos"])); recs.update(zip(g, zip([3] * len(g), [12] * len(g), [0] * len(g))))
    gk = np.array(sorted(recs), np.int64)
    by_tile, in_tile = {}, {}
    for ti in np.unique(pre_tile[gk]):
        by_tile[int(ti)] = gk[pre_tile[gk] == ti]
        in_tile[int(ti)] = np.where(pre_tile == ti)[0]
    _G.update(tiles=tiles, by_tile=by_tile, in_tile=in_tile, pre_type=np.asarray(pre_type, np.int64),
              pre_reg=np.asarray(pre_reg, np.int64), name=t3.tile_name)
    with mp.get_context("fork").Pool(NPROC) as pool:
        loc = dict(pool.map(_locals, sorted(by_tile), chunksize=8))
    # per tile files
    adir = f"{out_dir}/annot"
    shutil.rmtree(adir, ignore_errors=True); os.makedirs(adir)
    mt, lt, xt = (pl[k] for k in ("mrna_tile", "link_tile", "extra_tile"))
    m_by, l_by, x_by = ({int(t): np.where(a == t)[0] for t in np.unique(a[a >= 0])} for a in (mt, lt, xt))
    done, nbytes = [], 0
    for ti in sorted(set(m_by) | set(l_by) | set(loc)):
        t = tiles[ti]
        r = np.load(f"{RAW}/{t3.tile_name(t)}.npz")
        n_raw_m = int((np.diff(r["mrna_off"]) >= 2).sum())        # packer mRNA threads (none since gene-specific mRNAs)
        xi = x_by.get(ti, np.zeros(0, np.int64))
        d = dict(m=[None] * n_raw_m + m_rec[m_by.get(ti, [])].tolist(), l=l_rec[l_by.get(ti, [])].tolist(),
                 x=x_rec[xi].tolist() if (x_rec[xi, 0] == 0).any() else [], i=[])
        if loc.get(ti) is not None:
            L_ = loc[ti][np.argsort(loc[ti][:, 0])]
            d["i"] = [[int(lc), *map(int, recs[int(g)])] for lc, g in L_]
        if not (d["m"] or d["l"] or d["i"] or d["x"]):
            continue
        if not d["x"]:
            del d["x"]
        b = gzip.compress(json.dumps(d, separators=(",", ":")).encode(), compresslevel=9, mtime=0)
        open(f"{adir}/{t[0]}_{t[1]}_{t[2]}.json", "wb").write(b)
        done.append(ti); nbytes += len(b)
    ftypes, feat = features()
    genome = dict(note="polymer_annot.py: genome features (SGD R64), gene layouts and transcripts for the viewer's hover; "
                       "per-tile identities in annot/<x>_<y>_<z>.json",
                  ftypes=ftypes, feat=feat, genes=genes,
                  tx=dict(g=tx_g.tolist(), r=tx_rib.tolist(), n=tx_len.tolist(), d=drawn.tolist()),
                  l_codes=L_CODES, i_codes=I_CODES, tiles=done,
                  i_types=sorted({int(pre_type[g]) for g in gk}))
    b = gzip.compress(json.dumps(genome, separators=(",", ":")).encode(), compresslevel=9, mtime=0)
    open(f"{out_dir}/library/genome.json", "wb").write(b)
    log(f"polymer identities: {len(done)} tiles ({nbytes / 1e6:.1f} MB gz), {len(m_rec)} mRNA pieces, {len(l_rec)} DNA pieces, "
        f"{len(gk)} instances; genome.json {len(b) / 1e6:.2f} MB gz")
    check(out_dir, tiles, done, pl, loc, pre_type, pre, rng_n=check_tiles, log=log)
    return done


def check(out_dir, tiles, done, pl, loc, pre_type, pre, rng_n=300, log=print):
    """decode encoded tiles and compare: chain counts per section, first points and lengths of the chains, and type and
    position of every annotated instance against its pre-placed source (raises on a mismatch)"""
    _, pp, _ = wc3d.pre_instances(pre)
    rng = np.random.default_rng(0)
    sel = rng.choice(done, min(rng_n, len(done)), replace=False)
    mt, lt = pl["mrna_tile"], pl["link_tile"]
    n_inst = 0
    for ti in sel:
        t = tiles[ti]
        blob = gzip.decompress(open(f"{out_dir}/tiles/{t[0]}_{t[1]}_{t[2]}.bin.gz", "rb").read())
        D = we.decode_tile(blob)
        h = np.frombuffer(blob, np.uint32, 16)
        d = json.loads(gzip.decompress(open(f"{out_dir}/annot/{t[0]}_{t[1]}_{t[2]}.json", "rb").read()))
        assert len(d["m"]) == h[6] and len(d["l"]) == h[8], (t, len(d["m"]), h[6], len(d["l"]), h[8])
        assert "x" not in d or len(d["x"]) == h[13]
        mi = np.where(mt == ti)[0]; li = np.where(lt == ti)[0]
        n_raw = len(d["m"]) - len(mi)
        for k, j in zip(range(n_raw, len(d["m"])), mi):
            assert np.abs(D["mrna"][k][0] - pl["mrna"][j][0]).max() < 0.01 and len(D["mrna"][k]) == len(pl["mrna"][j])
        for k, j in enumerate(li):
            assert np.abs(D["linker"][k][0] - pl["links"][j][0]).max() < 0.02 and len(D["linker"][k]) == len(pl["links"][j])
        types, pos = D["inst"][0], D["inst"][1]
        if loc.get(ti) is not None:
            lc, g = loc[ti][:, 0], loc[ti][:, 1]
            assert (types[lc] == np.asarray(pre_type)[g]).all(), (t, int((types[lc] != np.asarray(pre_type)[g]).sum()))
            dist = np.abs(pos[lc] - pp[g]).max()
            assert dist < 0.01, (t, float(dist))
            n_inst += len(lc)
    log(f"check: {len(sel)} tiles decoded; chain counts and first points match; {n_inst} annotated instances at "
        f"pre-placed positions of their types")
