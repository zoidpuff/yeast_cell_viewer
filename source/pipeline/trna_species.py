"""Nuclear tRNA species (IMPROVEMENTS P12): one library type per isoacceptor instead of one tRNA-Phe model (1EHZ)
for all N_TRNA (3.3 million) tRNAs.

Sequences: the mature tRNAs of S. cerevisiae S288C in GtRNAdb (tRNAscan-SE on sacCer3; Chan & Lowe 2016 Nucleic
Acids Res 44:D184; data/trna/sacCer3-mature-tRNAs.fa, downloaded 2026-11-12): introns removed, but without the 3' CCA,
which yeast nuclear tRNA genes do not encode; it is appended here (tRNA-Phe then has the 76 nt of 1EHZ). The 275
genes fall into 42 isoacceptors (anticodon families, the initiator Met apart from the elongator); two genes of
undetermined identity (Und-NNN) are left out. Each isoacceptor is drawn as its most common mature sequence.
Abundance: the share of N_TRNA of each isoacceptor is its share of the genes. tRNA levels in growing yeast follow
gene copy number (Percudani et al. 1997 J Mol Biol 268:322: correlation with the measured content of 21 tRNA species);
per-species measurements are not used.
Structures (predict): the isoacceptor's sequence (with CCA) as one RNA chain in Boltz-2 (Passaro et al. 2025
bioRxiv 10.1101/2025.06.14.659707; an AlphaFold3-like model; no MSA for RNA), best of 5 samples by confidence; the
models (pipeline/data/trna/models/<lib_id>.cif) keep base identities but no modified nucleosides. Class II tRNAs
(Leu, Ser; 85-87 nt with CCA) have the long variable arm that the single tRNA-Phe model lacked. Library entries:
one bead per nucleotide (P atom), as for every RNA in the library.
"""
import os, re, json, shutil, subprocess, collections
import numpy as np
import pandas as pd

from config import DATA, LIB

FASTA = f"{DATA}/trna/sacCer3-mature-tRNAs.fa"
MODELS = f"{DATA}/trna/models"
_TOP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..")      # next to yeast_cell_viewer/
BOLTZ = os.environ.get("YC_BOLTZ", os.path.join(_TOP, "boltz_env", "bin", "boltz"))
BOLTZ_CACHE = os.environ.get("YC_BOLTZ_CACHE", os.path.join(_TOP, "boltz_cache"))      # model weights
EXCLUDE = {"Und"}
SAMPLES = 5


def species():
    """one row per isoacceptor: lib_id, isotype, anticodon, genes (count), seq (most common mature sequence + CCA),
    n_seqs (distinct sequences in the family), share (of all included genes)"""
    recs, name, seq = [], None, []
    for line in open(FASTA):
        if line.startswith(">"):
            if name:
                recs.append((name, "".join(seq)))
            name = line[1:].split()[0]; seq = []
        else:
            seq.append(line.strip().upper())
    recs.append((name, "".join(seq)))
    fam = collections.defaultdict(list)
    for n, s in recs:
        m = re.match(r"Saccharomyces_cerevisiae_tRNA-(\w+)-(\w{3})-\d+-\d+", n)
        if m.group(1) not in EXCLUDE:
            fam[(m.group(1), m.group(2).replace("T", "U"))].append(s)
    rows = []
    for (aa, ac), ss in fam.items():
        c = collections.Counter(ss)
        rows.append(dict(lib_id=f"tRNA_{aa}_{ac}", isotype=aa, anticodon=ac, genes=len(ss), n_seqs=len(c),
                         seq=c.most_common(1)[0][0] + "CCA"))
    d = pd.DataFrame(rows).sort_values("lib_id").reset_index(drop=True)
    d["share"] = d.genes / d.genes.sum()
    return d


def label(r):
    aa = {"iMet": "initiator Met"}.get(r.isotype, r.isotype)
    return f"tRNA-{aa} (anticodon {r.anticodon}; {r.genes} gene{'s' if r.genes > 1 else ''})"


def predict(d=None, gpu="0", log=print):
    """Boltz-2 models of every isoacceptor not yet in MODELS; writes MODELS/<lib_id>.cif and <lib_id>.json"""
    d = species() if d is None else d
    os.makedirs(MODELS, exist_ok=True)
    todo = [r for r in d.itertuples() if not os.path.exists(f"{MODELS}/{r.lib_id}.cif")]
    if not todo:
        return
    work = f"{MODELS}/_boltz_in"
    shutil.rmtree(work, ignore_errors=True); os.makedirs(work)
    for r in todo:
        open(f"{work}/{r.lib_id}.yaml", "w").write(f"version: 1\nsequences:\n  - rna:\n      id: A\n      sequence: {r.seq}\n")
    out = f"{MODELS}/_boltz_out"
    env = dict(os.environ, CUDA_VISIBLE_DEVICES=gpu)
    subprocess.run([BOLTZ, "predict", work, "--out_dir", out, "--diffusion_samples", str(SAMPLES), "--output_format", "mmcif",
                    "--cache", BOLTZ_CACHE, "--override", "--no_kernels"], check=True, env=env)
    for r in todo:
        pdir = [p for p in [f"{out}/boltz_results__boltz_in/predictions/{r.lib_id}"] if os.path.isdir(p)][0]
        shutil.copy(f"{pdir}/{r.lib_id}_model_0.cif", f"{MODELS}/{r.lib_id}.cif")
        conf = json.load(open(f"{pdir}/confidence_{r.lib_id}_model_0.json"))
        json.dump(dict(seq=r.seq, model="Boltz-2", samples=SAMPLES, **{k: conf[k] for k in conf if not isinstance(conf[k], dict)}),
                  open(f"{MODELS}/{r.lib_id}.json", "w"), indent=1)
        log(f"  {r.lib_id}: {len(r.seq)} nt, confidence {conf.get('confidence_score', float('nan')):.2f}, "
            f"pLDDT {conf.get('complex_plddt', float('nan')):.2f}")


def build_library(d=None, log=print):
    """library entries LIB/<lib_id>.npz from the models (one bead per nucleotide, centred)"""
    import gemmi, structlib as sl, build_structs as bs
    d = species() if d is None else d
    for r in d.itertuples():
        st = gemmi.read_structure(f"{MODELS}/{r.lib_id}.cif")
        X, K, _ = sl.beads_from_structure(st)
        assert len(X) == len(r.seq), (r.lib_id, len(X), len(r.seq))
        n, ns, rmax = bs.save(r.lib_id, X - X.mean(0), K, extra=dict(src="pred", pdb="", predictor="Boltz-2"), per=15)   # as tRNA (1EHZ)
        log(f"  {r.lib_id}: {n} beads, {ns} spheres, radius {rmax:.1f} nm")


def expected_rows(n_total):
    """inventory rows (lib_id, ing_id, region, expected) splitting n_total tRNAs over the isoacceptors, or the single
    tRNA-Phe row where the library entries are missing"""
    d = species()
    if not all(os.path.exists(f"{LIB}/{lid}.npz") for lid in d.lib_id):
        return [("tRNA", "tRNA", "cytosol", n_total)]
    return [(r.lib_id, r.lib_id, "cytosol", n_total * r.share) for r in d.itertuples()]


if __name__ == "__main__":
    import sys
    d = species()
    print(d[["lib_id", "genes", "n_seqs", "share"]].assign(nt=d.seq.str.len()).to_string())
    if "predict" in sys.argv:
        predict(d, gpu=os.environ.get("CUDA_VISIBLE_DEVICES", "0"))
    if "library" in sys.argv:
        build_library(d)
