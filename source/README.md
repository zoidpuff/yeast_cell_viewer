# Source code

This folder holds the code that builds the whole-cell model and the code of the viewer. The scientific description of the model (data, assumptions, design decisions, validation and limitations) is in [MODEL.md](MODEL.md), which is the same file as the one at the repository root.

## Layout

| Path | Contents |
|---|---|
| `pipeline/` | Python modules that build the model. Paths and machine settings are in `config.py`. |
| `scripts/` | The pipeline stages (`s0_library.py` to `s5_encode.py`), `run_all.sh`, `build_isolated.sh` (a build from a frozen checkout into its own folder), `build_stamp.py`, and the deploy and export scripts. See `scripts/README.md`. |
| `web/app/` | The viewer: Vite + TypeScript, with all WGSL shaders in `src/shaders.ts`. |
| `env/requirements.txt` | Python dependencies (Python 3.11). |

Where to look for a given part of the model:

| Topic | Module(s) |
|---|---|
| Abundances, localization rules, compartment volumes | `pipeline/inventory.py` |
| Curated complexes, copy numbers of complexes | `pipeline/ingredients.py`, `pipeline/assemblies.py`, `pipeline/extend_complexes.py` |
| AlphaFold homo- and heterodimers | `pipeline/heterodimers.py`, `scripts/s1_inventory.py` |
| Residue beads, sphere trees, pLDDT trimming | `pipeline/structlib.py`, `pipeline/build_structs.py` |
| Membrane orientation | `pipeline/membrane_orient.py` |
| Cell geometry (signed-distance fields) | `pipeline/cell3d.py`, `pipeline/tiles3d.py` |
| Real cell, nucleus and vacuole shapes (segmented meshes → distance fields) | `pipeline/cellshape.py`, `pipeline/cell_meshes.npz`, `scripts/cell_meshes.py` |
| Packing (random sequential addition, collision grid) | `pipeline/pack3d.py`, `pipeline/wc3d.py` |
| Organisation layer of the packing (STRING network propagation, partner attraction) | `pipeline/organisation.py` |
| Mitoribosomes on the inner membrane | `pipeline/mitoribosomes.py` |
| Eisosome coats (Pil1/Lsp1 lattice in the plasma-membrane furrows) | `pipeline/eisosomes.py`; furrows and endocytic invaginations in `pipeline/cell3d.py` |
| Chromosome arms (Rabl layout), transcribing Pol II, mtDNA nucleoids | `pipeline/wc3d.py` (`rabl_segments`), `pipeline/pack3d.py` (`arm_walks`), `pipeline/pol2.py`, `pipeline/mtdna.py` |
| RAC on ribosomes, Cvt complex, bud scar | `pipeline/rac.py`, `pipeline/cvt.py`, `pipeline/cell3d.py` (`T_TORUS`) |
| Pre-placed structures | `polysome.py`, `mrna.py`, `er_ribosomes.py`, `microtubules.py`, `cytoskeleton.py`, `actin_patches.py`, `atp_dimers.py`, `supercomplexes.py`, `fg_mesh.py` in `pipeline/` |
| Validation | `pipeline/val3d.py`, `scripts/s4_validate.py` |
| Viewer data format | `pipeline/webenc.py` (format in its docstring) |
| Simulation of a box (Brownian dynamics on the GPU, MD export; MODEL.md section 14) | `web/app/src/sim/`: `extract.ts` (box to rigid bodies), `params.ts` (diffusion, hydrodynamics, time-step limits), `engine.ts` (WebGPU kernels), `controller.ts` (run, record, play), `export.ts` and `run_openmm.py` (MD files) |

## Rebuilding the model

### Environment

```
python3.11 -m venv .venv && .venv/bin/pip install -r env/requirements.txt
cd web/app && npm ci              # Node 22
```

### Inputs

These are not in the repository; all are public downloads.

| Input | Where it goes | Source |
|---|---|---|
| PaxDb abundances | `pipeline/data/paxdb_integrated.txt` | PaxDb, *S. cerevisiae* (4932) whole-organism integrated dataset |
| UniProt proteome table | `pipeline/data/uniprot_yeast.tsv` | Proteome UP000002311 (`pipeline/fetch_uniprot.py`) |
| UniProt topology | `pipeline/data/uniprot_topology.tsv` | Same proteome: topology, transmembrane and intramembrane columns |
| Complex Portal | `pipeline/data/complexportal_559292.tsv` | EMBL-EBI Complex Portal, `complextab/559292.tsv` |
| IntAct pairs | `pipeline/data/intact_pairs.tsv` | Built from IntAct's yeast export by `pipeline/intact_pairs.py` |
| Paralog pairs | `pipeline/data/alliance_yeast_paralog_pairs.csv` | Alliance of Genome Resources |
| AlphaFold DB proteome | `YC_AFDB_TAR` | `UP000002311_559292_YEAST_v6.tar` from the AlphaFold DB FTP site |
| AlphaFold DB predicted dimers | `YC_DIMERS` (`homodimers/`, `heterodimers/`) | AlphaFold DB yeast homo- and heterodimer models |
| Segmented meshes of the model cell | `pipeline/cell_meshes.npz` (in the repository) | Extracted by `scripts/cell_meshes.py` from segmented meshes of EMPIAR-11416 (`YC_ORGANELLE_MESHES`); only needed to re-extract |
| Organelle shape statistics | `pipeline/organelle_stats.json` (in the repository) | Derived by `scripts/organelle_stats.py` from segmented meshes of EMPIAR-11416 cells (`YC_ORGANELLE_MESHES`); only needed to re-derive the statistics |
| STRING v12 yeast links | `pipeline/data/string/4932.protein.links.v12.0.txt.gz` | stringdb-downloads.org, `protein.links.v12.0` |
| SGD chromosomal features | `pipeline/data/sgd/SGD_features.tab` | downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab (gene coordinates for transcribing Pol II) |
| Engaged Pol II per gene | `pipeline/data/pelechano2010/files/pone.0015442.s002.xls` | Pelechano et al. 2010 PLoS One 5:e15442, Table S1 (Europe PMC supplementary files of PMC2982843) |
| CYCLoPs LOC-scores (WT1-3) | `pipeline/data/cyclops/WT{1,2,3}_LOCscore.xls` | thecellvision.org/cyclops, Supplemental, Online Tables 3-1 to 3-3 (Koh et al. 2015 G3) |
| Submitochondrial localization | `pipeline/data/vogtle2017_submito.xlsx` | Vögtle et al. 2017 Nat Commun 8:290, Supplementary Data 5 |
| Mitochondrial proteome | `pipeline/data/morgenstern2017/TableS1_mmc2.xlsx` | Morgenstern et al. 2017 Cell Rep 19:2836, Table S1 |
| Pool pair scores (stage 3 organisation layer) | `YC_POOL_SCORES` | `summary_models.parquet` of the in-house AlphaFold 3 pool predictions; not public. Without it, run stage 3 with `--random` |

PDB assemblies are downloaded from RCSB on demand and cached.

### Structure library

The residue-bead library (`pipeline/lib/`, one `.npz` per species, about 63 MB) is not included. Stage 0 (`scripts/s0_library.py`, or `scripts/run_all.sh 0`) builds its base from the AlphaFold proteome, the curated PDB complexes and the curated PDB homo-oligomers; stage 1 then adds the complexes and dimers it selects, and the builders of stage 3 add their models. Until 2026-11-13 the base was built by cells of the notebook that drove the first version of the pipeline (removed since; it is in the git history up to commit `3a509bf`, with the eight modules only it used: the 2D column model and its figures).

### Running the pipeline

```
scripts/run_all.sh        # inventory -> geometry -> packing -> validation -> web encoding
scripts/deploy_viewer.sh  # build the viewer and install it next to the encoded data
```

The first run of stage 1 spends most of its time on RCSB queries; after that a full rebuild takes about 1.5 hours on a 12-core machine, most of it the squeeze relaxation in stage 3. Seeds are fixed throughout.
