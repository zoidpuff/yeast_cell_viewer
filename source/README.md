# Source code

This folder holds the code that builds the whole-cell model and the code of the viewer. The scientific description of the model (data, assumptions, design decisions, validation and limitations) is in [MODEL.md](MODEL.md), which is the same file as the one at the repository root.

## Layout

| Path | Contents |
|---|---|
| `pipeline/` | Python modules that build the model. Paths and machine settings are in `config.py`. |
| `scripts/` | The five pipeline stages (`s1_inventory.py` to `s5_encode.py`), `run_all.sh`, and the deploy and export scripts. See `scripts/README.md`. |
| `web/app/` | The viewer: Vite + TypeScript, with all WGSL shaders in `src/shaders.ts`. |
| `notebook/` | The notebook that drove the first version of the pipeline. It is a historical record, but cells 3–4 are still how the structure library is first built (see below). |
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
| Packing (random sequential addition, collision grid) | `pipeline/pack3d.py`, `pipeline/wc3d.py` |
| Pre-placed structures | `polysome.py`, `mrna.py`, `er_ribosomes.py`, `microtubules.py`, `cytoskeleton.py`, `actin_patches.py`, `atp_dimers.py`, `supercomplexes.py`, `fg_mesh.py` in `pipeline/` |
| Validation | `pipeline/val3d.py`, `scripts/s4_validate.py` |
| Viewer data format | `pipeline/webenc.py` (format in its docstring) |

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

PDB assemblies are downloaded from RCSB on demand and cached.

### Structure library

The residue-bead library (`pipeline/lib/`, one `.npz` per species, about 63 MB) is not included. It is built once from the AlphaFold proteome and the curated PDB complexes by `build_structs.process_all_af` and `build_structs.process_complex` (notebook cells 3–4). Stage 1 then adds the complexes and dimers it selects.

### Running the pipeline

```
scripts/run_all.sh        # inventory -> geometry -> packing -> validation -> web encoding
scripts/deploy_viewer.sh  # build the viewer and install it next to the encoded data
```

The first run of stage 1 spends most of its time on RCSB queries; after that a full rebuild takes roughly 10–15 minutes on a recent multi-core machine. Seeds are fixed throughout.
