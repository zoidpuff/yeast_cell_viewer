# Pipeline scripts

Paths come from `pipeline/config.py`, which you can override with `YC_*` environment variables. All output goes to `work/` (`YC_WORK`). Times are from the log of build R on 64 cores (2026-11-12).

| Script | Output | Time |
|---|---|---|
| `s0_library.py` (not run by default: `run_all.sh 0`) | the base structure library in `pipeline/lib/` (`YC_LIB`): curated complexes, one AlphaFold model per protein, the curated PDB homo-oligomers, tRNA-Phe | 1.5 min |
| `s1_inventory.py` | `work/results/inventory/{ingredients.csv, oligo_map.json, complexes_*.csv, oligomers.csv}`; new library entries; then `pipeline/features.py` (per-protein feature table) | 4 min (about 1 h on a first run without `work/rcsb_cache.json`: RCSB queries) |
| `s2_geometry.py` (`--sphere`: analytic cell instead of the real shapes) | `work/yc3d/{geometry.json, tile_measure.parquet, tile_labvol.parquet, geom/}` | 1 min |
| `s3_pack.py` (`--random`: no organisation layer) | `work/yc3d/{tiles_raw/, tile_counts.parquet, prepass.npz, expected_totals.parquet, pack_summary.parquet, thinning.json}` | 47 min: 28 min serial builders, 19 min tile packing |
| `s4_validate.py` | `work/results/wholecell3d_validation.json`, `work/yc3d/libpack_meta.npz` | 1 min |
| `s5_encode.py` | `work/yc3d/web/data/` (needs the AFDB tar for per-residue pLDDT) | 3 min |
| `run_all.sh [first_stage]` | all of the above in order, after `build_stamp.py start` (writes `work/build.json`: commit, master seed `YC_SEED`, constants, inputs) | 71 min |
| `validate_all.sh` | every comparison with data, the placement checks, figures and the PDF of the validation document; ends with `build_stamp.py check` (lists result files older than the build) | 16 min (needs the GPUs) |
| `build_isolated.sh NAME [GIT_REF] [FIRST_STAGE]` | a build of one commit from a frozen git worktree in `../runs/NAME/` with its own library and output; `YC_RESUME=1` continues after a fix; `YC_LIB_FROM=empty` starts without a library | as above |
| `organisation_sweep.py` (optional) | — | prints how well personalized-PageRank profile correlations predict co-complex membership and shared localization, for STRING thresholds, damping values and profile variants; used to choose the parameters of `pipeline/organisation.py` | ~20 min |
| `cell_meshes.py` (optional) | — | `pipeline/cell_meshes.npz`: the segmented meshes of the model cell (`YC_ORGANELLE_MESHES`); the file is in the repository, so this only re-extracts it | 5 s |
| `organelle_stats.py` (optional) | — | `pipeline/organelle_stats.json` from segmented cryo-FIB/SEM meshes (`YC_ORGANELLE_MESHES`); the file is in the repository, so this only re-derives it | 5 s |
| `deploy_viewer.sh` | — | builds `web/app` (Node in `.node/`) and installs `index.html`, `assets/` and `data/` into the viewer folder; keeps the previous version as `*_prev` | 10 s |
| `check_inventory.py` | — | copy conservation of the inventory against the proteome | 5 s |

`s5_encode.py --annot-only` writes only `library/genome.json` and `annot/` (what each mRNA, DNA chain, nucleosome, polysome ribosome and Pol II is, for the viewer's hover; `pipeline/polymer_annot.py`) for tiles already encoded, and checks them against the decoded tiles (25 s). A full s5 run writes them too.

To run the whole pipeline: `scripts/run_all.sh` (or `scripts/run_all.sh 3` to restart from packing), then `scripts/deploy_viewer.sh`. It stops at the first failing stage.

**Reproducing a build.** The build is deterministic: the same commit, master seed (`YC_SEED`), inputs and structure library give the same cell, independent of the number of workers. Checked on 2026-11-13: `main` rebuilt by `build_isolated.sh` was identical to build R in all 5,912 tiles (24,484,185 instances), all pre-pass arrays and the expected totals. `work/build.json` records commit, seed, constants and inputs of a build. (`s3_pack.py --from-intermediates` and the `intermediates/` folder date from 2026-11-10 and no longer reproduce the current packing.)

**Inputs outside the bundle** (found by `pipeline/config.py` in `../inputs/` next to the bundle, or set with the environment variable):
- `YC_AFDB_TAR`: AlphaFold DB yeast proteome v6 tar (`inputs/UP000002311_559292_YEAST_v6.tar`).
- `YC_DIMERS`: AFDB predicted yeast dimers (`inputs/afdb_yeast_dimers/` with `homodimers/`, `heterodimers/` and `candidates.csv`).
- `YC_POOL_SCORES`: pair scores of the in-house AlphaFold 3 pool predictions (`inputs/summary_models.parquet`), for the organisation layer of stage 3.
