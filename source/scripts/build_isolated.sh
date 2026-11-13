#!/bin/bash
# Build the model from a frozen checkout of one commit into its own run folder, so the working tree can be edited
# while a build runs and two builds never share output or the structure library.
#
# Usage: scripts/build_isolated.sh NAME [GIT_REF] [FIRST_STAGE]
#   NAME         run name; the run lives in $YC_RUNS/NAME (default: <bundle>/../runs/NAME)
#   GIT_REF      commit or branch to build (default: HEAD)
#   FIRST_STAGE  passed to run_all.sh (default 1)
# Environment: YC_SEED (master seed, default 0), YC_NPROC, YC_VALIDATE=0 to skip validate_all.sh,
#              YC_LIB_FROM=empty to start from an empty structure library (needs stage 0)
#
# Layout of a run:
#   tree/   git worktree of GIT_REF (detached); pipeline/data and .venv are links to the bundle's
#   lib/    this run's structure library: a copy of the bundle's pipeline/lib taken at the start
#   work/   YC_WORK of the run; seeded with the download caches only (PDB files, RCSB answers, pool pairs)
#   build.log
set -euo pipefail
BUNDLE=$(cd "$(dirname "$0")/.." && pwd)
NAME=${1:?run name}; REF=${2:-HEAD}; FIRST=${3:-1}
RUNS=${YC_RUNS:-$BUNDLE/../runs}
RUN=$RUNS/$NAME
mkdir -p "$RUNS"
SHA=$(git -C "$BUNDLE" rev-parse "$REF")
if [ -e "$RUN/tree" ]; then
  # YC_RESUME=1: continue an existing run from FIRST_STAGE with the code of GIT_REF (after a fix); the run keeps its
  # library and output, and build.json records the restart and marks the build as mixed
  [ "${YC_RESUME:-0}" = 1 ] || { echo "$RUN exists; choose another name, remove it (git worktree remove $RUN/tree), or set YC_RESUME=1"; exit 1; }
  git -C "$RUN/tree" checkout -q --detach "$SHA"
  export YC_WORK="$RUN/work" YC_PDB="${YC_PDB:-$BUNDLE/work/pdb}" YC_INPUTS="${YC_INPUTS:-$(cd "$BUNDLE/.." && pwd)/inputs}"
  export YC_SEED=${YC_SEED:-0} YC_NPROC=${YC_NPROC:-$(nproc)} OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
  {
    echo "resume $NAME  commit $SHA  seed $YC_SEED  nproc $YC_NPROC  from stage $FIRST  $(date -Is)"
    cd "$RUN/tree"
    bash scripts/run_all.sh "$FIRST"; echo "run_all exit $?"
    if [ "${YC_VALIDATE:-1}" = 1 ]; then bash scripts/validate_all.sh; echo "validate exit $?"; fi
    echo "end $(date -Is)"
  } >> "$RUN/build.log" 2>&1
  exit 0
fi
mkdir -p "$RUN/work"
git -C "$BUNDLE" worktree add --detach "$RUN/tree" "$SHA" > /dev/null
ln -s "$BUNDLE/.venv" "$RUN/tree/.venv"
ln -s "$BUNDLE/pipeline/data" "$RUN/tree/pipeline/data"
[ -d "$BUNDLE/.chromelibs" ] && ln -s "$BUNDLE/.chromelibs" "$RUN/tree/.chromelibs"
if [ "${YC_LIB_FROM:-bundle}" = empty ]; then mkdir -p "$RUN/lib"; else cp -r "$BUNDLE/pipeline/lib" "$RUN/lib"; fi
ln -s "$RUN/lib" "$RUN/tree/pipeline/lib"
# download caches: without them stage 1 queries RCSB for about an hour and the result depends on the PDB of the day
[ -f "$BUNDLE/work/rcsb_cache.json" ] && cp "$BUNDLE/work/rcsb_cache.json" "$RUN/work/"
mkdir -p "$RUN/work/organisation"
[ -f "$BUNDLE/work/organisation/pool_pairs.parquet" ] && cp "$BUNDLE/work/organisation/pool_pairs.parquet" "$RUN/work/organisation/"
# an input that sits in work/ and has no generator in the repository (found by the first isolated build, 2026-11-13):
# the per-pair summary of the AlphaFold 3 pool predictions, read by features.py and condensates.py
mkdir -p "$RUN/work/interactions"
[ -f "$BUNDLE/work/interactions/pair_summary.parquet" ] && ln -s "$BUNDLE/work/interactions/pair_summary.parquet" "$RUN/work/interactions/"
# Brownian-dynamics crowding sweeps (simgpu/sweep.py): the calibration of the diffusion model's crowding term, read by
# validation/diffusion.py. They were run on boxes of an earlier build and are carried over, not recomputed per build
[ -d "$BUNDLE/work/sweeps" ] && cp -r "$BUNDLE/work/sweeps" "$RUN/work/sweeps"
# measurements that do not depend on the model (hours of image processing): reuse
mkdir -p "$RUN/work/validation"
for f in image_nc_raw.tsv proximity_sets.tsv; do [ -f "$BUNDLE/work/validation/$f" ] && cp "$BUNDLE/work/validation/$f" "$RUN/work/validation/"; done
export YC_WORK="$RUN/work" YC_PDB="${YC_PDB:-$BUNDLE/work/pdb}" YC_INPUTS="${YC_INPUTS:-$(cd "$BUNDLE/.." && pwd)/inputs}"
export YC_SEED=${YC_SEED:-0} YC_NPROC=${YC_NPROC:-$(nproc)} OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
{
  echo "run $NAME  commit $SHA  seed $YC_SEED  nproc $YC_NPROC  start $(date -Is)"
  cd "$RUN/tree"
  bash scripts/run_all.sh "$FIRST"; echo "run_all exit $?"
  if [ "${YC_VALIDATE:-1}" = 1 ]; then bash scripts/validate_all.sh; echo "validate exit $?"; fi
  echo "end $(date -Is)"
} > "$RUN/build.log" 2>&1
