#!/bin/sh
# Rebuild the whole model: inventory -> geometry -> packing -> validation -> web encoding.
# Usage: scripts/run_all.sh [first_stage]   e.g. "3" to restart from packing, "0" to rebuild the base structure library
# first (stage 0 is skipped by default: it rewrites pipeline/lib). Logs go to $YC_WORK/sN.log.
set -e
set -o pipefail          # a failing stage stops the run (its output goes through tee)
cd "$(dirname "$0")"
PY=../.venv/bin/python
WORK=$($PY -c "import sys; sys.path.insert(0, '../pipeline'); import config; print(config.WORK)")
first=${1:-1}
# numba's on-disk cache does not notice changes in functions called from other modules (e.g. cell3d.mem_sdfs inlined
# into pack3d / tiles3d): start every run from a clean cache
rm -f ../pipeline/__pycache__/*.nbi ../pipeline/__pycache__/*.nbc
# Python randomises string hashes per process; fix them so that set iteration order (and with it float sums and
# random draws downstream) cannot differ between runs
export PYTHONHASHSEED=0
# stamp the output folder: commit, master seed (YC_SEED), constants, inputs (scripts/build_stamp.py)
$PY build_stamp.py start "$first"
for s in s0_library s1_inventory s2_geometry s3_pack s4_validate s5_encode; do
  n=${s#s}; n=${n%%_*}
  [ "$n" -lt "$first" ] && continue
  echo "== $s"
  $PY "$s.py" 2>&1 | tee "$WORK/$s.log"
  # the per-protein feature table carries the inventory's copy numbers (validation baselines, diffusion and
  # condensate layers): refresh it whenever the inventory is rebuilt (was stale from 2026-11-11 to build R)
  # (an if, not "[ ] && ...": as the last command of the loop a false test made the script return 1 after stage 5)
  if [ "$s" = s1_inventory ]; then (cd ../pipeline && $PY features.py > "$WORK/features.log" 2>&1); fi
done
