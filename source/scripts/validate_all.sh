#!/bin/bash
# Re-run every comparison of the model with data after a rebuild (run_all.sh), regenerate the figures and the PDF of
# model_validation_and_tuning.md. Logs: $WORK/validate_all.log. Needs the GPUs (contact pass) and, for the PDF,
# the headless-browser libraries in .chromelibs (HANDOVER).
set -e
set -o pipefail
cd "$(dirname "$0")/.."
B=$(pwd)
PY="$B/.venv/bin/python"
WORK=$($PY -c "import sys; sys.path.insert(0, 'pipeline'); import config; print(config.WORK)")
export YC_NPROC=${YC_NPROC:-$(nproc)}
export PYTHONHASHSEED=0
log() { echo "$(date +%H:%M:%S) $*"; }
{
log "contact pass (decode tiles, GPU sphere-tree gaps)"
rm -f "$WORK/contacts/instances.npz"
(cd pipeline && $PY contacts.py)
log "pair statistics for kinases"
(cd pipeline && $PY encounters.py --kinases)
cd validation
log "localization benchmark (dev, test)"
for s in dev test; do $PY localization.py --split $s --out "$WORK/validation/latest" | grep -v Warn; done
log "proximity tests"
[ -f "$WORK/validation/proximity_sets.tsv" ] || $PY proximity_sets.py
$PY proximity_test.py | grep -A12 "gap 10"
$PY ppiseq_test.py | tail -3
$PY singh2020_test.py
log "localization of signalling proteins"
$PY kinase_localization.py | tail -1
log "ribosomes vs cryo-ET, literature numbers, copy numbers"
$PY ribosome_cryoet.py | grep -E "model_uM|nn_median"
(cd ../pipeline && $PY genome_map.py > /dev/null)
$PY chromatin_stats.py
$PY hic_test.py
$PY nucleosome_test.py > /dev/null
$PY mrna_stats.py > /dev/null
$PY mrna_localization.py > /dev/null
$PY trna_nucleolus.py | grep -E 'within_overlap'
$PY dnase_test.py 2>/dev/null | grep -E 'core_over_linker|protected'
$PY silent_chromatin_test.py > /dev/null 2>&1
[ -f ../pipeline/data/validation/microscopy/yeastrgb_dubreuil2019/raw/RAW-C_SWAT_GFP-ORF.zip ] && $PY image_nc.py 2>/dev/null | grep -E '"spearman"|auroc'
$PY er_ribosome_density.py | grep -E '"per_um2"'
$PY literature_checks.py
$PY copy_numbers.py > /dev/null
$PY diffusion.py > /dev/null
$PY anchors.py > /dev/null
log "placement checks (who is missing, tile seams, membrane shells and layers next to membranes)"
$PY placement_checks.py
log "negative control: no condensate in the log-phase cytosol (mean-field mixture, validation/cond_test.py)"
mkdir -p "$WORK/condensates"
$PY cond_test.py > "$WORK/condensates/cond_test.log" 2>&1 || echo "cond_test failed (see $WORK/condensates/cond_test.log)"
(grep -E "^log phase|^stress" "$WORK/condensates/cond_test.log" | cut -c1-160) || true
log "simulated microscopy"
$PY microscopy_gallery.py | grep -v Warn
$PY microscopy_idr0078.py | grep -v Warn | tail -1
log "figures"
$PY make_figures.py | grep -v Warn
cd ..
log "PDF"
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$B/.chromelibs FONTCONFIG_FILE=$B/.chromelibs/fonts/fonts.conf \
  $PY scripts/md_to_pdf.py model_validation_and_tuning.md docs/model_validation_and_tuning.pdf
log "result files against the build stamp"
$PY scripts/build_stamp.py check
log "done"
} 2>&1 | tee "$WORK/validate_all.log"
