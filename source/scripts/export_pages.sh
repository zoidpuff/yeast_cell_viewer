#!/bin/sh
# Build a self-contained static copy of the viewer for GitHub Pages (or any static host).
# Usage: scripts/export_pages.sh [out_dir]      (default: ~/Downloads/yeast_cell_pages)
# Re-running updates the folder in place and keeps its .git, so it can be committed and pushed again.
set -e
cd "$(dirname "$0")/.."
BUNDLE=$(pwd)
OUT=${1:-$HOME/Downloads/yeast_cell_pages}
export PATH="$BUNDLE/.node/bin:$PATH"
DATA=$(.venv/bin/python -c "import sys; sys.path.insert(0, 'pipeline'); import config; print(config.WEB_OUT)")
[ -f "$DATA/library/types.json" ] || { echo "no encoded data at $DATA (run s5_encode.py)"; exit 1; }

cd web/app
[ -d node_modules ] || npm ci
mkdir -p public && rm -f public/data
[ -d public/molstar ] || cp -r node_modules/molstar/build/viewer public/molstar
npx vite build

mkdir -p "$OUT"
cd "$OUT"
rm -rf index.html assets molstar data
cp "$BUNDLE/web/app/dist/index.html" index.html
cp -r "$BUNDLE/web/app/dist/assets" assets
cp -r "$BUNDLE/web/app/dist/molstar" molstar
cp -r "$DATA" data
cp "$BUNDLE/MODEL.md" MODEL.md       # full model description (also shown in the viewer under "Model details")
# GitHub Pages does not compress .bin files: pre-compress the uncompressed library files under their own names
# (the viewer recognises gzip data by its magic bytes and inflates it itself)
for f in data/library/beads.bin data/library/spheretree.bin data/library/types.json data/types_tiles.bin \
         data/geometry.json data/geometry_fields.bin data/bookmarks.json; do
  [ -f "$f" ] && gzip -9 -n -c "$f" > "$f.tmp" && mv "$f.tmp" "$f"
done
# build stamp: appended to every data URL, so browsers never mix cached files from two builds
printf '{"v": "%s"}\n' "$(date +%Y%m%d%H%M%S)" > data/version.json
# source code (pipeline, scripts, viewer source); inputs and the structure library are public downloads, see source/README.md
rm -rf source && mkdir -p source/pipeline source/scripts source/web/app source/env
cp "$BUNDLE"/pipeline/*.py source/pipeline/
cp "$BUNDLE"/scripts/*.py "$BUNDLE"/scripts/*.sh "$BUNDLE"/scripts/README.md source/scripts/
cp -r "$BUNDLE"/web/app/src "$BUNDLE"/web/app/index.html "$BUNDLE"/web/app/package.json "$BUNDLE"/web/app/package-lock.json \
      "$BUNDLE"/web/app/tsconfig.json "$BUNDLE"/web/app/vite.config.ts "$BUNDLE"/web/app/vite.single.config.ts source/web/app/
cp "$BUNDLE"/env/requirements.txt source/env/
cp "$BUNDLE/MODEL.md" source/MODEL.md     # the viewer build imports it from here
cp "$BUNDLE/README_source.md" source/README.md
# what MODEL.md and the READMEs refer to (review of 2026-11-12, comment 19): the validation and GPU-simulation code,
# the validation document, the backlog and progress log, and the two geometry inputs kept in the repository
mkdir -p source/validation source/simgpu source/docs
cp "$BUNDLE"/validation/*.py source/validation/
cp "$BUNDLE"/simgpu/*.py source/simgpu/
cp "$BUNDLE"/model_validation_and_tuning.md source/
cp -r "$BUNDLE"/docs/figures source/docs/
cp "$BUNDLE"/pipeline/cell_meshes.npz "$BUNDLE"/pipeline/organelle_stats.json source/pipeline/ 2>/dev/null || true
[ -f "$BUNDLE/../IMPROVEMENTS.md" ] && cp "$BUNDLE/../IMPROVEMENTS.md" source/IMPROVEMENTS.md
mkdir -p source/pipeline/data/validation
cp "$BUNDLE"/pipeline/data/validation/diffusion_measurements.tsv "$BUNDLE"/pipeline/data/validation/membrane_diffusion_measurements.tsv \
   source/pipeline/data/validation/ 2>/dev/null || true
[ -f "$BUNDLE/LICENSE" ] && cp "$BUNDLE/LICENSE" LICENSE
touch .nojekyll        # serve files as they are (no Jekyll processing)
cat > README.md <<'EOF'
# Whole yeast cell

An interactive, molecule-by-molecule model of a haploid *Saccharomyces cerevisiae* cell (tens of millions of
macromolecules, each a PDB or AlphaFold structure; MODEL.md gives the counts of the current build), rendered in the
browser with WebGPU.

Open `index.html` through a web server (e.g. GitHub Pages); it needs a WebGPU browser (recent Chrome or Edge on
desktop, recent Chrome on Android). On phones and tablets it has touch controls (pinch to zoom).

How the model was built (data sources, assumptions, design decisions, validation, densities and limitations)
is described in [MODEL.md](MODEL.md), which is also shown in the viewer under "Model details".

The code that builds the model and the viewer source are in [source/](source/); see source/README.md for how to rebuild.
EOF
echo "exported to $OUT: $(du -sh . | cut -f1), $(find . -type f -not -path './.git/*' | wc -l | tr -d ' ') files;" \
     "largest $(find . -type f -not -path './.git/*' -exec ls -l {} + | sort -k5 -n | tail -1 | awk '{print $5, $9}')"
