# Whole yeast cell

An interactive, molecule-by-molecule model of a haploid *Saccharomyces cerevisiae* cell (tens of millions of
macromolecules, each a PDB or AlphaFold structure; MODEL.md gives the counts of the current build), rendered in the
browser with WebGPU.

Open `index.html` through a web server (e.g. GitHub Pages); it needs a WebGPU browser (recent Chrome or Edge on
desktop, recent Chrome on Android). On phones and tablets it has touch controls (pinch to zoom).

How the model was built (data sources, assumptions, design decisions, validation, densities and limitations)
is described in [MODEL.md](MODEL.md), which is also shown in the viewer under "Model details".

The code that builds the model and the viewer source are in [source/](source/); see source/README.md for how to rebuild.
