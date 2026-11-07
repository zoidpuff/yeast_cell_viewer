export const ABOUT_HTML = `
<h2>About this model</h2>
<p>This is one complete haploid <i>Saccharomyces cerevisiae</i> cell (strain BY4741, log phase in YPD), built as a
structural whole-cell model in the spirit of Maritan <i>et al.</i> (J. Mol. Biol. 2022; whole <i>Mycoplasma genitalium</i> cell built with cellPACK)
and drawn in the style of David Goodsell. Every macromolecule is drawn from its experimental (PDB) or predicted
(AlphaFold DB) coordinates, at one bead per amino-acid residue or nucleotide when it is large enough on screen.
Further away, molecules switch to a sphere-tree stand-in and then to one sphere each. This only happens once
the detail is smaller than a pixel. Click any molecule to inspect it and to load its full-atom model in Mol*.</p>
<p>This page is a summary. The full description of the model (data sources, assumptions, design decisions, validation,
densities and limitations) is under <a href="#" id="aboutToModel" style="color:#8fc3ff">Model details</a>.</p>

<h3>What is in the cell</h3>
<ul>
<li><b>21.6 million macromolecule instances (most multi-chain assemblies), 5,638 molecular species.</b> Protein abundances are from the PaxDb integrated dataset rescaled to
4.2 × 10⁷ protein molecules per cell (Ho <i>et al.</i>), in a 42 µm³ cell. There are 169,000 80S ribosomes, 3.3 million tRNAs, 73,150 nucleosomes and 117 nuclear pore complexes.</li>
<li><b>Structures.</b> 55% of protein mass is drawn from PDB entries. These include the ribosome 4V88, nuclear pore 7N9F,
nucleosome 1ID3, 26S proteasome 4CR2, FAS 2UV8, Pol I/II/III, V-ATPase 6O7V, Pma1 7VH5, CCT, about 100 further Complex Portal
complexes with a yeast PDB assembly (e.g. phosphofructokinase, eIF2/eIF3, the mitochondrial ribosome, mediator, exosome, snoRNPs with
their RNA) and 133 homo-oligomers. Membrane complexes (ATP synthase, respiratory complexes II-IV, Sec translocon, OST, EMC, glucan
synthase, flippases, ...) are oriented in their membrane automatically from UniProt transmembrane and topology annotations.
14% are AlphaFold DB predicted complexes: 452 homodimers and 505 heterodimers (ipTM &ge; 0.5, partners sharing a
compartment; e.g. actin-cofilin, actin-profilin, RNR2-RNR4, SOD1-CCS1, eIF4G-eIF4A); 30% are AlphaFold monomers.
AlphaFold residues with pLDDT &lt; 50 are trimmed. Outer-membrane beta-barrels (TOM40, SAM50, MDM10) are oriented by their strands.</li>
<li><b>Nuclear pore FG mesh:</b> the disordered FG-repeat domains of NSP1, NUP49, NUP57 and NUP159 (128 per pore, 6.1 million
residues in all) are threaded as self-avoiding residue chains from where the 7N9F model ends, filling the pore channel and
reaching into the cytoplasm (checkbox "FG repeats").</li>
<li><b>mRNA and polysomes:</b> 30,000 full-length mRNAs (lognormal lengths, mean 1.5 kb including UTRs and poly(A);
40 million nucleotides, one bead per nucleotide) are threaded through the cytosol before the proteins are packed.
23,000 of them carry 91,000 of the 80S ribosomes (about 4 per mRNA): the mRNA runs through each ribosome's channel
(path from PDB 7ZW0), short linkers join neighbouring ribosomes, and the rest of the transcript leaves the first and
last ribosome as 3&prime; and 5&prime; tails. The other 7,000 mRNAs carry no ribosomes. Threads that run out of room end
early (86% of the target length is reached overall). Linker DNA joins consecutive nucleosomes.</li>
<li><b>ER-bound ribosomes:</b> 4,300 ribosomes sit on translocons in the ER and outer nuclear membrane: on the Ssh1 complex
(yeast ribosome-Ssh1 cryo-EM structure 2WW9) or on Sec61 (7KAH, placed through the shared Sss1 subunit). Their number is
limited by the free translocon pools (80% assumed ribosome-bound); the PM-facing side of the cortical ER is too narrow for ribosomes.</li>
<li><b>Microtubules:</b> 21 microtubules (11,700 &alpha;/&beta;-tubulin dimers) grow from &gamma;-tubulin rings on the spindle pole body,
opposite the nucleolus: 16 short kinetochore microtubules whose plus ends anchor the 16 chromosomes (Rabl configuration),
2 more nuclear ones and 3 cytoplasmic ones about 1 &micro;m long that curve around organelles. The 13-protofilament lattice
and the ring come from the native yeast &gamma;-TuRC-capped microtubule 8QV3; the lattice dimer is 5W3F (8.3 nm repeat).</li>
<li><b>Mitochondrial cristae:</b> all 11,800 ATP synthase copies form 5,900 dimers in rows along the rims of the lamellar
cristae, every 13 nm. Each dimer is two copies of the yeast monomer 8F39 placed in the dimer arrangement of 4B2Q
(cryo-ET), whose 84&deg; angle matches the curvature of the rim (crista lumen 18 nm wide). Complex IV is entirely in
1,560 respiratory supercomplexes III<sub>2</sub>IV<sub>2</sub> (6HU9) in the inner membrane, next to 1,400 free complex III dimers.</li>
<li><b>Actin cables:</b> 8 bundles of 8 filaments, 12 nm apart (33,000 actin subunits, 92 &micro;m), run 90 nm beneath the
plasma membrane, built from the yeast F-actin cryo-EM structure 9GO5 (167&deg; and 2.76 nm per subunit). Both strands of
every filament are coated head-to-tail by 6,600 tropomyosin TPM1 dimers (200 residues, one per 5 actins), following the
tropomyosin path of the mammalian actin-tropomyosin structure 5JLF wound as an ideal coiled coil.</li>
<li><b>Actin patches:</b> 22 endocytic patches under the plasma membrane, each a branched network of about 1,500 actin
subunits: short filaments branch at 73&deg; from Arp2/3 complexes (2,200 branches; geometry from the branch-junction
structure 8UXW) and are capped by capping protein (8F8Q), barbed ends toward the membrane. The yeast Arp2/3 complex
(all 7 subunits) and capping protein are AlphaFold models superposed on those templates. A further 4,800 Arp2/3
complexes are free in the cytosol. Assumption: 30% of actin in cables, 30% in patches, the rest monomeric
(much of it bound to profilin, cofilin or Srv2).</li>
<li><b>Membranes</b> (plasma membrane, cortical and perinuclear ER, nuclear envelope with pores, vacuole, a reticular
mitochondrion, Golgi/endosomes, peroxisomes, lipid droplets) are represented by signed-distance fields. The cell
outline, vacuoles and lipid droplets are those of a cell segmented from a cryo-FIB/SEM volume (EMPIAR-11416, Dumoux
<i>et al.</i> 2023), the nucleus that of a second cell; the shapes are smoothed, and the vacuoles' thickness along the
milling direction, which the segmentation loses, is reconstructed.
Near the camera every lipid is drawn individually. Class composition (Ejsing <i>et al.</i> 2009): PI 20%, PE 15%, PC 14%, PA 8%, PS 2%,
ergosterol 12%, 0.65 nm² per lipid. Lipids are removed wherever they would overlap a membrane protein.</li>
<li><b>Cell wall:</b> an inner β-1,3-glucan layer and an outer mannoprotein layer with mannan brushes. Bead counts follow a
wall mass of about 20% of dry mass.</li>
</ul>

<h3>How it was built</h3>
<p>Proteins were assigned to compartments from UniProt subcellular-location annotations, with GO as fallback. Membrane proteins were oriented by their
transmembrane segments. Everything was then packed into 200 nm tiles with collision detection on residue-level
sphere trees. Placement is not purely random: proteins of the same STRING interaction module are co-enriched in
domains a few hundred nm across, and about a fifth of all molecules were placed next to an interaction partner predicted
by AlphaFold 3 or linked in STRING (an assumed, tunable organisation; see Model details). Validation on the final model:</p>
<ul>
<li>Bead clash rate 1.9×10⁻⁶.</li>
<li>Containment violations 5.6×10⁻⁵.</li>
<li>Macromolecule densities (protein plus the RNA of ribosomes and tRNAs and the DNA of nucleosomes; lipids and small
molecules not counted): cytosol 93 mg/mL (70 protein + 23 RNA), nucleoplasm 56 mg/mL, mitochondrial matrix 111 mg/mL.</li>
<li>99.9% of requested copies were placed outside mitochondria.</li>
</ul>
<p>Positions are stored to 0.0055 nm and orientations to 0.004°.</p>

<h3>Honest caveats</h3>
<ul>
<li><b>Layout.</b> The organelle shapes and positions are a plausible illustration, not a tomogram of one specific cell.</li>
<li><b>Abundances.</b> Absolute abundances are uncertain to about 2-fold.</li>
<li><b>Mitochondria.</b> The mitochondrion saturated: only 65% of matrix and 86% of IMS proteins fit, so their copies are under-represented.</li>
<li><b>Chromatin</b> is a coarse Rabl-like arrangement of nucleosome chains in 59 fragments, each chromosome starting at its kinetochore microtubule. rDNA repeats are not modelled.
Linker DNA and mRNA paths are schematic.</li>
<li><b>Procedural content.</b> Lipids, β-glucan and mannan are generated procedurally in your browser from deterministic seeds. They are not stored
coordinates, and polysaccharides are only partly clash-checked.</li>
<li><b>Spindle pole body.</b> Only the &gamma;-tubulin rings are modelled, not the SPB plaques; PaxDb lists fewer SPC98 copies than 21 rings need.</li>
<li><b>Not included:</b> septins, fimbrin cross-links, ribosome-bound chaperones (NAC, RAC/Ssb), metabolites and water.
The FG nucleoporins absent from 7N9F (NUP100, NUP116, NUP145N, NUP42 and the basket NUP1, NUP60, NUP2) have no FG domains.</li>
<li><b>Dropped copies.</b> 9,946 transmembrane-protein copies assigned to soluble compartments were dropped, and 6,888 copies had no structure.</li>
</ul>
<p style="color:#a9a69e">Renderer: WebGPU compute culling and level-of-detail selection, sphere impostors with exact cut caps, analytic
compartment backgrounds and Goodsell-style outlines. Full-atom models are fetched live from RCSB PDB and AlphaFold DB
when you request them.</p>
<p style="color:#a9a69e">Data sources: RCSB Protein Data Bank (CC0) and the individual entries named above; AlphaFold Protein Structure
Database monomers and predicted complexes (CC-BY 4.0, Jumper <i>et al.</i> 2021; Varadi <i>et al.</i> 2024); PaxDb protein abundances;
Complex Portal; UniProt and Gene Ontology annotations; IntAct interactions; STRING v12 associations (CC-BY 4.0); EMPIAR-11416 cryo-FIB/SEM volume (CC0); Alliance of Genome Resources paralogs;
lipid composition from Ejsing <i>et al.</i> 2009; total protein from Ho <i>et al.</i> 2018.</p>
`;
