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
<li><b>35.8 million macromolecule instances (most multi-chain assemblies), 5,615 molecular species.</b> Protein abundances are from the PaxDb integrated dataset rescaled to
5.0 pg of protein per cell (measured 4–6 pg), 7.8 × 10⁷ protein molecules in a 42 µm³ cell. There are 301,000 80S ribosomes, 28,000 pre-ribosomes, 3.3 million tRNAs, 78,000 nucleosomes and 117 nuclear pore complexes.</li>
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
25,600 of them carry 121,000 of the 80S ribosomes (about 5 per mRNA): the mRNA runs through each ribosome's channel
(path from PDB 7ZW0), short linkers join neighbouring ribosomes, and the rest of the transcript leaves the first and
last ribosome as 3&prime; and 5&prime; tails. The other 4,400 mRNAs carry no ribosomes. Threads that run out of room end
early (86% of the target length is reached overall). Linker DNA joins consecutive nucleosomes.</li>
<li><b>ER-bound ribosomes:</b> 8,100 ribosomes sit on translocons in the ER and outer nuclear membrane: on the Ssh1 complex
(yeast ribosome-Ssh1 cryo-EM structure 2WW9) or on Sec61 (7KAH, placed through the shared Sss1 subunit). Their number is
limited by the free translocon pools (80% assumed ribosome-bound); the PM-facing side of the cortical ER is too narrow for ribosomes.
A further 15,000 free ribosomes sit next to the ER and outer nuclear membrane, bringing the density to the 550&ndash;900 per &micro;m&sup2;
seen in tomograms (West <i>et al.</i> 2011).</li>
<li><b>RAC on ribosomes:</b> all 91,000 ribosome-associated complexes (Hsp40 Zuo1 with Hsp70 Ssz1) sit on ribosomes at the
tunnel exit, posed as in the 80S&ndash;RAC structure 7Z3N (<i>C. thermophilum</i> RAC standing in for the yeast one).</li>
<li><b>Chromosomes:</b> the Rabl configuration: 16 centromeres at the spindle pole body, each of the 32 arms a nucleosome chain
running to its telomere at the nuclear envelope; telomeres gather in 7 foci, longer arms reaching farther from the pole
(telomere distance follows arm length, Spearman 0.71); chromosome XII runs through the nucleolus. 550 transcribing Pol II sit
on their genes (per-gene counts from run-on, Pelechano <i>et al.</i> 2010; nucleosomes evicted around each), DNA through
the cleft and nascent mRNA from the exit.</li>
<li><b>mtDNA:</b> 15 nucleoids in the matrix holding 35 copies of the 85.8 kb mitochondrial genome (3 Mb of DNA).</li>
<li><b>Cvt complex:</b> one cluster of precursor aminopeptidase I dodecamers, Ams1 and the receptor Atg19, 180 nm across, next
to the vacuole (the cytoplasm-to-vacuole targeting cargo).</li>
<li><b>Ribosome assembly:</b> 15,000 pre-ribosomes (about 9% of ribosomes) follow the assembly pathway from the nucleolus
to the cytoplasm: 90S particles (SSU processome, 5WLC) and nucleolar pre-60S (7R7A) in the nucleolus, Rix1&ndash;Rea1
pre-60S (6YLH) in the nucleoplasm, and Nmd3&ndash;Lsg1 pre-60S (5T62) and pre-40S (6FAI) in the cytoplasm, each with its
pre-rRNA. Their numbers follow the abundances of the assembly factors specific to each stage. In the nucleolus, 75 active
rRNA genes carry 50 RNA polymerase I each on nucleosome-free DNA, with nascent pre-rRNA growing longer along each gene and
90S particles on uncleaved transcripts (the "Christmas trees" of Miller spreads); the inactive repeats are nucleosomal.
The <i>rRNA genes</i> guided view shows them on their own.</li>
<li><b>Microtubules:</b> 21 microtubules (11,700 &alpha;/&beta;-tubulin dimers) grow from &gamma;-tubulin rings on the spindle pole body,
opposite the nucleolus: 16 short kinetochore microtubules whose plus ends anchor the 16 chromosomes (Rabl configuration),
2 more nuclear ones and 3 cytoplasmic ones about 1 &micro;m long that curve around organelles. The 13-protofilament lattice
and the ring come from the native yeast &gamma;-TuRC-capped microtubule 8QV3; the lattice dimer is 5W3F (8.3 nm repeat).</li>
<li><b>Mitochondrial cristae:</b> 11,900 of the 21,800 ATP synthase copies form 5,900 dimers in rows along the rims of the lamellar
cristae, every 13 nm. Each dimer is two copies of the yeast monomer 8F39 placed in the dimer arrangement of 4B2Q
(cryo-ET), whose 84&deg; angle matches the curvature of the rim (crista lumen 18 nm wide). Complex IV is entirely in
2,900 respiratory supercomplexes III<sub>2</sub>IV<sub>2</sub> (6HU9) in the inner membrane, next to 2,600 free complex III dimers.</li>
<li><b>Actin cables:</b> 15 bundles of filaments 12 nm apart (62,000 actin subunits, 171 &micro;m), run 90 nm beneath the
plasma membrane, built from the yeast F-actin cryo-EM structure 9GO5 (167&deg; and 2.76 nm per subunit). Both strands of
every filament are coated head-to-tail by 11,650 tropomyosin TPM1 dimers (200 residues, one per 5 actins), following the
tropomyosin path of the mammalian actin-tropomyosin structure 5JLF wound as an ideal coiled coil.</li>
<li><b>Actin patches:</b> 41 endocytic sites, each a plasma-membrane invagination (a tube about 40 nm wide, 0&ndash;100 nm deep
at random stages; Kukulski <i>et al.</i> 2012) inside a branched network of about 1,500 actin subunits: short filaments branch at 73&deg; from Arp2/3 complexes (4,200 branches; geometry from the branch-junction
structure 8UXW) and are capped by capping protein (8F8Q), barbed ends toward the membrane. The yeast Arp2/3 complex
(all 7 subunits) and capping protein are AlphaFold models superposed on those templates. A further 4,800 Arp2/3
complexes are free in the cytosol. Assumption: 30% of actin in cables, 30% in patches, the rest monomeric
(much of it bound to profilin, cofilin or Srv2).</li>
<li><b>Endocytic sites</b> also carry their coat at the tip (Sla2, Sla1, End3; numbers per site from Picco <i>et al.</i> 2015)
and amphiphysin Rvs161/167 on the tubes deeper than 50 nm.</li>
<li><b>Nucleus&ndash;vacuole junction:</b> 758 Vac8&ndash;Nvj1 complexes (5XJG) where a vacuole lobe lies against the nucleus.</li>
<li><b>Eisosomes:</b> 27 plasma-membrane furrows, about 300 nm long and 50 nm deep, whose half-pipe bottoms are lined by
4,700 Pil1 and Lsp1 BAR-domain dimers in their native lattice (8QB7, 8QB8; Kefauver <i>et al.</i> 2024), amphipathic
helices in the inner leaflet. The coats hold 6.7% of the cell's Pil1 and Lsp1; the rest stay free near the membrane.
The MCC markers Sur7 and Nce102 crowd the furrow walls.</li>
<li><b>Membranes</b> (plasma membrane, cortical and perinuclear ER, nuclear envelope with pores, vacuole, a reticular
mitochondrion, Golgi/endosomes, peroxisomes, lipid droplets) are represented by signed-distance fields. The cell
outline, vacuoles and lipid droplets are those of a cell segmented from a cryo-FIB/SEM volume (EMPIAR-11416, Dumoux
<i>et al.</i> 2023), the nucleus that of a second cell; the shapes are smoothed, and the vacuoles' thickness along the
milling direction, which the segmentation loses, is reconstructed.
Near the camera every lipid is drawn individually, with the composition of each membrane from isolated organelle
membranes (Zinser <i>et al.</i> 1991), cardiolipin in the inner mitochondrial membrane and PS/PE in the cytoplasmic leaflet
of the plasma membrane; 0.65 nm² per lipid. Lipids are removed wherever they would overlap a membrane protein.</li>
<li><b>Vacuolar polyphosphate:</b> the vacuole lumen holds little protein but is packed with short polyphosphate chains
(about 300 mM phosphate units, an assumption within the reported "hundreds of millimolar"), the cell's phosphate store,
drawn procedurally when you zoom in.</li>
<li><b>Cell wall:</b> an inner β-1,3-glucan layer and an outer mannoprotein layer with mannan brushes. Bead counts follow a
wall mass of about 20% of dry mass. The cell is drawn as a mother that has budded once: one bud scar, a chitin ring 1 &micro;m
across whose rim (half the wall's chitin, as microfibrils) stands about 50 nm above the surface.</li>
</ul>

<h3>How it was built</h3>
<p>Proteins were assigned to compartments from UniProt subcellular-location annotations, with GO as fallback, and measured shares from the yeast GFP collection (CYCLoPs) and mitochondrial proteomics (Vögtle 2017, Morgenstern 2017). Membrane proteins were oriented by their
transmembrane segments. Everything was then packed into 200 nm tiles with collision detection on residue-level
sphere trees. Placement is not purely random: proteins with similar STRING network neighbourhoods (personalized
PageRank profiles, no clustering) are co-enriched in regions a few hundred nm across, and about a quarter of all
molecules were placed next to an interaction partner predicted by AlphaFold 3 or a close network neighbour (an assumed,
tunable organisation within each compartment; see Model details). Random addition jams below cellular crowding, so
tiles that could not fit everything were relaxed: overlapping molecules pushed apart as rigid bodies and the missing copies
added in rounds. Validation on the final model:</p>
<ul>
<li>Bead clash rate 9×10⁻⁶.</li>
<li>Containment violations 1.4×10⁻⁴.</li>
<li>Macromolecule densities (protein plus RNA and DNA; lipids and small molecules not counted): cytosol 163 mg/mL
(127 protein + 36 RNA), nucleoplasm 103 mg/mL, nucleolus 137 mg/mL, mitochondrial matrix 107 mg/mL. With lipids and
small molecules this comes to roughly the measured dry mass of the cytoplasm.</li>
<li>98.2% of requested copies were placed outside mitochondria (94.7% overall).</li>
</ul>
<p>Positions are stored to 0.0055 nm and orientations to 0.004°.</p>

<h3>Simulating a piece of the cell</h3>
<p><b>Simulate a box</b> in the side panel cuts a 30–200 nm cube out of the packed cell and runs Brownian dynamics of
its molecules on your GPU, in WebGPU compute shaders written for this viewer. The rest of the cell stops rendering meanwhile.
Each molecule is a rigid body made of collision spheres at one of three resolutions: its sphere tree (<i>schematic</i>, about
30 residues per sphere, microseconds), clusters of about 4 residues (<i>mid</i>), or every exposed residue (<i>residue</i>,
nanoseconds). Diffusion comes from its shape (Kirkwood hydrodynamic radius and Stokes–Einstein), with mean-field
hydrodynamics for the crowding. Molecules repel on overlap and interact by:</p>
<ul>
<li><b>Screened electrostatics</b> between charged residues and nucleotides, at the cytosol's salt and pH.</li>
<li><b>Specific binding</b> at the predicted interfaces of 371 AlphaFold dimers.</li>
<li><b>mRNA links</b> between the ribosomes of a polysome.</li>
<li><b>An optional nonspecific stickiness.</b> Without it, proteins diffuse about 5× faster than measured in yeast.
Stickiness can close that gap, as in bacterial cytoplasm models, but here it does so by clustering the molecules into
networks. A higher viscosity, standing in for the many-body hydrodynamics the model lacks, is the cleaner option.</li>
</ul>
<p>Membranes act as walls, and membrane proteins diffuse within them. Filaments and chromatin stay fixed. The box is
periodic, or embedded in the cell by a tethered outer shell when it contains membranes.</p>
<ul>
<li><b>Running:</b> watch it live, or compute a stretch first and play it back at a chosen frame rate.</li>
<li><b>Tracking:</b> click molecules to track them; the others fade and the tracked paths are drawn.</li>
<li><b>Readout:</b> compares the simulated diffusion with yeast measurements.</li>
<li><b>Export MD files:</b> saves
  <ul>
  <li>a coarse-grained mmCIF;</li>
  <li>a manifest with each molecule's structure source, pose, charges, links and contacts;</li>
  <li>the trajectory, with a script that turns it into a DCD for VMD or ChimeraX;</li>
  <li>a ready-to-run OpenMM script with the same interactions.</li>
  </ul></li>
</ul>
<p>See Model details, section 14.</p>

<h3>Honest caveats</h3>
<ul>
<li><b>Layout.</b> The organelle shapes and positions are a plausible illustration, not a tomogram of one specific cell.</li>
<li><b>Abundances.</b> Absolute abundances are uncertain to about 2-fold.</li>
<li><b>Full compartments.</b> Small, crowded compartments saturated: only 40% of the mitochondrial matrix proteins, 80% of the IMS, 77% of the nucleolus and 81% of the Golgi/endosomes fit, so their copies are under-represented.</li>
<li><b>Histones</b> beyond the fixed number of nucleosomes are not removed: 80,000 surplus H2B are left out and about 50,000 free H4 and histone dimers sit in the nucleoplasm, although cells keep free histones low.</li>
<li><b>Chromatin</b> follows the Rabl layout at the level of arms (kb mapped linearly along each walk); there are no
sequence-level nucleosome positions or inter-chromosomal contacts. Linker DNA and mRNA paths are schematic.</li>
<li><b>Procedural content.</b> Lipids, β-glucan and mannan are generated procedurally in your browser from deterministic seeds. They are not stored
coordinates, and polysaccharides are only partly clash-checked.</li>
<li><b>Spindle pole body.</b> Only the &gamma;-tubulin rings are modelled, not the SPB plaques; PaxDb lists fewer SPC98 copies than 21 rings need.</li>
<li><b>Not included:</b> septins, fimbrin cross-links, the ribosome-bound chaperones NAC and Ssb, metabolites and water.
The FG nucleoporins absent from 7N9F (NUP100, NUP116, NUP145N, NUP42 and the basket NUP1, NUP60, NUP2) have no FG domains.</li>
<li><b>Dropped copies.</b> 15,982 transmembrane-protein copies assigned to soluble compartments were dropped, and 90,403 copies had no structure (80,000 of them the surplus histone H2B).</li>
</ul>
<p style="color:#a9a69e">Renderer: WebGPU compute culling and level-of-detail selection, sphere impostors with exact cut caps, analytic
compartment backgrounds and Goodsell-style outlines. Full-atom models are fetched live from RCSB PDB and AlphaFold DB
when you request them.</p>
<p style="color:#a9a69e">Data sources: RCSB Protein Data Bank (CC0) and the individual entries named above; AlphaFold Protein Structure
Database monomers and predicted complexes (CC-BY 4.0, Jumper <i>et al.</i> 2021; Varadi <i>et al.</i> 2024); PaxDb protein abundances;
Complex Portal; UniProt and Gene Ontology annotations; IntAct interactions; STRING v12 associations (CC-BY 4.0); EMPIAR-11416 cryo-FIB/SEM volume (CC0); Alliance of Genome Resources paralogs;
lipid composition from Ejsing <i>et al.</i> 2009; total protein from Ho <i>et al.</i> 2018.</p>
`;
