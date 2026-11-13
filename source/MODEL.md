# Whole yeast cell: model description

This document describes how the structural model shown in the viewer was built. It covers what the model contains, which data and structures it comes from, the assumptions and design decisions behind it, how it was validated, and where it is known to be wrong or incomplete. It is written for readers who want to judge the model scientifically or reuse its numbers.

All numbers refer to the current build (build S of 2026-11-13, candidate build 3 of round S: commit 9568da3, master seed 0, build id `20261113-155407-9568da3-s0`; 5,562 molecular species, 24.9 million packed instances, 25.3 million with the pre-placed structures) unless stated otherwise. Build S keeps the cell and the protein copy numbers of build R (3.5 pg of protein in 34 µm³ instead of 5.0 pg in 42 µm³; section 2.1); its inventory differs only in the nucleus/cytosol split, which now uses the compartment volumes as built (section 4.2), and in free H2B, which now has a structure (section 2.3). It changes the vacuole (one near-spherical body instead of 12 lobes), the mitochondria, the side of the nucleolus, the pores at the nucleus–vacuole junction and the packing algorithm (sections 3.3, 6.7, 6.10 and 8), and it is the first build made in an isolated run folder with a build stamp (section 13). Where a build R number is informative it is given as "(build R: x)"; numbers marked "build R" alone were not re-measured on build S. Numbers quoted for builds up to Q refer to the larger cell. Build S is one build at one seed: the spread between replicate builds is not yet measured (section 11).

## 1. What the model is

- **The cell.** One haploid *Saccharomyces cerevisiae* cell (strain BY4741 / S288C background) growing exponentially in rich glucose medium (YPD). The cell is unbudded, and the microtubules and chromatin are arranged as in G1. Since build R its size and protein content are also those of one G1 cell with one genome, not of the average cell of a growing culture (section 2.1).
- **What a molecule is.** Every macromolecule is an individual copy of an experimental (PDB) or predicted (AlphaFold DB) structure. It is drawn at one bead per amino-acid residue (Cα) or per nucleotide (P or C4′). There are 18.3 billion beads in all (build R: 17.9 billion).
- **How the cell was assembled.** Copy numbers come from proteomics. Compartment volumes and membrane areas come from imaging literature. The organelles are built as analytic geometry, and molecules are packed into them with residue-level collision checks.
- **What it is not:**
  - **Not a dynamic or physical simulation.** There are no forces, no diffusion and no binding equilibria. Except where stated below, molecules are placed at random positions and orientations by random sequential addition. A snapshot from a real cell would show interaction-driven clustering that the model lacks.
  - **Not a reconstruction of a particular cell.** The organelle volumes and areas are taken from measurements, but their shapes and positions are a plausible arrangement, not a tomogram.
- **Lineage.** The approach follows the cellPACK-style whole-cell models of *Mycoplasma genitalium* (Maritan *et al.* 2022). The visual style follows David Goodsell.

## 2. Composition

### 2.1 Proteins

- **Abundances:** PaxDb integrated whole-organism dataset for *S. cerevisiae* (taxon 4932), in ppm, matched to UniProt by ORF.
- **Absolute scale (2026-11-12, build R):** ppm values are scaled so that the proteins weigh 3.5 pg. With the abundance-weighted mean mass of 38.7 kDa this is 5.45 × 10⁷ protein molecules.
  - **Why not the bulk value:** bulk measurements give 4–5 pg (Carroll *et al.* 2011) and 6 pg (Sherman 2002) per haploid cell, and builds up to Q used 5.0 pg in the culture's mean cell volume of 42 µm³. Those are averages over a growing culture, which includes budded cells with up to two genomes. The model cell is one unbudded G1 cell with one genome, a fixed number of pores and one spindle pole body.
  - **Anchors (`validation/anchors.py`):** structures counted per cell say how much protein this cell holds. At 5.0 pg the core histones were 1.3 times what one genome's nucleosomes need (so 3.9 pg), and the scaffold and FG nucleoporins 1.5 times what the pores need (3.2 pg). 3.5 pg lies between the two.
  - **Check, not used:** the scale gives 211,500 ribosomes. Counts from rRNA are about 200,000 (Warner 1999) and 187,000 ± 56,000 (von der Haar 2008); a growth-law estimate for YPD is about 3 × 10⁵ (Kostinski & Reuveni 2021).
  - **At 3.5 pg:** core histones are 0.9 times the nucleosomes' need (median of four; H3 0.70, H2B 1.31), nucleoporins 1.36 times the need of 93 pores.
  - **The cell shrinks with it** (section 3.1), so that concentrations stay where they were: the whole-cell protein concentration is 103 mg/mL (119 before), near the 100 mg/mL commonly used (Ohsawa *et al.* 2026 and references there), above a recent Bradford-based 71 mg/mL (Kunzi *et al.* 2026, preprint).
- **Earlier total:** earlier versions used the unified total of Ho, Baryshnikova and Brown (2018), 4.2 × 10⁷ molecules. That comes to only 2.7 pg of protein and left the cytosol at about 90 mg/mL of macromolecules.
- **Uncertainty:**
  - PaxDb integrated values are not specific to BY4741 in YPD. Per-protein uncertainty is roughly 2-fold.
  - The absolute total is still the largest single uncertainty: the anchors span 3.2–3.9 pg, and bulk measurements of average cells 4–6 pg.
- **Proteins this cell does not have (2026-11-12; `inventory.NOT_IN_THIS_CELL`):** PaxDb mixes strains, media and stages. Curated lists set 73 proteins to zero before scaling: the genes deleted in BY4741 (His3, Leu2, Met17, Ura3), α-specific proteins in this MATa cell (8), meiosis and sporulation proteins (42), glucose-repressed proteins (11: Gal1, Pck1, Mls1, Icl1, Cyc7 and others) and proteins degraded before G1 (8: Mcd1, Cdc5, Sgo1, Pds1, mitotic cyclins). Together 0.27% of the abundance signal. The lists treat a symptom; a condition-matched abundance set would replace the last three.

### 2.2 Nucleic acids and fixed counts

| Component | Copies | Basis |
|---|---|---|
| 80S ribosomes | 211,500 | Ribosomal protein abundances (median over the r-protein families, with Asc1 and Stm1, section 6.5), less the r-protein content of the pre-ribosomes |
| Pre-ribosomes | 19,500 | Abundances of stage-specific assembly factors (section 5.2) |
| tRNA | 2.31 × 10⁶ | 3.3 × 10⁶ per average cell (assumption, literature order of magnitude) × 0.7, the size of this cell relative to the average cell: about 11 per ribosome; split over 42 isoacceptors by tRNA gene copies (Percudani *et al.* 1997; section 6.5) |
| mRNA | 18,100 | Per gene from Pelechano *et al.* 2010 (Table S1: molecules per cell compiled from Nagalakshmi *et al.* 2008 and Miura *et al.* 2008; 25,900 per average cell over 5,085 genes) × 0.7; each mRNA is a transcript of a named gene (section 6.5) |
| Nucleosomes | 69,714 | 64,487 measured dyads on the chromosome arms (section 6.3; before the measured map: genome ÷ 165 bp = 73,150), plus 5,227 on the nucleosomal rDNA (section 6.6) |
| Nuclear pore complexes | 93 | 10.8 per µm² of nuclear surface: the density of 117 pores (Winey *et al.* 1997, mean over the cell cycle; 65–182, rising from G1) on the earlier, larger nucleus. An assumption, checked 2026-11-13 against Table 1 of that study: 10 G1 nuclei have 86 ± 16 pores (65–113), 11.7 per µm², so 93 is inside the G1 range |

Because the ribosome count follows the ribosomal-protein abundances, it inherits the uncertainty in the proteome total (section 2.1 for the measurements).

### 2.3 Not included

- **Small molecules and ions:** metabolites, free amino acids, nucleotides, ions and water.
- **Storage carbohydrates:** glycogen and trehalose.
- **Polyphosphate.**
- **Proteins without a structure:** 7,200 copies (build R: 40,300; 33,100 of those were free histone H2B, whose library entries HTB1 and HTB2 were missing until stage 0 rebuilt the library for build S, section 13).
- **Transmembrane proteins annotated only to soluble compartments:** 11,200 copies (see 4.3).
- **Surplus nucleoporins:** 13,000 copies beyond a free pool of 25% of the copies in pores (section 6.1).
- **Proteins this cell does not have:** section 2.1.

Lipids and the cell-wall polysaccharides are drawn in the viewer but generated procedurally (section 7).

## 3. Geometry

### 3.1 Volumes

The cell is 34 µm³ with a 110 nm cell wall, the lower bound of 110–200 nm (Orlean 2012). Its outline and nucleus are real segmented shapes; since build S the vacuole is one near-spherical body with the shape statistics measured by soft X-ray tomography of G1 cells (section 3.3).

**Cell volume (2026-11-12, build R; was 42 µm³, the mean log-phase volume of Jorgensen *et al.* 2002).** With 3.5 pg of protein (section 2.1) and the fractions below, 34 µm³ is the volume at which the cytosol keeps the macromolecular density it had (about 170 mg/mL), so the crowding and the diffusion calibration of section 14.5 stay valid. It lies between the culture mean (42 µm³) and unbudded G1 cells of the same strain and medium by soft X-ray tomography (24.9 µm³; Chen *et al.* 2025). The measurements do not agree on it: the dry-mass density of the cytoplasm (section 10.3) alone would put the cell at 36–38 µm³, and the ribosome concentration by cryo-ET (23 µM; the model has 17 µM) at about 26 µm³. The volume is a choice inside that range, not a measurement.

Organelle volume fractions:
- **Nucleus, vacuole and lipid droplets:** set against segmented cryo-FIB/SEM cells (EMPIAR-11416; Dumoux *et al.* 2023; section 3.3) and the literature:
  - **Nucleus 7% (build R; was 8%).** The nucleus-to-cell ratio is constant across cell sizes at about 7% (fluorescence; Jorgensen *et al.* 2007). Soft X-ray tomography of unbudded G1 BY4741 in YPD gives 5.9 ± 1.3% (Chen *et al.* 2025) and 5.3–6.3% in control strains (Ohsawa *et al.* 2026, a preprint); the two segmented nuclei are 8.0% and 9.7% of their cells (both from cells cut by the edge of the imaged volume, which biases the fraction up); serial-section EM gives 10.5% (Yamaguchi *et al.* 2011). The nucleus is 2.0 µm³ inside the inner membrane (2.9 before). The nucleolus is 20% of the nucleus (Uchida *et al.* 2011).
  - **Vacuole 17.2% (build R; was 9%).** Chen *et al.* 2025 (Mol Biol Cell 36:ar132): 17.2 ± 6.0% in 287 unbudded G1 cells of BY4741 in YPD by soft X-ray tomography, the measurement matched in strain, medium and stage. The five segmented cells give 6–13% (a lower bound: the segmentation misses faint vacuoles), an older soft X-ray study 7% (Uchida *et al.* 2011), fluorescence 25% (Joyner *et al.* 2016).
  - **Lipid droplets 0.3%.** The 14 segmented droplets of the cell whose outline is used (the five cells range 0.3–0.8%, from the depositors' hand-corrected labels); soft X-ray tomography gave 1%.
- **Mitochondria 3% (2026-11-11; was 2%, and 1% before that).** Measurements span 1.7% (serial-section EM; Yamaguchi *et al.* 2011) to 4.77% (3D fluorescence in glucose; Di Bartolomeo *et al.* 2020 PNAS 117:7524; CEN.PK in minimal medium, and diffraction inflates 300 nm tubules); soft X-ray tomography gave 1% (Uchida *et al.* 2011). At 2% the matrix asked for about 270 mg/mL of protein and only 40% fitted. 3% is where the requested matrix protein reaches 1.1× the cytosol's density, the ratio of soft X-ray absorption of mitochondria to cytosol (Chen *et al.* 2025): a fit to that ratio, inside the measured range. In build S, with the rebuilt mitochondria (section 3.3), the placed matrix density is 1.05 times the cytosol's (build R: 0.89).
- **Golgi, endosomes and peroxisomes:** 1% together (assumption).
- **Cortical ER:** covers 48% of the plasma membrane (Wolf *et al.* 2012).

Measured on the final voxelised geometry (Monte Carlo volumes of stage 2, build S; `work/s2_geometry.log`):

| Compartment | Volume (µm³) | Share of cell | Build R (µm³) |
|---|---|---|---|
| Cytosol | 18.38 | 54.1% | 18.42 |
| Cell wall (glucan 3.66, mannan 1.82) | 5.48 | 16.1% | 5.48 |
| Vacuole (1 compartment) | 5.83 | 17.2% | 5.81 (12 compartments) |
| Nucleoplasm | 1.56 | 4.6% | 1.56 |
| Nucleolus | 0.47 | 1.4% | 0.47 |
| Mitochondrial matrix | 0.82 | 2.4% | 0.70 |
| Mitochondrial intermembrane space (incl. crista lumen) | 0.18 | 0.5% | 0.27 |
| ER lumen: cortical ER | 0.45 | 1.3% | 0.45 |
| ER lumen: perinuclear space | 0.29 | 0.9% | 0.29 |
| Lipid droplets (14) | 0.12 | 0.3% | 0.12 |
| Golgi / endosomes | 0.27 | 0.8% | 0.26 |
| Peroxisomes | 0.07 | 0.2% | 0.07 |
| **Total** (incl. 0.03 µm³ of pore channels) | **33.95** | | 33.95 |

The tile measures used for packing differ from these by up to 0.01 µm³ (cytosol 18.37, vacuole 5.84, matrix 0.82 µm³).

### 3.2 Membranes

Membrane mid-surface areas, measured on the final geometry:

| Membrane | Area, build S (µm²) | Build R (µm²) |
|---|---|---|
| Plasma membrane (with furrows and invaginations) | 47.4 | 47.4 |
| ER (cortical sheets) | 46.6 | 46.5 |
| Nuclear envelope (both membranes) | 17.7 | 17.7 |
| Mitochondrial inner membrane (boundary membrane and cristae) | 13.1 | 31.5 |
| Mitochondrial outer membrane | 11.7 | 14.3 |
| Vacuole | 16.0 | 40.1 (12 lobes) |
| Golgi / endosomes | 8.2 | 8.6 |
| Peroxisomes | 1.3 | 1.3 |

The internal membranes total 115 µm² in build S against 160 µm² in build R: one round vacuole has 40% of the membrane of the carved lobes, and the rebuilt mitochondria have short cristae (section 3.3).

**Plasma-membrane invaginations** (`cell3d.py`, primitive `T_SLOT`). The plasma membrane is otherwise a fixed depth under the outer surface; these pockets are added to it by a smooth union (rim radius about 8 nm), and their lumen counts as cell wall.
- **Eisosome furrows (MCC):** 23 in build S (22 in build R; of 30 tried; they keep 300 nm apart, and the smaller cell has less free plasma membrane), about 300 nm long (lengths drawn from 300 ± 50 nm; mean 268 nm in build S, 272 nm in build R) and 50 nm deep (Strádalová *et al.* 2009; Karotki *et al.* 2011: about 30 per cell). The cross-section is a U whose bottom is a half-cylinder with the membrane radius of the native Pil1/Lsp1 coat (10.2 nm midplane, 31 nm coat; in situ 30–60 nm, Bharat *et al.* 2018). The coat is in section 6.8.
- **Endocytic invaginations:** one per actin patch (29). Each is a tube with a rounded tip, 40 nm wide at the membrane midplane (the 6,400 nm² mean area of the vesicles they release is a 45 nm sphere; Kukulski *et al.* 2012, 2016), 0–100 nm deep, drawn uniformly because the patches are at random stages and invaginations reach about 100 nm at scission. Measured tips are narrower and invaginations deeper than this (section 11).
- **Sites:** where there is no cortical ER, the membrane normal is within 8° of radial, and the nucleus and vacuole are clear. Furrows run in random directions in the membrane plane.

**Bud scar** (`cell3d.py`, primitive `T_TORUS`). The model cell is taken to be a mother that has budded once (assumption), so it carries one bud scar: the chitin ring left on the mother's wall by cytokinesis.
- **Shape:** a ring 1.0 µm across (AFM: about 1 µm; Pillet *et al.* 2014), modelled as a torus whose tube (radius 50 nm) is centred on the outer surface, so the rim stands 50 nm proud of it (the height is an assumption). The torus is lowered by the surface's sag under the ring, and it is joined to the outer surface only by a smooth union (30 nm), so the plasma membrane and the wall layers under it keep their depths. The AFM profile's convex centre and the depression inside the rim are not modelled.
- **Chitin:** the rim is filled with procedural chitin microfibrils (section 7).
- **Birth scar:** not modelled. A daughter's birth scar has no chitin rim and is only a faint feature of the surface.

### 3.3 Organelle shapes

**A composite of real cells.** The cell outline, the nucleus and the lipid droplets are real shapes from a cryo-plasma FIB/SEM volume of plunge-frozen log-phase yeast (EMPIAR-11416; Dumoux *et al.* 2023), segmented into meshes by a separate pipeline (Cellpose plus manual curation). Up to build R the vacuoles were too; since build S the vacuole is one ellipsoid placed near them (below):
- **Outline and lipid droplets:** cell 08, the only segmented cell not cut by the edge of the imaged volume (44.2 µm³, scaled by 0.916 in length to 34 µm³; 0.983 to 42 µm³ before build R). Its segmented vacuoles decide between equally good positions of the vacuole.
- **Nucleus:** cell 05's. Cell 08's nucleus was not segmented; cell 05's is one of two real ones, and the other (cell 04's) has its short axis along the milling direction of a cut cell. It is scaled to 7% of the cell volume and placed where the segmentation pipeline put cell 08's stand-in nucleus (the cell's most open region, where its unsegmented nucleus probably is), oriented along that stand-in's principal axes, and moved 118 nm from there to stay 260 nm below the outer surface.
- **Not used:** the bundle's procedural mitochondria, cortical ER and peroxisomes.

**Cleaning the meshes** (`pipeline/cellshape.py`; the meshes are in `pipeline/cell_meshes.npz`, extracted by `scripts/cell_meshes.py`):
- **Smoothing:** each part is voxelised on a 20 nm grid, closed morphologically (radius 100 nm for the outline, 60 nm for the rest), which fills grooves a curtaining stripe of the FIB data cut into the outline and some vacuoles. It is then Gaussian-smoothed (σ 50 nm outline, 40 nm nucleus, 30 nm vacuoles), which removes the mesh facets and the 50 nm slice steps.
- **One vacuole (build S, 2026-11-13; `cellshape.VAC_MODE = "single"`).** Cryo soft X-ray tomography of unbudded G1 BY4741 in YPD (Chen *et al.* 2025, Supplementary file S8: 490 cells at 30 nm voxels) finds one vacuole object in 488 cells, with a corrected sphericity of 0.98 (median) and a best-fit ellipsoid of elongation 1.10 and flatness 1.07 (medians); the same study gives the 17.2% volume fraction used since build R. The vacuole is therefore an ellipsoid with these axis ratios and 17.2% of the cell volume. Of 81 orientations (80 random), it is put where it needs the least carving to stay 220 nm below the outer surface and 30 nm from the nucleus (ties: nearest to the volume-weighted centroid of cell 08's segmented vacuoles), then carved where it still violates either, grown back to the target volume and smoothed. Result: one compartment of 5.84 µm³ (semi-axes 1.22, 1.11 and 1.04 µm before carving; 1.8% of it carved by the nucleus and the cortex; centre 1.2 µm from the centroid of the segmented vacuoles), 16.0 µm² of membrane, sphericity 0.98. 30 nm voxels cannot tell one vacuole from lobes pressed flat against each other, so the tomography does not exclude lobes; the carved lobes of build R (sphericity 0.39, below) are far from what it measures.
- **The segmented lobes (builds up to R; `VAC_MODE = "lobes"`, kept in the code):** the segmentation loses the tops and bottoms of vacuoles (z extent over in-plane extent 0.55–0.61 across the 107 segmented vacuoles, against 1.02 for the cells), so each lobe of cell 08 was stretched along z to its smaller in-plane extent (stretch 1.0–3.0, median 2.0), smoothed more along z (σ 90 nm), blended 1:1 with its moment-matched ellipsoid, and one streak along the curtaining stripe was dropped. The 15 lobes were scaled to 17.2% of the cell (build R: their segmented size, in an outline shrunk by 0.916), carved 30 nm apart from each other and from the nucleus and 220 nm from the outer surface, and smoothed (σ 20 nm). In build R that gave 12 compartments (5.81 µm³) with 40.1 µm² of membrane and a sphericity of 0.39: at 17% of a 34 µm³ cell the lobes collided and were carved into boxy pieces with slits between them, and one lobe wrapped another (visual check of build R, 2026-11-13).
- **Lipid droplets:** the 14 segmented droplets are kept as volume-equivalent spheres at their positions (median radius 124 nm); 13 that overlap the vacuole or the nucleus are moved to the nearest free space (build R: 8). Spheres also remove the droplets' z-blur (z extent / in-plane extent 1.25).

**Distance fields.** Each cleaned surface is re-meshed (marching cubes) and the exact distance to it is computed at the grid nodes near the surface (closest triangle). Between nodes the distance is interpolated linearly, accurate to about 0.2 nm, so membranes, offsets (the 35 nm mannan layer, the 110 nm wall, the cortical-ER sheets) and protein orientations are as precise as with the analytic shapes. The fields are stored in sparse 8³-voxel bricks within a band around each surface, and the viewer uses the same fields (re-bricked, 7.5 MB compressed). `scripts/s2_geometry.py --sphere` builds the previous analytic cell instead.

Remaining shapes and fixed dimensions:

| Element | Shape and dimensions |
|---|---|
| Nuclear envelope | Perinuclear space 34 nm between the membrane mid-surfaces, following the real nucleus. 93 pores of 45 nm radius at the mid-plane, at least 150 nm apart, spread uniformly over the envelope area outside the nucleus–vacuole junction (section 6.10) and set along the local normal. The nucleolus is a crescent on the side facing the vacuole (20% of the nucleus; below). |
| Cortical ER | Sheets 27 nm beneath the plasma membrane, following its real shape, with a 20 nm lumen. |
| Mitochondria (rebuilt for build S, 2026-11-13) | Curved tubules of capsule segments 110 nm long, started 457–597 nm and kept 437–637 nm below the outer surface, 1–2 separate networks per cell with side branches (build S: 2 networks, 7 chains, 5 branch points, 9 tips; 10.4 µm of tubule). Outer diameter varies along a tubule, 250–414 nm in build S (median 338 nm; 339 ± 5 nm in glucose by live 4Pi microscopy, Egner *et al.* 2002; median 330 nm in cryo-tomograms of BY4741 in YPD, Chang *et al.* 2025). Intermembrane space 14 nm between the membrane mid-surfaces (13–16 nm in those tomograms). 169 short lamellar cristae (16 per µm of tubule) with an 18 nm lumen (rim-tube radius 9 nm, chosen to fit ATP synthase dimers; section 6.7), each reaching 50–200 nm into the tubule (gamma-distributed around 85 nm; never across it) from alternating sides that turn slowly along the tubule. Crista membrane about 2.6 µm², 0.24 times the outer membrane (0.28 measured on the YPD meshes of Chang *et al.* 2025; 0.72 in glycerol). 11 crista-free nucleoid sites (section 6.7). Until build R: a dozen straight tubules of 150 nm radius, an 8 nm intermembrane space and a disc across the tubule every 42 nm (239 cristae; crista membrane 1.4 times the outer membrane, the look of a respiring cell). |
| Other organelles | Golgi/endosome cisternae and peroxisomes are simple primitives. |

**Arrangement.** The outline and lipid droplets keep their arrangement in cell 08 (13 of the 14 droplets moved out of the vacuole's way), and the vacuole sits where it needs the least carving, nearest the centroid of cell 08's vacuoles among such positions (1.2 µm from it). The nucleus sits in that cell's open region. The nucleolus faces the vacuole, and the spindle pole body, opposite the nucleolus, faces open cytoplasm; the mitochondria are tubules along the cortex. The last three are procedural.
- **Nucleolus next to the vacuole (build S, 2026-11-13).** The nucleolus usually lies at the nuclear envelope next to the vacuole (Barbosa *et al.* 2015; Santana-Sosa *et al.* 2023), and the spindle pole body lies opposite the nucleolus (Yang *et al.* 1989). The crescent is a sphere centred beyond the envelope in the direction of the vacuole centroid, sized to 20% of the nucleus. Up to build R it was on the side away from the vacuoles; with one large vacuole that put the spindle pole body against the vacuole, and no cytoplasmic microtubule could be built (round S candidate builds 1 and 2; build R had 3).

**Mitochondrial sources.** Measurements of glucose-grown cells collected for the rebuild are in `pipeline/data/validation/literature/mitochondria/` (README of 2026-11-13). Crista shape: "most cristae were rather sheet-like", "almost none of them spanned the mitochondria", mostly on two opposing sides of the tubule (FIB-SEM, galactose; Stoldt *et al.* 2019); glucose cells have relatively few, poorly defined cristae (Klecker & Westermann 2021). The crista density and reach were measured on the YPD meshes of EMPIAR-12534 (Chang *et al.* 2025): 24 crista patches per µm² of outer membrane, median depth 85 nm, maximum about 200 nm. Number of separate networks: median 2 per haploid cell (Seel *et al.* 2023).

**Shape statistics** of all five segmented cells (`scripts/organelle_stats.py` → `pipeline/organelle_stats.json`) were used to choose the composite's parts and to check it. The cell 08 outline is mildly non-spherical: its volume-moment semi-axes are 2.41, 2.22 and 2.00 µm.

**Known geometric deviations from EM.** Serial-section EM (Yamaguchi *et al.* 2011) gives a nucleus of 10.5% of cell volume; this model has 7%. Mitochondria are 1.7% there and 2.9% here (matrix plus intermembrane space, 0.99 µm³ in build S; 0.98 in build R).

## 4. Localization

### 4.1 Source annotations

**UniProt subcellular-location text** is mapped to compartments by ordered regular-expression rules. Examples:
- "nucleolus"
- "mitochondrion inner membrane"
- "endoplasmic reticulum"
- "cell cortex" → plasma membrane

Annotations are filtered by their evidence:
- **GFP data preferred:** where a location token is backed by the genome-wide GFP study of Huh *et al.* (2003), the tokens kept are the GFP-backed ones plus those with experimental evidence from a dedicated study. A paper cited as location evidence in 40 or more UniProt entries counts as high-throughput. Tokens are kept in UniProt's order, giving the bulk steady-state localization without proteomics-only calls. Examples: GPD1 keeps its cytoplasm, so it is two-thirds cytosolic; PNC1 is mostly cytosolic and nuclear with 14% in peroxisomes.
- **Proteomics-only calls dropped:** tokens supported only by two high-throughput mitochondrial proteomics papers (Sickmann 2003, Reinders 2006) are discarded. Those studies assign many abundant cytosolic enzymes to mitochondria.

**Fallbacks:**
- **GO cellular component:** used when UniProt has no location. Complex terms are skipped, except pores, ribosomes and proteasomes. If GO lists the cytosol, only cytosolic and nuclear terms are kept.
- **Cytosol:** the default for proteins with no annotation at all (625 proteins).

**Curated fix:** Pil1 goes to the plasma membrane. Its UniProt location is only "lipid droplet", from an organelle proteomics screen, while GO and the eisosome structures place it in the eisosome coat (section 6.8).

### 4.2 Weights across compartments

**Measured shares (CYCLoPs).** For 2,930 proteins, the share of each compartment comes from imaging: the CYCLoPs LOC-scores of the three wild-type screens of the yeast GFP collection (Koh *et al.* 2015 G3 5:1223; fraction of cells showing the protein in each of 16 compartments, averaged over the screens).
- **Groups:** the scores set the share of each compartment group (nucleocytoplasm, ER with nuclear envelope, mitochondria, Golgi/endosomes, vacuole, peroxisomes, plasma membrane). Within a group the annotation-based rule below decides, because imaging is less reliable there: the classifier calls many abundant cytosolic enzymes "Nucleus" (ERG10, TKL1), and LOC-scores are fractions of cells, not concentrations.
- **Periphery:** "cell periphery" and bud calls count as plasma membrane for membrane proteins and proteins annotated there, and as cytosol otherwise (abundant cytosolic proteins look peripheral because the vacuole presses the cytoplasm against the cortex). Cortical patches count as plasma membrane.
- **Not used** when the scores and the annotation share no group (534 proteins; C-terminal GFP mislocalises some secreted, lumenal and imported proteins, e.g. Hsp60, cytochrome c), for membrane proteins without a membrane call, for secreted cell-wall proteins and for ribosomal proteins (GFP-tagged r-proteins accumulate in the nucleus unassembled).
- **New groups:** a group the annotation does not list is added only from a LOC-score of 0.5, never the nucleocytoplasm (where fusions that fail to import end up), and mitochondria only for proteins of the high-confidence mitochondrial proteome.
- **Nucleus and cytosol (2026-11-11).** Within the nucleocytoplasm group the annotation rule used to split the two by volume (equal concentration), so proteins that imaging shows enriched in the nucleus had 92% of their copies in the cytosol. The split now follows the CYCLoPs Nucleus − Cytoplasm LOC-score difference *d*: the nuclear/cytosolic concentration ratio is exp(5*d*) (12× at *d* = 0.5, 1/12 at *d* = −0.5). The factor 5 was fitted on the dev half of an independent comparison (section 9.1) and raised agreement on the held-out half with N-terminal SWAT calls from 64.0% to 66.8% (paired difference +2.9 points, 95% CI 2.1–3.7; 49 proteins fixed, none broken).
- **Effect:** 1.25 M protein copies change compartment (counted in the build before the proteome was rescaled to 5.0 pg, when it held 42 M copies; not recounted since). Examples: proteasome RPN1 and WTM1 are nuclear, COPI coat proteins Golgi, actin-patch proteins (ABP1, SAC6) cortical, profilin cytosolic.
- **Complexes inherit the split (2026-11-11):** the curated 26S proteasome took a fixed 50/50 nucleus/cytosol split; it now takes the copy-weighted split of its subunits' measured weights (81% of RPN1-containing copies placed in the nucleus). The proteasome is assembled in, and stays mostly in, the nucleus of growing cells (Laporte *et al.* 2008 J Cell Biol 181:737).

**Annotation-based rule** (proteins without usable CYCLoPs data, and the split within groups):
- **Rank order:** a protein's compartments get weights 1, ½, ¼, … in annotation order, with the first-listed location taken as primary.
- **Cytosol and nucleoplasm:** their shares are pooled and re-split by volume, i.e. equal concentration across the nuclear pore. Since build S (2026-11-13) the volumes are those of the built cell (nucleoplasm 1.65 and cytosol 18.44 µm³, against 1.56 and 18.42 as built; the rest is pore channels and the smoothing of the nuclear surface). Until build R the inventory used 1.90 µm³ of nucleoplasm, which includes the lumen of the nuclear envelope, so "equal concentration" came out 1.24 times enriched in the nucleus. The fitted nucleus rule below uses the same volumes.
- **Secondary mitochondrial locations:** a mitochondrial location listed after a cytosolic or nuclear primary one (e.g. SOD1, ADK1 in the intermembrane space) gets 0.1× weight.
- **High-confidence mitochondrial proteins:** for proteins in the high-confidence mitochondrial proteome whose first location is mitochondrial, any other location gets 0.1× weight. This is the reverse of the rule for cytosolic proteins with a minor mitochondrial pool. Example: PET9, the major ADP/ATP carrier, is 5% on peroxisomes rather than a third.

**Mitochondria.**
- **Sub-compartments, measured:** for 810 proteins the sub-compartment (outer membrane, inner membrane, intermembrane space, matrix) comes from SILAC-based submitochondrial fractionation (Vögtle *et al.* 2017 Nat Commun 8:290); it replaces annotated or generic locations for 252 inventory proteins. For proteins Vögtle does not cover, the submitochondrial classes of Morgenstern *et al.* 2017 (Cell Rep 19:2836; Matrix/IM, IMS/IM, OM) replace a generic "mitochondrion" (soluble proteins to the matrix or intermembrane space, membrane proteins to the inner membrane).
- **False mitochondrial assignments removed:** proteins outside Morgenstern's high-confidence mitochondrial proteome (901 proteins) that are depleted from pure mitochondria (their Class 3; 28 proteins, e.g. Glt1, Lys21, COPI coat) lose their mitochondrial share.
- **Remaining generic locations** (now 0.02 M copies, was 0.52 M): transmembrane proteins split 75/25 between the inner and outer membranes, soluble proteins 90/10 between matrix and intermembrane space.
- **Copy numbers:** Morgenstern's absolute copy numbers in glucose correlate with PaxDb (Spearman 0.86 for 759 mitochondrial proteins) and give the same mitochondrial share of protein mass (27% vs 25% over the 4,137 shared proteins); their total is 1.5× higher (63 M molecules), a different proteome scale. PaxDb is kept, so that the whole proteome stays on one scale.

### 4.3 Membrane proteins

- **Definition:** proteins with UniProt TRANSMEM segments. β-barrels (TOM40, SAM50, MDM10, MDM34) have no TRANSMEM annotation and are added explicitly.
- **Assigned to soluble compartments only:** these take their membrane from GO, reading "cell periphery" as plasma membrane. Copies still without a membrane are dropped (11,200).
- **ER membrane proteins:** 15% go to the nuclear envelope, in proportion to area, because the outer nuclear membrane is continuous with the ER.
- **Lipid-droplet transmembrane proteins:** placed in the ER membrane.

### 4.4 Peripheral proteins

Soluble proteins annotated to a membrane are placed in thin "peripheral" shells next to it:
- **Plasma membrane, nuclear envelope and lipid-droplet surface:** on the cytosolic side.
- **ER (2026-11-12):** a soluble ER-annotated protein goes into the lumen only with a signal peptide or a lumen annotation (Kar2, Pdi1, Ero1 ...); the rest sit in a 25 nm shell on the cytosolic face (COPII coat and Sar1, Cdc48, Ypt1, Get3, SRP subunits, Scp160, Acc1). Before this fix 181,000 copies of 72 such proteins were packed inside the ER lumen, the same error the rule below had removed from the other organelles.
- **Inner-membrane-associated mitochondrial proteins:** in the matrix.
- **Vacuole, Golgi/endosomes and peroxisomes (2026-11-11):** a soluble protein annotated to one of these organelles goes into its lumen only if it can get there:
  - vacuole and Golgi: proteins with a signal peptide (UniProt), annotated to the lumen, or Cvt cargo (Ape1, Ams1, Ape4, which reach the vacuole by a selective-autophagy route; Lynch-Day & Klionsky 2010; Yuga *et al.* 2011);
  - peroxisomes: the matrix stays the default (matrix proteins are imported without signal peptides); only proteins UniProt places on the cytoplasmic side (Pex1, Pex6, Pex4) move out;
  - all others go into a 20 nm shell on the cytosolic face of the organelle's membrane: free V-ATPase V1 subunits, the EGO complex, Tco89 and part of TORC1, Vac8, COPI, clathrin, COPII, Arf1/2 and Rab GTPases. Before this rule they were packed inside the organelles (680,000 copies; the Golgi lumen held mostly coat proteins).
  - Complexes follow the majority of their subunits.
- **No cell-wall share without a signal peptide (2026-11-11):** soluble proteins with a cell-wall annotation but no signal peptide (e.g. Ssa1/2 from surface proteomics, Fus3, Kss1) lose the wall share (132,000 Ssa1/2 copies).
- **GPI-anchored proteins (2026-11-12):** proteins with a GPI anchor (UniProt lipidation; Gas1, Ecm33, Yps1, Plb1 and others) are soluble by the transmembrane definition and annotated to the cell membrane, so the peripheral rule had put 64,000 copies on the cytosolic face of the plasma membrane. Their plasma-membrane and cytosol shares now go to the wall region; the model has no separate layer between the bilayer and the wall.
- **"Spindle pole" (2026-11-12):** the UniProt term used to send any protein to the nuclear envelope (43,000 calmodulin were spread over the whole envelope). It now does so only for the 16 components of the spindle pole body itself, which section 6.2 then draws; for all other proteins the term is ignored and their other annotations decide.
- **Capacity of the shells (build S, 2026-11-13; `wc3d.shell_capacity`):** a membrane-face shell is asked to hold at most the sphere volume fraction that the cytosol is asked to hold (0.41 in build S); beyond that, every species of the shell is reduced by the same factor and the rest of its copies go to the cytosol as the protein's soluble pool. In build S this affects two shells: the Golgi/endosome face asked for 0.67, so 0.61 of every species stays and 141,000 copies go to the cytosol; the nuclear-envelope face asked for 0.51, so 0.80 stays and 8,000 copies go to the cytosol. A peripheral protein exchanges between its membrane and a soluble pool, and the annotation gives no bound fraction; much of the Arf1 and coatomer pool is probably cytosolic at any moment, which the imaging-based Golgi weight cannot separate.
- **Result (build S):** the vacuole face holds 116,200 copies (99.2% of its request), the Golgi face 211,600 (97.7% of the capped request of 216,700). Build R, without the cap: the Golgi face was asked for 358,000, of which 84% fit in a 20 nm shell around 8.6 µm² of membrane, 17% of the placed copies ended outside the shell, and the vacuole face held 117,000 (99.9%).
- **In the viewer (fixed 2026-11-13):** the default region mask had 24 bits, so the vacuole, Golgi/endosome and peroxisome faces (regions 24–26, added later) were packed but not drawn, about 300,000 proteins in build R; this, not the packing, was the empty-looking layer around Golgi and endosome membranes in the published viewer of build R.

## 5. Structures

### 5.1 Bead representation and mass

- **Beads:** the Cα for amino acids, and P (or C4′) for nucleotides.
- **Mass:** for densities, a residue counts as 110 Da and a nucleotide as 330 Da.
- **What is missing from the beads:** residues missing from a structure (unresolved in the PDB entry, or trimmed from an AlphaFold model) are absent from the model and from its mass.

### 5.2 Choice of structure (in order of preference)

1. **Curated PDB complexes:**

   | Complex | PDB |
   |---|---|
   | 80S ribosome (with Asc1 and Stm1, section 6.5) | 4V88 |
   | Nuclear pore complex | 7N9F |
   | Nucleosome | 1ID3 |
   | 26S proteasome | 4CR2 |
   | Fatty acid synthase | 2UV8 |
   | CCT | 4V94 |
   | RNA polymerase I | 4C2M |
   | RNA polymerase II | 1WCM |
   | RNA polymerase III | 5FJ8 |
   | Pma1 hexamer | 7VH5 |
   | V-ATPase | 6O7V |

2. **Pre-ribosomes** (ribosome assembly intermediates; section 6.6), one yeast cryo-EM structure per stage, with all their RNA:

   | Stage | PDB | Compartment | Particles |
   |---|---|---|---|
   | 90S pre-ribosome (SSU processome) | 5WLC | nucleolus | 2,699 |
   | Nucleolar pre-60S (state E1) | 7R7A | nucleolus | 6,176 |
   | Nucleoplasmic pre-60S (Rix1–Rea1) | 6YLH | nucleoplasm | 1,886 |
   | Cytoplasmic pre-60S (Nmd3, Lsg1, Tif6) | 5T62 | cytosol | 6,126 |
   | Pre-40S | 6FAI | cytosol (85%), nucleoplasm (15%) | 2,654 |

   - **Copy number:** set by the stage's own assembly factors, as for complexes (below): the median over factors found only in that stage's structure, capped by all of its factors. The particles consume their factors before the Complex Portal step, so subcomplexes that are mostly particle-bound (UTP-B, PeBoW, Noc2–Noc3, Nog1–Tif6) are no longer placed free.
   - **Ribosomal proteins** in the particles are not taken from the pools; the mature 80S count drops by the particles' r-protein content instead (8,292 ribosome equivalents), so r-protein mass is conserved.
   - **Total:** 19,500 particles, 8.9% of ribosomes. At a production rate of about 2,000 ribosomes per minute for an average cell (Warner 1999), scaled to this cell, this is 2–4 minutes per stage.
   - **Placed (build S):** the nucleolus cannot hold its request, and its request is thinned by 0.907 before packing (section 8.2), so 4,801 of the 6,176 nucleolar pre-60S and 2,160 of the 2,699 90S particles are placed (329 of those on nascent transcripts); the other stages fit 99–100%. Build R, without thinning: 5,079 and 2,565 (349 on transcripts); the other stages 98–99%.
3. **Complex Portal complexes** (118 added) with a *S. cerevisiae* RCSB assembly. Rules:
   - **Acceptance:** the assembly must cover at least 80% of the subunit mass, and at least 60% of the residues must be modelled.
   - **RNA chains** are kept only where they belong to the complex: snoRNAs, snRNAs, telomerase RNA, and mitoribosomal rRNA.
   - **Copy number:** the median subunit pool, capped at 2× the limiting subunit. Complexes with four or more subunit types are capped at 2× the 25th-percentile subunit.
   - **Mitochondrially encoded subunits** (COB, COX1–3, ATP6, ATP8, OLI1, VAR1) are ignored when setting copy numbers, because mass spectrometry badly under-detects them.
   - **Subunit bookkeeping:** subunits are consumed from the monomer pool, so total protein copies are conserved.
   - **Location (2026-11-11):** the mean of the subunits' compartment weights, each subunit weighted by 1 / (number of Complex Portal complexes it belongs to): a subunit shared between complexes (Tor2 and Lst8 in TORC1 and TORC2, the Lsm ring of U6 snRNP and the cytoplasmic Lsm1-7) reports a mixture of their locations. Five of 118 complexes moved by more than 10% (U5 and U6 snRNPs more nuclear, CORVET towards endosomes, the G protein and Dnf1–Lem3 towards the plasma membrane).
4. **Homo-oligomers:**
   - **PDB:** 134 homo-oligomeric assemblies (112 found from the UniProt subunit text, 22 curated earlier).
   - **Curated from a PDB audit (2026-11-12; `assemblies.CURATED_HOMO`):** abundant proteins whose UniProt text does not name the oligomer and that were therefore monomers: Cdc48 hexamer (8U9P, substrate peptide left out), Tsa1 decamer (3SBC; a C47S mutant; the decamer exchanges with dimers), Gln1 decamer (3FKY; first 18 residues missing), Mmf1 and Hmf1 trimers (3QUW, 1JD1), Pst2 tetramer (5MP4), Ynk1 hexamer (3B54, an AlphaFold dimer before). Sah1 and Gdh1, also oligomers, have no *S. cerevisiae* entry in the PDB and stay monomers.
   - **AlphaFold:** 443 AlphaFold-predicted homodimers with ipTM ≥ 0.5.
5. **Heterodimers:** 476 AlphaFold-predicted heterodimers, from 1,250 models (build R: 475). Requirements:
   - ipTM ≥ 0.5.
   - The two partners share a compartment. A membrane protein may pair with a soluble partner on its cytosolic or matrix side.
   - Not both partners are subunits of an already-modelled complex.
   - Neither partner forms a larger homo-oligomer.
   - A pair with a histone is a native histone pair (H2A–H2B, H3–H4, H2A.Z–H2B, Cse4–H4): histone folds pair promiscuously in AlphaFold (H4–H2B, H2A.Z–TAF9 and others were accepted before 2026-11-12).

   **Copy allocation:**
   - Each protein offers its free copies to its accepted partners in proportion to the partners' copies, and a pair forms the smaller of the two offers (2026-11-12). Before, copies were split evenly over the partners: the rare paralog Btt1 took half of Egd2, and the obligate NAC heterodimer Egd1–Egd2 formed 61,000 dimers beside 123,000 free subunits (in the 5 pg cell). Now it forms 85,500 dimers with 1,550 free subunits.
   - Pairs supported by Complex Portal or a direct IntAct interaction bind all available copies; the others bind 50%. This is an assumption: without affinity data, stable and transient interactions cannot be told apart.
6. **AlphaFold monomers** for everything else. Residues with pLDDT < 50 are trimmed, which removes 3.4% of placed protein mass overall. The loss is concentrated in nuclear and nucleolar proteins: 15–16% of mass there, measured on an earlier build.

**Result by placed mass (build S):** PDB 52.0%, AlphaFold complexes 15.4%, AlphaFold monomers 30.1%, tRNA models predicted with Boltz-2 2.6% (build R: 52.4, 15.0, 30.1 and 2.6%).

**A systematic audit of the PDB (2026-11-12)** screened all 7,067 *S. cerevisiae* entries against the library. Of 869 distinct multi-protein subunit sets, 496 have at least 80% of their genes inside one PDB type of the model. What it found and fixed is in this section and in 6.5; the open candidates (a translating ribosome with tRNAs, eEF3, Ssb and eIF5A on the ribosome, 43S/48S complexes, remodelers on nucleosomes, eEF1A–eEF1Bα and others) are listed in the project's `pdb_audit_2026-11-12/README.md`.

### 5.3 Membrane orientation

- **Method:** UniProt TRANSMEM residues are mapped onto the structure, through RCSB's UniProt alignment for PDB chains. A search finds the membrane normal and offset that put the most TM residues, and the fewest other residues, inside a ±1.6 nm slab.
- **Side:** from UniProt topology domains: the N side (cytoplasm; matrix for the inner mitochondrial membrane; nucleoplasm for the nuclear envelope) faces +z, the P side (lumen, intermembrane space, vacuole interior, outside) −z. Without topology domains, the side carrying more soluble mass is taken as cytoplasmic.
  - **Fix (2026-11):** single-chain AlphaFold models had their side set only from "Cytoplasmic" domains, so 66 proteins annotated with other sides (mostly inner-membrane proteins: "Mitochondrial matrix" / "Mitochondrial intermembrane") were oriented by chance; 39 of them, e.g. TIM50, TIM21 and Mgm1, had their intermembrane-space domains in the matrix. They now use the same side rule as everything else.
- **β-barrels:** oriented by their β-strands, found from Cα geometry.
- **Library entries older than the side rule (found 2026-11-13):** when stage 0 rebuilt the structure library from the pipeline (section 13), 82 transmembrane proteins came out rotated by 180° against the library on disk, whose entries predated the side rule above; 61 of them were packed in build R (42,800 copies; Kre2, Ssh1, Ato3, Erv29, Ctr1 and others) upside down. Build S uses the rebuilt entries.
- **Check:** the method reproduces hand-oriented Pma1 and V-ATPase (z-coordinate correlation 1.000 and 0.9996).
- **Result:** in the packed cell, the median distance of TM residues from the membrane mid-surface is 0.72 nm, and 97.4% lie within 2 nm (build S; build R 0.73 nm and 97.3%).

## 6. Supramolecular organisation built explicitly

These structures are built before random packing. They then act as fixed obstacles for everything packed afterwards.

### 6.1 Nuclear pores and FG mesh

- **Pores:** 93 copies of 7N9F (the in-situ yeast pore), with the pore axis along the local nuclear-envelope normal.
- **None in the nucleus–vacuole junction (build S, 2026-11-13):** pores are absent from the envelope where it is closely apposed to the vacuole (Severs *et al.* 1976 and Pan *et al.* 2000, as cited by Kohler & Büttner 2021 and Otto & Thumm 2020). No pore centre is placed where the vacuole is closer than 92 nm (half the perinuclear space, the 45 nm junction contact distance of section 6.10 and a 30 nm margin; a pore's membrane torus has a 45 nm radius). Up to build R pores were spread over the whole envelope, including the junction.
- **Nucleoporin bookkeeping (2026-11-12):** each pore takes the measured copies of every nucleoporin: 8 to 48 per pore, about 550 proteins and 52 MDa in all (Kim *et al.* 2018 Nature 555:475, Supplementary Table 2A). Until then the subunit list came from Complex Portal, which gives no stoichiometry, so one copy of each nucleoporin per pore was booked and about 100,000 nucleoporin copies were placed as free proteins, mostly at the nuclear envelope.
  - **Free pool:** PaxDb still holds more scaffold and FG nucleoporins than the pores need (median 1.36 times at 3.5 pg). These exist almost only in pores, so their free pool is capped at 25% of the copies in pores (an assumption) and the surplus, 13,000 copies, is not placed. Sec13, Seh1, Dyn2 and the basket proteins Mlp1/2, which have other homes, are not capped.
- **FG mesh:** the disordered FG-repeat domains of the FG-Nups modelled in 7N9F are threaded as self-avoiding residue chains (0.38 nm per residue, low persistence) from the first modelled residue of each chain:

  | Nup | Residues | Copies per pore |
  |---|---|---|
  | NSP1 | 1–636 | 48 |
  | NUP57 | 1–286 | 32 |
  | NUP49 | 1–269 | 32 |
  | NUP159 | 388–1102 | 16 |

  - **Confinement:** the chains stay within 60 nm of their anchor, in the pore channel, cytoplasm or nucleoplasm.
  - **Result:** 4.9 of 5.5 million residues placed (89%).
  - **Missing:** FG-Nups absent from 7N9F (NUP100, NUP116, NUP145N, NUP42, NUP1, NUP60, NUP2), so the nuclear basket side is sparse.
  - **Caveat:** the chain conformations are excluded-volume random walks, not a physical model of FG cohesion.
- **Transport receptors in the channel (2026-11-12; `pipeline/pore_receptors.py`):** 20 receptors per pore (the Jurkat cell model's number; unchecked for yeast), drawn from the free pools of Kap95, Srp1, Kap123, Pse1, the Crm1 export complex, Mex67–Mtr2 and Ntf2 in proportion to their copies, are placed in the channel (within 22 nm of the pore axis and 25 nm of the envelope midplane) before the FG domains are threaded around them. Before, the packer never put a receptor inside a pore. All 1,860 fit (builds R and S).

### 6.2 Microtubules and the spindle pole body

- **Number:** 20 microtubules, made of 12,324 αβ-tubulin dimers (53% of the tubulin pool; build R: 12,831, 55%):
  - 15 kinetochore microtubules, 262 nm mean length (16 are asked for; the 16th did not fit at the plaque in builds R and S; its kinetochore sits 40 nm from the plaque in build R, not re-measured);
  - 2 other nuclear microtubules (4 asked for, about 400 nm): mean length 208 nm in build S (357 nm in build R);
  - 3 cytoplasmic microtubules, 1.2 µm (mean 1,173 nm).
- **Nucleation sites (build S, 2026-11-13):** up to 400 minus-end sites on the 90 nm plaque are tried per microtubule (80 before), since the rings nearly fill the plaque. With the nucleolus moved to the vacuole side (section 3.3) the spindle pole body faces open cytoplasm again, and the 3 cytoplasmic microtubules are built (round S candidate builds 1 and 2, with the spindle pole body against the vacuole, had none).
- **Lattice:** an ideal 13-protofilament, 3-start lattice with a seam. The lateral geometry comes from 8QV3 (native yeast γ-TuRC-capped microtubule, cryo-ET). The dimer is 5W3F, with an 8.3 nm axial repeat.
- **Caps:** each minus end has a γ-tubulin ring of 7 γ-TuSCs from 8QV3.
- **Growth:** from 90 nm plaques at the spindle pole body. Microtubules bend away from membranes, with a minimum radius of curvature of 300 nm.
- **Assumptions:** the lengths (G1-like), and the 70% cap on the use of the tubulin pool.
- **Spindle pole body plaques (2026-11-12; `pipeline/spb.py`):** the SPB is a layered cylinder in the envelope (cryo-EM of isolated SPBs, Bullitt *et al.* 1997 Cell 89:1077; integrative model of the core, Viswanath *et al.* 2017 Mol Biol Cell 28:3298). Its proteins had been spread through the nucleoplasm, so a simulated Spc42-GFP image showed no spot. Now each SPB protein goes into its layer, a disc of the plaque radius (90 nm) on the SPB axis: outer plaque (Spc72, Nud1; 45–75 nm above the outer membrane), Cnm67 (20–50 nm), central plaque (Spc42, Spc29, calmodulin; 3–25 nm), inner plaque (Spc110; 3–60 nm below the inner membrane), half-bridge (Sfi1, Kar1; beside the plaques).
  - **Core stoichiometry (2026-11-12):** Spc42 : Spc29 : Spc110 : Cmd1 : Cnm67 = 2 : 1 : 1 : 1 : 1 (Muller *et al.* 2005 Mol Biol Cell 16:3341, as used by Viswanath *et al.* 2017). The inner plaque is filled first; the Spc110 at the SPB is then that layer plus the 14 per nuclear γ-TuRC, calmodulin takes one per Spc110, and Spc42, Spc29 and Cnm67 are added in units of 2 : 1 : 1 until units stop fitting. Before, each protein was filled from its whole pool, which gave Spc42 : Spc29 : Cnm67 about 1 : 1 : 1 and calmodulin only on the inner-plaque Spc110.
  - **Result (build S):** Spc110 270 (238 in rings, 32 in the inner plaque), calmodulin 270, Spc42 297, Spc29 152, Cnm67 152, Nud1 536, Spc72 184, Sfi1 93, Kar1 38 (build R: Spc110 274, calmodulin 274, Spc42 239, Spc29 124, Cnm67 124, Nud1 413, Spc72 111, Sfi1 74, Kar1 38). The 90 nm disc runs out of room at about 1.1 Spc42 per Spc110 instead of 2 (build R: 0.9), so the core is short of Spc42, Spc29 and Cnm67 relative to Spc110. The absolute number per SPB is not measured here. Layer heights are approximate (unchecked). Simplification: the envelope has no fenestra, so the central plaque sits on its cytoplasmic face. The inventory copies exceed one SPB (573 Spc42), so what does not fit stays free; the simulated Spc42-GFP image shows one spot.

### 6.3 Chromatin

- **Rabl configuration:** each chromosome is two arms, both starting at its centromere, within 20 nm of its kinetochore-microtubule tip at the spindle pole body, and running to a telomere at the nuclear envelope. Arm lengths come from the S288C centromere positions (32 arms; chromosome XII-R is split where the rDNA array sits, so it runs into the nucleolus and out again, and the rDNA itself is in section 6.6).
- **Fibre (2026-11-11):** each arm is a worm-like backbone with the mass density and persistence length measured for interphase chromatin in living yeast (Bystricky *et al.* 2004 PNAS 101:16495: 110–150 bp/nm and 170–220 nm; the middle values 130 bp/nm and 195 nm are used): the backbone advances 1.27 nm per 165 bp nucleosome repeat. Nucleosomes sit at random within 60 nm of the backbone, each within reach (18 nm centre to centre, more across long linkers) of the previous one: dispersed, irregular clutches, not a regular 30 nm fibre, as in-situ cryo-ET of yeast nuclei suggests. The backbone drifts towards the arm's target as its remaining contour runs short, so arms wander freely at first and end at their telomeres; where nucleosomes repeatedly fail to fit, the backbone turns, and as a last resort the chain restarts nearby (a break).
  - **Why:** the earlier walk placed one nucleosome per 12.5 nm step (about 13 bp/nm), so short stretches were far too extended (rms 214 nm at 20 kb against 118–160 nm in vivo) while long stretches were compressed by the nucleus.
  - **Fibre radius (2026-11-11, late):** a first version with a 20 nm radius gave 122–146 nm at 14–20 kb (Bystricky's range) but failed the Hi-C test: 62% of contacts within chromosomes (Hi-C 37–41%), P(s) slope −1.6 (−1.09 to −1.27) and centromere clustering 13× (4–5×). Nine variants were built as pre-passes only (compaction 60–130 bp/nm, persistence 70–195 nm, radius 20–120 nm, and the earlier walk), and the 60 nm radius was chosen: it is the only variant whose R(s) lies inside the chain inferred by Arbona *et al.* 2017 (Genome Biol 18:81; 53–65 bp/nm, persistence 52–85 nm, from imaging and Hi-C together) at both 20 and 100 kb (197 and 411 nm against 163–223 and 393–553 nm), and inside Bystricky's at 100 kb (above it at 20 kb). The Hi-C test is therefore partly tuned: it was looked at in this choice.
  - **Telomere steering (fixed 2026-11-11, late):** the fibre walk turned the backbone in a fully random direction after every failed placement, which erased the drift towards the telomere focus, and the foci of short arms lay beyond their own fibre contour. Arms ended a median 790 nm from their targets (the earlier walk: 37 nm). Now a failed placement turns the backbone with a pull towards the target that grows as the contour runs short, the drift weight is capped at 5 (was 1.5), and each telomere target is moved towards the pole until it lies within 80% of its arm's contour of the start (`TEL_REACH`). In a pre-pass the median miss is 160 nm; the last nucleosome lies up to 60 nm off the backbone. The variant scan above was repeated after the fix (cis 0.47, slope −0.96 for the chosen radius) and the choice held.
- **Measured nucleosome positions (2026-11-11; `pipeline/nucleosome_map.py`):** the fibre is built from a measured nucleosome map: one nucleosome per measured dyad, in genomic order, so every nucleosome carries its genome position (`nuc_bp` in the pre-pass).
  - **Map:** dyads called from the H3 chemical-cleavage dyad maps (H3Q85C; W303 background, YPD) of Chereji *et al.* 2018 (Genome Biol 19:19; three replicates summed; sacCer3, the coordinates of the model's chromosomes). The signal is smoothed (Gaussian, 15 bp), and peaks are taken from the highest down, each excluding 135 bp on both sides, while the signal is above half its genome-wide median; below that the DNA counts as nucleosome-free. That gives 66,245 dyads. Dyads inside promoter nucleosome-depleted regions are then removed: the regions between the upstream stable nucleosome and the +1 nucleosome in Rossi *et al.* 2021 (Nature 592:309, Supplementary Data 1; 10,313 regions; dyads more than 73 bp inside both boundaries). That removes 1,758, leaving 64,487 nucleosomes on the arms (spacing median 169 bp, 4.3% of gaps over 300 bp). The chemical map of Brogaard *et al.* 2012 (Nature 486:496; 67,548 unique nucleosomes) is on the older sacCer2 coordinates (offsets of a few bp that grow along the arms), and Jiang & Pugh 2009 on an intermediate pre-sacCer3 release, so neither is used.
  - **Checks:** a nucleosome-free region 200–100 bp upstream of start codons (0.4 times the flanking density) and a phased array downstream. Confirmed origins sit in wider gaps (median 236 bp between the flanking dyads, against 174 bp at random positions), as expected from origin nucleosome positioning (Eaton *et al.* 2010 Genes Dev 24:748).
  - **Geometry:** the backbone advances by each nucleosome's spacing over 130 bp/nm. The reach from the previous nucleosome grows with the extra linker, 18 nm plus 0.34 nm per bp beyond 165 bp, at most 60 nm (across nucleosome-free regions).
  - **Count:** the inventory now takes the map's count (64,487 plus the rDNA) instead of the genome length over a 165 bp repeat (73,047). That is still 7–13% more than the 57,000–60,000 per cell measured by absolute occupancy (Oberbeckmann *et al.* 2019 Genome Res 29:1996; held out, not tuned to): a population map holds positions that a single cell does not all occupy.
  - **Held-out tests (`validation/nucleosome_test.py`):** the absolute occupancy at 13,807 restriction sites (Oberbeckmann 2019) is 0.75 at sites inside a placed model nucleosome and 0.57 outside (AUROC 0.75; 79% of sites covered against a mean occupancy of 0.72; build M). Of the MNase nucleosomes of Weiner *et al.* 2015 (unstressed), 79% have a model dyad within 40 bp. ATAC-accessible regions (Schep *et al.* 2015) hold 10% fewer model dyads per kb than the genome average; the merged regions are 1.3 kb wide, so the test is coarse.
  - **DNase I at base resolution (2026-11-12; `validation/dnase_test.py`; held out):** 83 million DNase I cuts in log-phase cells (Zhong *et al.* 2016 Genome Res 26:351; W303 in YPD, two replicates, lifted from sacCer2 to sacCer3). Around the model's dyads the cut density dips to 0.62 of the genome mean at the dyad and rises to 1.04 in the linkers 75–95 bp away; the core (±50 bp) is cut 0.66 times as often as the linkers, and 76% of nucleosomes are cut less in their core than in their own linkers. In 50 bp bins, stretches without a model nucleosome are cut 2.3 times the mean and fully covered ones 0.71 times. DNase-hypersensitive sites of BY4741 (Zelin *et al.* 2012 Mol Cell 48:459; 3,260 peaks) have a model occupancy of 0.62 against 0.76 genome-wide.
  - **Fit test:** in test pre-passes all 64,487 nucleosomes are placed (20 and 60 nm radius).
- **Telomeres:** gathered in 7 foci at the nuclear envelope (6–8 seen by immunofluorescence; Gotta *et al.* 1996), 25 nm inside the inner nuclear membrane. A focus's angle from the spindle-pole axis grows with the arm length of its telomeres, 35° for the shortest to 115° for arms of 700 kb or more. That angle range is an assumption; the trend is from Hi-C and imaging (Duan *et al.* 2010; Therizols *et al.* 2010). The foci stay out of the nucleolus.
- **Result (build S, in the nucleus of 2.0 µm³, with the nucleolus and the spindle pole body turned to the other side; section 3.3):** 62,882 of 64,487 nucleosomes placed (the rest evicted by Pol II and origins). All 32 telomeres are within 150 nm of the nuclear envelope (median 94 nm; build R 112 nm), a median 164 nm from their targets. Telomere-to-pole distance follows arm length (Spearman 0.92; build R 0.85). Against the contact maps (`validation/hic_test.py`), as provided and ICE-balanced: 44% / 50% of contacts within chromosomes (Hi-C 37–41% / 39–51%), P(s) slope −0.94 / −1.01 over 10–200 kb (−1.09 to −1.27 / −1.07 to −1.35), centromere clustering 6.5× / 5.9× (4.0–5.3× / 4.0–7.2×), telomere clustering 3.3× balanced (1.7–3.6×). Balanced, the cis share and both clustering measures are inside the range of the two maps; the contact decay is flatter than both. Build R had 51% cis, −1.04, 6.0× and 1.9× (balanced), build Q, in the 2.9 µm³ nucleus, 55%, −1.03, 7.1× and 2.3×. One chromatin realisation moves these measures by about this much (builds O and P differed as much with identical rules), so differences between builds R and S are not evidence for or against the new layout. Rms distances along the chromosome: 191, 289, 398 and 516 nm at 20, 50, 100 and 200 kb (build R: 195, 284, 376, 467; build Q: 193, 294, 403, 536). At 20, 50 and 100 kb that is inside the chain inferred by Arbona *et al.* 2017 (393–553 nm at 100 kb, inside which the fibre radius had been chosen; build R was 4% below it), at 200 kb below it (561–792 nm). Which chromosome pairs touch is not reproduced (balanced Spearman 0.14 with both maps), but the two maps also disagree with each other on that (Spearman −0.04 / 0.15).
- **Linker DNA:** joins consecutive nucleosomes, with a schematic path.
- **In the viewer (2026-11-12; `pipeline/polymer_annot.py`):** hovering or clicking chromatin gives its genome position (R64 coordinates):
  - **Nucleosomes:** the chromosome and the measured dyad, the gene or element there (SGD features: ORFs, RNA genes, centromeres, ARS, telomeres, Ty elements; otherwise the neighbouring genes and the distance upstream of a gene start), and the enriched histone marks.
  - **Linker DNA:** the bp range between the two wrapped cores.
  - **Pol II and the DNA through it:** the midpoint of its two flanking nucleosomes, and the gene there. Near gene ends that can be the neighbour of the transcribed gene. Nascent transcripts carry the same position.
  - **Without positions:** rDNA (nucleosomal or active) and mtDNA are named by kind only.
- **Transcribing Pol II:** per-gene numbers of engaged Pol II from genomic run-on (Pelechano *et al.* 2010 PLoS One 5:e15442, Table S1; 691 per cell over 4,670 genes, about 2% of all Pol II), at uniform positions along each gene (SGD coordinates).
  - **On the arms:** with the measured map, a polymerase replaces the nucleosome nearest its position in bp (earlier builds mapped kb linearly along each arm walk). The nucleosome at a polymerase's position and its two neighbours are removed (nucleosomes are evicted around elongating Pol II), and the polymerase takes the middle one's place. Its DNA axis comes from the elongation complex 1Y1W, superposed on the library 1WCM to 0.01 nm, and lies along the chain. DNA runs from the upstream nucleosome through the cleft to the downstream one.
  - **Nascent RNA:** threaded from the RNA exit, as long as the polymerase has transcribed, 2 nt per bead.
  - **Source pool:** engaged polymerases come out of the free Pol II pool.
- **Licensed replication origins (2026-11-11; `pipeline/origins.py`):** the model cell is G1-like, so each confirmed origin of OriDB (Siow *et al.* 2012 Nucleic Acids Res 40:D682; 410 'Confirmed', sacCer3) gets one MCM2–7 double hexamer (7W8G) and one ORC (7TJH).
  - **Placement:** the double hexamer goes in the gap between the two measured nucleosomes that flank the origin, its long axis (the DNA channel) along the line joining them. ORC sits beside it.
  - **Where it does not fit:** the fibre packs nucleosomes too densely for a 21 nm double hexamer between most flanking pairs. It is then moved up to 25 nm out of the fibre (the origin's free DNA looping out). Where that fails too, the two flanking nucleosomes are evicted (last resort: two on each side).
  - **Result (build S):** 365 of 410 origins licensed, with 356 ORC and 90 nucleosomes evicted (build R: 347, 340, 86; a test pre-pass before the fibre's final settings: 317, 296, 236).
  - **Limits:** one double hexamer per origin is a lower bound, since efficient origins load several (unchecked). The rest of the MCM (1,383 double hexamers in build R; not re-counted) and ORC (506 in build R) pools stay free in the nucleoplasm.
- **Kinetochores (2026-11-11; `pipeline/kinetochores.py`):** one per chromosome on its kinetochore-microtubule plus end: a CCAN bound to centromeric DNA with Cbf1 (8OVW) just beyond the plus end along the microtubule axis, and 6 MIND (Mis12) complexes (9S5N) in a ring around the plus end. Quantitative imaging counts per kinetochore in metaphase two Cse4, one or two copies of the inner proteins, 6–7 MIND, 8 Ndc80 complexes and 16–20 DAM/DASH (Joglekar *et al.* 2006 Nat Cell Biol 8:581, Table 1; checked 2026-11-12 against the paper, so the model's 6 MIND, 8 Ndc80 and 16 DASH are within the measured counts). The chromosome's two arms start beyond the CCAN. All 16 CCAN and 96 MIND fit (builds R and S).
  - **Outer kinetochore (2026-11-12):** per kinetochore 16 copies of each DASH/Dam1 subunit (Joglekar's count) in a ring 16–32 nm from the microtubule axis over its last 5–60 nm, and 8 copies of each Ndc80-complex subunit (Ndc80, Nuf2, Spc24, Spc25; 8 in Table 1 of the same study) in a band from the lattice to just beyond the tip; 3,050 of 3,072 subunits fit (build S; build R 2,908). Both complexes are placed as their AlphaFold subunits (no assembled library types), so the ring is approximate. Not built: the Cse4 nucleosome.
- **Chromatin proteins at their sites (2026-11-11; `pipeline/chromatin_sites.py`):** ChExMix peaks of the ChIP-exo survey of Rossi *et al.* 2021 (about 400 chromatin proteins, YPD, sacCer3; filtered peak files of the paper's repository, 188,489 peaks).
  - **Which proteins:** each protein that is a library type of its own in the nucleoplasm (about 280; subunits of complexes the inventory does not assemble, such as TFIIH, Mediator or SAGA, are placed one by one, so they gather at the same promoters).
  - **How many:** min(peaks, 30% of its nuclear copies) copies go on sites, drawn without replacement with probability proportional to the peak's tag count. The 30% bound share is an assumption (unchecked); the rest stay free.
  - **Geometry:** the DNA point of a site is interpolated between the two measured nucleosomes that flank it. The protein (random orientation) comes within 5 nm of it with one of its spheres, on the side away from the local nucleosome mass, and overlaps nothing else.
  - **Result (build S):** 20,226 of 50,931 wanted copies of 268 proteins fit (40%; build R 19,238 of 51,073 copies of 270 proteins, 38%); the fibre leaves little room next to most linkers. The peak calling excluded the rDNA, tRNA genes and telomeres, so Pol III at tRNAs and Sir/Rap1 at telomeres are under-represented.
- **tRNA genes, telomeres and silent chromatin (2026-11-12; `pipeline/trna_telomeres.py`):** the ChIP-exo peak files exclude tRNA genes and telomeres, so these sites are placed with the same geometry (`chromatin_sites.place`, one copy per site): one Pol III and one copy of each TFIIIC subunit in the library (Tfc1, Tfc7) on each of the 275 nuclear tRNA genes (SGD); 17 Rap1 at each of the 32 chromosome ends (TG1–3 tracts of about 300 bp with a Rap1 site every 18 bp on average: Marcand *et al.* 1997 Science 275:986; 15 before 2026-11-12); per Sir-bound nucleosome, two Sir3 and one each of Sir2 and Sir4 (saturating stoichiometries on reconstituted arrays, Swygert *et al.* 2014 Nat Commun 5:4751 and 2018 PNAS 115:12447, so upper bounds in cells; one each before). The chromatin is crowded, so a share does not fit (build S: 811 of 2,305 copies: Pol III on 102 of 275 tRNA genes, Rap1 89, Sir2 123, Sir3 82, Sir4 53; build R: 842 of 2,317, with Rap1 90, Sir2 135, Sir3 95, Sir4 57).
  - **Which nucleosomes are Sir-bound (2026-11-12; `pipeline/silent_chromatin.py`):** those where Sir3 ChIP-seq is enriched at least 5-fold over input in both replicates (anti-Sir3 ChIP of cycling W303 cells in YPD; Brothers & Rine 2022 eLife 11:e75653), within 30 kb of a chromosome end or in HML or HMR: 234 placed nucleosomes at 28 of the 32 ends in build S (109 within 3 kb, 70 at 3–10 kb, 55 at 10–30 kb; build R 237). Before, every nucleosome within 3 kb of an end and in HML/HMR was taken (529).
    - **Threshold and region, set before the test:** the same laboratory's untagged and *sir3Δ* controls reach about 1 and 4–5 at silent loci, and highly expressed genes are "hyper-ChIPable" (Teytelman *et al.* 2013 PNAS 110:18602): 33 nucleosomes at TDH3 and other such genes pass the threshold and are excluded by the 30 kb limit.
    - **Caveat:** the strains express Sir3 fused to a methyltransferase, not untagged Sir3.
    - **Held-out test (`validation/silent_chromatin_test.py`):** 77% of the Sir-bound nucleosomes are H4K16-hypoacetylated (below the genome's 10th percentile; Weiner *et al.* 2015; Sir2 deacetylates H4K16), against 28% under the earlier rule and 14% of the nucleosomes only the earlier rule took. Sir4 ChIP-seq (Goodnight & Rine 2020 eLife 9:e58910; with an untagged control at 0.9) is enriched 36-fold on them, against 4-fold under the earlier rule (build S; build R 78% and 37-fold).
  - **Check (conflicting evidence):** by FISH, tRNA genes cluster at the nucleolus (Thompson *et al.* 2003 Science 302:1399), an arrangement the model has no rule for; the balanced contact maps show no excess of trans contacts between tRNA-gene bins (1.05× and 1.00×; the model 1.14× in build S, 0.94× in build R), so no rule was added.
  - **Test (2026-11-12, `validation/trna_nucleolus.py`):** FISH finds the active SUP53 tRNA gene overlapping the nucleolus in 52% of 440 cells and the same gene with an inactive promoter in 13% of 715 (Thompson *et al.* 2003); in wild type 76% of the nuclei with a single tRNA-gene cluster have it at the nucleolus (Haeusler *et al.* 2008 Genes Dev 22:2204). In the model (build S), 23% of the 275 tRNA genes lie within 250 nm of the nucleolus (a FISH spot's radius, set before comparing), against 33% of nucleosomes, and SUP53 is 520 nm away: tRNA genes are no closer to the nucleolus than chromatin in general, so the test fails. (Build R: 35% against 43%, SUP53 740 nm; build Q, in the larger nucleus: 33% against 36%.)
- **A relaxed polymer instead of the fibre walk (2026-11-12, tested, not adopted; `pipeline/chromatin_polymer.py`, `validation/chromatin_polymer_scan.py`):** after the Jurkat cell model, every arm was simulated as a coarse polymer (2 kb beads; Arbona's 60 bp/nm and 70 nm persistence; soft excluded volume; confined to the nucleoplasm and clear of the nucleolus, pores and microtubules; centromeres tethered at their kinetochores, telomeres within 60 nm of the envelope, no foci imposed) with Langevin dynamics in OpenMM, and the measured nucleosomes were placed along the relaxed chain (`pack3d.guided_walks`). Its R(s) lies in the Arbona envelope. Scored on ICE-balanced Hi-C, it improved the contact decay (P(s) slope −1.24 to −1.42 against −1.07 and −1.35; the fibre walk −1.03) but made the chromosomes more territorial (cis share 0.63–0.77 against 0.39–0.51; fibre 0.57) and the centromeres more clustered (8–20× against 4.0–7.2×; fibre 6.4×), in all variants tried (longer runs, weaker repulsion, Bystricky's 130 bp/nm and 195 nm, looser telomere tethers; 4–8 replicas each). The fibre walk stays.
- **Not modelled:** inter-chromosomal contacts beyond the Rabl layout, histone variants and marks as structure (H2A.Z, H3K4me3, H3K36me3, H3K79me3, H3K56ac, H3K9ac and H4K16ac of the nearest unstressed MNase nucleosome within 40 bp in Weiner *et al.* 2015 are attached to 83% of nucleosomes in `work/chromatin/nucleosome_genome.parquet`, as annotation only; readers are not placed by them), cohesin and condensin.

### 6.4 Actin

The ACT1 pool (144,600 copies at 3.5 pg; 206,500 in builds J to Q) is split by assumption: 30% in cables, 30% in patches, and the rest monomeric or in profilin, cofilin and Srv2 heterodimers.

**Cables:**
- **Number:** 11 cables, 81 filaments, 43,400 subunits and 120 µm of filament in total (build S; build R 82 filaments).
- **Position:** 90 nm beneath the plasma membrane, along great circles.
- **Bundle:** filaments sit on a hexagonal lattice 12 nm apart (assumed fimbrin bundle spacing).
- **Filament structure:** yeast F-actin 9GO5, 167° and 2.76 nm per subunit.
- **Tropomyosin:** 8,269 TPM1 coiled-coil dimers (build S; 8,206 in build R), one per 5 actins per strand. The path is transferred from the mammalian actin–tropomyosin structure 5JLF, because no yeast structure exists. The chains are wound as an ideal coiled coil.

**Endocytic patches:**
- **Number:** 29 patches of about 1,500 actin each, one at each endocytic invagination (section 3.2; the number of sites follows the plasma-membrane area). The network grows in the cytosol around the tube.
- **Amphiphysin on the tubes:** Rvs161–Rvs167 heterodimers (AlphaFold model) on the tubular part of every invagination deeper than 50 nm, between the neck and the coated tip region (the last 35 nm). Their number rises linearly to 100 at 100 nm depth: Rvs assembles from about 50 nm and peaks at about 100 Rvs167 at scission (Picco *et al.* 2015 eLife 4:e04535). Each dimer lies around the tube with its concave face on the membrane. 422 of the 695 called for fit (build S; build R 386 of 748); random addition jams before the dense lattice that would cover the tube. The library entry of the heterodimer is made by the builder itself since build S (`amphiphysin.model`); before, the only copy was a leftover of 2026-11-07 that no stage wrote, found when a build started from an empty library.
- **Endocytic coat:** Sla2, Sla1 and End3 on the tip of each invagination (the last 35 nm; at flat sites, the membrane patch under the actin), with numbers per site from Picco *et al.* 2015. Sla2 is about 40 throughout. Sla1 (about 90) and End3 (about 60) are present at the start and disassemble during invagination, here linearly to none at 100 nm depth (assumption). Orientations are random, since there is no structure of the coat. Each molecule sits as close to the membrane as it fits.
- **Not built:** clathrin, Pan1, the epsins Ent1/2. Picco *et al.* also count about 3,000 actin per site at peak, twice the 1,500 used here.
- **Branching:** dendritic networks branch at 73° from Arp2/3 complexes (2,998 branches in build S; 3,071 in build R). The branch-junction geometry is 8UXW. The yeast Arp2/3 subunits are AlphaFold models superposed on the fission-yeast template, 32–74% identity.
- **Capping:** capping protein (from 8F8Q) caps filaments after 20 ± 6 subunits. Filaments that end against an obstacle stay uncapped.
- **Orientation:** barbed ends face the membrane.
- **Collision checks:** residue by residue inside each patch. Some overlaps remain inside the built structures: 0.25–0.5% of their beads lie within 0.3 nm of a bead of another built actin component (tropomyosin on its filament, filaments within a bundle, Arp2/3 at its branch), in every build (81 of 17,154 sampled beads in build S, 101 of 36,786 in build R). The overlap check of section 9 reports them separately.

Not built: fimbrin cross-links, septins, formins.

### 6.5 Ribosomes and mRNA

- **Gene-specific mRNAs (2026-11-11; `pipeline/mrna_genes.py`):** every mRNA is a transcript of a named gene, drawn from per-gene molecules per cell (Pelechano *et al.* 2010, Table S1, a population mean of 25,900; times 0.7 for the G1 cell like the proteome, section 2.1: 18,100 mRNAs of about 3,600 genes in one draw). Its length is 5′ UTR + CDS + 3′ UTR + poly(A), with the CDS the summed coding segments of the SGD ORF (introns left out; see the fix below); UTRs and poly(A) are per gene where measured: TIF-seq UTRs in YPD for 4,298 genes (Pelechano *et al.* 2013 Nature 497:127) and PAL-seq total-RNA tails for 3,454 genes (Subtelny *et al.* 2014 Nature 508:66); other genes take 50, 120 and 33 nt (near the data medians). Ribosome load is weighted by per-gene translational efficiency from the same ribosome profiling (Subtelny *et al.* 2014; clipped to 0.1–10 times the median). Tables and their derivation: `pipeline/data/mrna/README.md`.
- **Polysomes:** 70% of cytosolic ribosomes are on mRNAs (assumption). Each mRNA carries Poisson(CDS × TE / s) ribosomes (TE: translational efficiency, above), with one global s chosen so that the total matches (167 nt of CDS per ribosome at TE 1 in builds R and S); mRNAs drawing none are free. This is close to the average density from polysome profiling quoted in section 14.3 (one ribosome per 156 nt; Arava *et al.* 2003), which was not used; per-gene densities from that study can serve as a test.
  - **Chain geometry:** the mRNA between two ribosomes is a flexible single strand of (spacing − 28 nt footprint) nucleotides; the distance from one ribosome's exit to the next one's entry is drawn from the end-to-end distribution of a Gaussian chain of that contour (persistence 1.5 nm, 0.6 nm per nt; rms about 14 nm), within a 70° cone around the exit direction. Earlier builds used a fixed 3–12 nm linker in a 40° cone and gave up after 60 failed tries, which ended many chains early; the polysome builder now keeps 95% of the requested ribosomes in an otherwise empty cell (60% before).
  - **Spacing (2026-11-12):** the free nucleotides between neighbouring ribosomes are drawn from an exponential distribution around the mRNA's mean (ribosomes load independently; before, every linker of an mRNA held the mean), and the walk keeps a drawn linker while it tries up to 40 orientations of the next ribosome (before, a short linker that collided was replaced by a new, mostly longer one), with up to 1,200 failed tries per chain. In an empty-cell test this raised ribosomes in tight pairs (exit to entry < 8 nm) from 2.8% to 13.8% of chain ribosomes and placement from 90% to 94%.
  - **Not modelled:** the stacked relative orientations of ribosomes in tightly packed polysomes (cryo-ET).
  - **mRNA channel:** the path through each ribosome comes from 7ZW0, superposed on 4V88.
  - **Collided ribosomes (2026-11-12, `polysome.disome_pose`):** where the drawn free mRNA between two neighbours is at most 5 nt, the trailing ribosome takes its pose in the yeast collided disome 6I7O (Ikeuchi *et al.* 2019 EMBO J 38:e100276; both ribosomes superposed on the library 80S by sequence-aligned rRNA phosphates, median deviation 0.10 nm). The cut-off is from disome footprints (57–63 nt against 27–29 nt for monosomes; Meydan & Guydosh 2020 Mol Cell 79:588). The coarse sphere trees overlap in the structure's pose, so the trailing ribosome is moved 4.75 nm apart along the line between the centres (26.4 instead of 21.6 nm; exit-to-entry 8.2 instead of 4.8 nm), keeping the structure's relative orientation; the poly(A)-stalled disome 6T83 (Tesina *et al.* 2020) lies within 15° and 2.9 nm of 6I7O. How often it happens follows from independent loading alone. If the pose clashes, the ribosome is placed as before with a 5-nt linker.
  - **Result (build S):** 116,723 ribosomes on 17,566 mRNAs (mean 6.6 per mRNA; 7.6 requested). Chains stop where they run out of space, so 88% of the 132,700 ribosomes meant for polysomes are in them (62% of cytosolic ribosomes, against the assumed 70%; build R 115,139 and 61%; the generic mRNAs before 2026-11-11 reached 44%); the rest are free. 3,704 ribosomes were placed collided behind their neighbour, so 3.8% of cytosolic ribosomes are in collided pairs, against 5.8% by disome profiling (Zhao *et al.* 2021 Genome Biol 22:16; test; build O had none). Exit-to-entry distances: median 13.7 nm (5–95% 7.8–31.9 nm); with the collided pairs, 10.2% of cytosolic ribosomes are in tight pairs (under 8 nm), against 14.5–18.6% by cryo-ET (build R 10.0%, build O 7.2%, build M 1.4%). These shares are the same as in builds Q and R: cytosolic crowding was kept when the cell was made smaller.
- **Full-length mRNA:** all mRNAs are threaded at full length through the cytosol before protein packing.
  - **Lengths:** gene-specific (above). The 5′ tail is the 5′ UTR plus the part of the CDS upstream of the ribosomes (assumed centred on the CDS); the 3′ tail the rest. Each thread records its gene (`mrna_gene` in the pre-pass).
  - **Representation:** one bead per nucleotide, 0.6 nm apart.
  - **Polysome tails:** the untranslated remainder leaves the first and last ribosome as 3′ and 5′ tails. The 5′ side gets 10–40% of it.
  - **Result (build S):** 85% of the target nucleotides were placed (as in build R). Threads that run out of room end early. In all, 19.2 million nucleotides of mRNA are drawn, including the ribosome channels (27.7 million in build P, with 1.4 times the mRNAs).
  - **Proteins on the ends (2026-11-11; `pipeline/mrnp.py`):** Pab1 on each 3′ end, max(1, poly(A) // 27 nt) copies up to 3 (a footprint of about 27 nt, unchecked; per-gene tails, so mostly one), and one eIF4E (Cdc33) on each 5′ cap. Each is an AlphaFold monomer at a random orientation with one of its spheres touching the end bead. Copies come out of the free pools (Pab1 100,100, Cdc33 28,800). In build S, 15,299 Pab1 (84% of the 18,158 wanted) and 17,529 eIF4E (97% of 5′ ends) fit (build R: 15,046 and 17,425). eIF4G, eIF4A and the closed loop are not modelled.
- **ER-bound ribosomes:** 5,675 ribosomes are docked on translocons in the ER and the outer nuclear membrane (5,586 fit in build S; 5,507 in build R).
  - **Templates:** 4,449 on Ssh1 (2WW9, yeast 80S–Ssh1) and 1,226 on Sec61 (7KAH, placed through the shared Sss1 subunit).
  - **Assumption:** 80% of free translocons are ribosome-bound.
  - **ER-associated ribosomes:** tomograms show 550–900 ribosomes per µm² on the cytoplasmic face of the cortical ER (West *et al.* 2011 J Cell Biol 193:333). The translocon pools fall far short of that, so more free ribosomes are packed in the band next to the ER and the outer nuclear membrane, up to 725 per µm² with the docked ones (the outer nuclear membrane is assumed to match). 11,000 and 5,400 are requested; the envelope face is capped since build S (section 4.4), which leaves 4,300 for it. In build S, 6,924 of 10,953 fit next to the ER and 3,378 of 4,341 next to the outer nuclear membrane (build R: 6,121 of 10,952 and 3,263 of 5,394); the density reached with the tethered polysomes is measured below.
  - **Not modelled:** ER-bound polysomes and their mRNA. The docked ribosomes are single ribosomes, not members of the chains of secretory mRNAs.
- **mRNAs at the mitochondrial surface (2026-11-11; `pipeline/mrna_localization.py`):** polysomes of genes whose translating ribosomes are enriched next to the mitochondrial outer membrane are tethered to it. The genes are the 223 −CHX calls of proximity-specific ribosome profiling in glucose (Williams *et al.* 2014 Science 346:748; 147–196 of them are expressed in one draw, about 1.5% of mRNAs).
  - **Tether:** the first ribosome starts in the cytosol 16–30 nm from the outer membrane, and every ribosome of the chain stays within 45 nm of it. After 300 failed starts a chain is placed untethered. Free (untranslated) mRNAs of these genes are not tethered.
  - **Test pre-pass:** the tethered ribosomes sit a median 30 nm from the outer membrane (99% within 50 nm, against 3% of the others), but only 69% of the requested ribosomes fit on the tethered chains.
  - **Held-out test:** the mitochondrial-fraction mRNA shares of Saint-Georges *et al.* 2008 (PLoS One 3:e2293; galactose) in `validation/mrna_localization.py`. Build S: 275 polysomes of 147 genes tethered, Spearman 0.45 and AUROC 0.66 over 561 genes. Build R, with the same mRNAs and tethers: 0.35 and 0.61. Builds N to Q, with 1.4 times the mRNAs: 0.41–0.44 and 0.66–0.68. With fewer mRNAs per gene, the share of a gene's mRNAs near mitochondria is a coarser number (often one or two molecules), and fewer class I genes are expressed in the draw; the rule is unchanged. Builds R and S share the mRNA draw and the tethers and differ in geometry and packing, so the difference between them (0.35 against 0.45) is a measure of how much one realisation moves this test. Build M: 386 polysomes of 170 genes tethered; over 610 genes Spearman 0.44 between each gene's share of mRNAs with a ribosome within 50 nm of the outer membrane and its mitochondrial-fraction share, AUROC 0.66 for class I against classes III–IV (untethered build: −0.05 and 0.47). The model is all-or-none per gene, while fractionation finds a median 43% of class I mRNAs on mitochondria.
- **mRNAs translated on the ER (2026-11-12, IMPROVEMENTS M5):** polysomes of the secretome genes (signal sequence or transmembrane domain) whose translating ribosomes are enriched next to the translocon are tethered to the cortical ER the same way. The genes come from ER proximity-specific ribosome profiling in glucose (Jan, Williams & Weissman 2014 Science 346:1257521, Table S5: BirA–Sec63 with cycloheximide, enrichment above the 95th percentile of non-secretome genes, their definition; 873 genes). Seeds lie 16–30 nm from the ER membrane and at least 45 nm from the plasma membrane: the cortical ER's face towards the plasma membrane is a ~30 nm gap without room for ribosomes. The translocon-docked single ribosomes stay as before; joining them into these chains is not done.
  - **Held-out test:** membrane-polysome enrichment per gene by fractionation in YPD (Chartron, Hunt & Frydman 2016 Nature 536:224, Supplementary Table 1) in `validation/mrna_localization.py`. Build M, without the ER tether: Spearman −0.08 between the share of a gene's mRNAs with a ribosome within 50 nm of the ER and its log2 enrichment (AUROC 0.49 for the top decile). Build S: 1,741 polysomes of 556 genes tethered; Spearman 0.35, AUROC 0.89 (build R the same; builds N to P 0.35 and 0.90–0.91). Both methods ultimately reflect signal sequences and transmembrane domains, so a sequence rule would also score well.
  - **ER ribosome density (`validation/er_ribosome_density.py`):** with the tethered polysomes, 19,300 ribosome centres lie within 30 nm of the cortical ER, 827 per µm² of its cytoplasmic face (23.3 µm²; West *et al.* 2011: 550–900; build R 796, build M about 470). Builds N to R agree within 2%; build S is 4% above build R; round N quoted about 790 from a one-off area estimate. Partly tuned: the ER tether was chosen with this number in view.
- **RAC on ribosomes:** the ribosome-associated complex (Zuo1–Ssz1; 63,600 copies) sits on the 60S subunit at the tunnel exit and is ribosome-bound in vivo, so every copy rides on a ribosome. In build S that is 39,200 on polysome ribosomes, where its own sphere tree must clear the neighbouring ribosomes, and 24,400 on free ribosomes (build R: 38,700 and 25,000).
  - **Pose:** from the *C. thermophilum* 80S–RAC structure 7Z3N (RAC-1). Its rRNA is superposed on the library 4V88 by 4,558 sequence-aligned phosphates (median deviation 0.27 nm).
  - **Stand-in:** its Zuo1 and Ssz1 (58% and 47% identical to the yeast proteins) stand in for the yeast RAC. The yeast RAC structure 7X3K was deposited without its ribosome and, superposed on 7Z3N, would float off it (2.4 nm deviation; another state).
- **Asc1 and Stm1 (2026-11-12, found by the PDB audit of section 5.2):** both are chains of the 80S structure 4V88 and were also placed as free proteins (316,000 and 188,000 copies in build P). They are now booked to the ribosome (`ingredients.RIBO_EXTRA`): Stm1 has fewer copies than there are ribosomes, so none is free; 21,300 Asc1 remain free. Stm1 sits in the mRNA channel of idle ribosomes (4V88 is from starved cells), so translating ribosomes are drawn with a protein they do not carry; a structure of a translating 80S is an open item of the audit.
- **Not placed on ribosomes:** NAC (Egd1–Egd2), Ssb, eEF3, eIF5A and eEF1A are free in the cytosol, and no tRNA is on a ribosome. Yeast structures exist for all but NAC (PDB audit, `pdb_audit_2026-11-12/`); not built yet.
- **tRNAs (2026-11-12; `pipeline/trna_species.py`):** one type per isoacceptor instead of one tRNA-Phe model (1EHZ) for all of them.
  - **Sequences:** the 275 nuclear tRNA genes of S288C in GtRNAdb (mature, introns removed; Chan & Lowe 2016 Nucleic Acids Res 44:D184) form 42 isoacceptors (anticodon families, initiator and elongator Met apart; two genes of undetermined identity are left out). Each is drawn as its most common sequence, with the 3′ CCA that yeast nuclear tRNA genes do not encode: 74–78 nt, and 85–87 nt for the class II tRNAs (Leu, Ser) with their long variable arm.
  - **Abundance:** the 2.31 million tRNAs (3.3 million for the population-mean cell, an assumption, times 0.7; section 2.2) are split by gene copies (tRNA levels follow gene copy number in growing yeast; Percudani *et al.* 1997 J Mol Biol 268:322). Build S places 2,309,616, all in the cytosol (build R 2,309,466).
  - **Structures:** predicted with Boltz-2 (Passaro *et al.* 2025, bioRxiv; an AlphaFold3-like model; best of 5 samples), confidence 0.73–0.89. The predicted tRNA-Phe is within 1.2 Å P-atom RMSD of the crystal structure 1EHZ (not an independent test: 1EHZ is surely in Boltz-2's training data). Modified nucleosides are not modelled. The viewer opens each model with "Load atomic model" (`data/models/`).
- **Introns (found and fixed 2026-11-12, build O):** through build N, the CDS length of each gene was its SGD ORF span, which includes introns.
  - **Affected:** the 208 intron-containing genes, mostly ribosomal-protein genes. RPL25, for example, got 843 nt instead of 429.
  - **Share of build N:** 18% of mRNA molecules and 17.5% of polysome ribosomes. Intron sequence was 5.2% of all mRNA nucleotides. These mRNAs were too long and, since ribosome load scales with CDS length, carried too many ribosomes.
  - **Fix:** the CDS is now the sum of the ORF's SGD CDS segments. TIF-seq 5′ UTRs are genomic distances, so a 5′ UTR intron (SGD; 21 of the genes here, mostly ribosomal-protein genes) is subtracted where the measured UTR is longer than the intron. RPS8B's 5′ UTR, for example, goes from 389 to 29 nt.
  - **Effect in build O:** 5.8% fewer mRNA nucleotides; ribosomes on mRNAs 168,017 instead of 171,872 (median 5 per mRNA instead of 6); the mRNA tests are unchanged within noise (above).
- **In the viewer (2026-11-12; `pipeline/polymer_annot.py`):** hovering or clicking an mRNA names its gene and the part it belongs to: the 5′ tail, the stretch through ribosome *k* of *n* (counted from the 5′ end), the 3′ tail, or a free mRNA.
  - **Positions:** nucleotides are counted from the end of the tail, the cap or the 3′ end, where eIF4E and Pab1 sit. A thread that ran out of room ends early, and its drawn end counts as its end. The regions (5′ UTR, CDS, 3′ UTR, poly(A)) follow from the gene's layout.
  - **Between ribosomes:** the mRNA there is drawn shorter than the nucleotides it stands for, since each linker is its end-to-end distance. Those pieces are named by their ribosome only.
  - **Ribosomes:** polysome ribosomes name the mRNA they translate and their rank on it.

### 6.6 Ribosome assembly: pre-ribosomes and the rDNA

- **Pre-ribosomes:** 19,500 assembly intermediates in five stages, from the nucleolus to the cytoplasm (structures and counts in section 5.2). They are packed like other complexes in their compartments, except the 90S particles still attached to nascent transcripts (below).
- **rDNA:** about 150 tandem repeats of 9.1 kb on chromosome XII. Half of them (75) are active (psoralen accessibility; Dammann *et al.* 1995).
  - **Active genes:** nucleosome-depleted (Merz *et al.* 2008; later ChIP studies argue for unphased nucleosomes instead), with 50 RNA polymerase I complexes on the 6.7 kb 35S transcription unit, one every 134 bp (Miller spreads: 50 ± 20; French *et al.* 2003). At 60 nt/s this finishes about 2,000 transcripts per minute, the ribosome production rate.
  - **Geometry:** each gene is a persistent self-avoiding walk of polymerase sites, 45.6 nm of DNA apart, inside the nucleolus. The DNA (2 bp per bead) is a smooth spline through the sites and passes through each polymerase's cleft: Pol I is oriented by the DNA axis of the yeast Pol I elongation complex 5M3F, superposed on the library's 4C2M. The walks are not connected into the real tandem array.
  - **Nascent pre-rRNA:** one transcript per polymerase from its RNA exit, as long as the polymerase has transcribed (0–6.7 kb), threaded self-avoiding at 2 nt per bead (0.6 nm per bead, about the rise of base-paired RNA). About 70% of transcripts are cleaved at site A2 while still being made (Kos & Tollervey 2010): those keep only the part beyond A2. Uncleaved transcripts beyond A2 carry a 90S particle on their 5′ end (5′ETS and 18S inside it).
  - **Inactive repeats and the spacers of active ones:** nucleosomal chromatin in the nucleolus, 5,227 nucleosomes in persistent walks.
- **Result (build S):** 73 of the 75 genes fit, with 3,650 Pol I (of 10,000; the rest are free in the nucleolus), 329 transcripts with a 90S (of 666 that should have one), and 57% of the nascent RNA placed (build R: 3,601 Pol I, 349 of 649, 57%). The number of active genes and their polymerases are per genome and do not scale with the cell, so the smaller nucleolus is more crowded: 170 µm of transcribed DNA in 0.47 µm³ (0.67 before build R) puts gene axes about 50 nm apart, and threads that run out of room end early. Real nascent transcripts are compact RNPs, so they need less room than these threads.

### 6.7 Mitochondrial inner membrane

- **ATP synthase dimers:** 3,880 of the 15,300 ATP synthase copies form 1,940 dimers (build S), which fill the rim sites that fit (72% of 2,703 tried); the rest stay monomers in the flat inner membrane. Build R, with a crista across the tubule every 42 nm, had 6,257 dimers (52% of 11,935 sites). The short cristae of glucose-grown cells (section 3.3) have far less rim, so most ATP synthase is now monomeric in the inner boundary membrane; in glucose cells ATP synthase is 2–3 times less abundant per mitochondrial mass than in respiring cells (Klecker & Westermann 2021, citing Di Bartolomeo *et al.* 2020); the model's copy number is PaxDb's (section 2.1). Each dimer is two copies of the yeast monomer 8F39 arranged as in the cryo-ET dimer 4B2Q (dimer angle 83.5°).
  - **Placement:** in rows along the crista rims, every 13 nm (assumed spacing). The residue test of each rim site uses only the geometry near its crista and is compiled (build S; same tests and the same result on build R's geometry).
  - **Missing subunits:** 8F39 lacks the dimer-specific subunits e and k.
- **Respiratory supercomplexes:** complex IV is entirely in 2,031 III₂IV₂ supercomplexes (6HU9). The remaining 1,850 complex III dimers are free (builds R and S).
- **A crowded inner membrane (build S):** the inner membrane has 13.1 µm² instead of 31.5 for 309,800 requested membrane-protein copies (build R: 301,200; more ATP synthase stays monomeric now), and holds 87.4% of them (build R: 98.2%). Carriers lose most: Yhm2 40%, Ggc1 42%, Odc1 and Odc2 37–38% placed. Membrane regions are not thinned (section 8.2).
- **Mitoribosomes on the membrane:** in yeast the vast majority of mitoribosomes are tethered to the inner membrane, exit tunnel towards it, through Mba1 and the 21S rRNA expansion segment 96-ES1 (cryo-ET; Pfeffer *et al.* 2015 Nat Commun 6:6019).
  - **Structure:** the 74S monosome of 9HM0 (yeast mitoribosome with mRNA; the translational activator Aep3 left out).
  - **Frame:** the axis from the large subunit's centre to the exit-site proteins (uL23m MRP20, uL24m MRPL40, uL29m MRPL4) is the membrane normal, pointing into the matrix; the monosome is lowered until its nearest bead is 1 nm above the matrix-side leaflet, which puts the exit site 5.1 nm above it (the gap Mba1 and Oxa1 bridge). No residue enters the bilayer.
  - **Copies:** 90% of the smaller subunit pool (assumption from "vast majority"): 1,367 monosomes, of which 1,364 fit in build S (1,311 in build R), along the inner boundary membrane and at crista tips. In build R, with cristae 24 nm apart along the whole tubule, the 37 nm tall monosome had no room between them. The remaining 152 large and 169 small subunits are free in the matrix.
- **mtDNA nucleoids (revised for build S, 2026-11-13):** a haploid cell in glucose has 1.2–1.35 nucleoids per µm of tubule, about 0.8 µm apart, with one within 0.5 µm of 69% of the tubule tips (Osman *et al.* 2015), and about 21 ± 4 copies of the 85.8 kb genome per haploid cell in YPD (Göke *et al.* 2020), about 0.45 per fL of cell (Seel *et al.* 2023; both W303). The geometry chooses the sites (1.27 per µm of tubule, near tips first, at least 550 nm apart and 250 nm from branch points, off the axis towards the inner boundary membrane), and no crista starts within 110 nm of one: a nucleoid is a compact body in crista-free matrix (size from Göke *et al.* 2020; position by analogy with mammalian cryo-ET, Kukat *et al.* 2015; no yeast EM places nucleoids relative to cristae). Build S: 11 nucleoids holding 15 genome copies (0.45 per fL × 34 fL; 1–3 per nucleoid), 1.27 Mb of DNA drawn at 2 bp per bead. Up to build R: 15 nucleoids anywhere in the matrix with 35 copies (haploid cells: 10–20 nucleoids, 25–50 genome copies, as summarised by Zelenaya-Troitskaya *et al.* 1998 Genetics 148:1763), 3.0 Mb.
  - **Walks:** each nucleoid is a cluster of 5 kb loops, self-avoiding walks (10 bp per 3.4 nm step) that start near the nucleoid centre and are drawn back to it. Ninety percent of each nucleoid's DNA lies within 53 nm (median) of its centre (build R, threaded between cristae every 42 nm: 63 nm).
  - **Abf2:** the packaging protein (35,600 copies) is packed free in the matrix.
- **Not modelled:** III₂IV₁ supercomplexes, confinement of supercomplexes to the flat crista membranes, and mitochondrial polysomes.

### 6.8 Eisosomes

- **Structure:** the near-native Pil1 and Lsp1 lattices of Kefauver *et al.* 2024 (8QB7, 8QB8). These were solved on eisosome tubules isolated from yeast plasma membranes, with the membrane lipids still bound.
- **Lattice:** the seven dimers of each entry are related by screw operations about the tube axis (0.03 nm rms). Two of them span the lattice: 13.6° with −4.0 nm along the axis, and −37.7° with −3.1 nm. That gives 44.8 nm² per dimer at the coat's mid-radius.
- **Membrane:** the amphipathic helix (residues 39–48) sits 11.5 nm from the tube axis. Taking it as 1.3 nm above the bilayer midplane (assumption: an interfacial helix in the cytoplasmic leaflet) puts the midplane at 10.2 nm, which is the furrow-bottom radius in section 3.2.
- **Coat:** whole dimers on the half-cylinder bottom of each furrow, within ±60° of its lowest line, so that no dimer runs onto the flat walls. The lattice phase is random per furrow. Each site holds a Pil1 or an Lsp1 dimer in proportion to the pools; the two lattices are nearly identical.
- **Result (build S):** 3,757 dimers on 23 furrows (2,013 Pil1, 1,744 Lsp1; 163 per furrow), 7.6% of the cell's Pil1 and Lsp1 (build R: 3,660 on 22 furrows, 7.4%). The rest stays near the membrane as free monomers and dimers. At about 2,500 molecules per eisosome (unchecked), the measured copy numbers imply coats much larger than one half-pipe 300 nm long: longer furrows, coated walls, or a cytosolic pool.
- **MCC membrane proteins:** the furrows are the membrane compartment of Can1 (MCC). Its tetraspanning markers Sur7 and Nce102 are concentrated there (Malínská *et al.* 2003; Grossmann *et al.* 2008). 80% of their copies (assumption) go onto the flat furrow walls, oriented like other plasma-membrane proteins, as many as fit. The rest stay in the rest of the plasma membrane.
- **Not built:** Can1 and other transporters that move in and out of the MCC with conditions; PI(4,5)P₂ or sterol enrichment in the furrows; the stretched state of the lattice.

### 6.9 Cvt complex

- **What it is:** the cytoplasm-to-vacuole targeting cargo of growing cells. Precursor aminopeptidase I (prApe1) dodecamers assemble through their propeptides into one large complex with α-mannosidase Ams1, coated by the receptor Atg19. It waits at the vacuole to be wrapped into a Cvt vesicle of about 150 nm (Baba *et al.* 1997; Lynch-Day & Klionsky 2010).
- **Built:** one cluster in the cytosol, 40 nm from the vacuole. It holds half of the Ape1 dodecamers (the library 12-mer 5JGF; the propeptides are not resolved) and half of the Ams1 dimers. Both halves are an assumption for the share still in precursor form. They are packed at random to 35% volume fraction, with all Atg19 on the surface.
- **Result (build S):** radius 84 nm (168 nm across, close to a 150 nm Cvt vesicle; build R 162 nm, 180 nm at 5.0 pg), 229 Ape1 dodecamers, 137 Ams1 dimers, 254 Atg19. The other halves of Ape1 and Ams1 stay in the vacuole as mature enzymes.

### 6.10 Nucleus–vacuole junction

- **What it is:** where the vacuole lies against the nucleus, the outer-nuclear-membrane protein Nvj1 binds Vac8 on the vacuolar membrane (Pan *et al.* 2000 Mol Biol Cell 11:2445).
- **Built:** all 532 Nvj1, each with one Vac8, as the Vac8–Nvj1 crystal structure 5XJG (Jeong *et al.* 2017 PNAS 114:E4539). They sit on the cytoplasmic face of the vacuolar membrane wherever it is within 45 nm of the outer nuclear membrane, which in build S is 1.01 µm² (build R, with 12 lobes: 0.86 µm²; 0.46 µm² in the larger cell of the builds up to Q).
- **No pores (build S):** the envelope of the junction carries no nuclear pore (section 6.1). Vac8 is against the membrane (its palmitoylated N-terminus is not resolved), with the Nvj1 fragment towards the nucleus.
- **Not drawn:** Nvj1's transmembrane anchor and linker, which are not resolved.

## 7. Procedural content (generated in the browser)

These are generated deterministically from fixed seeds when you zoom in. They are not stored coordinates and do not count in the densities.

- **Lipids:**
  - **Packing:** 0.65 nm² per lipid per leaflet.
  - **Composition per membrane** (Zinser *et al.* 1991 J Bacteriol 173:2026, isolated organelle membranes): phospholipid classes (PC, PE, PI, PS, PA, cardiolipin, other) from their Table 5 and ergosterol from the ergosterol : phospholipid ratios of their Table 4. ER: light microsomes; Golgi/endosomes: secretory vesicles; nuclear envelope: nuclei; lipid-droplet monolayer: ER phospholipids without sterol. Examples: inner mitochondrial membrane 12% cardiolipin, outer 6%; plasma membrane 33% of phospholipid PS; ER 41% PC and 21% ergosterol.
  - **Plasma membrane:** ergosterol 40 mol% (assumption; the measured ratio of 3.3 would mean 77 mol%, more than a bilayer holds). Leaflets (consensus): PS in the cytoplasmic leaflet, PE 80% there, sphingolipids (in "other") in the exoplasmic leaflet. Other membranes are drawn symmetric.
  - **Caveats:** mitochondria, microsomes and peroxisomes were isolated from lactate- or oleate-grown cells, the others from YPD; fractions carry some cross-contamination.
  - **Masking:** lipids within 0.95 nm of a protein bead are removed.
  - **In the viewer:** clicking a membrane lists the composition of that membrane (and leaflet). In cut views the zoom level at which lipids are drawn is tested at the cut slab for every tile (2026-11-13); tested per tile, it drew lipids only on the tiles nearest the eye, a wedge on membranes seen face-on.
- **Cell wall:**
  - **Layers:** an inner β-1,3-glucan layer (semi-flexible chains, 0.018 beads/nm³) and an outer mannan layer (branched brushes on the mannoproteins).
  - **Amount:** bead counts are set so the wall is about 20% of cell dry mass.
  - **Chitin:** 1.5% of the wall (Klis *et al.* 2002). Half of it (assumption: chitin is mostly in bud scars) fills the bud-scar rim as microfibrils running around the ring, one bead per 1.04 nm of an ~18-chain α-chitin microfibril (36 GlcNAc per bead, like the glucan beads).
  - **Collisions:** only partly checked.
- **Vacuolar polyphosphate:** the vacuole is the cell's phosphate store. The VTC complex in its membrane synthesises polyphosphate (polyP) straight into the lumen.
  - **Amount:** reported as "hundreds of millimolar" of phosphate units (Bru *et al.* 2025 eLife 14:RP108181). Here it is 300 mM in the lumen (assumption), which averages to about 27 mM over the cell. Whole-cell measurements range from about 23 to 250 mM; the location is disputed between older fractionations and the current vacuolar view.
  - **Chains:** short and linear, lognormal lengths around 24 units (the mean measured in polyP-rich yeast, Christ *et al.* 2020).
  - **Drawing:** one bead per 4 phosphate units (1.2 nm of chain), generated only inside the cut slab, kept 3 nm inside the vacuolar membrane and clear of vacuolar proteins.
  - **Not drawn:** counter-ions (Mg²⁺, K⁺, basic amino acids). It does not count in the protein densities.
- **Motion (Display → Motion; off by default):** illustrative animation, not molecular dynamics. Nothing is simulated and neighbouring molecules are not checked for collisions while they move.
  - **Wobble and flexing:** a small rigid-body wobble whose amplitude falls with molecular size (about 0.35 nm translation for small proteins). Membrane proteins move only within the bilayer and rotate only about the membrane normal. Residues with low pLDDT flex by up to 0.45 nm. Free polymers (mRNA, linker DNA, wall polysaccharides, FG repeats) wiggle along their length.
  - **Static:** actin filaments, microtubules, nuclear pores and nucleosomes stay still, so that their assemblies stay intact.

## 8. Packing

### 8.1 Spatial decomposition

- **Tiles:** the cell is cut into 200 nm cubic tiles; 5,912 of them intersect the cell (builds R and S; 6,416 in the larger cell).
- **Copies per tile:** each species' expected copy number per compartment is split over tiles in proportion to the compartment's volume (or membrane area) in each tile, modulated by the large-scale organisation of section 8.4. A multinomial draw with stochastic rounding conserves the cell totals.
- **Concentrations:** total macromolecular crowding is uniform within each compartment; the composition varies (section 8.4).

### 8.2 Algorithm

Random sequential addition in each tile, in rounds by size class (since build S).

**Before packing: what each region is asked to hold (build S, 2026-11-13; `pipeline/wc3d.py`).** The measure is the requested sphere volume fraction of a region: the summed sphere volume of the packing shapes of its copies per region volume (the spheres of one molecule overlap, so this exceeds the fraction they really fill). The cytosol is asked for 0.41.
- **Membrane-face shells** are capped at the cytosol's fraction; what a shell is asked for beyond it goes to the cytosol, every species by the same factor (section 4.4: Golgi/endosome face ×0.61, nuclear-envelope face ×0.80).
- **Over-full volume regions are thinned uniformly.** Regions asking for more than 1.1 times the cytosol's fraction are packed first in a pilot: 16 tiles per region, drawn in proportion to the region's volume in them, not adjacent, without neighbours, one round. A region whose pilot tiles held less than 98% of their requested sphere volume has the request of every one of its species multiplied by 0.97 times the share held (`work/yc3d/thinning.json`). In build S the nucleolus asked for 0.64 and its pilot tiles held 0.935, so every nucleolar species is packed at 0.907 of its inventory copies; the matrix asked for 0.56 and held 0.991, so it is not thinned. This replaces a shortfall that fell on mid-sized species (build R: Hsp60 49%, Sod2 37% and Ssc1 96% placed in the matrix) by one reported factor. The inventory itself is not changed, and the validation reports placed copies against the inventory as well as against the request. Membrane regions are not thinned. (Candidate build 1 drew the pilot tiles with the most of the region, its roomiest parts; they overstated what the matrix holds.)

**Placement:**
- **Rounds (build S):** every tile is visited in 6 rounds. Round 0 places the membrane proteins and the largest soluble species; later rounds place smaller size classes, split by bounding radius so that each class holds about the same requested sphere volume (thresholds 15.6, 8.1, 7.1, 5.6 and 5.0 nm in build S). Each round packs the tiles in the 8 checkerboard passes below, in the reverse pass order in every second round, so that each tile is as often ahead of a neighbour as behind it. Within a round, molecules go in by decreasing size. Up to build R each tile was packed once, membrane proteins first and then the rest by decreasing size: first-pass tiles over-filled the 5 nm at their faces (molecules protrude into still empty neighbours) and last-pass tiles under-filled it (section 9).
- **Soluble molecules:** a uniformly random position in the target compartment and a uniformly random orientation, or, for a share of copies, a position next to an interaction partner (section 8.4).
- **Membrane proteins:** projected onto the exact membrane surface, with their membrane normal along the local surface normal and a random spin about it.
- **Retries (2026-11-12, after the Jurkat cell model):** up to 2,000 trials per copy, and a job stops after 1,200 copies in a row fail (was 200 and 25). With the old budget tiles stopped early: on test tiles the Golgi face rose from 57% to 66% placed, the envelope periphery from 43% to 51%, the vacuole face from 57% to 63%, the cytosol of crowded tiles from 92% to 100%, the mitochondrial matrix only from 47% to 51% (space-limited).
- **Flat alignment:** in thin lumens and membrane shells (ER lumen, matrix, intermembrane space, Golgi lumen and the peripheral shells) half of the trials align the molecule's flattest principal axis with the normal of the bounding membrane; a small further gain in shells (vacuole face 63% to 67% on a test tile).
- **Membrane-face shells (build S):** a membrane-face protein proposed next to a partner (section 8.4), or moved by the squeeze below, must still have its centre within its shell width plus 4 nm of its membrane. In build R 4–17% of the placed copies of these shells had ended outside them (Golgi face 17%, envelope 11%, lipid droplets 7%, vacuole 7%, ER 5%, plasma membrane 4%).
- **Squeeze relaxation:** random sequential addition jams at roughly 35–40% sphere-tree volume fraction, short of cellular crowding. In the last round (since build S; before, after the single pass), in a tile where copies did not fit, those copies are inserted at random positions in their compartment, in rounds of at most 15% of the tile's soluble molecules (4 rounds).
  - **Relaxation:** after each insertion round, all soluble molecules of the tile are relaxed as rigid bodies for 60 iterations. Overlapping spheres of different molecules push each other apart (translation and a small rotation per molecule, at most 0.6 nm and 0.06 rad per iteration). Spheres crossing a membrane or leaving their compartment are pushed back along the membrane's distance gradient, and centres stay in the tile.
  - **Fixed:** membrane proteins, pre-placed structures and finished neighbour tiles do not move.
  - **Acceptance:** the configuration is then made exact again. Molecules are accepted one by one, previously placed ones first, only if they overlap nothing accepted so far and lie in their compartment, as in RSA. A round that does not increase the count is undone.
  - **Effect** on a test tile of the mitochondrial matrix: 72% → 89% placed at the earlier protein total, and 49% → 56% at the current one (with the library sphere trees; not re-measured with the packing shapes below). Cytosol and nucleoplasm tiles fit everything by RSA alone.
  - Copies that still don't fit are left out.

**Collision checks:**
- **Packing shapes (build S, 2026-11-13; `wc3d.PACK_PER = 15`):** for packing, each structure is approximated by k-means clusters of about 15 residues (or nucleotides) per sphere, made from the beads and cached in `work/yc3d/pack_trees_15.npz`; each sphere is padded by a 0.42 nm bead radius, and spheres may interpenetrate by up to 0.3 nm. The library's sphere trees of about 30 residues per sphere (40 for the curated complexes) are still what the viewer's simulation and the contact statistics use. As packing shapes they were fat: their union is 2.4 times the union of the residue beads for the most abundant matrix proteins, 1.9 times in the cytosol and 3.8 times in the nucleolus. On six matrix-rich test tiles of build R, 87.9% of the matrix request was placed with the library trees, 93.1% with 15 residues per sphere (2.6 times the CPU time) and 94.9% with 8 (7 times); species under 80% placed: 13, 3 and 9.
- **Compartments:** a copy must lie entirely within its compartment and clear of membranes, which are treated as 2 nm half-thickness slabs.

**Parallelism:**
- **Phases:** in each round, tiles are packed in 8 checkerboard passes, so neighbouring tiles are never packed simultaneously. A tile packed in a later round reads back what earlier rounds placed in it and its neighbours.
- **Halos:** each tile sees the already placed molecules of its neighbours through a halo.
- **Obstacles:** the pre-placed structures of section 6 are fixed obstacles.

### 8.3 Consequences

- **The density is not an equilibrium value.** Random sequential addition jams at lower density than real cytoplasm, and the outcome depends on insertion order; the squeeze relaxation (section 8.2) pushes past the jamming limit but is not an equilibrium either. Where the requested copies exceed what fits, the shortfall is reported, not hidden; where a compartment is thinned before packing (the nucleolus, by 0.907), the factor is reported.
- **Tile seams are gone (build S).** In 918 fully cytosolic tiles the sphere volume fraction within 5 nm of a tile face is 0.41 in first-pass and 0.40 in last-pass tiles, against 0.41 in the interior (15–35 nm from the faces). Build R, packed in one visit per tile: 0.49 and 0.20 against 0.40.
- **Membrane-face shells (build S):** none of the sampled membrane-face copies lies outside its shell (build R: 4–17%). Cytosolic proteins are still sparse in the 4–20 nm next to the Golgi/endosome membranes, where the shell is (0.16 of their density further out; vacuole 0.56, nuclear envelope 0.68), while all molecules together are at 1.25 times the density further out there (`placement_checks.json`, halos).
- **The shortfall is smaller and spread more evenly (build S).** Placed copies are 99.3% of the packing request and 99.1% of the inventory, 99.5% outside mitochondria, and 98.1% of the inventory by sphere volume (build R: 97.5% and 98.9% of the request, 96.0% by volume; builds N to Q 96.5–96.7% and 98.5–98.6%; build M 95.9 / 97.9% before the larger attempt budget). The cytosol (99.7%) and nucleoplasm (99.8%) fit nearly everything. The matrix holds 98.5% of its copies at 178 mg/mL (build R: 75% at 150 mg/mL; 132 in build P, 107 with mitochondria at 2% of the cell), with the finer packing shapes and a larger matrix (0.82 against 0.70 µm³; round S candidate build 1, with the old mitochondria and the library trees, held 73.6% of its inventory). The nucleolus holds 92.7% of its thinned request, 84.1% of its inventory (build R 80%). Shortfalls that remain: the mitochondrial inner membrane (87.4%; section 6.7), the cytosolic faces of the lipid droplets (95.0%), the ER (95.3%) and the nuclear envelope (96.8%); the Golgi face holds 97.7% of a capped request (build R 84% of the full one). Species placed at under half of their inventory copies: 10 of 4,469 (35,000 requested copies; build R 76); under 80%: 60 (build R 238; `validation/placement_checks.py`). Counted per region, 12 species fall below half of the region's share: six inner-membrane carriers (Yhm2, Ggc1, Odc1, Odc2, Mrs3, Yea6), three nucleolar proteins (Tif3, Esf1, Rrp12), two Golgi-face complexes (AP-3, Ric1–Rgp1) and the ER membrane protein Frt2.
- **Interaction-driven organisation is a hypothesis layer.** Section 8.4 places partners near each other and co-enriches network neighbours, but with assumed strengths. There are no explicit condensates, and the layer is not fitted to any measured protein positions.

### 8.4 Organisation from association data

Uniform random placement is itself an assumption: that proteins have no spatial organisation beyond their compartment. The model replaces it with a weak, data-driven organisation at two scales. It never moves a protein out of its compartment: it only changes where in its compartments a protein goes, and next to whom. Neither scale is calibrated against measured protein positions, which hardly exist at this resolution. The strengths are tuned so that the packing stays random-looking but is measurably organised, and total crowding stays uniform. Code: `pipeline/organisation.py`; `scripts/s3_pack.py --random` builds the same cell without this layer.

**Data:**
- **STRING v12** (Szklarczyk *et al.* 2023, Nucleic Acids Res. 51:D638): 271,240 links with combined score ≥ 400 (medium confidence) between 5,881 yeast proteins, weighted by score. STRING includes indirect associations (shared pathways, co-expression, text mining) that a structure prediction cannot see.
- **AlphaFold 3 pool predictions** (in-house, unpublished): every protein was predicted in many random multi-protein pools, giving a chain-pair score for about 10.3 million protein pairs. The best (lowest) chain-pair PAE min-min per pair is used.
  - **Calibration against known interactions:** pairs at ≤ 2 Å are 73% Complex Portal pairs and 72% have IntAct evidence; pairs at 2–4 Å, 43% and 48%; pairs at 6–8 Å are still 50 and 10 times above the background of 0.06% and 0.6%.
  - **Used:** 15,262 pairs at ≤ 8 Å after filtering. A pair is dropped when a paralog of one partner scores better with the other, which removes homologous-interface artefacts such as Pol II subunits scored against their Pol III paralogs.

**Network neighbourhoods (no clustering).** An earlier version split STRING into modules by Louvain clustering; the module resolution and minimum size were arbitrary and changed the result. The current version uses network propagation instead:
- **Profiles:** each protein gets its personalized PageRank vector on the weighted STRING graph: the stationary distribution of a random walk that continues along an edge with probability 0.5 (the damping) and otherwise restarts at the protein. A low damping keeps the profile local, concentrated on the protein's near neighbourhood. The protein's own entry is dropped and the vector renormalised. A complex takes the mean profile of its subunits; 5,452 of the 5,562 species of build S have a profile (the rest are prebuilt assemblies, tRNA and proteins without STRING links; 5,542 of 5,605 in an earlier build).
- **Similarity:** the Pearson correlation *r* of two profiles. It is high for proteins that share a network neighbourhood, whether or not they are linked directly. Typical neighbourhoods are small: 1.0% of species pairs have *r* ≥ 0.4. Examples: TDH3's closest neighbours are TDH1 (0.92), TDH2 (0.91), PGK1 (0.78), TPI1, ENO1, GPM1, PGI1 and FBA1; Sec61's are Sec63, Sec62, the SRP receptor, Ssh1 and OST.
- **Choice of parameters** (`scripts/organisation_sweep.py`): STRING thresholds 400 and 700, damping 0.3–0.85 and four profile variants were scored by how well *r* predicts co-complex membership (Complex Portal co-members against random pairs) and shared localization (same coarse compartment, for proteins with a dominant location). Every variant separates co-complex pairs almost perfectly (AUROC 0.993–0.998); shared localization is predicted only moderately (AUROC 0.56–0.61), much as for pairs that are not directly linked. Log-transformed profiles scored slightly higher on localization but are not local (a typical protein has thousands of partners at *r* > 0.3). Among the local variants, STRING ≥ 400 with damping 0.5 scored best (localization 0.594, co-complex 0.996).

**Large scale (hundreds of nm): network-informed sampling of where each species goes.**
- **Correlated fields:** the species correlation matrix *R* is used as the covariance of the species' log-density fields over the cell. It is approximated as *R* ≈ *E Eᵀ* from its 128 leading eigenvectors (54% of the trace). Each of the 128 components gets an independent smooth Gaussian random field *F*ₖ over the cell (correlation length 500 nm, unit variance). A species' copies are distributed over tiles in proportion to region volume × exp(0.8 × Σₖ *E*ᵢₖ *F*ₖ).
- **Interpretation:** the log-densities of two species co-vary across the cell by 0.8² × *r*: network neighbours are co-enriched in the same places, unrelated species vary independently, and species without a profile are spread uniformly. Nothing is clustered and no module sizes are chosen.
- **Crowding:** an iterative correction keeps each tile's total macromolecular volume proportional to its compartment volume, so composition varies between regions but crowding does not.
- **Effect** (abundant species in fully cytosolic tiles): the correlation of two species' tile densities follows their network similarity (Spearman 0.58 over 86 species, 0.71 over the 24 most abundant; the module version gave 0.44 on the same measure). Pairs with *r* ≥ 0.5 have tile densities correlated at +0.56; unrelated abundant pairs are slightly anti-correlated (−0.18), because the crowding correction makes abundant species compete for volume within a tile. Total crowding varies as little as with random placement (CV 0.040 vs 0.041).

**Meso and small scale (nm to tens of nm): partners placed next to each other.**
- **Contact partners:** pool pairs. Strength falls from 1 at PAE ≤ 2 Å to 0 at 8 Å, is multiplied by (1 + *r*⁺)/2 so that pairs that also share a network neighbourhood pull fully and network-unsupported pairs at half strength, and is divided by √(degree/20) for proteins with more than 20 pool partners, so hubs do not dominate.
- **Loose partners:** the species most correlated with it (*r* ≥ 0.4) that are not contact partners. Strength rises from 0.2 at *r* = 0.4 to 1 at *r* = 1. Each species keeps up to 16 partners, contact partners first.
- **Attraction:** when a copy is placed, each partner p contributes a_p = strength × A × n_p / (8 + n_p). Here n_p is the number of copies of p already placed in the tile and its halo, and A is 0.6 for contact and 0.25 for loose partners. The pull saturates with the partner's local copy number, like mass action.
  - With probability min(0.6, Σ a_p), the copy is proposed next to a randomly chosen placed copy of a partner (chosen in proportion to a_p): in contact for contact partners, with a 2–15 nm gap for loose partners.
  - The copy gets up to 40 such proposals; if none passes the collision and compartment checks, it is placed at random as before.
- **Effect:** 7.1 M copies drew a partner and 5.1 M of them (about 24% of all placed molecules) were placed next to one. A molecule's contact partner lies within contact distance 8.4% of the time, against 1.3% for a control species of the same local abundance, and a loose partner within 15 nm 19.9% of the time (control 14.0%). Placement completeness is unchanged (97.8%). Since build S a proposal next to a partner must, for a membrane-face protein, also lie in its shell (section 8.2). These shares were measured on an earlier build; in build S, 8.7 M copies tried for a partner and 4.9 M (20% of placed molecules) were placed next to one (build R: 9.5 M and 4.7 M, 19%).

**What it is not:** strengths and bound fractions are assumptions, as no affinities are used; the partners are not bound in a specific orientation (pre-docked pairs are a separate future step); and smaller species are pulled towards larger ones that are already placed, since species are placed in order of decreasing size. The same profiles could inform any later packing method (for example a denser, relaxed packing): the layer only supplies a covariance for where species go and a list of partners.

**In the viewer:** Display → colour "Network neighbourhood (STRING)" colours each species from a two-dimensional t-SNE map of the profiles (distance 1 − *r*), so network neighbours get similar colours. Clicking a molecule and then "Network neighbours" highlights it with its 30 most correlated species (*r* ≥ 0.2) and lists them.

## 9. Validation

| Check | Result |
|---|---|
| Residue-bead clashes (beads of different molecules closer than 0.3 nm) | Packed molecules: 3.4 × 10⁻⁶ of 5.3 M sampled beads in 64 test cubes (18 beads, none between tiles; build S; build R 6.4 × 10⁻⁶, builds N to Q 3.7–9.3 × 10⁻⁶), under the 10⁻⁵ acceptance limit; 8 of the 18 are Pab1 or eIF4E on the mRNA end it binds, which they are allowed to overlap. Including the built actin structures: 1.87 × 10⁻⁵ (99 beads; build R 2.65 × 10⁻⁵, 133 beads), above the limit. The other 81 beads are overlaps inside those structures (tropomyosin on its filament 33, filament on filament 18; section 6.4), 0.25–0.5% of their beads in every build; the smaller cell has the same cortex in fewer cubes' worth of volume, so the sample holds more of them. The two rates are reported separately since build R (`clash_rate_without_actin_structures`; four other cube samples of build R gave 3.8–6.6 × 10⁻⁶). The squeeze relaxation presses molecules to the 0.3 nm sphere-tree tolerance far more often than random addition alone (6.3 × 10⁻⁷ before it) |
| Containment (beads outside their compartment) | 6.7 × 10⁻⁵ of 19 M beads; 0.047% of 30,000 sampled instances have a bead outside (build S; build R 4.1 × 10⁻⁵ and 0.04%; round S candidate build 1 9.1 × 10⁻⁵). Relaxed molecules sit against compartment walls, and the packer tests packing shapes, not residue beads. Why the rate rose from build R is not yet analysed (section 11) |
| Membrane-protein orientation vs local surface normal | median 0.0°, all of 10,000 sampled within 1° |
| TM residues within 2 nm of the membrane mid-surface | 97.4% (median 0.72 nm; build R 97.3%, 0.73 nm) |
| Telomeres (Rabl layout) | 32 of 32 within 150 nm of the nuclear envelope, median 164 nm from their targets; telomere–pole distance vs arm length Spearman 0.92 (build S; build R 0.85, earlier builds 0.86–0.90; section 6.3) |
| Eisosome coat: membrane-inserted residues within 2 nm of the furrow mid-surface | all sampled (median 1.36 nm) and all other residues on the cytoplasmic side (builds R and S; an earlier build: 99.2% and 98.9%) |
| ATP synthase dimer TM residues within 2 nm of the mid-surface | 98.8% (builds R and S) |
| Placed / requested copies | 99.3% of the request, 99.1% of the inventory; 99.5% outside mitochondria, above the 99% acceptance target (build S). Lowest: mitochondrial inner membrane 87.4%, nucleolus 92.7% of its request (84.1% of its inventory after thinning by 0.907), LD cytosolic face 95.0%, ER cytosolic face 95.3%, nuclear-envelope periphery 96.8%, Golgi cytosolic face 97.7% (of a capped request); matrix 98.5%, intermembrane space 99.0%. Build R: 97.5% overall, 98.9% outside mitochondria (nucleolus 80%, Golgi face 84%, matrix 75%, intermembrane space 87%) |
| Placement checks (`validation/placement_checks.py`, since build S) | placed 99.1% of the inventory by copies and 98.1% by sphere volume (build R 97.5% and 96.0%); species under half placed 10, under 80% 60 (build R 76 and 238); tile seams: sphere volume fraction within 5 nm of a tile face 0.41 / 0.40 (first / last pass) against 0.41 inside (build R 0.49 / 0.20 against 0.40); membrane-face copies outside their shell 0% (build R 4–17%) |
| Inventory protein mass vs proteome | within 0.4% before packing |
| Protein scale against counted structures (`validation/anchors.py`; check) | At 3.5 pg: core histones 0.90 of what the nucleosomes need (H3 0.70, H2A 0.77, H4 1.04, H2B 1.31); scaffold and FG nucleoporins a median 1.36 times what 93 pores need (10th–90th percentile 0.99–2.15); ribosomes 1.06–1.13 times the 187,000–200,000 counted from rRNA. Each anchor alone would put the protein total at 3.9, 2.6 and 3.1–3.3 pg. The outer kinetochore (Ndc80 complex 6–10 times its 8 per kinetochore) is not an anchor: these proteins have soluble pools |
| Organisation layer (section 8.4) | contact partners in contact 8.4% vs 1.3% (abundance-matched control); loose partners within 15 nm 19.9% vs 14.0%; tile-density correlation vs network similarity, Spearman 0.58; total crowding CV 0.040 (random repack 0.041) (an earlier build; not re-measured) |

### 9.1 Validation against independent data (2026-11-11)

These checks use data that never entered the model, with held-out splits fixed before any tuning (`validation/splits.py`: proteins by a hash of the ORF name, half dev and half test; diffusion measurements by a hash of the paper's DOI). Changes were tuned on dev only; test numbers were computed once per finished change. Code and outputs: `validation/`, `work/validation/`.

**Localization** (`validation/localization.py`). The model's compartment shares are compared with genome-wide localization that does not come from the C-terminal GFP collection behind CYCLoPs and Huh *et al.* (2003): the N-terminally tagged SWAT library under its native promoter and under the NOP1 promoter (Weill *et al.* 2018) and organelle spatial proteomics (hyperLOPIT; Nightingale *et al.* 2019). Classes: cytosol, nucleus, nucleolus, ER and nuclear envelope, mitochondrion, Golgi/endosome, vacuole, cell periphery, peroxisome, lipid droplet; SWAT "punctate" counts as Golgi/endosome, peroxisome or lipid droplet. A protein agrees when the class holding most of its copies is among the reference calls (top-1).

| Reference (test half) | Proteins | Top-1 before | Top-1 after the nucleus rule (4.2) | Macro-F1 after |
|---|---|---|---|---|
| SWAT, native promoter | 1,715 | 64.0% | 66.8% | 0.48 |
| SWAT, NOP1 promoter | 2,614 | 72.6% | 73.8% | 0.50 |
| hyperLOPIT (assigned proteins only) | 468 | 84.6% | 84.6% | 0.75 |

- **Placed copies.** The table scores the inventory's compartment shares. Scored on the copies actually placed in the rebuilt cell, test top-1 is 66.6% (SWAT native), 74.2% (SWAT NOP1) and 84.8% (hyperLOPIT) in build S (build R 67.4, 74.6 and 84.8%; 67.0, 74.1 and 84.5% when the table was made): packing shortfalls, the thinning of the nucleolus and the shell caps (section 8.2) hardly change the comparison.
- **Ceiling.** The references disagree with each other: SWAT (native) and CYCLoPs share a class for only 74% of the proteins both cover. Agreement with CYCLoPs, an input, is a fit (78% top-1), not a test.
- **Where the errors are** (dev half): proteins with measured shares agree far more often (73%) than those placed from annotation alone (55%), and the 625 proteins without any annotation, defaulted to the cytosol, agree least (40%). Peroxisome and lipid-droplet calls are rarely matched (F1 0.10 and 0.03), Golgi/endosome poorly (0.28).
- **Concentration instead of copies.** Scoring by the compartment of highest concentration (shares divided by volumes) agreed worse with every reference (it over-calls the small compartments), so copies are the primary measure.
- **Not adopted:** membrane association of soluble proteins from lipid anchors (myristoyl, palmitoyl, prenyl; UniProt) or membrane-binding domains (PH, PX, FYVE, C2, BAR, ENTH; InterPro, and lipid-array hits of Gallego *et al.* 2010). A bound fraction for domain proteins made dev agreement worse: most of them image as cytosolic. Lipid anchors helped dev by 0.1 point and not test, so the rule is switched off (`membrane_anchor.py`, `BOUND_F = 0`).

**Nucleus against cytoplasm from raw images** (2026-11-12; `validation/image_nc.py`; the first test on pixels rather than on published calls). The C-terminal SWAT library re-imaged for YeastRGB (Dubreuil *et al.* 2019; library: Meurer *et al.* 2018): mNeonGreen fusions under the native promoter, every strain with the envelope marker Nup49-mCherry, 6,145 strains with up to 225 single-cell crops each (9.7 GB of 16-bit images). Per cell the nucleus is segmented from the marker and the nuclear/cytoplasmic intensity ratio measured; the segmentation and thresholds were fixed on four strains before the comparison. Blur compresses the measured ratios, so ranks are compared.
- **Result (test half, 531 proteins, build S):** Spearman 0.81 between the model's and the measured nuclear/cytoplasmic concentration ratio; proteins the model puts in the nucleus are separated from those it puts in the cytosol with AUROC 0.98. Build R: 0.80 and 0.98 (538 proteins); build Q: 0.81 and 0.98 (562 proteins).
- **How independent:** the images are not an input. The nucleus rule (section 4.2) was tuned with CYCLoPs scores of the other half of the proteins, also C-terminal fusions, so tag artefacts shared by C-terminal fusions would not show.
- **Abundance (check):** the mean cellular intensity of 4,886 strains against the model's copies per cell, Spearman 0.74 (0.81 over the brightest 70%).

**Copy numbers** (`validation/copy_numbers.py`; check; the same in builds R and S, which share the copy numbers). At 3.5 pg the model's copies are a median 1.23 times those of the mass-spectrometric proteome of Das *et al.* 2026 in glucose (Spearman 0.81; 46% within 2-fold; 1.75 times at 5.0 pg) and 0.37 times the unified data set of Ho *et al.* 2018 (Spearman 0.93; 0.53 at 5.0 pg). The two data sets differ from each other by a median factor of 3.3, and both are culture averages. The check reads the per-protein feature table, which had not been regenerated since 2026-11-11; the build script now refreshes it (`scripts/run_all.sh`).

**Diffusion** (`validation/diffusion.py`): a literature table of 93 measurements from 25 papers, every value quoted from the paper (`pipeline/data/validation/diffusion_measurements.tsv`); 35 in-vivo values of soluble molecules in unstressed cells form the benchmark (free GFP and tandem GFPs, tagged cytoplasmic and nuclear proteins, 40 nm GEM particles). Free GFP alone ranges from 8.7 to 21 µm²/s between laboratories and methods, so the benchmark has absolute tracer anchors and, separately, within-paper ratios of a molecule to the same paper's free tracer. Section 14.5 describes the model; on the held-out papers:

| Model | Ratio RMSE (log10) | Absolute RMSE (log10) |
|---|---|---|
| Stokes–Einstein × one constant | 0.20 | 0.48 |
| + size-dependent crowding from Brownian dynamics of packed boxes | **0.09** | **0.21** |
| + binding propensity (RNA-interactome hits, pool partners) | 0.17 (build R 0.18) | 0.25 |

Build S, 15 absolute values and 5 ratios in the test split, which since 2026-11-12 includes the 40 nm GEM values (section 14.5; they are why the size-blind null is worse than the 0.31 quoted before). Build S gives the same values to two decimals as build R (ratio 0.09 against 0.09, absolute 0.21 against 0.21). The crowding term comes from Brownian dynamics of boxes cut from an earlier build; the cytosol's density is unchanged (170 mg/mL in build S, 168 in build R, 170 in the build the boxes came from), so it was not re-run. Two free-histone rows (H2B) use the H2B–H2A.Z dimer of the library as their species since build R, because the non-native H4–H2B pair they used before is no longer built (section 5.2) and the library has no free H2A–H2B dimer.

**Condensates** (`pipeline/condensates.py`, IMPROVEMENTS H5). Can membrane-less compartments emerge from generic interaction data at the model's concentrations? Interaction energies came only from the AlphaFold 3 pool scores (PAE converted to the likelihood ratio of a known IntAct interaction, including weak pairs at PAE 8–25 Å, which are still 1.5–10× enriched; pairs within one stable complex excluded), from RNA-interactome capture (Beckmann *et al.* 2015; Matia-González *et al.* 2015; Mitchell *et al.* 2013) and from IDR sticker residues. A multicomponent second-virial free energy over all species gives the stability of the mixture and, by linear response, the enrichment of each protein in an RNA-rich phase. The weights were chosen on the nucleus only (nucleolar members, CD-CODE: AUROC 0.68) and then frozen. In the cytosol, CD-CODE P-body members were enriched (AUROC 0.73, against 0.57 for permuted RNA-binding labels), but:
- the pool network added nothing (the same 0.73 with RNA binding alone and with a permuted network);
- the negative controls failed: stress-granule-only members (0.72) and the hypoxia-induced glycolytic G-body enzymes (0.91) were enriched just as much in log phase;
- the cytosol was stable at the calibrated interaction scale, in log phase and after an *in-silico* stress that releases all mRNA.

Re-run on build S as part of the validation run (`validate_all.sh`): P-body members AUROC 0.74, stress-granule-only members 0.70, G-body enzymes 0.89, the same with RNA binding alone (0.74) and with a permuted pool network (0.74), 0.59–0.64 with permuted RNA-binding labels; the cytosol stays stable in log phase and after the in-silico stress, and the stress changes none of the enrichments. So this shows that RNA-binding proteins co-enrich with RNA, not that P-bodies emerge. The model places no condensates; specific P-body assembly needs information that generic pair scores do not carry (valence, the decapping/Lsm network, regulation).

## 10. Densities, mass and volume fractions

### 10.1 Macromolecular densities

Densities are the mass of placed macromolecules per compartment volume. They count protein plus the RNA of ribosomes, pre-ribosomes, tRNAs and complexes, the DNA of nucleosomes and active rRNA genes, nascent pre-rRNA, and mRNA. Lipids, polysaccharides and small molecules are not counted.

| Compartment | Total, build S (mg/mL) | Build R (mg/mL) | Build P, 5.0 pg in 42 µm³ (mg/mL) | Protein / RNA+DNA split (earlier build, mg/mL) |
|---|---|---|---|---|
| Cytosol | 170 | 168 | 170 | 129 / 36 |
| Nucleoplasm | 111 | 110 | 106 | 99 / 7 |
| Nucleolus | 167 | 164 | 166 (136 before the larger attempt budget) | 108 / 28 |
| Mitochondrial matrix | 178 | 150 | 132 | 101 / 6 |
| Mitochondrial intermembrane space | 25 | 14 | 16 | 21 / 0 |
| ER lumen | 21 | 21 | 25 (44 before the ER topology fix moved 181,000 cytosolic-face copies out) | 44 / 0 |
| Golgi / endosome lumen | 2 | 2 | 3 | 2 / 0 |
| Peroxisomes | 97 | 97 | 120 | 120 / 0 |
| Vacuole lumen | 1 | 1 | 3 | 3 / 0 |

Builds R and S kept the densities of the large compartments by design: the cell volume was chosen so that the cytosol stays near 170 mg/mL (section 3.1). The vacuole lumen is more dilute than before build R because the vacuole is 17% of the cell instead of 9%. In build S the matrix and the intermembrane space are denser: the matrix now holds 98.5% of its proteins (finer packing shapes, section 8.2), and the intermembrane space is smaller (0.18 against 0.27 µm³) and holds more of its proteins (99.0% against 87%).

The protein / nucleic-acid split and the dry volume fractions (cytosol 11.4%, nucleoplasm 7.6%) were computed on an earlier build and are not recomputed by s4.

- **Organelle lumens:** since the 2026-11-11 rule (section 4.4) the vacuole and Golgi/endosome lumens hold only proteins that can reach them (signal peptide or Cvt cargo); their coat and peripheral proteins are counted in the cytosol, whose shells on the organelle faces they occupy. The low lumenal densities do not include glycans, cargo in transit or vacuolar polyphosphate.
- **Dry volume fraction:** computed with partial specific volumes of 0.73 mL/g for protein and 0.55 mL/g for nucleic acid.
- **Hydrated volume fraction:** with a hydration shell of about 0.3–0.4 g water per g, the hydrated volume fraction is roughly 1.5× the dry one: about 17% in the cytosol.
- **Protein and RNA/DNA split:** protein is the placed species' mass minus their nucleotides (330 Da each); RNA/DNA is the rest of the validated total, which also holds the mRNA and other threads. Where volume rounding made the remainder slightly negative it is shown as 0.
- **Excluded volume:** the volume excluded to a diffusing particle is larger still, and grows with the particle's size. It can be computed from the coordinates as accessible volume versus probe radius; this has not been done yet.

### 10.2 Whole-cell totals

| Quantity | Value |
|---|---|
| Protein in the inventory | 3.5 pg (5.0 pg through build Q) |
| Protein placed | about 2.9 pg (build R; not re-measured; 4.2 of 5.0 pg before) |
| RNA and DNA placed | about 0.8 pg (build R; not re-measured; 1.0 pg before): rRNA of ribosomes 0.64, tRNA 0.10, pre-ribosomes about 0.05, DNA 0.013 |
| mRNA | 0.011 pg (19.2 million nucleotides; builds R and S) |

The build R values of placed protein and RNA are estimates from the compartment densities and volumes above and the placed counts of section 2.2, not a separate tally. Build S places more of the matrix (section 8.3) and fewer nucleolar copies (thinning) and fewer mtDNA (1.27 instead of 3.0 Mb); the totals were not re-estimated.

Placed protein is lower than the inventory for four reasons:
- residues not resolved in PDB entries;
- AlphaFold residues trimmed for pLDDT < 50 (3.4% of mass);
- copies that did not fit (mostly mitochondria, nucleolus and Golgi);
- dropped and structure-less copies (including surplus histone H2B and nucleoporins beyond the pores' need, sections 6.1 and 11).

With full-length sequences, cytosolic protein would be about 150 mg/mL.

### 10.3 Comparison with dry-mass measurements

Measured cytoplasmic dry-mass densities are around 200 mg/mL. These are rough estimates of all components, averaged over the cytoplasm (cytosol plus organelles, excluding the vacuole):

| Component | Estimate (mg/mL) | Basis |
|---|---|---|
| Macromolecules as placed | about 160 | This model (build S: cytosol 170, nucleoplasm 111, nucleolus 167, matrix 178; build R 168, 110, 164, 150; volume-weighted over the cytoplasm without the vacuole, 158 in build S and 155 in build R) |
| Full-length protein correction | about +20 | Section 10.2 |
| Lipids of internal membranes and lipid droplets | about 30–35 (build R; not re-estimated) | Build R: about 160 µm² of internal membrane plus 0.12 µm³ of lipid droplets, about 0.8 pg of lipid in 22.7 µm³ of cytoplasm (20–25 mg/mL in the larger cell, which had nearly the same membrane area). Build S has 115 µm² of internal membrane (section 3.2), so this term is smaller |
| Small molecules and ions | probably 30–40 | Free amino acids, nucleotides, glutathione, sugar phosphates, phosphate and K⁺, a few hundred mM in total |
| Glycogen and trehalose | low in log phase on glucose | Can exceed 10% of dry weight after the diauxic shift |

Together these give about 240–255 mg/mL in build R, above the measurements (225–240 mg/mL in the larger cell of builds J to Q); in build S the lipid term is smaller, which was not re-estimated. A newer measurement by refractive-index tomography of wild-type cells in glucose gives about 205–225 mg/mL with the vacuole excluded (Kunzi *et al.* 2026, a preprint; read from their Fig. 1B; check added 2026-11-12). The same study's whole-cell protein concentration, 71 mg/mL by Bradford assay normalised to the Coulter volume, is closer than before but still lower than the model's 103 mg/mL (3.5 pg in 34 µm³; 119 mg/mL at 5.0 pg in 42 µm³). **An unresolved conflict:** the dry-mass measurements would be matched by the same 3.5 pg in a cell of about 36–38 µm³, while the ribosome concentration measured by cryo-ET (23 µM; the model has 17 µM in the cytosol) asks for a smaller cytosol, about 26 µm³ of cell at this composition. Builds R and S sit between the two (build S: 17.1 µM of cytosolic ribosomes); the cell volume was set by keeping the cytosolic density, not by either measurement. The small-molecule estimate is the least certain term.

## 11. Known limitations

| Area | Limitation |
|---|---|
| Physics | Static packing (random sequential addition plus an overlap-removal squeeze): no physical interactions or dynamics and random orientations. The optional motion is illustrative (section 7). |
| Replicates | Build S is one build at master seed 0. The build is deterministic (section 13), but replicate builds at other seeds have not been made, so no number in this document has a build-to-build spread yet. Builds R and S share the copy numbers and the mRNA draw and still differ, for example, by 0.35 against 0.45 in the mitochondrial mRNA test and 1.9× against 3.3× in balanced telomere clustering (sections 6.3, 6.5). |
| Organisation | Partner proximity and network co-enrichment (section 8.4) follow assumed strengths, not measurements. Partners are near each other, not bound in specific orientations. |
| Abundances | About 2-fold per-protein uncertainty. The proteome total is set to 3.5 pg for one unbudded G1 cell, between the stoichiometric anchors (section 2.1; bulk measurements of cultures give 4–6 pg); tRNA and mRNA counts are culture means times 0.7, an assumption. 73 proteins of other genotypes, mating type or conditions are set to zero by a curated list, which is incomplete. tRNA species are split by gene copies, and their structures are predictions without modified nucleosides. Relative abundances do not match counted structures exactly: histone H3 and H2A are 0.70 and 0.77 of what the nucleosomes need (the nucleosome count is fixed, so the missing copies are simply not booked), H2B is 1.3 times (since build S the surplus is packed as free H2B; in build R it had no library entry and was left out), and nucleoporins beyond 1.25 times the pores' need are not placed. |
| Ribosomal proteins | Ribosomal proteins are 15.7% of the model's protein mass, against 29–32% in three proteomics datasets of fast-growing yeast (Metzl-Raz *et al.* 2017, Björkeroth *et al.* 2020, Xia *et al.* 2022; none is BY4741 in YPD), 15–19% in the mass-spectrometry studies collected by Ho *et al.* 2018 and 2–7% in its GFP and TAP studies. At 3.5 pg of protein, 30% would be about 460,000 ribosomes, which no count supports. The counts themselves ask for somewhat more than the model's 211,500 mature ribosomes (11,500 per µm³ of cytosol): about 14,000 per µm³ of cytoplasm by cryo-ET (Delarue *et al.* 2018) and 14,517 ± 3,193 per µm³ of cytosol (Goetz *et al.* 2025, preprint), and about 220,000 in a G1 cell by serial-section EM (Yamaguchi *et al.* 2011, abstract only). The copy-number round (IMPROVEMENTS plan U) sets the ribosomes from the cryo-ET density. |
| Geometry | A composite of two segmented cells, scaled to 34 µm³. The vacuole is one ellipsoid with the measured shape statistics of G1 cells, not a segmented shape. Organelle volume fractions are from different studies and strains (section 3.1). The cell volume follows from the protein total and an assumed cytosolic density; dry-mass measurements favour a somewhat larger cell and the cryo-ET ribosome concentration a smaller one (section 10.3). Procedural mitochondria (two tubular networks in build S) built to measured dimensions of glucose-grown cells but not to any segmented cell; nucleoid positions relative to cristae are by analogy with mammalian cells. No bud, so no cell-cycle-dependent structures. |
| Mitochondria and small compartments | Mitochondria are 3% of the cell, against 4.77% measured in glucose (Di Bartolomeo *et al.* 2020, another strain). Matrix 98.5% filled (178 mg/mL; build R 75%), intermembrane space 99.0% (build R 87%). The inner membrane, 13.1 µm² with short cristae, holds 87.4% of its proteins (carriers 37–42%; section 6.7), and only 1,940 ATP synthase dimers fit on the crista rims. The nucleolus is thinned to 0.907 of its inventory before packing and holds 84.1% of it (build R 80%). The Golgi cytosolic face keeps 0.61 of its requested copies, the rest is packed in the cytosol (section 4.4). Mitochondrially encoded subunits are inferred from their complexes. |
| Packing artefacts | Containment: 6.7 × 10⁻⁵ of sampled beads lie outside their compartment in build S, against 4.1 × 10⁻⁵ in build R (round S candidate build 1: 9.1 × 10⁻⁵); the rise has not been analysed. Cytosolic proteins are depleted to 0.16 of their density in the 4–20 nm next to Golgi and endosome membranes, where the membrane-face shell is (section 8.3). |
| Disorder | AlphaFold residues with pLDDT < 50 are deleted rather than represented, except the FG domains of 4 Nups. The simulation restores their effect on the hydrodynamic radius only (section 14.2). |
| Condensates | None are placed. A test of whether P-bodies emerge from generic interaction data was negative for specificity (section 9.1). |
| Localization | Agrees with independent N-terminal SWAT calls for 67–74% of proteins and with spatial proteomics for 85% (section 9.1); the references agree with each other for about 74%. Bulk steady state. Measured shares (CYCLoPs) only between compartment groups and only where imaging and annotation agree; within groups, and for proteins without imaging data, dual localizations are weighted by rank. 625 proteins default to the cytosol. |
| Nucleolus | Active rRNA genes are self-avoiding walks, not the real tandem array; nascent transcripts are free threads rather than compact RNPs and fit only partly (57%). The number of active genes is per genome, so the smaller nucleolus of builds R and S (0.47 µm³) is more crowded. One structure per pre-ribosome stage. No 5S rDNA (Pol III) or snoRNA genes. tRNA genes are not gathered at the nucleolus (section 6.3; the FISH test fails, 23% against 52% in build S). |
| Cytoskeleton | No septins, fimbrin or formins. The spindle-pole-body core has its measured stoichiometry but no measured copy number (section 6.2). Actin partition and filament lengths are assumptions. 0.25–0.5% of the beads of the built actin structures overlap a neighbouring built component (section 9). Only 2 of 4 intended non-kinetochore nuclear microtubules and 15 of 16 kinetochore microtubules are built (section 6.2). |
| Plasma membrane | Eisosome coats hold only 7.6% of Pil1/Lsp1; only Sur7 and Nce102 are concentrated in the furrows. The endocytic coat is partial (no clathrin) and randomly oriented, and invagination depths are assumed uniform. |
| Endocytic sites | The sites lack the ribosome exclusion zone around the actin network (181 ± 41 nm across, reaching about 100 nm beyond the tip; Picco *et al.* 2018, Buser & Drubin 2013). Invaginations are shallow: depths are drawn up to 100 nm, and the deepest in round S candidate build 1 was 92 nm, against about 100 nm on average and up to 140 nm at scission by tomography (Kukulski *et al.* 2012). The tube tips are too wide: 40 nm across, against a measured tip radius of 11.6 ± 2.6 nm, about 23 nm across (Kukulski *et al.* 2012 data in Picco *et al.* 2018). Actin per site is half the measured peak (section 6.4). |
| Nuclear pore basket | Not built. The only deposited basket model (8ZZC, integrative) resolves about a quarter of each Mlp1/Mlp2 as coarse beads, and their AlphaFold models are compact rather than the long coiled coils of the basket filaments. Mlp1, Mlp2 and Pml39 are packed free next to the nuclear envelope. |
| Contact sites | Only the nucleus–vacuole junction. The ER here is cortical and perinuclear, so ER–mitochondria (ERMES) and other cytoplasmic-ER contacts have no place. |
| Lipids and wall | Lipid compositions from isolated organelle membranes (1991), partly from non-glucose cultures; only the plasma membrane has leaflet asymmetry. Procedural cell wall. One bud scar (rim height and chitin share assumed), no birth scar; lateral-wall chitin not drawn. |
| mRNA | The mRNA between ribosomes is drawn shorter than its nucleotides, and the tails carry the rest (section 6.5). Ribosomes carry Stm1 (the idle-ribosome structure 4V88) and no tRNAs or elongation factors. Too few mRNAs are free of ribosomes: 570 of 18,136 (3%), while polysome profiling finds a median of 72% of a gene's mRNA on ribosomes (Arava *et al.* 2003, Supplementary Table 3; compared on build R, 2026-11-13). With 18,100 mRNAs a typical gene has few copies, so per-gene localization tests are noisy (mitochondrial mRNA test 0.45 / 0.66 in build S, 0.35 / 0.61 in build R with the same mRNAs; section 6.5). In a simulated box, releasing polysome ribosomes speeds 40 nm particles by 2% (measured 40–50%; section 14). |
| Missing classes | Metabolites, ions, water, glycogen, trehalose; polyphosphate only procedurally (amount assumed). |

## 12. Data sources and references

### 12.1 Data

| Source | Use | Terms |
|---|---|---|
| RCSB Protein Data Bank | Experimental structures | CC0 |
| AlphaFold Protein Structure Database | Monomers and predicted homo- and heterodimers (Jumper *et al.* 2021, Nature 596:583; Varadi *et al.* 2024, Nucleic Acids Res. 52:D368) | CC-BY 4.0 |
| PaxDb | Integrated whole-organism protein abundances for *S. cerevisiae* | |
| UniProt | Proteome UP000002311: sequences, subcellular location, transmembrane segments, topology, subunit structure | |
| Gene Ontology | Cellular-component annotations | |
| Complex Portal (EMBL-EBI) | Curated complexes for taxon 559292 | |
| IntAct | Physical interactions, to decide heterodimer binding | |
| Alliance of Genome Resources | Paralog pairs | |
| AlphaFold 3 pool predictions (in-house, unpublished) | Pair scores for partner placement | not public |
| STRING v12 | Functional and physical associations, for network-propagation profiles and loose partners | CC-BY 4.0 |
| CYCLoPs (Koh *et al.* 2015) | LOC-scores of the three wild-type GFP-collection screens (thecellvision.org/cyclops), for compartment shares | |
| Vögtle *et al.* 2017 | Submitochondrial localization of 810 proteins (Supplementary Data 5) | CC-BY 4.0 |
| Morgenstern *et al.* 2017 | High-confidence mitochondrial proteome, enrichment classes, submitochondrial classes (Table S1) | |
| Zinser *et al.* 1991 | Phospholipid classes and ergosterol of isolated organelle membranes (Tables 4, 5) | |
| EMPIAR-11416 (cryo-plasma FIB/SEM of yeast), segmented into cell, vacuole, lipid-droplet and nucleus meshes | The cell outline and lipid droplets (cell 08), the position of the vacuole (cell 08's vacuoles; their shapes up to build R) and the nucleus (cell 05) of the model cell (section 3.3) | CC0 |
| Chen *et al.* 2025, Supplementary file S8 (soft X-ray tomography of unbudded G1 BY4741 in YPD) | Number, sphericity and ellipsoid axis ratios of the vacuole (section 3.3) | |
| EMPIAR-12534 (in-cell cryo-ET of BY4741 in YPD; Chang *et al.* 2025), 50 segmented membrane meshes | Outer-membrane diameter, membrane spacing, crista area, density and depth, measured here on the YPD meshes (`pipeline/data/validation/literature/mitochondria/`; section 3.3) | |

### 12.2 Literature behind the parameters

| Reference | Used for |
|---|---|
| Carroll KM *et al.* (2011). *Mol. Cell. Proteomics* 10:M111.007633; Sherman F (2002). *Methods Enzymol.* 350:3–41 | Protein mass per cell (4–5 pg; 6 pg) |
| Ho B, Baryshnikova A, Brown GW (2018). *Cell Systems* 6:192–205 | Earlier protein total (molecules per cell) |
| Jorgensen P *et al.* (2002). *Science* 297:395–400 | Cell volume |
| Uchida M *et al.* (2011). *Yeast* 28:227–236 | Organelle volume fractions (soft X-ray tomography) |
| Dumoux M *et al.* (2023). *eLife* 12:e83623; EMPIAR-11416 | Cryo-plasma FIB/SEM volume of yeast; segmented cells used for nucleus, vacuole and lipid-droplet characteristics |
| Jorgensen P *et al.* (2007). *Mol. Biol. Cell* 18:3523–3532 | Nuclear-to-cell volume ratio |
| Yamaguchi M *et al.* (2011). *Microscopy* 60:321–335, doi:10.1093/jmicro/dfr052 | Serial-section EM comparison |
| Delarue M *et al.* (2018). *Cell* 174:338–349.e20, PMID 29937223 | Ribosome concentration by cryo-ET (about 14,000 per µm³ of cytoplasm) |
| Goetz SK *et al.* (2025). bioRxiv preprint, doi:10.64898/2025.12.22.696049 (Europe PMC PPR1230379; not peer reviewed) | Ribosome concentration by cryo-ET (14,517 ± 3,193 per µm³ of cytosol) |
| Orlean P (2012). *Genetics* 192:775–818 | Cell-wall thickness and composition |
| Wolf W *et al.* (2012). *PLoS ONE* 7:e39703 | Cortical ER coverage |
| Winey M *et al.* (1997). *Mol. Biol. Cell* 8:2119–2132 | NPC number |
| Ejsing CS *et al.* (2009). *PNAS* 106:2136–2141 | Lipid composition |
| Huh WK *et al.* (2003). *Nature* 425:686–691 | GFP localization |
| Kefauver JM *et al.* (2024). *Nature* 632:664–671 | Pil1/Lsp1 eisosome lattice on native membranes (8QB7, 8QB8) |
| Lynch-Day MA, Klionsky DJ (2010). *FEBS Lett.* 584:1359–1366; Yuga M *et al.* (2011). *J. Biol. Chem.* 286:13704–13713 | Cvt pathway; Cvt cargo (Ape1, Ams1, Ape4) reaches the vacuole lumen without a signal peptide |
| Laporte D *et al.* (2008). *J. Cell Biol.* 181:737–745 | 26S proteasome nuclear in growing cells |
| Bystricky K *et al.* (2004). *PNAS* 101:16495–16500 | Chromatin compaction and persistence in living yeast (fibre) |
| Arbona JM *et al.* (2017). *Genome Biol.* 18:81 | Chromatin compaction and persistence from imaging and Hi-C (fibre radius choice) |
| Pelechano V, Chávez S, Pérez-Ortín JE (2010). *PLoS One* 5:e15442 | mRNA molecules per gene (Table S1); engaged Pol II |
| Arava Y *et al.* (2003). *PNAS* 100:3889–3894 | Ribosome association of yeast mRNAs (test, not used; ribosome-free mRNAs, section 11) |
| Chereji RV, Ramachandran S, Bryson TD, Henikoff S (2018). *Genome Biol.* 19:19 | H3 chemical-cleavage dyad maps: measured nucleosome positions of the fibre |
| Rossi MJ *et al.* (2021). *Nature* 592:309–314 | Promoter nucleosome-depleted regions (NDR filter); ChIP-exo sites of chromatin proteins |
| Oberbeckmann E *et al.* (2019). *Genome Res.* 29:1996–2009 | Absolute nucleosome occupancy and count (held-out test) |
| Weiner A *et al.* (2015). *Mol. Cell* 58:371–386 | MNase nucleosome positions, unstressed (held-out test) |
| Schep AN *et al.* (2015). *Genome Res.* 25:1757–1770 | ATAC-accessible regions (held-out test) |
| Brogaard K, Xi L, Wang JP, Widom J (2012). *Nature* 486:496–501 | Chemical nucleosome map (count cross-check; sacCer2, not used) |
| Eaton ML *et al.* (2010). *Genes Dev.* 24:748–753 | Nucleosome-free regions at origins (check of the dyad calls) |
| Siow CC *et al.* (2012). *Nucleic Acids Res.* 40:D682–D686 | OriDB confirmed replication origins (MCM double hexamers, ORC) |
| Joglekar AP *et al.* (2006). *Nat. Cell Biol.* 8:581–585 | Kinetochore protein copy numbers |
| Marcand S, Gilson E, Shore D (1997). *Science* 275:986–990 | Rap1 sites every 18 bp on ~300 bp telomeric tracts (17 Rap1 per end) |
| Swygert SG *et al.* (2014). *Nat. Commun.* 5:4751; Swygert SG *et al.* (2018). *PNAS* 115:12447–12452 | Two Sir3 and one to two Sir2–Sir4 per nucleosome on reconstituted arrays (upper bounds in cells) |
| Muller EG *et al.* (2005). *Mol. Biol. Cell* 16:3341–3352 | SPB core stoichiometry Spc42:Spc29:Spc110:Cmd1:Cnm67 = 2:1:1:1:1 (as used by Viswanath *et al.* 2017) |
| Haeusler RA *et al.* (2008). *Genes Dev.* 22:2204–2214 | tRNA-gene clusters at the nucleolus in wild type (test) |
| Ikeuchi K *et al.* (2019). *EMBO J.* 38:e100276 (PDB 6I7O); Tesina P *et al.* (2020). *EMBO J.* 39:e103365 (PDB 6T83) | Collided-disome pose (6I7O used; 6T83 within 15° and 2.9 nm) |
| Meydan S, Guydosh NR (2020). *Mol. Cell* 79:588–602 | Disome footprints of 57–63 nt: at most about 5 nt between collided ribosomes |
| Zhao T *et al.* (2021). *Genome Biol.* 22:16 | 5.8% of ribosomes in collided disomes in exponential yeast (test) |
| Chan PP, Lowe TM (2016). *Nucleic Acids Res.* 44:D184–D189 | GtRNAdb mature tRNA sequences (42 isoacceptors) |
| Percudani R, Pavesi A, Ottonello S (1997). *J. Mol. Biol.* 268:322–330 | tRNA levels follow tRNA gene copy number in growing yeast |
| Passaro S *et al.* (2025). bioRxiv 10.1101/2025.06.14.659707 | Boltz-2, used to predict the tRNA structures |
| Kim SJ *et al.* (2018). *Nature* 555:475–482 | Copies of each nucleoporin per pore (Supplementary Table 2A), for the nucleoporin bookkeeping and the pore anchor |
| von der Haar T (2008). *BMC Syst. Biol.* 2:87 | 187,000 ± 56,000 ribosomes per cell (check of the protein scale) |
| Brothers M, Rine J (2022). *eLife* 11:e75653 | Sir3 ChIP-seq of cycling cells: which nucleosomes are Sir-bound |
| Goodnight D, Rine J (2020). *eLife* 9:e58910 | Sir4 ChIP-seq with untagged control (held-out test of the Sir-bound nucleosomes) |
| Teytelman L, Thurtle DM, Rine J, van Oudenaarden A (2013). *PNAS* 110:18602–18607 | Highly expressed genes give false ChIP signal (why Sir-bound nucleosomes are limited to chromosome ends) |
| Zhong J *et al.* (2016). *Genome Res.* 26:351–364 | Per-base DNase I cuts in log-phase yeast (held-out test of nucleosome positions) |
| Zelin E *et al.* (2012). *Mol. Cell* 48:459–470 | DNase-hypersensitive sites of wild-type BY4741 (held-out test) |
| Dubreuil B *et al.* (2019). *Nucleic Acids Res.* 47:D1245–D1249; Meurer M *et al.* (2018). *Nat. Methods* 15:598–600 | YeastRGB images of the C-SWAT mNeonGreen library with a nuclear-envelope marker (nucleus/cytoplasm test on raw images) |
| Das H *et al.* (2026). *Histochem. Cell Biol.* 164:8 | Absolute proteome of glucose-grown yeast (copy-number check) |
| Gade VR *et al.* (2025). *Cell Rep.* 44:116204 (preprint 10.1101/2024.11.14.623620); Kunzi M *et al.* (2026), bioRxiv 10.64898/2026.08.29.747980 | 40 nm GEM diffusion in budding yeast (test); cytoplasmic dry mass and protein concentration (check) |
| Ohsawa S *et al.* (2026). bioRxiv 10.64898/2026.09.03.749057 | Nucleus 5.3–6.3% of the cell in control strains (check) |
| Pelechano V, Wei W, Steinmetz LM (2013). *Nature* 497:127–131 | Per-gene UTR lengths (TIF-seq, YPD) |
| Subtelny AO *et al.* (2014). *Nature* 508:66–71 | Per-gene poly(A) tail lengths (PAL-seq) and translational efficiency |
| Williams CC, Jan CH, Weissman JS (2014). *Science* 346:748–751 | mRNAs translated at the mitochondrial surface (−CHX calls) |
| Jan CH, Williams CC, Weissman JS (2014). *Science* 346:1257521 | ER-enriched secretome genes (Sec63 proximity ribosome profiling) |
| Chartron JW, Hunt KCL, Frydman J (2016). *Nature* 536:224–228 | Membrane-polysome enrichment per gene (held-out test of ER mRNAs) |
| Bullitt E, Rout MP, Kilmartin JV, Akey CW (1997). *Cell* 89:1077–1086; Viswanath S *et al.* (2017). *Mol. Biol. Cell* 28:3298–3314 | Spindle pole body layers (cryo-EM; integrative model of the core) |
| Thompson M, Haeusler RA, Good PD, Engelke DR (2003). *Science* 302:1399–1401 | tRNA genes cluster at the nucleolus (check; Hi-C does not show it) |
| Saint-Georges Y *et al.* (2008). *PLoS One* 3:e2293 | Mitochondrion-localized mRNAs (held-out test) |
| Chen J *et al.* (2025). *Mol. Biol. Cell* 36:ar132; Di Bartolomeo F *et al.* (2020). *PNAS* 117:7524–7535 | Organelle volumes by soft X-ray tomography (vacuole fraction; since build S also the single near-spherical vacuole, Supplementary file S8); mitochondrial volume in glucose (checks) |
| Barbosa AD *et al.* (2015). *Mol. Biol. Cell* 26:3641–3657; Santana-Sosa S, Matos-Perdomo E, Ayra-Plasencia J, Machín F (2023). *Int. J. Mol. Sci.* 24:9829; Yang CH, Lambie EJ, Hardin J, Craft J, Snyder M (1989). *Chromosoma* 98:123–128 | Nucleolus at the envelope next to the vacuole; nucleolus opposite the spindle pole body |
| Severs NJ, Jordan EG, Williamson DH (1976). *J. Ultrastruct. Res.* 54:374–387; Kohler V, Büttner S (2021). *Contact* 4:25152564211016608; Otto FB, Thumm M (2020). *Int. J. Mol. Sci.* 21:4506 | No nuclear pores in the nucleus–vacuole junction (Severs 1976 and Pan 2000, as cited by the two reviews) |
| Egner A, Jakobs S, Hell SW (2002). *PNAS* 99:3370–3375 | Mitochondrial tubule diameter in glucose (live 4Pi microscopy) |
| Chang YT *et al.* (2025). *J. Cell Biol.* 224:e202407110 (EMPIAR-12534) | Mitochondrial membrane spacing and cristae in BY4741 in YPD (cryo-ET; meshes measured here) |
| Stoldt S *et al.* (2019). *PNAS* 116:9853–9858; Klecker T, Westermann B (2021). *Open Biol.* 11:210238 | Crista shape and placement (short, sheet-like, on two opposing sides); few cristae in glucose; ATP synthase per mitochondrial mass in glucose |
| Osman C *et al.* (2015). *PNAS* 112:E947–E956; Göke A *et al.* (2020). *Mol. Biol. Cell* 31:527–545; Seel A *et al.* (2023). *Nat. Struct. Mol. Biol.* 30:1549–1560; Kukat C *et al.* (2015). *PNAS* 112:11288–11293 | Nucleoids per µm of tubule and their spacing; mtDNA copies per cell and per fL; separate networks per cell; nucleoids in crista-free matrix (mammalian) |
| Metzl-Raz E *et al.* (2017). *eLife* 6:e28034; Björkeroth J *et al.* (2020). *PNAS* 117:21804–21812; Xia J *et al.* (2022). *Nat. Commun.* 13:2819 | Ribosomal-protein share of the proteome in fast-growing yeast (limitation, section 11) |
| Picco A *et al.* (2018). *Mol. Biol. Cell* 29:1346–1358; Buser C, Drubin DG (2013). *Microsc. Microanal.* 19:381–392 | Endocytic tip radius and ribosome exclusion zone (limitation, section 11) |
| Strádalová V *et al.* (2009). *J. Cell Sci.* 122:2887–2894; Karotki L *et al.* (2011). *J. Cell Biol.* 195:889–902 | Eisosome furrow dimensions and number |
| Bharat TAM, Hoffmann PC, Kukulski W (2018). *Structure* 26:879–886 | In-situ eisosome coat diameter |
| Kukulski W *et al.* (2012). *Cell* 150:508–520; Kukulski W *et al.* (2016). *eLife* 5:e16036 | Endocytic invagination depth at scission, vesicle size, tip radius |
| Pillet F *et al.* (2014). *BMC Biol.* 12:6 | Bud-scar diameter and profile (AFM) |
| Gotta M *et al.* (1996). *J. Cell Biol.* 134:1349; Duan Z *et al.* (2010). *Nature* 465:363; Therizols P *et al.* (2010). *PNAS* 107:2025 | Telomere foci; Rabl layout and arm-length dependence |
| Pelechano V, Chávez S, Pérez-Ortín JE (2010). *PLoS One* 5:e15442 | Engaged Pol II per gene |
| West M *et al.* (2011). *J. Cell Biol.* 193:333 | Ribosome density on the cortical ER |
| Picco A *et al.* (2015). *eLife* 4:e04535 | Molecules per endocytic site; Rvs161/167 assembly on the invagination |
| Pan X *et al.* (2000). *Mol. Biol. Cell* 11:2445; Jeong H *et al.* (2017). *PNAS* 114:E4539 | Nucleus–vacuole junction (no pores there, section 6.1); Vac8–Nvj1 structure (5XJG) |
| Kisonaite M *et al.* (2023). *Nat. Struct. Mol. Biol.* 30:670 | RAC on the 80S ribosome (7Z3N) |
| Baba M *et al.* (1997). *J. Cell Biol.* 139:1687 | Cvt vesicle size |
| Bru S *et al.* (2025). *eLife* 14:RP108181; Christ JJ *et al.* (2020). *Biotechnol. Bioeng.* 117:2089 | Vacuolar polyphosphate store; polyP chain length |
| Klis FM *et al.* (2002). *FEMS Microbiol. Rev.* 26:239–256 | Cell-wall composition (chitin 1–2%) |
| Sickmann A *et al.* (2003), *PNAS*; Reinders J *et al.* (2006), *J. Proteome Res.* | Excluded high-throughput organelle calls |
| Dolgalev G *et al.* (2023). *Int. J. Mol. Sci.*, doi:10.3390/ijms24032081 | Protein density comparison |
| Kostinski S, Reuveni S (2021). *Phys. Rev. Research*, doi:10.1103/physrevresearch.3.013020 | Ribosome count comparison |
| Warner JR (1999). *Trends Biochem. Sci.* 24:437–440 | Ribosome production rate; about 200,000 ribosomes per cell |
| French SL, Osheim YN, Cioci F, Nomura M, Beyer AL (2003). *Mol. Cell. Biol.* 23:1558–1568 | Active rRNA genes, Pol I per gene, elongation rate (Miller spreads) |
| Dammann R *et al.* (1995). *Mol. Cell. Biol.* 15:5294–5303; Merz K *et al.* (2008). *Genes Dev.* 22:1190–1204 | Active rDNA is nucleosome-depleted (psoralen accessibility; Hmo1 chromatin) |
| Kos M, Tollervey D (2010). *Mol. Cell* 37:809–820 | Co-transcriptional pre-rRNA cleavage |
| Barandun J *et al.* (2017), Cruz VE *et al.* (2022), Kater L *et al.* (2020), Malyutin AG *et al.* (2017), Scaiola A *et al.* (2018); Neyer S *et al.* (2016), *Nature* 540:607 | Pre-ribosome structures; Pol I elongation complex (5M3F) |
| Maritan M *et al.* (2022). *J. Mol. Biol.* 434:167351 | Whole-cell modelling approach |
| Weill U *et al.* (2018). *Nat. Methods* 15:617–622 | N-terminal SWAT GFP localization (independent validation, 9.1) |
| Nightingale DJH, Oliver SG, Lilley KS (2019). *Curr. Opin. Chem. Biol.* 48:86–95 | hyperLOPIT spatial proteome of yeast (validation, 9.1) |
| Gallego O *et al.* (2010). *Mol. Syst. Biol.* 6:430 | Protein–lipid arrays (membrane-anchor test, 9.1) |
| Beckmann BM *et al.* (2015). *Nat. Commun.* 6:10127; Matia-González AM *et al.* (2015). *Nat. Struct. Mol. Biol.* 22:1027–1033; Mitchell SF *et al.* (2013). *Nat. Struct. Mol. Biol.* 20:127–133 | RNA-interactome capture (condensate and diffusion tests, 9.1) |
| Rostam N *et al.* (2023). *Nat. Methods* 20:673–676 | CD-CODE condensate membership (condensate test, 9.1) |

### 12.3 PDB entries used for assemblies and templates

| Assembly or template | PDB entries |
|---|---|
| Ribosome and polysomes | 4V88, 7ZW0 |
| Pre-ribosomes and rDNA transcription | 5WLC, 7R7A, 6YLH, 5T62, 6FAI; Pol I elongation complex 5M3F |
| ER translocons | 2WW9, 7KAH |
| Nuclear pore and nucleosome | 7N9F, 1ID3 |
| Microtubules | 8QV3, 5W3F |
| Actin and actin-binding proteins | 9GO5, 5JLF, 8UXW, 8F8Q |
| Mitochondrial inner membrane | 8F39, 4B2Q, 6HU9; mitoribosome 9HM0 |
| Eisosomes | 8QB7, 8QB8 |
| RAC on the ribosome; Pol II elongation complex; Cvt cargo; nucleus–vacuole junction | 7Z3N; 1Y1W; 5JGF (Ape1 dodecamer); 5XJG |
| Other curated complexes | 4CR2, 2UV8, 4V94, 4C2M, 1WCM, 5FJ8, 7VH5, 6O7V |

The other complexes and homo-oligomers are listed per species in the viewer's inspect panel: click any molecule.

## 13. Reproducibility and data access

- **Pipeline.** The model is built by a Python pipeline (numba-accelerated) in six stages: the base structure library (stage 0), inventory, geometry, packing, validation and encoding (`scripts/run_all.sh`; `scripts/README.md` lists the stages). Build S took 87 minutes from stage 0 to the encoded tiles on 40 workers of a 64-core cluster node, 55 of them for the tile packing (all six rounds, the last with the overlap-removal squeeze); the comparisons with data (`scripts/validate_all.sh`, which needs the GPUs) took another 17 minutes. The real shapes come from `pipeline/cell_meshes.npz` (in the repository). The organisation layer of the packing (section 8.4) needs STRING (public) and the pool pair scores (not public; the derived partner lists are reproducible only with them). `work/interactions/pair_summary.parquet` (the pool pairs read by `features.py`) is an input kept in `work/` without a generator in the repository.
- **Stage 0, the base structure library (since 2026-11-13; `scripts/s0_library.py`).** The curated complexes from their RCSB assemblies (membrane complexes oriented), one AlphaFold DB model per protein of the inventory (residues with pLDDT < 50 left out, transmembrane proteins oriented from the UniProt topology), the curated PDB homo-oligomers and tRNA-Phe (1EHZ). Until build R these entries were made only by cells of a historical notebook, so the library on disk was state that no script could rebuild. Rebuilt and compared entry by entry, 5,561 of 5,647 shared entries were identical; the differences were the 82 transmembrane proteins oriented before the side rule (section 5.3) and the missing H2B entries (section 2.3). Build S was built with stage 0 (7,091 library entries); the amphiphysin heterodimer is now made by its builder (section 6.4). The library path can be set with `YC_LIB`.
- **Isolated builds (`scripts/build_isolated.sh NAME [GIT_REF] [FIRST_STAGE]`).** A build of one commit runs from a frozen git worktree of that commit in its own run folder (`../runs/NAME/`: `tree/` the worktree, `lib/` its own copy of the structure library, `work/` its own output, seeded only with download caches, and `build.log`), so code can be edited while a build runs and builds never share output or library. `YC_RESUME=1` continues a run from a later stage after a fix (the build is then marked as mixed); `YC_LIB_FROM=empty` starts from an empty library. Build S is the run `runs/S3_seed0` of commit 9568da3.
- **Master seed (`YC_SEED`, `config.seed`).** Every builder has its own random stream, numbered by its historical seed (tiles: 1000 plus the tile index). At master seed 0 the stream number itself is the seed, which reproduces the reference build; any other `YC_SEED` shifts all streams together, giving an independent replicate of the cell. Streams that make the structure library, the encoding or validation samples do not depend on it. Builds run with `PYTHONHASHSEED=0`.
- **Build stamp (`scripts/build_stamp.py`, `work/build.json`).** At the start of a build: build id, commit, whether code under `pipeline/`, `scripts/` or `validation/` was uncommitted, seed, start time, host, workers, every module-level constant of `pipeline/*.py`, and size and modification time of the inputs. After the validation run, `build_stamp.py check` lists result files older than the build and writes `work/validation/MANIFEST.json` (file, build id, age); for build S all 30 result files are of this build. With several builds a day, numbers of an earlier build had reached this document twice (2026-11-12).
- **Determinism.** The build is deterministic: the same commit, master seed, inputs and structure library give the same cell, independent of the number of workers and of the hash seed. Checked on 2026-11-13: `main` (build R) rebuilt by `build_isolated.sh` with 56 instead of 64 workers was identical to build R in all 5,912 tiles (24,484,185 instances at the same positions), all 128 pre-pass arrays and the expected totals; the seed plumbing at master seed 0 reproduced build R in all 5,912 tiles as well. Replicates at other seeds have not been built yet (section 11).
- **Source code.** The pipeline and the viewer source are in `source/`; `source/README.md` explains how to rebuild. The historical notebook and the eight modules only it used (the 2D column model and its figures) were removed in round S.
- **This repository** holds the built viewer and the encoded data.
  - The tiles under `data/tiles/` store every instance as a species index, a position quantised to about 0.0055 nm and an orientation quantised to about 0.004°.
  - The residue-bead library is in `data/library/`.
- **Full-atom coordinates.** The viewer fetches full-atom models live from RCSB and AlphaFold DB for display in Mol*.

## 14. Simulating a box of the cell (experimental)

The viewer can cut a cube of 30–200 nm out of the packed cell and run Brownian dynamics (BD) of its molecules on the GPU. While the simulation runs, the whole-cell rendering and tile streaming pause. Only the box and the membranes inside it are drawn. The aim is to watch the packed cell move on physically grounded diffusion. It is not to predict specific biochemistry. The model has excluded volume, a nonspecific attraction, screened electrostatics, specific binding at the interfaces of predicted complexes, polysome mRNA links and mean-field hydrodynamics, but no flexibility and no reactions.

The approach follows BD models of bacterial cytoplasm built from structures (McGuffee & Elcock 2010; Ando & Skolnick 2010) and of concentrated protein solutions (Mereghetti & Wade 2012; Balbo *et al.* 2013). All-atom MD of a cytoplasm (Yu *et al.* 2016) needs supercomputers for 20 ns. Rigid-body BD reaches microseconds of a 100 nm box in a browser.

Approaches considered:

| Approach | Resolution and time scale | Fit for a box of this cell in a browser |
|---|---|---|
| All-atom MD with explicit water (Yu *et al.* 2016) | atoms; tens of ns for ~100 nm | far beyond a browser GPU |
| Martini 3 coarse-grained MD (Souza *et al.* 2021); a whole minimal cell has been built in Martini (Stevens *et al.* 2023) | about 4 heavy atoms per bead, explicit CG water; ns–µs | the export (section 14.7) is the starting point for it; too heavy to run interactively |
| Elastic-network CG with implicit solvent (OpenMM custom forces) | one bead per residue; ns–µs | used for the exported OpenMM script |
| Rigid-body Brownian dynamics with implicit solvent and mean-field hydrodynamics (McGuffee & Elcock 2010; Mereghetti & Wade 2012) | residues or sphere trees; µs | **the model in the viewer** |
| Particle-based reaction–diffusion: Smoldyn (Andrews 2017), ReaDDy 2 (Hoffmann *et al.* 2019) | point particles; ms–s | no molecular shape, so it loses what the packed model provides |
| Lattice reaction–diffusion master equation: Lattice Microbes (Roberts *et al.* 2013), used for a whole minimal cell (Thornburg *et al.* 2022) | ~10 nm voxels; minutes of cell time | a natural next level for the whole cell; needs kinetic parameters this model does not have |

### 14.1 What is simulated

- **Molecules.** Every packed instance whose centre lies in the box becomes a rigid body with its packed position and orientation.
- **Schematic model.** Contacts use the library sphere tree of each species (about 30 residues per sphere). The default step is 0.2 ns, which reaches microseconds.
- **Mid model.** Contacts use clusters of about 4 residues (2 nucleotides): k-means of each species' beads, per kind, seeded by runs of the Morton-ordered beads. A cluster's radius is its bead radius plus 0.8 times the rms distance of its members from the centre. Measured along rays toward the centres of 13 species, the cluster surface then lies 0.13–0.22 nm outside the residue surface on average (0.3 nm for the ribosome), so molecules touch with a small gap instead of visibly overlapping as in the schematic model. Only clusters holding an exposed residue are kept, plus core spheres (below). The default step is 0.2 ns.
- **Residue model.** Contacts use one sphere per residue bead: radius 0.36 nm for protein, 0.6 nm for nucleotide. Only solvent-exposed beads are used, about 60% of them. A bead is exposed if a probe bead in contact along one of 26 directions overlaps none of the molecule's other beads; buried beads cannot meet another molecule before its surface beads do. The default step is 0.1 ns, which reaches hundreds of nanoseconds.
- **Core spheres (mid and residue models).** With surface spheres only, molecules are hollow shells. In a test, a strong attraction (mid model, 1.5 kT) pulled them through each other: grid cells filled far beyond their capacity, neighbours were missed, and molecules moved faster than in water. Each molecule therefore also carries core spheres: its sphere tree shrunk by 1 nm, about 1 nm inside the residue surface. Cores meet only other cores, with overlap repulsion and nothing else (no charge, attraction or contacts). They stop molecules passing through each other without adding to the cost of surface contacts. Keeping every buried bead instead (solid bodies) did the same job at about 1.5× the cost.
- **Drawing.** All models draw the residue beads.
- **Held fixed, acting only as obstacles:**
  - structures built explicitly before packing and anchored to filaments, membranes or DNA (section 6: filaments, lattices, assemblies; instance flag bit 1 from the encoder). Polysome ribosomes are the exception: they move;
  - nucleosomes and nuclear pores.
- **Membrane proteins** (membrane regions) diffuse laterally in their membrane (section 14.2), or stay fixed if that option is off.
- **Box faces.** Three boundary conditions:

  | Boundary | Behaviour | Use |
  |---|---|---|
  | Periodic | Every molecule moves. Centres are wrapped into the box and contacts use the nearest periodic image, as in standard MD of bulk solutions. | Default for boxes without membranes, such as cytoplasm or nucleoplasm. |
  | Embedded | The box is closed, but molecules in its outer 8 nm keep moving while tied to their starting positions by a harmonic tether (1 kT/nm², about 1 nm rms per axis). Molecules near the faces still see crowded surroundings, and the cell outside holds them in place. | Default for boxes containing membranes, which do not continue across the faces. |
  | Closed walls | Molecule centres reflect at the faces; parts of a molecule may stick out. | Isolated boxes. Beware the depletion layer next to the walls. |

  The first version froze an 8 nm shell completely, which immobilised about 60% of the molecules of a 60 nm box.
- **Membranes.** Membranes are walls. A grid over the box stores the signed distance to the nearest bilayer midplane, evaluated once from the same analytic membrane distance functions as the renderer. Soluble molecules are kept at least 2 nm (half a bilayer) plus their sphere radius from every bilayer, on their own side.
- **Polysomes.** Ribosomes of one polysome stay linked by their mRNA (section 14.3). The packer placed each next ribosome so that its mRNA entry lies 3–12 nm beyond the previous one's exit, within 0.7 rad of the exit direction, with its entry channel pointing back along the linker. Pairs of ribosomes in the box that meet all three conditions to the tile quantisation are linked. Links to ribosomes outside the box are dropped. A 200 nm cytoplasm box holds 125 ribosomes with 12 links joining 21 of them.
- **Not included as bodies.** mRNA threads, chromatin fibres beyond the nucleosomes, cell-wall polysaccharides and lipids are absent.

### 14.2 Diffusion coefficients

**Membrane proteins** diffuse in the tangent plane of their membrane. The plane comes from the gradient of the membrane distance grid at the centre. The protein is returned to its starting depth after every step, and it rotates only about the membrane normal. Its coefficients follow Saffman & Delbrück (1975) for a cylinder of radius *a* in a membrane of viscosity μ_m and thickness h between water of viscosity μ:

```
D = kT / (4π μ_m h) · (ln(μ_m h / (μ a)) − 0.5772)
D_r = kT / (4π μ_m h a²)
```

- **Parameters:** μ_m = 0.1 Pa·s for a fluid bilayer and h = 4 nm (assumptions). *a* is half the hydrodynamic radius, kept within 1–5 nm, as a rough size of the membrane-embedded part. D depends on *a* only logarithmically: about 4 µm²/s for *a* = 1.5 nm.
- **Calibration per membrane (2026-11-11; `pipeline/membrane_diffusion.py` → `library/membrane_diffusion.json`, viewer option "calibrated per membrane", on by default).** The fluid-bilayer value is multiplied by a factor per membrane: the median of measured D ÷ the Saffman–Delbrück D of the same protein (`pipeline/data/validation/membrane_diffusion_measurements.tsv`; FRAP of transmembrane proteins in unperturbed cells):
  - plasma membrane 5.9 × 10⁻⁴ (GFP-Sso1 and GFP-Snc1, 0.0025–0.0029 µm²/s; Valdez-Taubas & Pelham 2003);
  - vacuole 0.016 (an Sso1 anchor retargeted to the vacuole, 0.068 µm²/s; same paper);
  - ER 0.035 (Sec61 and Heh2-derived reporters, 0.15–0.18 µm²/s; Meinema *et al.* 2013);
  - nuclear envelope 0.079 (0.09–0.92 µm²/s; same paper);
  - mitochondrial membranes, Golgi/endosomes and peroxisomes: no yeast data; the ER factor is used (assumption).
- **Caveats of the calibration:** FRAP measures micrometre-scale mobility over seconds, including domains, obstacles and large immobile fractions (40–80% on the plasma membrane; Spira *et al.* 2012), and the plasma-membrane values are approximate in the source (D = r²/4t½). Applying them as one factor to the short-time motion of a 30–200 nm box makes membrane proteins as slow on nanosecond scales as on micrometre scales, which is the same compromise as the cytoplasm factor (section 14.5). Immobile fractions are not modelled.
- **Not modelled:** membrane domains as such, the cell wall pressing on the plasma membrane, and hydrodynamic coupling within the membrane.
- **Measured motion:** in a 100 nm cortical-ER box over 5 µs, membrane proteins moved a median 2.7 nm and free soluble molecules 18 nm (fluid bilayer, before the calibration). In an 80 nm cortical-ER box the median membrane-protein D is 4.3 µm²/s uncalibrated and 0.15 µm²/s calibrated (mostly ER proteins).

**Hydrodynamic radius.** R_h of each species comes from its residue beads. The Kirkwood double sum 1/R_K = ⟨1/r_ij⟩ runs over bead pairs, subsampled to 700 beads. It is scaled by 6/5, the value for a uniformly filled sphere, where R_K = 5R/6, and a 0.3 nm hydration shell is added. Elongated species get a larger R_h than globular species of the same mass.

**Trimmed disordered regions** (2026-11-11). The library leaves out AlphaFold residues with pLDDT < 50, which makes IDR-rich proteins too small hydrodynamically. For each AlphaFold monomer species, `features.py` regrows the trimmed IDRs (runs of ≥ 10 residues) as self-avoiding C-alpha coils from their anchor residues (virtual bond angle 2.15 rad, chosen so that free coils of 50–400 residues follow the IDP scaling laws R_g ≈ 0.254 *N*^0.522 nm and R_h ≈ 0.249 *N*^0.509 nm within about 10%; Bernadó & Blackledge 2009, Marsh & Forman-Kay 2010; the exponents are (unchecked) against the papers' text). The ratio of the full chain's Kirkwood R_h to the folded part's (`library/rh_scale.json`; median 1.1, up to 3) multiplies the species' R_h in the simulation.

**Diffusion in water** at 30 °C (η = 0.797 mPa·s) follows Stokes–Einstein–Debye:

```
D_t = kT / (6π η R_h)
D_r = kT / (8π η R_h³)
```

The Viscosity setting divides both coefficients. It stands for a fluid phase more viscous than water.

**Mean-field hydrodynamic interactions** (on by default) follow Mereghetti & Wade (2012) and Balbo *et al.* (2013). Every molecule's coefficients are multiplied by the short-time self-diffusion factors of hard spheres at the box's macromolecular volume fraction φ, using second-order expansions that include lubrication (Cichocki, Ekiel-Jeżewska & Wajnryb 1999):

```
D_t / D_t0 = 1 − 1.832 φ − 0.219 φ²
D_r / D_r0 = 1 − 0.631 φ − 0.726 φ²
```

- **Volume fraction.** φ is the hydrated volume of all molecules in the box (partial specific volume 0.73 mL/g plus 0.35 g/g bound water) divided by the box volume. Membranes and the lumen of other compartments count as free volume, so boxes that cut through organelles underestimate φ.
- **Typical values.** Cytoplasm boxes give φ ≈ 0.2, which makes D_t 0.63 times its value in water.
- **What it leaves out.** Explicit many-body hydrodynamics. Ando & Skolnick (2010) used Stokesian dynamics with lubrication on equivalent spheres at 250–350 mg/mL (sphere volume fractions 0.43–0.60) and found GFP at D_L/D₀ ≈ 0.08 at 300 mg/mL. That density is far above this model's (cytoplasm boxes 160–185 mg/mL, R_h-sphere volume fraction 0.32), so their number does not carry over.
- **Far-field many-body hydrodynamics tested** (2026-11-11, `simgpu/`): the Ewald-summed Rotne–Prager–Yamakawa mobility of all molecules of a 120 nm cytoplasm box (polydisperse, with the overlap correction of Zuk *et al.* 2014; Rotne & Prager 1969; periodic sum after Hasimoto 1959) slowed GFP-sized molecules to 0.56 of their speed in water, against 0.62 with no hydrodynamic interactions and 0.36 with the mean-field factors. Pairwise far-field hydrodynamics leaves each molecule's own mobility at its free value; the short-time slowdown comes from the near field (lubrication), which is ill-defined for the overlapping hydrodynamic spheres of shaped molecules (0.6 overlapping neighbours per molecule in the packed box). The mean-field factors stay the default.

### 14.3 Interactions and integration

**Overlap.** A harmonic penalty applies between spheres of different bodies: U = ½ k (r_i + r_j − d)² for d < r_i + r_j, with k = 20, 30 and 50 kT/nm² in the schematic, mid and residue models.

**Stickiness.** A short-range nonspecific attraction of depth ε acts per sphere pair:

```
U = −ε cos²(π (d − r_i − r_j) / (2 · 0.6 nm))   for 0 ≤ d − r_i − r_j ≤ 0.6 nm
```

McGuffee & Elcock (2010) needed such an attraction to reproduce GFP diffusion in *E. coli*. It is applied per sphere pair, so the effective well between two molecules is several ε, and its effect differs between models (section 14.5).

**Electrostatics** (on by default). Charged spheres interact by the screened Coulomb (Debye–Hückel) potential of two charged spheres of radii a_i, a_j:

```
U = l_B q_i q_j exp(−κ (d − a_i − a_j)) / ((1 + κ a_i)(1 + κ a_j) d)      (kT units; d ≥ a_i + a_j, else its contact value)
```

- **Bjerrum length** l_B = 0.72 nm (water at 30 °C).
- **Screening.** κ = 3.30 √I nm⁻¹ for ionic strength I in M. The default I = 0.3 M comes from the in vivo-like yeast assay medium of van Eunen *et al.* (2010): 300 mM K⁺, 245 mM glutamate, 50 mM phosphate, 20 mM Na⁺, 2 mM free Mg²⁺, about 0.33 M ionic strength. That gives a Debye length of 0.55 nm, so electrostatics act over about a nanometre. The Salt setting changes it.
- **Charges:**
  - Asp and Glu −1, Lys and Arg +1, His neutral at the cytosolic pH of 7.2 (Orij *et al.* 2009).
  - Each nucleotide −0.25 e. On a double helix counterion condensation leaves the fraction 1/ξ = b/l_B of the phosphate charge (Manning 1969; b = 0.17 nm).
  - A sphere carries the summed charge of the residues it represents. In the residue model, buried residues give their charge to the nearest exposed bead.
- **Cutoff.** Three Debye lengths beyond contact, with the force shifted to zero there.
- **Where the charges come from.** The structure library keeps only bead coordinates, so residue identities are recovered from the source structures by `pipeline/simdata.py`:
  - **Matching.** A library species' protein beads are its source's C-alpha atoms after a rigid motion (centring, membrane orientation) and the removal of some residues (low pLDDT, unused chains). Windows of four consecutive C-alphas are matched by their six internal distances. A window pose that brings 90% of the beads within 0.03 nm of a source C-alpha is accepted, and every bead takes the residue of its nearest one.
  - **Sources:** AlphaFold DB monomer models, the AlphaFold DB dimers and PDB assemblies, fetched from RCSB when not cached.
  - **Composites:** structures such as a ribosome with its translocon take the part found in the structure of the species they extend.
  - **Coverage:** 98.6% of the 3.25 M protein beads are matched. The 11 species that are not are custom-built scaffolds (actin patch, tropomyosin path, γ-TuRC, ATP synthase dimers) and four small proteins; they stay neutral.
  - **Checks:** net charges match the sequence (actin −10 in the library, about −12 from its sequence). The ribosome carries −4,017 e of bare charge (5,209 nucleotides, 2,191 basic and 999 acidic protein residues). That shrinks to about −110 e after condensation.
  - **Output:** `library/charges.bin`, one byte per bead.

**Specific binding** (on by default). Proteins that form a confident AlphaFold DB dimer (ipTM ≥ 0.5; the same predicted heterodimers and homodimers as the packing, section 5) attract each other at their predicted interface when both are free in the box:
- **Contacts.** Native contacts are the inter-chain residue pairs with C-alpha atoms within 0.8 nm in the dimer model, mapped onto the library beads of the two monomer species.
- **Coverage.** 371 dimer models qualify, out of 1,848 confident ones; the others lack a monomer species of one of the partners in the library. The median interface has 61 contacts.
- **Energy.** Every contact pair between two molecules of the partner species gets a well of the same cos² shape as the stickiness, of depth E/n. E is the energy of a complete interface, 20 kT by default, and n the number of distinct contact pairs at the model's resolution (residues, clusters or sphere-tree spheres), so the total does not depend on the model.
- **Orientation.** The attraction is orientation-specific, since only the predicted interface residues attract: a Gō-like interface (Levy, Wolynes & Onuchic 2004).
- **One energy for all.** Affinities are not known per pair, so E is the same for every complex. Complexes already packed as one rigid body (section 5) are unaffected.
- **Output:** `library/contacts.json` and `contacts.bin`.

**Polysome mRNA links** (on by default). Linked ribosomes (section 14.1) are joined between the mRNA exit of one and the entry of the next by a worm-like chain (Marko & Siggia 1995):

```
F = (kT / l_p) · (1 / (4 (1 − x/L_c)²) − 1/4 + x/L_c)        (x capped at 0.95 L_c)
```

- **Contour length.** In growing yeast, ribosomes sit one per 156 nt of ORF on average (0.64 per 100 nt; Arava *et al.* 2003). The ribosome covers the 29 nt of its channel, so a linker holds about 127 nt. At 0.49 nm per nucleotide (the contour length fitted to rU40 by Chen *et al.* 2012), L_c = 62 nm.
- **Stiffness.** Chen *et al.* (2012) measured a persistence length l_p of 1.0 nm at high salt and 2.5 nm without added salt for single-stranded rU40. 1.5 nm is used.
- **Result.** The rms end-to-end distance is about 14 nm, so a polysome stays together as a loose cluster. Every linker gets the mean length.
- **Drawing.** The links are drawn as dark brown lines between the anchors.

**Integration.** Each step:
1. **Binning.** The spheres go into two uniform grids. Grid A holds spheres up to the 90th-percentile radius, with cells of that diameter plus the attraction range. Grid B, with larger cells, holds the bigger spheres (schematic and mid) and the cores, so the common small spheres need not scan far for the rare big ones. Per-cell capacities come from the fullest cell of the starting configuration ×1.5; an earlier fixed capacity of 24 overflowed and gave asymmetric forces.
2. **Forces.** Each sphere scans the cells within its interaction range of both grids and skips:
   - cells whose box is farther than that range;
   - cells holding only spheres of its own molecule (each cell records the lowest and highest body id binned into it);
   - in the shell between contact range and the electrostatic cutoff, uncharged spheres (flagged in the cell items);
   - core spheres, unless the scanning sphere is a core itself.

   Cells can still fill beyond capacity when molecules cluster tightly. The panel reads the grid occupancy back every two seconds, and every eight command buffers while computing ahead. On any overflow it doubles the capacities before the next submission.
3. **Links.** Bond forces are computed per body.
4. **Reduction.** One workgroup per body sums its force and torque over its spheres (parallel reduction).
5. **Move.** An overdamped Euler–Maruyama step:
   ```
   dx = D_t F dt + √(2 D_t dt) ξ
   dθ = D_r τ dt + √(2 D_r dt) ξ        (kT = 1)
   ```
   The deterministic displacement is capped at 0.15 nm per step, both for the body centre and for the body surface under rotation.

All steps run in WebGPU compute shaders, up to 512 steps in one compute pass. Computing ahead, each command buffer holds about 0.1 s of GPU work, so that no submission trips the GPU watchdog.

**Speed.** Before the two-level grid, the cell skipping and the parallel per-body sums, a 60 nm residue box (surface beads, no electrostatics) ran at 279 steps/s. Half the time went to the per-body sums, which one thread per body had done serially. After: 2,145 steps/s (7.7×). A 120 nm schematic box went from 643 to 1,546 steps/s (2.4×). The rewrite left diffusion unchanged within noise (section 14.5). Steps per second with the cores, Apple-silicon laptop GPU, cytoplasm boxes:

| Model | Box | Collision spheres | Electrostatics on | Off |
|---|---|---|---|---|
| Residue | 60 nm (185 molecules) | 96,432 | 634 (0.063 µs/s) | 1,631 (0.16 µs/s) |
| Mid | 60 nm | 39,900 | 1,185 (0.24 µs/s) | 2,097 (0.42 µs/s) |
| Schematic | 120 nm (1,476 molecules) | 41,847 | 944 (0.19 µs/s) | 1,441 (0.29 µs/s) |

About half of the residue model's surface beads are charged, and each scans 1.65 nm further, so electrostatics cost 2.6× there; switching them off is the quickest speed-up.

### 14.4 Time step and convergence

Long steps overestimate diffusion. The table gives the diffusion of the GFP-sized class (R_h < 3 nm) in a 120 nm cytoplasm box, mean-field hydrodynamics on, time-averaged MSD at a 0.7 µs lag, in µm²/s:

| Stickiness ε | dt 0.02 ns | 0.05 ns | 0.1 ns | 0.2 ns | 0.5 ns |
|---|---|---|---|---|---|
| 0 | | 61.9 | 60.5 | 60.1 | 71.8 |
| 1 kT | 42.4 | 40.0 | 47.3 | 58.4 | |
| 2 kT | 13.1 | 16.0 | 24.8 | | |

- **Without attraction,** dt ≤ 0.2 ns is converged.
- **With attraction,** the 0.6 nm well needs about dt ≤ 0.05/ε ns. Contact forces summed over many sphere pairs are stiff, and a longer step lets molecules escape the well too easily. A broader 1.2 nm well did not help: more pairs fall inside it.
- **Mid model.**
  - **Without attraction,** the GFP-sized class gives 51.7, 51.6 and 51.0 µm²/s at dt 0.05, 0.1 and 0.2 ns in an 80 nm box: converged at 0.2 ns.
  - **With attraction,** the schematic rule is not enough. At 1 kT the class gave 17.9 µm²/s at dt 0.05 ns and 11.8 at 0.01 ns. Its clusters meet in many pairs per contact, which makes the summed contact forces stiff. The panel uses dt ≤ 0.01/ε ns; smaller steps were not tested.
  - **Unstable at 1.5 kT.** Even dt 0.01 ns failed there: large complexes moved faster than in water (28.6 against 23.6 µm²/s). The Stickiness slider is therefore capped at 1 kT in the mid model and 0.5 kT in the residue model (2 kT in the schematic model).
- **The panel enforces these limits.** The residue model is assumed to scale the same way (dt ≤ 0.0125/ε ns, at most 0.1 ns); this has not been measured.
- **Native contacts** make wells deeper than ε for the few bound pairs (E/n per pair: about 0.3 kT per residue pair, a few kT per sphere-tree pair at E = 20 kT). The time step is not lowered for them, so bound complexes in the schematic model may come apart somewhat too easily.

### 14.5 Calibration against yeast measurements

The readout reports, per size class, the time-averaged mean squared displacement MSD(τ) of the moving molecules over the recording, as D = MSD/(6τ). It uses the longest lag whose rms displacement stays below a quarter of the free interior of the box; longer lags feel the confinement by the box faces.

**Before electrostatics** (schematic model, 120 nm cytoplasm box of the current packing, φ = 0.198, periodic boundary, all 1,476 molecules moving, converged steps, mean-field hydrodynamics on unless noted):

| Stickiness ε | GFP-sized (R_h < 3 nm) | R_h 3–6 nm | R_h > 6 nm (ribosomes, large complexes) |
|---|---|---|---|
| 0, hydrodynamics off | 69 µm²/s (0.62 of water) | 42 | 7.6 |
| 0 | 43 (0.39) | 29 | 4.1 |
| 1 kT | 24 (0.21) | 12 | 2.4 |
| 1.5 kT | 10.0 (0.090) | 6.9 | 1.9 |
| 2 kT | 8.0 (0.072) | 6.8 | 2.4 |

The earlier frozen-shell boundary gave the same numbers within noise (62, 40, 23, 9.0 and 6.3 µm²/s for the GFP-sized class). Before the proteome-mass correction (L1; φ = 0.106, 1,051 molecules) the same box gave 77, 62, 40–42 and 13 µm²/s for the GFP-sized class. The large class holds only 15–25 molecules, so its values scatter. After the engine rewrite (section 14.3) the same box gave 41.3 and 10.4 µm²/s at ε = 0 and 1.5 kT.

**With electrostatics** (0.3 M; GFP-sized class, µm²/s; schematic in the 120 nm box, mid and residue in an 80 nm box of 446 molecules, 113 of them GFP-sized; 1 µs runs; mid and residue with core spheres):

| Stickiness ε | Schematic | Mid | Residue |
|---|---|---|---|
| 0 | 42.2 | 54.2 | 50.9 |
| 0.25 kT | | | 56.5 |
| 0.5 kT | | | 31.9 |
| 1 kT | 23.9 | 17.9 (11.8 at dt 0.01 ns) | |
| 1.25 kT | | 22.7 (unstable) | |
| 1.5 kT | 10.0 | 24.6 (unstable) | |
| 2 kT | 6.6 | | |

- **Electrostatics change little at cytoplasmic salt.** Without attraction they slowed the GFP-sized class by 0–4% (schematic 41.3 → 42.2, mid about 51 → 49). The Debye length is 0.55 nm, and most charged groups never come that close.
- **Without attraction the three models agree** within the noise of these runs: 42–54 µm²/s. The schematic model is slowest, since its sphere trees are bulkier than the residues.

**Stickiness slows molecules by clustering them.** Clusters were counted in 80 nm boxes after 1 µs. Two molecules count as in contact when any two of their surface spheres lie within 0.3 nm of touching.

| Model, ε | Largest clusters (of 446 molecules) | Single molecules | Contacts per molecule | GFP-sized D (µm²/s) |
|---|---|---|---|---|
| Schematic, 0 | 20, 18, 14 | 192 | 0.9 | 46 |
| Schematic, 1 kT | 399 | 16 | 3.8 | 21 |
| Schematic, 1.5 kT | 247, 123 | 7 | 4.4 | 10.4 |
| Schematic, 2 kT | 429 | 6 | 4.5 | 6.4 |
| Mid, 0 | 15, 10, 9 | 255 | 0.6 | 47 |
| Mid, 1 kT | 184, 109, 90 | 7 | 4.6 | 15 |
| Mid, 1.5 kT (unstable steps) | 116, 53, 51 | 4 | 12.0 | 27 |
| Residue, 0.5 kT | 285, 91, 14 | 12 | 9.3 | 25 |

- **Schematic model.** The stickiness that matches the GFP measurement (1.5–2 kT) joins most molecules into one network within a microsecond, with voids between its strands. The earlier calibration matched the number this way.
- **Mid model.** Its clusters of ~4 residues meet in many more pairs per contact than sphere-tree spheres, so the same ε binds much more strongly. At 1 kT most molecules aggregate. The rise in D at 1.25–1.5 kT is an integration artefact, not physics: those steps were too long for the stiff many-pair contacts (section 14.4).
- **Residue model.** Most molecules aggregate already at 0.5 kT.
- **What the gap really is (corrected 2026-11-11).** This section used to say that full many-body hydrodynamics closed the gap in *E. coli* (Ando & Skolnick 2010) and suggested a Viscosity of 5–6× in its place. That result was obtained at a much higher density (section 14.2), and adding far-field many-body hydrodynamics to this model's boxes changed little. The remaining gap is now calibrated on measurements instead (below).
- **Defaults.** Stickiness defaults to 0. Viscosity defaults to **2.75×**, the calibration below. The panel's readout says what each choice does for the model in use.

**Calibration on a literature benchmark** (2026-11-11; `validation/diffusion.py`, section 9.1). Diffusion of a molecule is modelled as D = D₀(R_h) × crowding(R_h) × *s*: D₀ from Stokes–Einstein with the R_h of the measured construct (the species that diffuses, with its trimmed IDRs restored and its fluorescent tag added as a bead model), crowding(R_h) from Brownian dynamics of packed boxes with mean-field hydrodynamics on 4 GPUs (`simgpu/sweep.py`: ln(D/D₀) = −1.04 − 0.099 R_h/nm in a 120 nm cytoplasm box, −1.55 − 0.142 R_h/nm in a 100 nm nucleoplasm box with its fixed chromatin), and *s* a remaining compartment factor fitted on the free tracers of the training papers: 0.365 for the cytoplasm (hence Viscosity 2.75×), 0.26 for the nucleoplasm (single-molecule tracking only). On the held-out papers the size dependence from the packed boxes beats size alone (within-paper ratios RMSE 0.11 vs 0.21 log10 units). Directly simulated 40 nm GEM particles move at D/D₀ = 0.04–0.15 in a 200 nm box, consistent with the fit; with *s* this predicts 0.25 µm²/s against 0.3 measured in both yeasts (Delarue *et al.* 2018; Garner *et al.* 2023; Molines *et al.* 2022). A binding term fixed in advance (RNA-interactome hits plus abundance-weighted pool partners) fitted the training ratios slightly better but made the held-out ones worse, so it is not used: what slows Srv2 or Ede1 (0.48 and 0.59 µm²/s, about 10× and 2× below the size prediction) is not captured by generic scores.

**Polysomes and large particles: a test that fails** (2026-11-12; `simgpu/gade_test.py`, links in `simgpu/engine.py`; nothing was fitted). Dissolving polysomes with hippuristanol speeds 40 nm GEM particles 1.4–1.5 times and 20 nm ones 1.3 times, while removing 20–25% of ribosomes in *stm1Δ* changes nothing (Gade *et al.* 2025). The server-side engine got the polysome links the browser has (neighbouring ribosomes of one mRNA tied by entropic springs), and six 200 nm boxes of the cytosol of build Q (about 85 ribosomes and 20–44 links each) were run for 10 µs with and without the links, and with 22% of the ribosomes removed; 72 tracers per condition, mean-field hydrodynamics.
- **Result:** without links the 40 nm particles are 1.02 times as fast (95% interval 0.98–1.06 at a lag of 300 ns), the 20 nm particles 0.98 times; with fewer ribosomes 1.05 times (measured 1.0).
- **What it means:** tying ribosomes together does not slow large particles in this model; ribosome crowding itself does, slightly. The measured effect needs something the simulation lacks: the mRNA as a physical chain between and beyond the ribosomes (it is not a body here), a mesh spanning several polysomes, or times beyond the 10 µs reached (tracking measures at 10–100 ms).
- **Not re-run on builds R and S:** the cytosol's density and polysome statistics are the same (sections 6.5 and 10.1).

**Measured values:**
- GFP in the fission-yeast cytoplasm: 8.6 ± 2.2 µm²/s, and 5.6 ± 2.8 µm²/s in the nucleus, by FRAP at room temperature (Delivani *et al.* 2013, checked in the text). This is roughly 0.1 of GFP in water.
- Free GFP in the budding-yeast cytoplasm by FCS: 8.7–21 µm²/s depending on the laboratory (Slaughter *et al.* 2007; Boeke *et al.* 2014; Pack *et al.* 2017). The benchmark table lists all values with quotes.
- 40 nm GEM particles in fission yeast: 0.4–0.6 µm²/s by particle tracking (Tan *et al.* 2025).
- 40 nm GEM particles in budding yeast (added 2026-11-12, held out: read from figures and placed in the test split whatever their paper hash, since they were seen after the fit): about 0.135 and 0.21 µm²/s in two experiments of Gade *et al.* 2025 (Fig. 2 of the preprint) and 0.20 in Kunzi *et al.* 2026 (preprint, Fig. 1D). The model predicts 0.25, so at this size it is 1.2–1.8× too fluid; Delarue *et al.* 2018's ~0.3, which is in the training split, sits on the other side.
- In budding yeast, ribosome crowding controls the diffusion of particles ≥ 20 nm but not of molecules ≤ 5 nm (Delarue *et al.* 2018).

**Interpretation:**
- Excluded volume, electrostatics and mean-field hydrodynamics slow GFP-sized molecules to about 0.3–0.5 of their speed in water. That is still about 2.7× faster than the benchmark's free tracers (and about 5× faster than the slowest GFP measurement), as found for bacterial models.
- A uniform nonspecific attraction can close the gap numerically, as in McGuffee & Elcock (2010). Here it does so by clustering the molecules into networks, so it is a poor stand-in for what the model lacks (many-body hydrodynamics, transient specific binding, metabolites, the mRNA and other polymers left out).
- Crowding itself matters: doubling the macromolecular density (L1) slowed GFP-sized molecules from 0.55 to 0.36 of their speed in water with no attraction.

### 14.6 Running, recording and playback

- **Placing the box.** "Place box" shows the outline at the centre of the view and frames it. While placing, the box follows the view centre: pan, orbit and zoom as usual, or click a point to bring it to the centre. The outline is drawn on top of everything, including the parts outside the cut slab.
- **Changing settings.** Changing the model, boundary, run mode, box size or the membrane option restarts a running simulation in place, with the same box.
- **Live.** Each drawn frame computes the chosen number of steps, and the molecules move as they are computed.
- **Compute first, then play.** The chosen length is computed with command buffers submitted back to back (up to eight in flight, no waiting for the display) behind a progress bar. The recording then plays in a loop. The Playback slider sets the rate, 1–60 recorded frames per second (24 by default).
- **Tracking.** Clicking a molecule in a running simulation tracks it; up to eight can be tracked, and clicking again drops one. The other molecules are drawn faded: one pixel in four, lightened, still clickable, without outlines. Each tracked molecule's centre path over the recording, up to the frame shown, is drawn as a coloured line on top of everything. A segment that crosses a periodic face is left out, since the molecule re-enters on the other side. The panel lists each tracked molecule's displacement since the start of the recording and its own diffusion estimate from successive frames.
- **Recording.** About 300 frames are recorded per run. In live mode recording continues indefinitely; past the frame budget (about 60 MB of poses) every other frame is dropped and the interval doubles, so the spacing stays uniform.
- **Throughput:** see the speed table in section 14.3. With stickiness, the shorter step lowers it in proportion. Live mode is slower, since it waits for the display.

### 14.7 Export of MD-ready files

"Export MD files" saves five files:

| File | Contents |
|---|---|
| `yeast_box.cif` | Coarse-grained mmCIF. One pseudo-atom per residue (CA) or nucleotide (P) bead, one chain per molecule. Chain ZZZZ holds membrane surface beads from the wall grid. Occupancy is 1 for moving molecules, 0 for fixed ones. |
| `yeast_box.json` | Manifest: per molecule the species, genes, UniProt, structure source (PDB entry or AlphaFold DB model), region, pose (position, quaternion), R_h, D_t, D_r, mobility; plus the settings, φ and the hydrodynamic factors. |
| `yeast_box_trajectory.json` | The recorded rigid-body trajectory, up to 100 frames. |
| `run_openmm.py` | Runs the box in OpenMM (Eastman *et al.* 2017, 2024) as a coarse-grained model, with the same boundary (periodic boundary conditions for periodic boxes, walls otherwise, plus the shell tethers of embedded boxes) and the same interactions where enabled. |
| `trajectory_to_dcd.py` | Expands the browser trajectory into a bead-level DCD that VMD, ChimeraX, MDAnalysis or mdtraj open together with the mmCIF. |

The manifest also lists the interactions of the browser run, at residue level: the charge of every charged bead, the mRNA links (between the beads nearest the channel anchors, with contour and persistence lengths), and the native contact pairs between every two molecules of partner species, with their share of the interface energy.

The OpenMM model has:
- an elastic network within each moving molecule (bead pairs closer than 1.2 nm, 500 kJ/mol/nm²);
- Weeks–Chandler–Andersen excluded volume between molecules, plus the same cos² attraction;
- Debye–Hückel electrostatics between the charged beads (cut off three Debye lengths beyond contact), worm-like-chain mRNA links and cos² native-contact wells, as exported;
- massless (fixed) beads for membranes, filaments and chromatin, with fixed–fixed pairs switched off;
- soft box walls;
- Langevin dynamics at 303 K with a 20 fs step after energy minimisation; DCD output.

Tested: a 60 nm plasma-membrane box gave 92,606 beads, 2,896 of them membrane beads, 18,037 moving and 266,635 network bonds. It runs on the OpenCL platform of OpenMM 8.6. On the cluster (2026-11-11) a 60 nm cytoplasm box (175,045 beads, 2.4 M network bonds) runs on CUDA at about 700–1,100 ns/day per RTX 4090. It first blew up within 80–7,200 steps: ribosome models carry about 20 bead pairs 0.001–0.1 nm apart (overlapping protein and RNA beads), which got network bonds of near-zero length. Pairs closer than 0.1 nm now get no bond, and the excluded volume continues linearly below 0.8 σ; the box then ran 50,000 steps (1 ns).

For atomistic or Martini systems, take each molecule's structure from the manifest's source and apply its pose to the structure's bead frame (CA/P positions centred as in the library).

### 14.8 Limitations

- Rigid bodies: no internal flexibility, and disordered tails stay in their AlphaFold pose.
- No reactions.
- Electrostatics are Debye–Hückel between point charges at sphere centres: no dielectric boundary, no ion correlation, no explicit Mg²⁺. In the schematic and mid models a sphere's summed charge sits at its centre.
- Specific binding exists only for pairs with a confident predicted dimer whose partners are both monomer species, with one interface energy for every pair. No other complexes, no measured affinities.
- mRNA links use the mean linker length. Free mRNA threads, the mRNA inside polysomes and links leaving the box are absent.
- Hydrodynamics only as a mean-field rescaling with one volume fraction for the whole box. The remaining slowdown is a calibrated, uniform viscosity (2.75×, section 14.5): it does not distinguish what causes it (fluid-phase viscosity, missing polymers such as mRNA and IDRs, transient binding), and specific slow-downs of individual proteins (Srv2 about 10× below its size prediction) are not reproduced.
- Membrane proteins: Saffman–Delbrück with a fluid bilayer gives about 1–4 µm²/s, while yeast plasma-membrane proteins measured by FRAP and tracking move at 0.0004–0.05 µm²/s (Bianchi *et al.* 2018; Slaughter *et al.* 2013). Since 2026-11-11 a factor per membrane, fitted on FRAP data, closes this gap on average (section 14.2), but nothing mechanistic (domains, the cell wall, immobile fractions) is modelled, and three of the eight membranes borrow the ER factor.
- Filaments and chromatin are fixed.
- No mRNA, lipids or cell-wall polymers as bodies.
- In embedded and closed boxes the faces confine the moving molecules, so diffusion over distances near the box size is not meaningful; the readout's lag limit avoids it. Periodic boxes let molecules meet their own images across a small box.
- Results depend on the time step beyond the limits in 14.4.

### 14.9 References for the simulation

| Reference | Used for |
|---|---|
| Saffman PG, Delbrück M (1975). *PNAS* 72:3111–3113 | Lateral and rotational diffusion of membrane proteins |
| Valdez-Taubas J, Pelham HRB (2003). *Curr. Biol.* 13:1636–1640; Meinema AC *et al.* (2013). *Traffic* 14:487–501; Spira F *et al.* (2012). *Nat. Cell Biol.* 14:640–648 | Membrane-protein diffusion in yeast (calibration per membrane) |
| van Eunen K *et al.* (2010). *FEBS J.* 277:749–760 | In vivo-like ionic composition of the yeast cytosol (ionic strength) |
| Orij R, Postmus J, Ter Beek A, Brul S, Smits GJ (2009). *Microbiology* 155:268–278 | Cytosolic pH 7.2 in growing yeast (residue charge states) |
| Manning GS (1969). *J. Chem. Phys.* 51:924–933 | Counterion condensation on polyelectrolytes (effective nucleotide charge) |
| Levy Y, Wolynes PG, Onuchic JN (2004). *PNAS* 101:511–516 | Gō-like models of protein binding through native interface contacts |
| Arava Y *et al.* (2003). *PNAS* 100:3889–3894 | Ribosome density on yeast mRNAs: 0.64 per 100 nt (one per 156 nt) |
| Chen H *et al.* (2012). *PNAS* 109:799–804 | Persistence length (1.0–2.5 nm) and contour length of single-stranded RNA |
| Marko JF, Siggia ED (1995). *Macromolecules* 28:8759–8770 | Worm-like-chain force–extension interpolation |
| McGuffee SR, Elcock AH (2010). *PLoS Comput. Biol.* 6:e1000694 | BD of a structure-based cytoplasm; attraction calibrated to GFP diffusion |
| Ando T, Skolnick J (2010). *PNAS* 107:18457–18462 | Crowding plus hydrodynamic interactions account for in-vivo diffusion; GFP D/D0 0.31–0.46 without and 0.08 with hydrodynamics |
| Mereghetti P, Wade RC (2012). *J. Phys. Chem. B* 116:8523–8533; Balbo J *et al.* (2013). *Biophys. J.* 104:1576–1584 | Mean-field treatment of hydrodynamic interactions in BD |
| Cichocki B, Ekiel-Jeżewska ML, Wajnryb E (1999). *J. Chem. Phys.* 111:3265–3273 | Short-time translational and rotational self-diffusion vs volume fraction |
| Yu I *et al.* (2016). *eLife* 5:e19274 | All-atom MD of a bacterial cytoplasm (scale comparison) |
| Delivani P, Chacón JM, Schroth-Diez B, Tolić-Nørrelykke IM (2013). *Bio-protocol* 3:e941 | GFP diffusion in fission-yeast cytoplasm and nucleus (FRAP) |
| Slaughter BD *et al.* (2007). *PNAS* 104:20320–20325; Boeke D *et al.* (2014). *Mol. Syst. Biol.* 10:756; Pack CG *et al.* (2017). *BMB Rep.* 50:478–483 | Free GFP and tagged proteins in budding-yeast cytoplasm (FCS; benchmark) |
| Garner RM *et al.* (2023). *Biophys. J.* 122:767–783; Molines AT *et al.* (2022). *Dev. Cell* 57:466–479 | 40 nm GEM diffusion in fission yeast (benchmark) |
| Slaughter BD *et al.* (2013). *Nat. Commun.* 4:1380; Bianchi F *et al.* (2018). *Nat. Commun.* 9:501 | Plasma-membrane protein diffusion in yeast |
| Rotne J, Prager S (1969). *J. Chem. Phys.* 50:4831–4837; Zuk PJ *et al.* (2014). *J. Fluid Mech.* 741:R5; Hasimoto H (1959). *J. Fluid Mech.* 5:317–328 | Rotne–Prager–Yamakawa mobility, unequal and overlapping spheres, periodic (Ewald) sum (`simgpu/hydro.py`) |
| Bernadó P, Blackledge M (2009). *Biophys. J.* 97:2839–2845; Marsh JA, Forman-Kay JD (2010). *Biophys. J.* 98:2383–2390 | Size scaling of disordered chains (IDR hydrodynamic radii) |
| Tan C, Lanz MC, Swaffer M, Skotheim J, Chang F (2025). *Mol. Biol. Cell* 36:ar51 | 40 nm GEM diffusion in fission yeast |
| Delarue M *et al.* (2018). *Cell* 174:338–349 | Ribosome crowding sets mesoscale diffusion in budding yeast |
| Eastman P *et al.* (2017). *PLoS Comput. Biol.* 13:e1005659; Eastman P *et al.* (2024). *J. Phys. Chem. B* 128:109–116 | OpenMM |
| Souza PCT *et al.* (2021). *Nat. Methods* 18:382–388; Stevens JA *et al.* (2023). *Front. Chem.* 11:1106495 | Martini 3; Martini model of an entire minimal cell |
| Andrews SS (2017). *Bioinformatics* 33:710–717; Hoffmann M, Fröhner C, Noé F (2019). *PLoS Comput. Biol.* 15:e1006830 | Particle-based reaction–diffusion (Smoldyn, ReaDDy 2) |
| Roberts E, Stone JE, Luthey-Schulten Z (2013). *J. Comput. Chem.* 34:245–255; Thornburg ZR *et al.* (2022). *Cell* 185:345–360 | Lattice Microbes; whole minimal-cell simulation |
