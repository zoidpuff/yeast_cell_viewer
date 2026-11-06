# Whole yeast cell: model description

This document describes how the structural model shown in the viewer was built. It covers what the model contains, which data and structures it comes from, the assumptions and design decisions behind it, how it was validated, and where it is known to be wrong or incomplete. It is written for readers who want to judge the model scientifically or reuse its numbers.

All numbers refer to the current build (5,638 molecular species, 21.6 million placed instances) unless stated otherwise.

## 1. What the model is

- **The cell.** One haploid *Saccharomyces cerevisiae* cell (strain BY4741 / S288C background) growing exponentially in rich glucose medium (YPD). The cell is unbudded, and the microtubules and chromatin are arranged as in G1.
- **What a molecule is.** Every macromolecule is an individual copy of an experimental (PDB) or predicted (AlphaFold DB) structure. It is drawn at one bead per amino-acid residue (Cα) or per nucleotide (P or C4′). There are 14.2 billion beads in all.
- **How the cell was assembled.** Copy numbers come from proteomics. Compartment volumes and membrane areas come from imaging literature. The organelles are built as analytic geometry, and molecules are packed into them with residue-level collision checks.
- **What it is not:**
  - **Not a dynamic or physical simulation.** There are no forces, no diffusion and no binding equilibria. Except where stated below, molecules are placed at random positions and orientations by random sequential addition. A snapshot from a real cell would show interaction-driven clustering that the model lacks.
  - **Not a reconstruction of a particular cell.** The organelle volumes and areas are taken from measurements, but their shapes and positions are a plausible arrangement, not a tomogram.
- **Lineage.** The approach follows the cellPACK-style whole-cell models of *Mycoplasma genitalium* (Maritan *et al.* 2022). The visual style follows David Goodsell.

## 2. Composition

### 2.1 Proteins

- **Abundances:** PaxDb integrated whole-organism dataset for *S. cerevisiae* (taxon 4932), in ppm, matched to UniProt by ORF.
- **Absolute scale:** ppm values are rescaled to 4.2 × 10⁷ protein molecules per cell, the unified total of Ho, Baryshnikova and Brown (2018), who used the same 42 µm³ cell volume. This gives 2.7 pg of protein per cell (abundance-weighted mean mass 39 kDa).
- **Uncertainty:**
  - PaxDb integrated values are not specific to BY4741 in YPD. Per-protein uncertainty is roughly 2-fold.
  - The absolute total is the largest single uncertainty in the model. 4.2 × 10⁷ molecules in 42 µm³ is 1.0 × 10⁶ per µm³, at the low end of published normalised estimates (1–3 × 10⁶ per µm³; Dolgalev *et al.* 2023). If the true total is higher, every protein density in this document scales up with it.

### 2.2 Nucleic acids and fixed counts

| Component | Copies | Basis |
|---|---|---|
| 80S ribosomes | 169,000 | Ribosomal protein abundances (median over the r-protein families) |
| tRNA | 3.3 × 10⁶ | Assumption, literature order of magnitude |
| mRNA | 30,000 | Assumption; mean length 1,500 nt including UTRs and poly(A) |
| Nucleosomes | 73,150 | 12.07 Mb nuclear genome (rDNA repeats not expanded) ÷ 165 bp repeat length |
| Nuclear pore complexes | 117 | Winey *et al.* 1997 |

The ribosome count is probably low. A growth-law estimate for YPD is about 3 × 10⁵ (Kostinski & Reuveni 2021). Because the count follows the ribosomal-protein abundances, it inherits the uncertainty in the proteome total.

### 2.3 Not included

- **Small molecules and ions:** metabolites, free amino acids, nucleotides, ions and water.
- **Storage carbohydrates:** glycogen and trehalose.
- **Polyphosphate.**
- **rDNA and pre-rRNA:** the rDNA array and nascent pre-rRNA in the nucleolus, and most pre-ribosomal particles.
- **Proteins without a structure:** 6,888 copies.
- **Transmembrane proteins annotated only to soluble compartments:** 9,946 copies (see 4.3).

Lipids and the cell-wall polysaccharides are drawn in the viewer but generated procedurally (section 7).

## 3. Geometry

### 3.1 Volumes

The cell is a sphere of 42 µm³ (mean log-phase volume; Jorgensen *et al.* 2002) with a 110 nm cell wall, the lower bound of 110–200 nm (Orlean 2012).

Organelle volume fractions:
- **Nucleus, vacuole and lipid droplets:** soft X-ray tomography (Uchida *et al.* 2011): nucleus 7%, vacuole 7%, lipid droplets 1%. The nucleolus is 20% of the nucleus.
- **Mitochondria:** raised from the measured 1% to 2%. The assumption is that X-ray tomography under-segments thin tubules. At 1%, the mitochondrial proteome would imply about 480 mg/mL of protein.
- **Golgi, endosomes and peroxisomes:** 1% together (assumption).
- **Cortical ER:** covers 48% of the plasma membrane (Wolf *et al.* 2012).

Measured on the final voxelised geometry (4 nm grid):

| Compartment | Volume (µm³) | Share of cell |
|---|---|---|
| Cytosol | 27.81 | 66.2% |
| Cell wall (glucan 4.09, mannan 2.01) | 6.10 | 14.5% |
| Vacuole | 2.94 | 7.0% |
| Nucleoplasm | 2.02 | 4.8% |
| Nucleolus | 0.59 | 1.4% |
| Mitochondrial matrix | 0.61 | 1.5% |
| Mitochondrial intermembrane space (incl. crista lumen) | 0.24 | 0.6% |
| ER lumen: cortical ER | 0.50 | 1.2% |
| ER lumen: perinuclear space | 0.29 | 0.7% |
| Lipid droplets | 0.42 | 1.0% |
| Golgi / endosomes | 0.34 | 0.8% |
| Peroxisomes | 0.08 | 0.2% |
| **Total** | **41.98** | |

### 3.2 Membranes

Membrane mid-surface areas, measured on the final geometry:

| Membrane | Area (µm²) |
|---|---|
| Plasma membrane | 52.6 |
| ER (cortical sheets and perinuclear) | 51.8 |
| Nuclear envelope | 18.3 |
| Mitochondrial inner membrane | 28.4 |
| Mitochondrial outer membrane | 12.3 |
| Vacuole | 9.9 |
| Golgi / endosomes | 8.9 |
| Peroxisomes | 1.5 |

### 3.3 Organelle shapes

All membranes are analytic signed-distance fields (spheres, capsules, pills and cortical-ER patches), so containment and membrane distance are exact everywhere. Shapes and fixed dimensions:

| Element | Shape and dimensions |
|---|---|
| Nucleus | A sphere, with the nucleolus as a crescent on the side away from the vacuole. Nuclear envelope: perinuclear space 34 nm between the membrane mid-surfaces; 117 pores of 45 nm radius at the mid-plane. |
| Cortical ER | Sheets 27 nm beneath the plasma membrane, with a 20 nm lumen. |
| Mitochondria | One reticulum of tubules with 150 nm radius. Intermembrane space 8 nm. Lamellar cristae with an 18 nm lumen (rim-tube radius 9 nm, chosen to fit ATP synthase dimers; section 5.6). |
| Other organelles | Vacuole, lipid droplets, Golgi/endosome cisternae and peroxisomes are simple primitives. |

**Arrangement.** Organelle positions follow a plausible layout:
- The nucleus is off-centre.
- The spindle pole body faces away from the nucleolus.
- The mitochondria form branched tubules along the cortex.

This layout is illustrative; it is not a reconstruction of a particular cell.

**Known geometric deviations from EM.** Serial-section EM (Yamaguchi *et al.* 2011) gives a nucleus of 10.5% of cell volume; this model has about 7%. Mitochondria are 1.7% there and about 2% here.

## 4. Localization

### 4.1 Source annotations

**UniProt subcellular-location text** is mapped to compartments by ordered regular-expression rules. Examples:
- "nucleolus"
- "mitochondrion inner membrane"
- "endoplasmic reticulum"
- "cell cortex" → plasma membrane

Annotations are filtered by their evidence:
- **GFP data preferred:** where a location token is backed by the genome-wide GFP study of Huh *et al.* (2003), only GFP-backed tokens are kept. This gives the bulk steady-state localization.
- **Proteomics-only calls dropped:** tokens supported only by two high-throughput mitochondrial proteomics papers (Sickmann 2003, Reinders 2006) are discarded. Those studies assign many abundant cytosolic enzymes to mitochondria.

**Fallbacks:**
- **GO cellular component:** used when UniProt has no location. Complex terms are skipped, except pores, ribosomes and proteasomes. If GO lists the cytosol, only cytosolic and nuclear terms are kept.
- **Cytosol:** the default for proteins with no annotation at all (625 proteins).

### 4.2 Weights across compartments

- **Rank order:** a protein's compartments get weights 1, ½, ¼, … in annotation order, with the first-listed location taken as primary.
- **Cytosol and nucleoplasm:** their shares are pooled and re-split by volume, i.e. equal concentration across the nuclear pore.
- **Secondary mitochondrial locations:** a mitochondrial location listed after a cytosolic or nuclear primary one (e.g. SOD1, ADK1 in the intermembrane space) gets 0.1× weight.
- **Mitochondrial sub-compartments:** taken from UniProt or GO where named (inner membrane, outer membrane, intermembrane space, matrix). Otherwise transmembrane proteins are split 75/25 between the inner and outer membranes, and soluble proteins 90/10 between matrix and intermembrane space.

### 4.3 Membrane proteins

- **Definition:** proteins with UniProt TRANSMEM segments. β-barrels (TOM40, SAM50, MDM10, MDM34) have no TRANSMEM annotation and are added explicitly.
- **Assigned to soluble compartments only:** these take their membrane from GO, reading "cell periphery" as plasma membrane. Copies still without a membrane are dropped (9,946).
- **ER membrane proteins:** 15% go to the nuclear envelope, in proportion to area, because the outer nuclear membrane is continuous with the ER.
- **Lipid-droplet transmembrane proteins:** placed in the ER membrane.

### 4.4 Peripheral proteins

Soluble proteins annotated to a membrane are placed in thin "peripheral" shells next to it:
- **Plasma membrane, nuclear envelope, ER and lipid-droplet surface:** on the cytosolic side.
- **Inner-membrane-associated mitochondrial proteins:** in the matrix.

## 5. Structures

### 5.1 Bead representation and mass

- **Beads:** the Cα for amino acids, and P (or C4′) for nucleotides.
- **Mass:** for densities, a residue counts as 110 Da and a nucleotide as 330 Da.
- **What is missing from the beads:** residues missing from a structure (unresolved in the PDB entry, or trimmed from an AlphaFold model) are absent from the model and from its mass.

### 5.2 Choice of structure (in order of preference)

1. **Curated PDB complexes:**

   | Complex | PDB |
   |---|---|
   | 80S ribosome | 4V88 |
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

2. **Complex Portal complexes** (126 added, 22 of them membrane complexes) with a *S. cerevisiae* RCSB assembly. Rules:
   - **Acceptance:** the assembly must cover at least 80% of the subunit mass, and at least 60% of the residues must be modelled.
   - **RNA chains** are kept only where they belong to the complex: snoRNAs, snRNAs, telomerase RNA, and mitoribosomal rRNA.
   - **Copy number:** the median subunit pool, capped at 2× the limiting subunit. Complexes with four or more subunit types are capped at 2× the 25th-percentile subunit.
   - **Mitochondrially encoded subunits** (COB, COX1–3, ATP6, ATP8, OLI1, VAR1) are ignored when setting copy numbers, because mass spectrometry badly under-detects them.
   - **Subunit bookkeeping:** subunits are consumed from the monomer pool, so total protein copies are conserved.
3. **Homo-oligomers:**
   - **PDB:** about 130 homo-oligomeric assemblies.
   - **AlphaFold:** 452 AlphaFold-predicted homodimers with ipTM ≥ 0.5.
4. **Heterodimers:** 505 AlphaFold-predicted heterodimers, from 1,250 models. Requirements:
   - ipTM ≥ 0.5.
   - The two partners share a compartment. A membrane protein may pair with a soluble partner on its cytosolic or matrix side.
   - Not both partners are subunits of an already-modelled complex.
   - Neither partner forms a larger homo-oligomer.

   **Copy allocation:**
   - Each protein's free copies are split evenly among its accepted partners.
   - Pairs supported by Complex Portal or a direct IntAct interaction bind all available copies; the others bind 50%. This is an assumption: without affinity data, stable and transient interactions cannot be told apart.
5. **AlphaFold monomers** for everything else. Residues with pLDDT < 50 are trimmed, which removes 3.4% of placed protein mass overall. The loss is concentrated in nuclear and nucleolar proteins: 15–16% of mass there, measured on an earlier build.

**Result by mass of placed protein:** PDB 55.3%, AlphaFold complexes 14.3%, AlphaFold monomers 30.4%.

### 5.3 Membrane orientation

- **Method:** UniProt TRANSMEM residues are mapped onto the structure, through RCSB's UniProt alignment for PDB chains. A search finds the membrane normal and offset that put the most TM residues, and the fewest other residues, inside a ±1.6 nm slab.
- **Side:** from UniProt topology domains. Without them, the side carrying more soluble mass is taken as cytoplasmic.
- **β-barrels:** oriented by their β-strands, found from Cα geometry.
- **Check:** the method reproduces hand-oriented Pma1 and V-ATPase (z-coordinate correlation 1.000 and 0.9996).
- **Result:** in the packed cell, the median distance of TM residues from the membrane mid-surface is 0.73 nm, and 97.2% lie within 2 nm.

## 6. Supramolecular organisation built explicitly

These structures are built before random packing. They then act as fixed obstacles for everything packed afterwards.

### 6.1 Nuclear pores and FG mesh

- **Pores:** 117 copies of 7N9F, with the pore axis along the local nuclear-envelope normal.
- **FG mesh:** the disordered FG-repeat domains of the FG-Nups modelled in 7N9F are threaded as self-avoiding residue chains (0.38 nm per residue, low persistence) from the first modelled residue of each chain:

  | Nup | Residues | Copies per pore |
  |---|---|---|
  | NSP1 | 1–636 | 48 |
  | NUP57 | 1–286 | 32 |
  | NUP49 | 1–269 | 32 |
  | NUP159 | 388–1102 | 16 |

  - **Confinement:** the chains stay within 60 nm of their anchor, in the pore channel, cytoplasm or nucleoplasm.
  - **Result:** 6.1 of 6.9 million residues placed (89%).
  - **Missing:** FG-Nups absent from 7N9F (NUP100, NUP116, NUP145N, NUP42, NUP1, NUP60, NUP2), so the nuclear basket side is sparse.
  - **Caveat:** the chain conformations are excluded-volume random walks, not a physical model of FG cohesion.

### 6.2 Microtubules and the spindle pole body

- **Number:** 21 microtubules, made of 11,505 αβ-tubulin dimers:
  - 16 kinetochore microtubules, 229 nm mean length;
  - 2 other nuclear microtubules, 311 nm;
  - 3 cytoplasmic microtubules, 1.0 µm.
- **Lattice:** an ideal 13-protofilament, 3-start lattice with a seam. The lateral geometry comes from 8QV3 (native yeast γ-TuRC-capped microtubule, cryo-ET). The dimer is 5W3F, with an 8.3 nm axial repeat.
- **Caps:** each minus end has a γ-tubulin ring of 7 γ-TuSCs from 8QV3.
- **Growth:** from 90 nm plaques at the spindle pole body. Microtubules bend away from membranes, with a minimum radius of curvature of 300 nm.
- **Assumptions:** the lengths (G1-like), and the 70% cap on the use of the tubulin pool.
- **Not built:** the spindle pole body's own plaques (Spc42, Spc72 layers).

### 6.3 Chromatin

- **Walks:** 16 chromosomes are self-avoiding, persistent nucleosome walks (12.5 nm step), with lengths proportional to chromosome size.
- **Rabl-like layout:** each chromosome starts within 20 nm of its kinetochore microtubule tip.
- **Result:** 73,150 nucleosomes in 59 fragments; walks that get trapped restart elsewhere.
- **Linker DNA:** joins consecutive nucleosomes, with a schematic path.
- **Not modelled:** chromosome territories, telomere anchoring, the rDNA array, nucleosome-depleted regions and transcribing Pol II.

### 6.4 Actin

The ACT1 pool (111,000) is split by assumption: 30% in cables, 30% in patches, and the rest monomeric or in profilin, cofilin and Srv2 heterodimers.

**Cables:**
- **Number:** 8 cables of 8 filaments, 33,300 subunits and 92 µm of filament in total.
- **Position:** 90 nm beneath the plasma membrane, along great circles.
- **Bundle:** filaments sit on a hexagonal lattice 12 nm apart (assumed fimbrin bundle spacing).
- **Filament structure:** yeast F-actin 9GO5, 167° and 2.76 nm per subunit.
- **Tropomyosin:** 6,610 TPM1 coiled-coil dimers, one per 5 actins per strand. The path is transferred from the mammalian actin–tropomyosin structure 5JLF, because no yeast structure exists. The chains are wound as an ideal coiled coil.

**Endocytic patches:**
- **Number:** 22 patches of about 1,500 actin each.
- **Branching:** dendritic networks branch at 73° from Arp2/3 complexes (2,169 branches). The branch-junction geometry is 8UXW. The yeast Arp2/3 subunits are AlphaFold models superposed on the fission-yeast template, 32–74% identity.
- **Capping:** capping protein (from 8F8Q) caps filaments after 20 ± 6 subunits. Filaments that end against an obstacle stay uncapped.
- **Orientation:** barbed ends face the membrane.
- **Collision checks:** residue by residue inside each patch.

Not built: fimbrin cross-links, septins, formins.

### 6.5 Ribosomes and mRNA

- **Polysomes:** 70% of cytosolic ribosomes are in polysomes (assumption), with sizes 1 + Poisson(4).
  - **Chain geometry:** each ribosome's mRNA entry lies 3–12 nm beyond the previous one's exit, within a 40° cone. Neighbouring ribosomes end up 25–31 nm apart, as in cryo-ET polysomes.
  - **mRNA channel:** the path through each ribosome comes from 7ZW0, superposed on 4V88.
  - **Result:** 93,000 ribosomes on 23,000 mRNAs.
- **Full-length mRNA:** all 30,000 mRNAs are threaded at full length through the cytosol before protein packing.
  - **Lengths:** lognormal, mean 1,500 nt, σ = 0.5 in log space, minimum 300 nt.
  - **Representation:** one bead per nucleotide, 0.6 nm apart.
  - **Polysome tails:** the untranslated remainder leaves the first and last ribosome as 3′ and 5′ tails. The 5′ side gets 10–40% of it.
  - **Result:** 86% of the target nucleotides were placed. Threads that run out of room end early.
- **ER-bound ribosomes:** 4,264 ribosomes are docked on translocons in the ER and the outer nuclear membrane.
  - **Templates:** 3,349 on Ssh1 (2WW9, yeast 80S–Ssh1) and 915 on Sec61 (7KAH, placed through the shared Sss1 subunit).
  - **Assumption:** 80% of free translocons are ribosome-bound.
  - **Not modelled:** ER-bound polysomes and their mRNA.
- **Not placed on ribosomes:** the ribosome-associated factors NAC, RAC/Ssb and eEF1A are free in the cytosol.

### 6.6 Mitochondrial inner membrane

- **ATP synthase dimers:** all 11,762 ATP synthase copies form 5,881 dimers. Each dimer is two copies of the yeast monomer 8F39 arranged as in the cryo-ET dimer 4B2Q (dimer angle 83.5°).
  - **Placement:** in rows along the crista rims, every 13 nm (assumed spacing).
  - **Missing subunits:** 8F39 lacks the dimer-specific subunits e and k.
- **Respiratory supercomplexes:** complex IV is entirely in 1,561 III₂IV₂ supercomplexes (6HU9). The remaining 1,422 complex III dimers are free.
- **Not modelled:** III₂IV₁ supercomplexes, and confinement of supercomplexes to the flat crista membranes.

## 7. Procedural content (generated in the browser)

These are generated deterministically from fixed seeds when you zoom in. They are not stored coordinates and do not count in the densities.

- **Lipids:**
  - **Packing:** 0.65 nm² per lipid per leaflet.
  - **Composition:** whole-cell, from Ejsing *et al.* (2009): PI 20.3%, PE 14.9%, PC 14.3%, PA 8.3%, PS 1.8%, ergosterol 12%, other 28.4%. The same composition is used for every membrane, with no leaflet asymmetry and no organelle-specific lipids (e.g. cardiolipin).
  - **Masking:** lipids within 0.95 nm of a protein bead are removed.
- **Cell wall:**
  - **Layers:** an inner β-1,3-glucan layer (semi-flexible chains, 0.018 beads/nm³) and an outer mannan layer (branched brushes on the mannoproteins).
  - **Amount:** bead counts are set so the wall is about 20% of cell dry mass.
  - **Collisions:** only partly checked.

## 8. Packing

### 8.1 Spatial decomposition

- **Tiles:** the cell is cut into 200 nm cubic tiles; 6,416 of them intersect the cell.
- **Copies per tile:** each species' expected copy number per compartment is split over tiles in proportion to the compartment's volume (or membrane area) in each tile. A multinomial draw with stochastic rounding conserves the cell totals.
- **Concentrations:** uniform within each compartment.

### 8.2 Algorithm

Random sequential addition, run independently in each tile.

**Placement:**
- **Order:** membrane proteins first, then the rest by decreasing size.
- **Soluble molecules:** a uniformly random position in the target compartment and a uniformly random orientation.
- **Membrane proteins:** projected onto the exact membrane surface, with their membrane normal along the local surface normal and a random spin about it.
- **Retries:** a bounded number of trials per copy. Copies that don't fit after that are left out.

**Collision checks:**
- **Sphere trees:** each structure is approximated by k-means clusters of about 30 residues, each sphere padded by a 0.42 nm bead radius. Spheres may interpenetrate by up to 0.3 nm.
- **Compartments:** a copy must lie entirely within its compartment and clear of membranes, which are treated as 2 nm half-thickness slabs.

**Parallelism:**
- **Phases:** tiles are packed in 8 checkerboard phases, so neighbouring tiles are never packed simultaneously.
- **Halos:** each tile sees the already placed molecules of its neighbours through a halo.
- **Obstacles:** the pre-placed structures of section 6 are fixed obstacles.

### 8.3 Consequences

- **The density is not an equilibrium value.** Random sequential addition jams at lower density than real cytoplasm, and the outcome depends on insertion order. Where the requested copies exceed what random sequential addition can fit, the shortfall is reported, not hidden.
- **The shortfall is concentrated in mitochondria.** Placed copies are 97.9% of requested overall and 99.8% outside mitochondria. In the mitochondrial matrix only 65% fit, so matrix protein density is under-represented.
- **No interaction-driven organisation.** There are no metabolons and no condensates; enzymes are not co-localised beyond what the explicit structures in section 6 provide.

## 9. Validation

| Check | Result |
|---|---|
| Residue-bead clashes (beads of different molecules closer than 0.3 nm) | 1.4 × 10⁻⁶ of 2.9 M sampled beads in 64 test cubes; 0 between tiles |
| Containment (beads outside their compartment) | 4.3 × 10⁻⁵ of 17 M beads |
| Membrane-protein orientation vs local surface normal | median 0.0°, 100% within 1° |
| TM residues within 2 nm of the membrane mid-surface | 97.2% (median 0.73 nm) |
| ATP synthase dimer TM residues within 2 nm of the mid-surface | 98.8% |
| Placed / requested copies | 97.9% overall; 99.8% outside mitochondria; matrix 65%; nucleolus 95% |
| Inventory protein mass vs proteome | within 0.4% before packing |

## 10. Densities, mass and volume fractions

### 10.1 Macromolecular densities

Densities are the mass of placed macromolecules per compartment volume. They count protein plus the RNA of ribosomes, tRNAs and complexes, the DNA of nucleosomes, and mRNA. Lipids, polysaccharides and small molecules are not counted.

| Compartment | Protein (mg/mL) | RNA/DNA (mg/mL) | Total (mg/mL) | Dry volume fraction |
|---|---|---|---|---|
| Cytosol | 67 | 22 | 89 | 6.1% |
| Nucleoplasm | 56 | 6 | 62 | 4.4% |
| Nucleolus | 87 | 6 | 94 | 6.7% |
| Mitochondrial matrix | 108 | 4 | 113 | 8.1% |
| Mitochondrial intermembrane space | 27 | 0 | 27 | 1.9% |
| ER lumen | 32 | 0 | 32 | 2.3% |
| Golgi / endosomes | 64 | 0 | 64 | 4.7% |
| Peroxisomes | 83 | 0 | 83 | 6.1% |
| Vacuole | 5 | 0 | 5 | 0.4% |

- **Dry volume fraction:** computed with partial specific volumes of 0.73 mL/g for protein and 0.55 mL/g for nucleic acid.
- **Hydrated volume fraction:** with a hydration shell of about 0.3–0.4 g water per g, the hydrated volume fraction is roughly 1.5× the dry one: about 9–10% in the cytosol.
- **Excluded volume:** the volume excluded to a diffusing particle is larger still, and grows with the particle's size. It can be computed from the coordinates as accessible volume versus probe radius; this has not been done yet.

### 10.2 Whole-cell totals

| Quantity | Value |
|---|---|
| Protein in the inventory | 2.7 pg |
| Protein placed | 2.2 pg |
| RNA and DNA placed | 0.64 pg |
| mRNA | 0.02 pg |

Placed protein is lower than the inventory for four reasons:
- residues not resolved in PDB entries;
- AlphaFold residues trimmed for pLDDT < 50 (3.4% of mass);
- mitochondrial copies that did not fit;
- dropped and structure-less copies.

With full-length sequences, cytosolic protein would be about 80–85 mg/mL.

### 10.3 Comparison with dry-mass measurements

Measured cytoplasmic dry-mass densities are around 200 mg/mL, more than twice the macromolecular density here. These are rough estimates of the missing components, averaged over the cytoplasm (cytosol plus organelles, excluding the vacuole):

| Component | Estimate (mg/mL) | Basis |
|---|---|---|
| Macromolecules as placed | about 90 | This model |
| Full-length protein correction | about +15 | Section 10.2 |
| Lipids of internal membranes and lipid droplets | about 20–25 | About 130 µm² of internal membrane plus 0.4 µm³ of lipid droplets in this model, about 0.8 pg of lipid |
| Small molecules and ions | probably 30–40 | Free amino acids, nucleotides, glutathione, sugar phosphates, phosphate and K⁺, a few hundred mM in total |
| Glycogen and trehalose | low in log phase on glucose | Can exceed 10% of dry weight after the diauxic shift |

Together these give about 150–160 mg/mL. The remainder is within the uncertainty of the absolute proteome size (section 2.1). A higher protein total than 4.2 × 10⁷ molecules would close most of the gap.

## 11. Known limitations

| Area | Limitation |
|---|---|
| Physics | Static random sequential addition: no interactions, no dynamics, random orientations, and densities below the jamming limit of real cytoplasm. |
| Abundances | About 2-fold per-protein uncertainty. The proteome total is probably low. Ribosomes are probably low (169,000 vs about 3 × 10⁵ estimated). Fixed counts for tRNA and mRNA are assumptions. |
| Geometry | Illustrative layout. Nucleus smaller than in EM. A single mitochondrial reticulum. No bud, so no cell-cycle-dependent structures. |
| Mitochondria | Matrix only 65% filled; mitochondrially encoded subunits are inferred from their complexes. |
| Disorder | AlphaFold residues with pLDDT < 50 are deleted rather than represented, except the FG domains of 4 Nups. |
| Localization | Bulk steady state from annotations. 625 proteins default to the cytosol. Dual localizations are weighted by rank, not measured. |
| Nucleolus | No rDNA, Pol I transcription units or pre-rRNA; few pre-ribosomes. |
| Cytoskeleton | No septins, fimbrin, formins or spindle-pole-body plaques. Actin partition and filament lengths are assumptions. |
| Lipids and wall | One lipid composition for all membranes. Procedural cell wall. |
| Missing classes | Metabolites, ions, water, glycogen, trehalose, polyphosphate. |

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

### 12.2 Literature behind the parameters

| Reference | Used for |
|---|---|
| Ho B, Baryshnikova A, Brown GW (2018). *Cell Systems* 6:192–205 | Total protein molecules per cell |
| Jorgensen P *et al.* (2002). *Science* 297:395–400 | Cell volume |
| Uchida M *et al.* (2011). *Yeast* 28:227–236 | Organelle volume fractions (soft X-ray tomography) |
| Yamaguchi M *et al.* (2011). *Microscopy* 60:321–335, doi:10.1093/jmicro/dfr052 | Serial-section EM comparison |
| Orlean P (2012). *Genetics* 192:775–818 | Cell-wall thickness and composition |
| Wolf W *et al.* (2012). *PLoS ONE* 7:e39703 | Cortical ER coverage |
| Winey M *et al.* (1997). *Mol. Biol. Cell* 8:2119–2132 | NPC number |
| Ejsing CS *et al.* (2009). *PNAS* 106:2136–2141 | Lipid composition |
| Huh WK *et al.* (2003). *Nature* 425:686–691 | GFP localization |
| Sickmann A *et al.* (2003), *PNAS*; Reinders J *et al.* (2006), *J. Proteome Res.* | Excluded high-throughput organelle calls |
| Dolgalev G *et al.* (2023). *Int. J. Mol. Sci.*, doi:10.3390/ijms24032081 | Protein density comparison |
| Kostinski S, Reuveni S (2021). *Phys. Rev. Research*, doi:10.1103/physrevresearch.3.013020 | Ribosome count comparison |
| Maritan M *et al.* (2022). *J. Mol. Biol.* 434:167351 | Whole-cell modelling approach |

### 12.3 PDB entries used for assemblies and templates

| Assembly or template | PDB entries |
|---|---|
| Ribosome and polysomes | 4V88, 7ZW0 |
| ER translocons | 2WW9, 7KAH |
| Nuclear pore and nucleosome | 7N9F, 1ID3 |
| Microtubules | 8QV3, 5W3F |
| Actin and actin-binding proteins | 9GO5, 5JLF, 8UXW, 8F8Q |
| Mitochondrial inner membrane | 8F39, 4B2Q, 6HU9 |
| Other curated complexes | 4CR2, 2UV8, 4V94, 4C2M, 1WCM, 5FJ8, 7VH5, 6O7V |

The other complexes and homo-oligomers are listed per species in the viewer's inspect panel: click any molecule.

## 13. Reproducibility and data access

- **Pipeline.** The model is built by a Python pipeline (numba-accelerated) in five stages: inventory, geometry, packing, validation and encoding. A full rebuild takes roughly 10–15 minutes on a recent multi-core desktop or laptop, with fixed seeds throughout. The pipeline code is not part of this repository.
- **This repository** holds the built viewer and the encoded data.
  - The tiles under `data/tiles/` store every instance as a species index, a position quantised to about 0.0055 nm and an orientation quantised to about 0.004°.
  - The residue-bead library is in `data/library/`.
- **Full-atom coordinates.** The viewer fetches full-atom models live from RCSB and AlphaFold DB for display in Mol*.
