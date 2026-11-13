# Model accuracy: improvement backlog

Measured 2026-11-03 against the current build: 25.2 M instances, 5,483 types. The scripts and raw output are in `yeast_cell_dev_bundle/analysis/`. Rerun `analyze.py` and `rcsb_coverage.py` with `yeast_cell_dev_bundle/.venv/bin/python` to refresh the numbers. "Mass" means the share of all placed protein mass (copies × mass).

## H. Data-driven realism: system-level properties against real data (plan of 2026-11-11)

**Direction (the user and a collaborator, 2026-11-11).** Match system-level properties of the cell to the structural model, bottom-up and against real data:
- recover known organisation emergently (P-bodies and other membrane-less organelles, metabolons, the nucleolus, membrane-adjacent proteins) from weak interactions, scaffolds, RNA binding and PTMs, using annotations, prediction datasets and the AlphaFold 3 pool scores;
- compare protein localization systematically with microscopy, find systematic errors, fix them, and keep cases out as a test;
- calibrate the simulation's diffusion on literature measurements by compartment, with per-protein parameters (size, charge, modifications, interactions) and held-out cases.

**Rules fixed before any tuning (pre-registration).**
- **Circularity:** CYCLoPs (Koh *et al.* 2015) and the Huh *et al.* 2003 GFP calls are inputs of the localization. Comparing against them is a fit, never a test. Tests use data the model has not seen: N-terminally tagged SWAT libraries and organelle spatial proteomics (hyperLOPIT).
- **Held-out localization split:** proteins are split by `sha256(ORF) mod 10`: 0–4 dev, 5–9 test. Error analysis and rule changes use only dev. Test numbers are computed once per finished change and reported next to dev.
- **Held-out diffusion split:** measurements are split by paper (`sha256(DOI) mod 10 < 4` → test, stratified so that every compartment keeps a training paper). Global parameters are fitted on train only.
- **Held-out condensates:** the nucleolus is the only condensate that may be used to set the interaction scale. CD-CODE P-body membership is a pure test set and is not looked at until the sampler is frozen. Stress granules (and hypoxia G-bodies, starvation filaments) must *not* form in the log-phase cell: negative controls. Their appearance after an *in-silico* stress (polysome run-off) is a prediction.
- **Emergence vs built in:** no condensate or localization label enters the energies. Every positive result is compared with degree-preserving scrambles of the interaction network and with ablations (no RNA term, no IDR term, no pool term).
- **Model errors kept apart from packing shortfalls:** the localization benchmark is computed both on the inventory's assigned shares and on the placed copies (tiles), so the known L1 regressions (nucleolus 77%, Golgi 81%, matrix 40% placed; free histones; 44% polysomes) show up separately from localization errors.

**Items** (ranked; each finished one gets a Progress entry):

| # | Item | Question and data | Method | Held-out test and success metric | Compute here | Risks |
|---|---|---|---|---|---|---|
| H1 | Localization benchmark | Where does the model put each protein, against independent microscopy (SWAT N-GFP, Weill *et al.* 2018 / Yofe *et al.* 2016) and spatial proteomics (hyperLOPIT, Nightingale *et al.* 2019)? | Coarse classes (cytosol, nucleus, nucleolus, ER/NE, mitochondrion, Golgi/endosome, vacuole, PM/periphery, peroxisome, lipid droplet); model shares from the inventory and from placed copies | Test half of the ORF split. Top-1 agreement, per-class precision/recall, macro-F1, with bootstrap CIs; CYCLoPs agreement reported as the training fit | minutes, CPU | reference datasets disagree with each other (agreement between them is the ceiling) |
| H2 | Diffusion benchmark | Measured D of proteins, complexes and particles in yeast by compartment and method (FRAP, FCS, SPT) | Literature table with construct, tag, mass, method, condition; every value quoted from the paper | Test papers. RMSE of log10 D, Spearman within compartment; nulls: current model, Stokes–Einstein × one constant | none | few S. cerevisiae values per compartment; method biases (FRAP vs FCS vs SPT) |
| H3 | GPU Brownian dynamics with many-body hydrodynamics | Does full hydrodynamics close the 5× diffusion gap, as Ando & Skolnick 2010 found for *E. coli*? | `simgpu/`: CuPy rigid-body BD of packed boxes; Ewald-summed polydisperse RPY; sweeps over 4 GPUs | Engine checks (Hasimoto law, browser agreement); then H2's test set | ~1 GPU-min per µs of a 120 nm box | far-field RPY lacks lubrication; near-field hydrodynamics on shaped bodies is ill-defined |
| H4 | Per-protein diffusion model | Why do some proteins diffuse slower than their size predicts? | D_i = D0(R_h incl. IDRs) × crowding(R_h) from BD × hydrodynamic factor × 1/(1 + binding to slow partners: ribosomes/mRNA via charge and RNA-binding evidence, pool-score partners); 2–4 global parameters | H2 test papers; must beat both nulls | CPU + BD sweeps | over-fitting with few measurements; tag effects |
| H5 | Emergent condensates | Do P-bodies (and the nucleolus) emerge from generic, data-derived interactions at real concentrations, while stress granules do not? | Interaction energies from pool PAE/ipTM (including weak pairs), RNA-interactome capture, IDR stickers (aromatic/Arg, RG/RGG, Q/N) and charge; multicomponent Flory–Huggins stability analysis over all species, then a GPU particle sampler for the unstable modes to get partition coefficients and droplets | Calibrate the scale on the nucleolus only. Test: P-body AUROC of partition coefficients vs other cytosolic proteins > 0.7 and above scrambles; SG/G-body AUROC ≈ 0.5 in log phase | GPU hours | mean-field misses valence; pool scores are noisy for weak pairs; P-bodies are small in log phase |
| H6 | Membrane-adjacent proteins from lipidation and lipid-binding domains | Are N-myristoylated, palmitoylated, prenylated and PH/PX/FYVE/C2/BAR proteins at membranes in the model? | Membrane affinity from lipid anchors and domains; peripheral-shell shares replaced by a fitted bound fraction | Test half of the ORF split: periphery/punctate calls of SWAT and hyperLOPIT PM/organelle classes | minutes; s1–s3 rebuild | membrane identity of the anchor (PM vs Golgi vs vacuole) needs rules |
| H7 | Localization fixes | Systematic errors found on dev (H1): default-cytosol proteins, GO fallbacks, rank weights, packing shortfalls | Rule changes in `inventory.py`, checked on dev, then test | H1 metrics on test must improve | s1–s5 rebuild (~30 min on 64 cores, to be measured) | fixes that help dev but not test are reverted |

**How the 4 GPUs are used.** The browser engine runs headless on a 4090 but is dispatch-bound (about laptop speed). Sweeps use `simgpu` (CuPy, one process per GPU: `CUDA_VISIBLE_DEVICES=k`); OpenMM CUDA runs exported boxes (fixed 2026-11-11). The 4090s are weak in double precision, so everything runs in float32 except the mobility assembly in tests.

**What can realistically emerge.** The nucleolus (rRNA plus GAR/RGG proteins at high concentration): likely, but it is the calibration case. P-bodies: plausible if the decapping/Lsm/Dhh1/Pat1 network has strong pool scores and RNA binding, but they are small and few in log phase, so the metric is enrichment, not droplet count. Membrane-adjacent proteins: yes, by construction of the energies, so the test is *which* proteins, measured against held-out microscopy. Metabolons: glycolytic G-bodies are hypoxia-induced and should not form; a weak co-enrichment of glycolytic enzymes would be the most to expect.

## I. Imaging, proximity and pair statistics; fixes from a second literature round (plan of 2026-11-11, afternoon)

**Direction (the user, 2026-11-11).** A new literature round on yeast cell biology and microscopy data; simulate fluorescence microscopy from the model and compare with fluorescent-protein images; compare with proximity data (BioID/TurboID, LUX-MS and similar); extract pair statistics (how often A meets B) for kinetic models such as kinase signalling; a validation document with figures; a list of data that could not be retrieved.

**Reviews** (three parallel agents, all citations checked in Europe PMC/Crossref): `literature_reviews_2026-11-11/` (cell biology, microscopy data, proximity data). Downloads in `pipeline/data/validation/{proximity,microscopy,literature,biogrid}/` with READMEs. Unretrievable data: `data_requests.md`.

**Rules** as in section H: nothing tuned on test data; every comparison labelled test, check or fit; negative results reported.

| # | Item | Data | Status |
|---|---|---|---|
| I1 | Literature and data reviews; data request list | Europe PMC, BioGRID, IDR, CZ cryoET portal, OSF | done |
| I2 | Simulated fluorescence microscopy (`pipeline/microscopy.py`): Debye PSF, widefield/confocal, z-stacks, spot and object counts | gallery; idr0078 screen (Mattiazzi Usaj 2020) | done; next: nuclear/cytoplasmic ratio and intensity vs held-out SWAT images (YeastRGB), texture vs DeepLoc crops |
| I3 | Pair-contact atlas of the packed cell (`pipeline/contacts.py`, GPU sphere-tree gaps 2/5/10 nm) | — | done |
| I4 | Proximity benchmarks (`validation/proximity_test.py`, `ppiseq_test.py`) | BioGRID PCA, proximity labelling, XL-MS, kinase–substrate; PPiSeq with negatives | done; next: TurboID with controls (Favre 2026), cross-link distances on placed copies |
| I5 | Pair statistics for kinetic models (`pipeline/encounters.py`) and BD check of encounter rates (`simgpu/encounter_bd.py`) | — | done |
| I6 | Validation document `model_validation_and_tuning.md` (+ PDF, `scripts/md_to_pdf.py`) | all of the above and earlier rounds | done; regenerate after each round |
| I7 | Ribosomes vs in-situ cryo-ET (`validation/ribosome_cryoet.py`) | Delarue 2018; CZ cryoET 10007 | done: 17 vs 23 µM, spacing right, polysomes too regular |
| I8 | **Fix: cytosolic faces of vacuole, Golgi and peroxisomes** (soluble proteins without signal peptide were packed in the lumen: V1, EGO, TORC1, COPI, clathrin, Arf) and no wall share without a signal peptide (Ssa1/2) | UniProt signal peptides and locations | done (rebuilt); Golgi face 73% placed |
| I9 | **Fix: 26S proteasome nucleus/cytosol split from its subunits** (was 50/50; mostly nuclear, Laporte 2008) | CYCLoPs via subunits | done (81% of Rpn1 nuclear) |
| I10 | Eisosome coats hold 7% of Pil1/Lsp1 (microscopy: Pil1 should be punctate) | Pil1/Lsp1 copies; Kefauver 2024 lattice | open: larger or more coats, or the surplus held next to furrows |
| I11 | Spindle pole body plaques (Spc42 diffuse in the model, one spot in cells) | SPB composition | open |
| I12 | Kinase localization test set: TORC1 and Sch9 on the vacuole, Bcy1 and Pho80–Pho85 nuclear, Slm1/2 and Pkh1/2 at eisosomes | review Part D | open (partly addressed by I8) |
| I13 | Polysome compaction and ribosome concentration (cryo-ET: 23 µM, tighter and more variable spacing; ~14.5% tightly packed) | Delarue 2018; Goetz 2025 (preprint); Xie 2024 | open |
| I14 | Membrane-protein diffusion: one factor per membrane | Valdez-Taubas 2003; Meinema 2013 | done (calibrated; fit) |
| I15 | Nucleoplasm too dilute relative to the cytosol (soft X-ray absorption 0.90 vs model density 0.64) | Chen 2025 | open |
| I16 | Endocytic sites twice too many in simulated images (polarization to the bud missing) | idr0078 | open |
| I17 | ER band ribosomes: docked ones match cryo-ET (7,715 per cell; Wu 2026), the extra free ones are doubtful | Wu 2026; Jan 2015 | open |
| I18 | ER–mitochondria contacts (ERMES), nucleus–vacuole gap 21 nm (model builds within 45 nm) | Collado 2019; Wozny 2023 | open |

## K. Chromatin and the genome (plan of 2026-11-11, evening)

**Direction (the user, 2026-11-11).** More accurate chromatin: map chromosome positions to genes, and use accessibility and other epigenetic readouts to build the structure.

**Now** (MODEL.md 6.3): 32 arms as persistent self-avoiding nucleosome walks (12.5 nm steps) from the centromere at the SPB to a telomere focus; kb positions mapped linearly along each walk; 552 engaged Pol II on their genes (Pelechano 2010); active rDNA in the nucleolus. Not modelled: nucleosome positions on the sequence, promoter nucleosome-depleted regions, histone variants and marks, folding measured by Hi-C, DNA-bound proteins other than Pol I/II, kinetochores, origins, condensin/cohesin, tRNA gene clustering.

**First measurement** (`validation/chromatin_stats.py`): the walks are too extended at short range and too compact at long range against in-vivo imaging (Bystricky et al. 2004: worm-like chain, persistence 170–220 nm, 110–150 bp/nm): rms distance at 20 kb 214 nm (imaging 118–160), at 100 kb 347 nm (412–552). Contacts within 30 nm fall as s^-0.89 over 10–200 kb; 66% of contacts are within one arm.

**Ideas** (ranked; data survey running, `literature_reviews_2026-11-11/review_chromatin.md`):

| # | Idea | Data | Use | Test (held out) |
|---|---|---|---|---|
| K1 | Genome coordinates on every nucleosome, gene and Pol II (sacCer3); viewer search for a gene → its locus; colour chromatin by gene, expression, mark | SGD coordinates (have), nucleosome maps | annotation, viewer | — |
| K2 | Nucleosomes at measured dyad positions; promoter nucleosome-depleted regions; +1/−1 nucleosomes with H2A.Z | chemical map (Brogaard 2012), MNase consensus (Jiang & Pugh 2009), H2A.Z ChIP | input | ATAC / MNase accessibility, in-situ cryo-ET spacing |
| K3 | Fibre compaction from physics: irregular nucleosome clutches with the measured mass density and persistence length, instead of a 12.5 nm walk | Bystricky 2004 (fit); in-situ cryo-ET of yeast nuclei | input (compaction) | Hi-C / Micro-C P(s), locus distances |
| K4 | 3D folding: Rabl + confinement + rDNA/nucleolus + telomere anchoring as a polymer ensemble, compared with Hi-C (kept as test); optionally Hi-C-restrained models with held-out chromosomes | Duan 2010 Hi-C; Micro-C (Hsieh 2015; Ohno 2019) | test (or input with held-out chromosomes) | contact maps, centromere/telomere clustering by arm length |
| K5 | DNA-bound proteins at their sites: transcription factors, remodelers, GTFs from ChIP-exo of ~400 proteins; bound vs free fractions | Rossi et al. 2021 | input | nuclear diffusion benchmark (bound fractions, SMT), proximity in the nucleus |
| K6 | Landmarks: 16 kinetochores on centromeres; MCM double hexamers on licensed origins (the cell is in G1); ORC; condensin at the rDNA; tRNA genes and Pol III near the nucleolus; Sir proteins at silent loci and telomeres | OriDB, ORC/MCM ChIP, kinetochore structures, tRNA coordinates | input | locus imaging (tRNA clustering, telomere foci) |
| K7 | Histone marks as annotations that place readers: H3K36me3 (Rpd3S) over gene bodies, H3K4me3 at 5′ ends, acetylation (SAGA/NuA4) at promoters | genome-wide ChIP (e.g. Weiner 2015) | input (placement), viewer colour | proximity of readers, ChIP of readers |
| K8 | Gene positioning: highly expressed and induced genes at nuclear pores; silent loci at the periphery | locus imaging datasets | input or test | locus–NE distances |

**Status** (2026-11-11, late): K1 done (`genome_map.py`, `library/gene_loci.json`, viewer gene search). K3 done (`pack3d.fiber_walks`); the first 20 nm fibre failed Hi-C (cis 0.62–0.66 vs 0.37–0.41, P(s) slope −1.6), so the radius is 60 nm, inside the Arbona 2017 R(s) envelope (Hi-C now partly tuned). K4 done as a test (`validation/hic_test.py`: Lazar-Stefanita 2017 G1 and Swygert 2019 log-phase Hi-C). K2 done: the fibre is built from dyads called on the Chereji 2018 H3 chemical-cleavage maps (`nucleosome_map.py`; Brogaard 2012 is on sacCer2, Jiang & Pugh 2009 on a pre-sacCer3 release), with per-nucleosome backbone advance and linker reach; H2A.Z not yet. K6 partly done: MCM double hexamers and ORC on the 410 confirmed OriDB origins (`origins.py`); kinetochores (CCAN + 6 MIND on each kMT plus end, `kinetochores.py`); condensin, cohesin, tRNA genes and Sir proteins open. K5 done (`chromatin_sites.py`: ~280 proteins at Rossi 2021 ChExMix peaks, 30% bound share assumed; 36% of the wanted copies fit in build M). K7 partly done: seven marks and H2A.Z per nucleosome as annotation (`genome_map.marks`, Weiner 2015 t = 0); reader placement and viewer colouring open. K8 open. Held-out nucleosome test `validation/nucleosome_test.py` (Oberbeckmann 2019 occupancy and count, Weiner 2015 MNase, Schep 2015 ATAC); the map is NDR-filtered with Rossi 2021 promoter NDRs.

## M. mRNA (plan of 2026-11-11, evening)

**Direction (the user, 2026-11-11).** Make the mRNA more accurate. Survey: `literature_reviews_2026-11-11/review_mrna.md` (agent; downloads in `pipeline/data/validation/mrna/`).

| # | Item | Data | Use | Test (held out) | Status |
|---|---|---|---|---|---|
| M1 | Gene-specific mRNAs: copies per gene, length from CDS + UTRs + poly(A) | Pelechano 2010 (copies), SGD | input | smFISH counts per gene | done (`mrna_genes.py`) |
| M2 | Per-gene UTRs, poly(A) tails and translational efficiency | TIF-seq (Pelechano 2013), PAL-seq and ribosome profiling (Subtelny 2014) | input | Arava 2003 ribosome density per gene | done; Arava data not yet found |
| M3 | Polysome linkers as flexible single strands (Gaussian chain of the free nucleotides) | — | geometry | cryo-ET tight fraction 14.5–18.6% | done; tight fraction 7.0% after random spacing (round N); stacked orientations missing |
| M4 | mRNAs translated at the mitochondrial surface | Williams 2014 (−CHX calls, glucose) | input | Saint-Georges 2008 MLR % (fractionation, galactose) | done (`mrna_localization.py`, `validation/mrna_localization.py`) |
| M5 | ER-translated (secretory) mRNAs as membrane-bound polysomes that include the translocon-docked ribosomes | Jan 2014 (Sec63/Ssh1 proximity profiling) | input | Chartron 2016 SRP-selective profiling | done in round N (tethered to the cortical ER; Chartron test Spearman 0.35, AUROC 0.91); the docked ribosomes are still single |
| M6 | mRNP proteins on the threads: Pab1 on poly(A), eIF4F at the cap (closed loop), compaction | Pab1 footprint, smFISH 5′–3′ distances | input | smFISH end-to-end distances | Pab1 and eIF4E on the ends done (`mrnp.py`); eIF4G/4A, closed loop and compaction open |

## N. Ideas from the Jurkat cell model, and the open items of K and M (plan of 2026-11-11, night)

**Direction (the user, 2026-11-11).** Continue the open items, and take what is valuable from the Jurkat T-cell model (OpenOrganelle jrc_jurkat-1, model v7; built on a fork of this pipeline). Report reviewed: `Jurkat Cell 1 Model.zip` (project folder). Not applicable to yeast: EM geometry, lamina, speckles, CTCF, H1.

| # | Item | Jurkat precedent | Use / test | Status |
|---|---|---|---|---|
| N1 | Attempt budget in thin lumens and shells (matrix 62%, ER periphery 48%, Golgi face 71%, nucleolus 78%) | v7: 10x tries in thin lumens; matrix 71% to 86%, ER lumen 66% to 82% | placement; acceptance 99% outside mitochondria | done: 2,000 tries / 1,200 failures for all regions (thin-only stole room from the cytosol), flat alignment in shells (`wc3d.PACK_TRIES`, `FLAT_MEM`) |
| N2 | Chromatin from a coarse polymer (OpenMM; Rabl tethers, nucleolus, confinement, excluded volume), measured nucleosomes back-mapped along it | 10 kb polymer, back-mapped; P(s) slope from flat to -1.29 | Hi-C (held out apart from the fibre radius choice): cis share, P(s), chromosome pairs (Tjong 2012: tethering + exclusion predict them) | tested, not adopted (`chromatin_polymer.py`, `pack3d.guided_walks`, `validation/chromatin_polymer_scan.py`): P(s) better, cis share and CEN clustering worse in all variants; the fibre walk stays. Found on the way: the two Hi-C maps disagree with each other on chromosome pairs (Spearman -0.04 as provided, 0.15 balanced), and the data measures shift with normalisation; hic_test now also scores ICE-balanced maps |
| N3 | Locus positions against live-cell imaging | territory order by gene density (held out) | Berger 2008 gene territories (relative to nucleolus and envelope) | open: data not in hand (data_requests.md) |
| N4 | Distance-defined regions for structure-bound proteins: SPB, microtubule plus ends, lattice band, kinetochore pools | pericentriolar sphere, +TIP and MAP bands | simulated foci (Spc42 diffuse: I11), localization | done: SPB plaques (`spb.py`), outer kinetochore DASH ring and Ndc80 (`kinetochores.outer`); plus-end and lattice bands not done |
| N5 | Transport receptors in the pore channel | 20 per pore within 25 nm of the axis | NPC composition | done (`pore_receptors.py`, before the FG mesh) |
| N6 | Localization scored on shares, mass-weighted; image N/C ratio on a held-out half | Tier 1 and 2 scores | localization benchmark | share agreement (plain and mass-weighted) done; image N/C test open |
| N7 | Per-gene ledger of placed shares across builds | provenance ledger, mass flow | documentation | done (`validation/ledger.py`; snapshots from build M on) |
| N8 | Diffusion with binding to partners, weak stickiness, chromatin binding; train / held-out split | held-out within 3x: 24% to 65% | diffusion benchmark | not redone: a partner-binding term was already tried and did not help on the held-out rows (2026-11-11); the benchmark (36 rows) is too small to fit seven terms |
| N9 | Open items of K and M: tRNA genes with Pol III/TFIIIC, condensin at rDNA and tRNA genes, Sir/Rap1 at telomeres; tight polysomes; ER polysomes (M5); marks in the viewer; nuclear BD rerun | — | — | done: tRNA genes (Pol III, TFIIIC), Rap1 and Sir2/3/4 (`trna_telomeres.py`), ER polysomes (M5), polysome spacing and orientation (tight pairs), ER topology fix, nuclear BD rerun; open: condensin, tRNA-gene clustering at the nucleolus (fails), marks in the viewer |

## P. Round-2 data and the ETH IBC condensate/cytoplasm labs (plan of 2026-11-12, afternoon)

**Direction (the user, 2026-11-12).** Integrate the round-2 data (`data_requests_round2_2026-11-12.zip`) and the ideas from the Michaels, Weis and Neurohr labs (ETH IBC) as seen fit, at a good pace. Working mode for this round: tests that only read the model run on the current build (O) without a rebuild; model changes are batched, each checked by a quick targeted test, then one rebuild, one validation and one docs pass.

| # | Item | Data | Kind | Rebuild | Status |
|---|---|---|---|---|---|
| P1 | Rap1 per chromosome end 15 → 17 | Marcand 1997 Science (a Rap1 site every 18 bp on ~300 bp of TG1–3 repeat) | input | yes | done: 17 per end; 94 of 544 fit (crowded chromatin) |
| P2 | Sir3 2 per silent nucleosome; Sir2–Sir4 1 | Swygert 2014 Nat Commun (2 Sir3 per nucleosome), Swygert 2018 PNAS (1–2 Sir2–4 per nucleosome); both in vitro, saturating | input | yes | done: Sir3 173, Sir2 310, Sir4 119 placed |
| P3 | SPB core in its stoichiometry Spc42:Spc29:Spc110:Cmd1:Cnm67 = 2:1:1:1:1 (now 204:190:54:54:221) | Viswanath 2017 MBoC (citing Muller 2005); no copies per SPB | input | yes | done: Spc110 316 (rings + plaque), Cmd1 316, Spc42:Spc29:Cnm67 = 260:136:136; the 90 nm disc caps Spc42 at ~0.8 per Spc110 |
| P4 | Kinetochore counts | Joglekar 2006 Table 1: Ndc80 complex 8, DAM–DASH 16–20, MIND 6–7 (the model's numbers) | check | no | done (doc): within the measured counts |
| P5 | Collided ribosomes | Zhao 2021 Genome Biol: monosome:disome 16.3:1 after RNase I, 5.8% of ribosomes in disomes (exponential yeast) | test | no | done: 3.9% vs 5.8% (partial; with P11) |
| P6 | tRNA genes at the nucleolus | Thompson 2003 Science (active SUP53 overlaps the nucleolus in 52% of 440 cells, inactive 13%); Haeusler 2008 Genes Dev (single cluster in 83% of WT nuclei, 76% of them nucleolar) | test (a rule only if it can be kept separate from these data) | no | done: fails (33% of tRNA genes within 250 nm vs 52%); no rule added |
| P7 | Budding-yeast 40 nm GEM diffusion as held-out benchmark rows (forced to test: seen after the fit) | Gade 2025 Cell Rep / preprint (≈0.135 and ≈0.21 µm²/s, figure); Kunzi 2026 bioRxiv (≈0.20, figure) | test | no | done: model 0.25 vs 0.135-0.21 (1.2-1.8x too fluid) |
| P8 | Dry mass and protein concentration | Kunzi 2026: dry mass ≈205–225 mg/mL (vacuole excluded), whole-cell protein 71 mg/mL, buoyant density 1.083 g/mL (WT only; PKA inhibition as the negative control) | check | no | done: dry mass consistent; protein 71 mg/mL conflicts with 5 pg (left) |
| P9 | Nucleus / cell volume | Ohsawa 2026 bioRxiv (5.3–6.3% in control strains) | check | no | done: model 8% vs 5.3-6.3% (upper end) |
| P10 | Polysome connectivity in the GPU BD engine; the Gade perturbations (polysome disassembly, fewer idle ribosomes, less mRNA) | Gade 2025 Fig. 2 and 4 | test | no (simulation) | open: needs links between polysome ribosomes in simgpu |
| — | Locus maps (Berger 2008), receptors per pore | figures only; no yeast count | — | — | open |
| P11 | Collided ribosomes in the disome pose (yeast 6I7O) where the drawn linker is <= 5 nt | Ikeuchi 2019 (6I7O), Meydan & Guydosh 2020 (footprints); tested by P5 | input (geometry) | yes | done |
| P12 | tRNA: one type per isoacceptor (42), structures predicted (the user's request) | GtRNAdb sequences; gene copies (Percudani 1997); Boltz-2 | input | yes | done: 3.30 M placed; Mol* opens the models |

## Q. Pixels, accessibility, silent chromatin, polysome links, PDB audit (plan of 2026-11-12, late)

**Direction (the user, 2026-11-12).** Start the next update from the gaps named in the session (no test on raw pixels; no DNase data; silent-chromatin extent a fixed 3 kb; histone marks unused; the Gade perturbations untestable) and use much more data. During the round the user asked for a systematic audit of the PDB against the model (the collided disome 6I7O had been overlooked), by an Opus subagent.

| # | Item | Data | Kind | Rebuild | Status |
|---|---|---|---|---|---|
| Q1 | Nucleus/cytoplasm intensity ratios from raw images | YeastRGB C-SWAT raw grids (9.7 GB; 6,145 strains, Nup49-mCherry marker); Dubreuil 2019 | test (held-out half) | no | done: test half Spearman 0.81, AUROC 0.98 (562 proteins) |
| Q2 | Sir-bound nucleosomes from Sir3 ChIP-seq instead of a fixed 3 kb | Brothers & Rine 2022 (input); H4K16ac (Weiner 2015) and Sir4 ChIP (Goodnight & Rine 2020) held out | input + test | yes | done: 238 Sir-bound nucleosomes at 28 of 32 ends; 78% H4K16-hypoacetylated (old rule 28%); Sir4 ChIP 36x (old rule 4x) |
| Q3 | DNase I accessibility at base resolution | Zhong 2016 per-base cuts (lifted from sacCer2); Zelin 2012 DHS peaks | test | no | done: cores cut 0.66x as often as linkers; 76% of nucleosomes protected |
| Q4 | Links between polysome ribosomes in the GPU Brownian dynamics; the Gade perturbations (was P10) | Gade 2025 Fig. 2 | test | no | done: fails (1.02x, CI 0.98-1.06, against 1.4-1.5x measured; fewer idle ribosomes 1.05x against 1.0x) |
| Q5 | Image intensity against copy numbers | YeastRGB per-ORF table | check | no | done: Spearman 0.74 (0.81 for the brightest 70%) |
| Q6 | PDB audit: assemblies, poses and states in the PDB that the model does not use | RCSB (Opus subagent) | survey | — | done: 7,067 entries screened; Asc1/Stm1 double count fixed; 7 homo-oligomers added; candidates in `pdb_audit_2026-11-12/` |

## R. A consistent G1 cell: protein total, cell volume and organelle fractions chosen together (plan of 2026-11-12, night)

**Direction (the user, 2026-11-12).** "Okay that sounds like a sensible solution. I'm sure it will require a lot of tuning and validation that nothing breaks. You are welcome to implement it." Problem (review comment 9, `validation/anchors.py`): at 5.0 pg of protein, nucleoporins are 1.5x what 117 pores need, core histones 1.3x what one genome's nucleosomes need, ribosomes 1.5-1.6x the ~2 x 10^5 of the literature. Per-cell counts say less protein; the density measurements say the cytoplasm is about right. Both hold in a smaller cell.

| # | Item | Basis | Kind |
|---|---|---|---|
| R1 | Protein total 3.5 pg (was 5.0) | histones per nucleosome (3.9 pg) and nucleoporins per pore (3.2 pg); not the ribosome count, which stays a check | input |
| R2 | tRNA and mRNA counts scaled by the same factor (0.7) | both were per-cell numbers of an average cell; keeps tRNA per ribosome and the ribosome spacing on mRNA (164 nt, Arava 156) | input |
| R3 | Cell volume 34 fL (was 42) | the volume at which the cytosol keeps its macromolecular density (171 mg/mL), so the diffusion calibration and crowding stay valid; dry mass would prefer ~36-38 fL, the cryo-ET ribosome concentration ~26 fL | choice, documented |
| R4 | Vacuole 17.2% (was 9%), nucleus 7% (was 8%) | Chen 2025 (BY4741, YPD, unbudded G1, soft X-ray): vacuole 17.2 +- 6.0%, nucleus 5.9 +- 1.3%; nucleus kept at the size-invariant ratio of Jorgensen 2007 (7%), between the X-ray values (5.3-6.3%) and the FIB-SEM cells (8-9.7%); imaging of loci assumes a nuclear radius near 0.9-1 um | input |
| R5 | Pores: areal density kept (117 on the old envelope = 10.8 per um2) | Winey 1997: 65-182 per nucleus, rising from G1; G1 density not readable here (unchecked) | assumption |
| R6 | Volumes read from the built geometry everywhere (no hard-coded 26.62 um3) | validation scripts, condensates | fix |
| R7 | Regression: every test of validate_all against build Q; chromatin (R(s), Hi-C, telomeres) in the smaller nucleus; packing completeness; viewer | | test |

## S. Stack review: reproducible builds, placement artefacts, vacuoles (plan of 2026-11-13)

Asked for by the user: a review of the technical stack (two subagents: code audit, survey of other methods; reports and
the scripts behind every measured claim in `stack_review_2026-11-13/`), then "You can start work on the codebase".
Branch `stack-cleanup`. Verdict of the review: keep the pipeline (a cellPACK-style greedy construction; no published
method does the whole job better), fix its measured artefacts, and make builds reproducible and comparable.

**Part 1, plumbing (no change to the model at master seed 0):**
- S1. `scripts/build_isolated.sh`: build one commit from a frozen git worktree into its own run folder (own output, own
  copy of the library), so code can be edited while a build runs and replicates can run side by side.
- S2. Master seed `YC_SEED` (`config.seed`), `PYTHONHASHSEED=0`.
- S3. Build stamp (`scripts/build_stamp.py`): `build.json` with commit, seed, all module constants, inputs; stale
  result files are listed after every validation run.
- S4. Stage 0 (`scripts/s0_library.py`): the base structure library built by the pipeline, not by the old notebook.
- S5. `validation/placement_checks.py`: who is dropped, tile seams, membrane shells, layers next to membranes.
- S6. Baseline: rebuild `main` in isolation; is the build reproducible run to run? Then replicates at other seeds:
  every validation number gets a build-to-build spread.
- S7. Validation run completed: stress-condensate negative control, copy ledger, polysome diffusion test on the
  current build (or stated as without power at reachable lags).
- S8. Delete dead code (notebook-only modules, unreachable in-tile mRNA code, switched-off branches) once S4 is in a build.

**Part 2, changes to the model (one build, then replicates; every change inspected visually):**
- S9. Vacuoles: lobes that the scaling pushes together are fused, not carved 30 nm apart (boxy lobes, slits and a
  wrapped lobe in build R: the lobes of cell 08 are at native size inside an outline shrunk by 8%).
- S10. Membrane-face shells: a capacity per shell from its volume; copies beyond it go to the cytosol as the soluble
  pool of the protein (the Golgi face asks for about 200 mg/mL, is 84% placed, 17% of its copies end outside the shell
  and cytosolic proteins are 10 times depleted next to Golgi and endosome membranes); shell test in the attraction
  path and the squeeze.
- S11. Over-full compartments (matrix 75%, nucleolus 80%): thin all species by one reported factor before packing
  instead of dropping mid-sized species (Hsp60 49%, Sod2 37% placed in build R). First test whether finer collision
  shapes fit more (survey: Bentopy reaches 240-300 g/L by plain sequential addition on 0.5 nm voxels).
- S12. Tile seams: interleave size classes across the 8 passes; anchors for partner attraction only where the
  proposal can land in the tile.
- S13. Inventory uses the built compartment volumes (nucleoplasm 1.90 analytic against 1.56 um3 built).
- S14. Library: the 82 transmembrane proteins whose library entry predates the topology side rule, and free H2B
  (no library entry: 33,100 copies never packed), come right with S4.
- S15. Speed: numba kernels for `atp_dimers.place` and `mrnp.place` (18 of 47 min of stage 3).
- Mitochondrial cristae crossing each other in the guided view of build R: check the geometry.

**Later (prototypes, from the survey):** chromatin as a draw from an equilibrated restrained-polymer ensemble (R(s) at
100-200 kb, tRNA genes at the nucleolus); ribosome-scale Brownian dynamics in a periodic box of >= 1,000 ribosomes for
the polysome test; relaxation with growing radii in place of the squeeze; conformer sets for disordered regions.
**Not worth it (survey):** generative models of whole-cell configurations (no training data), wave function collapse
or Wang tiles for molecules, a Martini-scale yeast cell, porting to cellPACK/Mesoscope.

**Copy numbers (user's suggestion, after Part 2):** non-uniform scaling from the average cell to the G1 cell; ribosome
number (cryo-ET 1.4-1.7 times the model). Data being fetched by a subagent (requests to the user only for what cannot
be downloaded).

## T. Sub-projects that use the structural model: simulating processes (plan of 2026-11-13, afternoon)

Asked for by the user: "start some sub-projects/branches related to using the structural model ... run simulations of
some processes using perhaps some coarse grained or simplified versions of our structural model", first endocytic
vesicle release (the user chose this reading of "endosome release"), with visuals and videos. Later the same day:
the axisymmetric membrane approach is fine; simulate "some form of budding"; and which simulations could teach
something new.

**T1. Endocytic vesicle release** (branch `sim-endocytosis`, git worktree `subprojects/endocytosis`, code in
`processes/`; results in `subprojects/endocytosis_results/`; parameter sheet with 523 rows and verified quotes in
`pipeline/data/validation/literature/endocytosis/`).
- Stage 0, membrane: invagination and scission as an axisymmetric elastic membrane under turgor (`processes/membrane/axisym.py`, `processes/endocytosis/scission.py`). The same solver does outward buds.
- Stage 1, release: the freed vesicle with its coat and actin shell among explicit ribosomes and large complexes (GPU Brownian dynamics, `processes/common/mesobd.py`, `release.py`); small proteins are in the diffusion coefficients (`calibrate.py`).
- Stage 2, transit: the bare vesicle as a walker in the whole cell's geometry to its acceptor compartment (`transit.py`).
- Membranes in stages 1 and 2 are rigid walls; scission and fusion are events, not simulated membrane dynamics.

**T3. Membrane-bending proteins, and the three methods feeding each other** (user's request of 2026-11-13, evening):
- Each bending mechanism becomes a term of the axisymmetric membrane: coat adaptors with amphipathic helices (Ent1/Ent2 ENTH helix 0) as a spontaneous curvature from the inserted area, BAR scaffolds (Rvs161/167 N-BAR, Bzz1 and Syp1 F-BAR) as an anisotropic curvature that prefers a tube, crowded disordered domains as a lateral pressure, the actin network as the pulling force. Copy numbers and timing per site come from the measurements (Picco 2015 and others), footprints and helix areas from the structures in the model's library.
- Knock-outs in silico (no helix, no Rvs, no Bzz1, no Syp1, no crowding) against the measured phenotypes of the deletion strains: the test that each term is right.
- Feedback: the solved membrane shapes become the invagination profiles of the structural model (fixes the 40 nm tips and the shallow tubes of build S), and the release and transit runs start from that geometry.
- Parameters: `pipeline/data/validation/literature/membrane_bending/` (being collected).

**T2. Candidate questions where the model's copy numbers and geometry decide the answer** (given to the user, none started):
tRNA supply to ribosomes per codon against ribosome-profiling dwell times; whether SRP can find every secretory nascent
chain in time; whether crowding alone explains the stiffening of the cytoplasm on volume loss (reuses the T1 engine);
the largest mobile particle per compartment (pore-size analysis of the packed structure); Cdc42 polarity on the real
cell surface (the cell's own budding); import and pore-traffic budgets. Each needs a novelty check in the literature
before it is started.

## Progress

**2026-11-03: inventory v3 (items 0.1, A1, A2 soluble).** Pipeline scripts are in `yeast_cell_dev_bundle/scripts/`; see the README there.

| Metric | Before (v2) | After (v3) |
|---|---|---|
| Protein mass in multi-chain assemblies | 34.7% | 52.7% |
| Multi-protein complex types | 11 curated + 7 from Complex Portal | 11 curated + 102 from Complex Portal (RCSB assemblies) |
| Homo-oligomers | 23 PDB | 130 PDB + 313 AF-predicted dimers |
| Mass by source: PDB / AF complex / AF monomer | 51.2 / 0 / 48.8% | 54.8 / 8.1 / 37.1% |
| Placed instances | 25.2 M | 22.7 M (fewer, larger objects) |
| Mitochondrial matrix / IMS filled | 55% / 42% | 63% / 47% |
| Clash rate | 4.6e-6 | 3.2e-6 |

Other changes:
- The complex rule is now: subunits present ≥ 80% of mass and ≥ 60% of residues modeled.
- Complex RNA is kept only where it belongs to the complex: snoRNA, snRNA, telomerase RNA, and mitochondrial rRNA for the mitoribosome subunits.
- Complex copy numbers are capped at 2× the limiting subunit. This fixed a 140× inflation of MNL1, plus smaller ones of NTH1 and CRM1. Inventory mass is now within 0.1% of the proteome.
- Bug fixed along the way: the mRNA thread count is now reported correctly in the validation JSON.

New open issues:
- **Nucleolus placement** dropped to 94% (was about 99.8%) because of large pre-ribosomal complexes.
- **Mitochondrial ribosome:** only 218 of 1,168 large subunits fit into the crowded matrix (E4).
- **TSA1 and other UniProt "dimers"** whose PDB assembly is a larger oligomer stay monomers.

**2026-11-03 (second round): membranes (items B1, A2 membrane, A1 membrane, B2 for mitochondria).**

| Metric | v3 first round | Now |
|---|---|---|
| Protein mass in assemblies | 52.7% | 54.0% |
| Complex Portal complexes | 102 | 125, including 20 membrane complexes |
| Membrane homo-oligomers | 0 | 40 (37 AF dimers, 3 PDB) |
| Molecules placed | 96.8% | 98.0% |
| TM copies dropped (annotated soluble) | 63k | 43k |
| Mitochondrial matrix / IMS filled | 63% / 47% | 68% / 68% |
| Clash rate / containment | 3.2e-6 / 2.3e-4 | 0 in sample / 5.2e-5 |
| TM residues of oriented assemblies within 2 nm of the membrane midplane | — | 97.6% (median 0.73 nm) |

**Automatic membrane orientation** (`pipeline/membrane_orient.py`):
- **Method:** UniProt TRANSMEM residues are mapped onto PDB chains through RCSB's UniProt alignment. A search then finds the membrane normal that maximizes TM residues inside a ±1.6 nm slab and minimizes other residues there. The side comes from UniProt topology domains; without them, the side carrying more soluble mass is taken as cytoplasmic.
- **Validation:** it reproduces the hand-coded orientations of Pma1 (z correlation 1.000) and V-ATPase (0.9996).

**Membrane complexes placed include:** ATP synthase (8F39, 11.7k copies), respiratory complexes II, III and IV, the Sec translocon, OST, EMC, PMT, VTC, β-1,3-glucan synthase, SPOTS, the Drs2 and Dnf1 flippases, ceramide synthase, IPC synthase, SAM and GPI transamidase.

**Mitochondrial sub-compartments** are now parsed from UniProt: inner membrane, outer membrane, intermembrane space and matrix. This replaces the fixed 75/25 and 90/10 splits. Rules:
- A secondary mitochondrial location behind a cytosolic primary one (SOD1, ADK1) gets 0.1× weight.
- A TM protein listed in both soluble and membrane compartments goes entirely to the membranes.

**Copy numbers:**
- mtDNA-encoded subunits (COB, COX1–3, ATP6, ATP8, OLI1, VAR1) are ignored when estimating complex copy numbers, because mass spectrometry badly under-detects them (COB 40, ATP8 0).
- Complexes with 4 or more subunit types are capped at 2× their 25th-percentile subunit.

**Still open:**
- **β-barrels without TRANSMEM annotations:** TOM40, and with it the TOM complex, can't be oriented yet.
- **Remaining dropped TM copies:** 43k TM copies are annotated only as soluble.
- **ER vs nuclear envelope, Golgi vs endosome:** not yet split (rest of B2).
- **A3 (heterodimers):** done in the third round, see below.

**2026-11-03 (third round): heterodimers (item A3)** in `pipeline/heterodimers.py`.

234 of the 1,250 AFDB heterodimer models are now packing ingredients. Protein mass in assemblies went from 54.0% to 55.3%; AF complexes are 9.6% of mass.

Filters (counts are pairs removed):
- **442:** a partner has no free copies left after complexes were assigned.
- **259:** no IntAct physical-interaction evidence. `data/intact_pairs.tsv` was built from IntActYeast.txt by `pipeline/intact_pairs.py`.
- **133:** ipTM < 0.7.
- **70:** no shared compartment. A membrane protein may pair with a soluble partner on its cytosolic or matrix side.
- **70:** same-family isoform pairs with fewer than 3 publications and no direct-interaction record, e.g. SSA1–SSA2, TDH1–TDH3, HXK1–HXK2, PDC1–PDC5. Complex Portal co-membership does not count for these, because Complex Portal lists interchangeable isoforms in one entry.
- **30:** both partners are subunits of a complex that is already modeled, so the pair would only be proteome surplus.
- **12:** a partner forms a larger homo-oligomer (e.g. SNZ/SNO, TKL), unless Complex Portal lists the pair together.

Copy allocation:
- Each protein's free copies are split evenly among its accepted partners, so no protein is over-used.
- Pairs with Complex Portal support (non-isoform) or a direct-interaction record bind all available copies. Others bind 50%, an assumption: without affinity data, stable and transient interactions can't be told apart.

Largest new pairs: RNR2–RNR4, actin–cofilin, SAM1–SAM2, PDA1–PDB1, SOD1–CCS1, actin–profilin, BMH1–BMH2, eIF4G–eIF4A, SEC13–SEC31, COPI subcomplexes, TPS1–TPS2, RVS161–RVS167.

Validation after this round: placed 98.0%, clash rate 1.4e-6, containment 6.4e-5, 97.7% of oriented TM residues within 2 nm of the bilayer midplane, inventory mass within 0.4% of the proteome. On the Mac's GPU: 10 ms/frame for the cytoplasm view, 31 ms/frame for the whole cell, no console errors.

**2026-11-03 (fourth round): relaxed dimer filters, membrane completeness, polysomes, actin cables.**

| Metric | Third round | Now |
|---|---|---|
| Heterodimers / AF homodimers | 234 / 349 | 512 / 453 |
| Protein mass in assemblies | 55.3% | 61.1% |
| Mass by source: PDB / AF complex / AF monomer | 55.2 / 9.6 / 35.2% | 55.0 / 14.6 / 30.4% |
| TM copies dropped | 43.8k | 9.9k |
| Membrane proteins on the nuclear envelope | 19k | 92k |
| Ribosomes in polysomes | 0 | 92.7k on 23.0k mRNAs (4.0 per mRNA) |
| Actin in filaments | 0 | 33.3k subunits, 74 filaments in 10 cables (92 µm) |
| Clash rate / containment | 1.4e-6 / 6.4e-5 | 0 in sample / 1.3e-5 |

**Dimer filters relaxed** at the user's request:
- **Heterodimers:** ipTM ≥ 0.5. IntAct evidence is no longer required; it now only decides full binding versus 50% binding. The isoform rule is removed, and paralog pairs are flagged in `heterodimers.csv` instead.
- **Homodimers:** ipTM ≥ 0.5.
- **Kept:** the shared-compartment rule and the two consistency rules (not both partners in the same modeled complex; no partner that is modeled as a larger homo-oligomer).

**Membranes:**
- **β-barrels** (TOM40, SAM50, MDM10, MDM34) have no TRANSMEM annotation. They are now flagged as integral, and their β-strands (from Cα geometry) define the membrane segment. The TOM core complex (9J99) is now placed.
- **TM proteins annotated only to soluble compartments** take their membranes from GO, with "cell periphery" read as plasma membrane.
- **ER membrane proteins** are shared with the nuclear envelope by area (15%), because the outer nuclear membrane is ER.

**Polysomes** (`pipeline/polysome.py`):
- **Channel geometry:** the mRNA channel comes from 7ZW0 superposed on the 4V88 ribosome (median deviation 0.16 nm).
- **Assembly:** each ribosome's mRNA entry sits 3–12 nm beyond the previous one's exit, within a 40° cone and at random spin; neighbouring ribosomes end up 25–31 nm apart.
- **Copy numbers:** 70% of cytosolic ribosomes go into polysomes, with sizes 1 + Poisson(4). These are assumptions. Free mRNA threads are reduced by one per polysome.

**Actin cables** (`pipeline/cytoskeleton.py`):
- **Geometry:** the F-actin subunit and helix come from 9GO5 (167.2° and 2.76 nm per subunit; reproduced exactly).
- **Layout:** cables are 8-filament bundles 90 nm beneath the plasma membrane, along random great circles.
- **Assumption:** 30% of ACT1 is in cables.
- **Viewer:** a new guided view, "Actin cable".

**Still open:**
- **Actin patches** (Arp2/3 branched networks), tropomyosin and fimbrin on cables, microtubules and septins.
- **ER-bound polysomes** on Sec61.
- **mRNA outside polysomes:** UTRs and untranslated stretches of polysome mRNAs are not drawn.
- **Golgi vs endosome:** not split (the geometry has one compartment).
- **Frame time:** close-up views got slower on the GPU (cytoplasm 10 → 12–25 ms, cortex 66 ms) because there are more residues per object.

**2026-11-05 (fifth round): microtubules and the spindle pole body, ER-bound ribosomes.**

| Metric | Fourth round | Now |
|---|---|---|
| Microtubules | 0 | 21 (16 kinetochore, 2 other nuclear, 3 cytoplasmic); 11.5k tubulin dimers |
| γ-tubulin rings | 0 | 21, one on each minus end |
| Ribosomes docked on ER translocons | 0 (3% floated near the ER) | 4,285 placed (3,365 on Ssh1, 920 on Sec61) |
| Chromosome walks | start near the SPB | start at their kinetochore microtubule tips |
| Placed / clash rate / containment | 98.0% / 0 in sample / 1.3e-5 | 98.0% / 6.8e-7 / 6.4e-5 |
| TM residues within 2 nm of midplane | 97.5% | 97.4% |

**Microtubules** (`pipeline/microtubules.py`):
- **Lattice geometry:** taken from 8QV3, the native yeast γ-TuRC capping a spindle microtubule (cryo-ET, 8.2 Å). Two neighbouring protofilaments give the lateral step (measured 28.1° and 1.02 nm). The lattice is built as the ideal 13_3 lattice, with its seam between protofilaments 12 and 0.
- **Lattice dimer:** 5W3F, the yeast microtubule structure at 3.7 Å, superposed on 8QV3's first dimer (0.16 nm deviation).
  - 8QV3's own dimer is about 4% longer along the axis (α–β spacing 4.35 vs 4.18 nm), so its dimers would interpenetrate at a normal repeat.
  - The axis follows the dimer's α→β direction, which differs from the two-protofilament fit by 1.2°.
  - The repeat is 8.3 nm, the shortest at which stacked dimers don't interpenetrate; only the α N-terminal Met comes within 0.3 nm.
- **γ-TuRC:** 7 γ-TuSCs from 8QV3, repeated two protofilaments at a time.
  - The 7th γ-TuSC's overlapping half (γ-tubulin 14 and its SPC98) is left out, because on the ideal helix it would sit inside protofilament 0.
  - Nuclear rings keep the SPC110 fragments.
- **Layout:**
  - **Placement:** the SPB sits opposite the nucleolus. Minus ends are spread over 90 nm plaques and fan out up to 30° in the nucleus and 70° in the cytoplasm.
  - **Growth:** microtubules grow layer by layer and bend away from membranes, with a curvature radius of at least 300 nm.
  - **Why bending is needed:** the SPB faces the vacuole, so straight cytoplasmic microtubules stopped at about 250 nm.
  - **Lengths (G1 assumptions):** kinetochore microtubules 250 ± 60 nm, other nuclear ones 400 nm, cytoplasmic ones 1.2 µm. They are capped at 70% of the αβ-tubulin pool, which scales them by 0.87.
- **Chromatin:** microtubules are placed before chromatin, and each chromosome walk starts within 20 nm of its kinetochore microtubule tip.
- **Ring subunits:** the γ-TuSC complexes (62) are used first, then free TUB4/SPC97 and SPC110 (from its CMD1 heterodimer). PaxDb lists fewer SPC98 copies than 21 rings need (64 short).

**ER-bound ribosomes** (`pipeline/er_ribosomes.py`):
- **Docking geometry:** 2WW9 (yeast 80S–Ssh1, 8.6 Å) places the translocon on the 4V88 ribosome through 601 tunnel-exit residues (0.19 nm deviation).
- **Two types:**
  - **Ssh1:** SSH1, SSS1 and SBH2 from 2WW9.
  - **Sec61:** SEC61, SSS1 and SBH1 from 7KAH, placed through the shared SSS1. Its channel lies on 2WW9's Ssh1 channel within 0.29 nm (median).
- **Orientation:** both are oriented from the translocon TM segments (score 0.93), with the ribosome forced to the cytosolic side. 2WW9's Ssh1 model doesn't follow Ssh1's residue numbering, so its topology-based sign came out inverted.
- **Copy numbers:** 80% of the free channel pools is taken as ribosome-bound (assumption). The shared SSS1 pool limits both channels together and is split in proportion to the two pools.
- **Result:** 2.6% of ribosomes on the ER and outer nuclear membrane. These replace the ribosome ingredient's former 3% "ER" share.
- **SBH1:** all of it is already in the SEC complex, so the Sec61-type β chains are drawn but not counted against the pool.

**Viewer:**
- **Guided views:** two new ones, "Spindle pole body" and "Cytoplasmic microtubule". The legend has entries for microtubule and γ-tubulin ring.
- **GPU frame times:** spindle pole body 12 ms, cytoplasmic microtubule 5 ms, cortical ER 29 ms, cytoplasm 22 ms. The only console messages were occasional tile-fetch warnings from the local Python server; those tiles are re-requested automatically.

**Pipeline:** `config.py` now falls back to `~/Downloads/UP000002311_559292_YEAST_v6.tar` and `~/Projects/bioemu_docker/afdb_yeast_dimers` when the `YC_*` variables are unset.

**Still open:**
- **Cytoskeleton:** actin patches (Arp2/3 is present only as AF heterodimers), tropomyosin and fimbrin on cables, and septins.
- **SPB plaques:** the Spc42 layer and the Spc72 outer plaque are not built.
- **ER-bound polysomes:** the docked ribosomes are single, with no mRNA.
- **Abundances:** the ribosome-associated factors NAC, RAC/Ssb and eEF1A are still free in the cytosol.
- **Remaining items** from the fourth round: UTRs, the Golgi/endosome split, and frame time.

**2026-11-05 (sixth round): Arp2/3 complex and actin patches.** The viewer is also published on GitHub Pages (repository `yeast_cell_viewer`); `scripts/export_pages.sh` rebuilds the upload folder.

| Metric | Fifth round | Now |
|---|---|---|
| Arp2/3 | split into 5 pairwise AF heterodimers | one 7-subunit complex: 6,972 copies (2,169 at patch branches, the rest free) |
| Actin patches | 0 | 22 patches, 33.0k actin, 2,169 branches, 1,251 capping proteins |
| ACT1 allocation (of 111k) | cables 33k; the rest monomeric or in heterodimers | cables 33k, patches 33k; the remaining 44k is split between free actin and its cofilin, profilin and Srv2 heterodimers |
| Placed / clash rate / containment | 98.0% / 6.8e-7 / 6.4e-5 | 98.0% / 6.7e-7 / 6.6e-5 |

**Arp2/3 complex and capping protein** (`pipeline/actin_patches.py`):
- **Arp2/3:** the AlphaFold models of ARP2, ARP3, ARC40, ARC35, ARC18, ARC19 and ARC15 (pLDDT ≥ 50) are superposed on their fission-yeast homologues in the branch junction 8UXW (2.7 Å). Sequence identity is 32–74%, and each subunit fits to 0.05–0.22 nm.
- **Capping protein:** CAP1/CAP2 AlphaFold models superposed on CapZ at the capped barbed end 8F8Q (38% and 52% identity).
- **Copies:** Arp2/3 copies follow the usual complex rule (median subunit pool, capped at 2× the limiting one), so ARC40 and ARC19 run about 1.5× over their pools. The complex is built in s1 before the heterodimer step, so pairs of its own subunits are no longer formed.
- **Branch geometry:** the library F-actin subunit (9GO5) fits every actin in both templates to 0.04 nm. 8UXW gives the Arp2/3 pose and the first daughter actin relative to the mother subunit they sit on (branch angle 73°). 8F8Q gives the cap pose relative to the terminal actin. +z of the actin frame points to the barbed end.

**Patches:**
- **Sites:** cytosolic points directly under the plasma membrane (no cortical ER in between).
- **Growth:**
  - 3 primer filaments, barbed ends toward the membrane.
  - Dendritic branching from random mother subunits; daughters may not point away from the membrane.
  - Filaments are capped after 20 ± 6 subunits.
  - Each patch is confined to 130 nm laterally and 220 nm deep.
- **Collisions:**
  - Sphere trees check against everything else.
  - Inside a patch, the Arp2/3 complex, the first two daughter actins and the caps are also checked residue by residue against everything except their own mother filament.
  - Result: no inter-instance contacts beyond the native junction contacts, about 3 Cα pairs within 0.3 nm per branch, as in the templates.
- **Capping:** about half the barbed ends are capped. Filaments that stop against the membrane or another filament have no room for a cap.
- **Assumptions:** 30% of ACT1 in patches (and 30% in cables), 1,500 actin per patch, cap length 20 ± 6.

**Not done: RAC/Ssb on ribosomes.** The yeast RNC–RAC entries hold only a Zuo1 fragment with 94 nt of rRNA (7X34) or RAC alone at 6 Å (7X3K). No entry has the ribosome model needed to dock it reliably, and there is no yeast NAC–ribosome structure.

**Viewer: Show panel.**
- **Own checkboxes:** Actin (filaments, Arp2/3, capping protein, tropomyosin) and Microtubules (tubulin, γ-TuRC). Ribosomes now also covers the ER-docked ribosomes.
- **Compartment filters:** groups with their own checkbox (ribosomes, tRNA, nucleosomes, pores, actin, microtubules) are no longer hidden by them. The shader checks a "class" type-flag bit.
- **Shortcuts:** an "only" link on each row, plus Hide all and Show all buttons.

**Viewer: fine cut-plane steps.** Cut plane panel has a "Fine step" row: −/+ buttons (repeat while held) and a step size of 0.1–50 nm, default 0.5 nm, remembered per browser. Alt+wheel (one step per notch) and Alt+, / Alt+. use the same step. The camera moves with the plane, as with the Position slider.

**Real cell shapes** (`pipeline/cellshape.py`; MODEL.md 3.3). The model cell is a composite of two segmented cryo-FIB/SEM cells: the outline, vacuoles and lipid droplets of cell 08 (the only uncut one) and the nucleus of cell 05.
- **Cleaning:** voxelised at 20 nm, closed (fills the curtaining-stripe groove) and Gaussian-smoothed.
- **Vacuoles:** the segmentation flattens them, so their z extent is reconstructed: a per-lobe z stretch up to 3×, stronger z smoothing, and a 1:1 blend with the moment-matched ellipsoid. One streak artefact is dropped; the lobes are rescaled to 9% and carved 30 nm apart from each other and the nucleus.
- **Distance fields:** exact distance (closest triangle) on the grid, stored as sparse bricks, appended to `G` after `cell3d.G_N` (`cell3d.field_sdf`). Builder code uses dense grids (`cellshape.Grid`); `ray_depth` gives cortex points.
- **Viewer:** the same fields (`geometry_fields.bin`, 7.5 MB gz) in the fill shader, the section map and the worker's wall generator. The tile-state grid is widened to 26³ (±2.6 µm).
- **Result:** placement and validation are unchanged (97.8% placed, clash 1.9e-6).
- **Gotchas:**
  - the bricked field is clamped (cell 260 nm inside), so deeper queries need the dense grid;
  - marching cubes emits zero-area triangles that broke the distance kernel;
  - the old ±2.2 µm tile key wrapped around for tiles beyond it.

**Geometry matched to segmented cryo-FIB/SEM cells** (EMPIAR-11416; `scripts/organelle_stats.py` → `pipeline/organelle_stats.json`; MODEL.md 3.1, 3.3):
- **Nucleus:** an ellipsoid (1.25 : 1.10 : 1) of 8% of the cell, closer to the centre (0.3 of the radius). The pores sit along the local envelope normal. One first-order ellipsoid distance is used in the pipeline, the WGSL fill shader and the section map.
- **Vacuole:** 23 lobes resampled from the 107 segmented vacuoles, 9% of the cell, grown as touching clusters. Lobe spread and surface-to-volume match the cells (7.5 vs 6.8–9.9 µm⁻¹).
- **Lipid droplets:** 16 droplets resampled from the 83 segmented ones, 0.5%.
- **Effect:** cytosol 93 mg/mL (was 89), nucleoplasm 54 (was 62); placement completeness unchanged.
- **Gotcha:** numba `cache=True` keeps stale inlined copies of `cell3d` functions in modules that were not edited (`tiles3d`, `mrna`, ...). After changing `cell3d`, delete the `*.nbi`/`*.nbc` files in `pipeline/__pycache__`.
- **Superseded by the real shapes below** (still available with `s2_geometry.py --sphere`).

**Organisation layer of the packing** (`pipeline/organisation.py`; MODEL.md section 8.4). Replaces uniform random placement with a weak, assumed organisation:
- **Large scale:** 31 STRING modules (Louvain, score ≥ 700) with smooth random fields (500 nm) that modulate how each species is split over tiles, with a crowding correction so only the composition varies.
- **Meso/small scale:** partners from the AlphaFold 3 pool scores (contact, PAE ≤ 8 Å, paralog-filtered) and STRING (loose, 2–15 nm), proposed next to placed partner copies during packing with a mass-action-like saturating probability.
- **Result:** about 20% of molecules placed next to a partner; contact partners in contact 8 times as often as random. Completeness and clash rate unchanged.
- **Viewer:** "Interaction module (STRING)" colouring with a clickable legend that highlights one module.
- **Knobs:** `BETA` (domain contrast), `FIELD_L` (domain size), `A_CONTACT`/`A_LOOSE`/`P_MAX` (partner pull).
- **Next:** pre-docked pairs and triplets from the pool structures as sampling units, placed by the same rules; compartment assignment of the 625 default-cytosol proteins from their partners.

**Ribosome assembly: pre-ribosomes and the rDNA** (`pipeline/preribosomes.py`, `pipeline/rdna.py`; MODEL.md 5.2, 6.6; items C2 rDNA and C3 pre-ribosomes):
- **Pre-ribosomes:** five stages with all their RNA: 90S 5WLC (2,075) and nucleolar pre-60S 7R7A (4,748) in the nucleolus, Rix1–Rea1 pre-60S 6YLH (1,450) in the nucleoplasm, Nmd3–Lsg1 pre-60S 5T62 (4,709) and pre-40S 6FAI (2,040) in the cytoplasm. Counts from the stage-specific assembly factors, which the particles consume; the 80S count drops by the particles' r-protein content (6,375 ribosome equivalents). 15,000 particles, 8.9% of ribosomes: about a minute per stage at 2,000 ribosomes per minute.
- **rDNA:** 75 active genes (half of ~150 repeats) with 50 Pol I each on a smooth DNA path (Pol I oriented by the 5M3F elongation complex), nascent pre-rRNA per polymerase (70% cleaved at A2), 90S on uncleaved 5' ends; 5,227 nucleosomes on inactive repeats and spacers. 385 of 678 tethered 90S fit, 58% of the nascent RNA: the nucleolus is crowded, and threads that run out of room end early.
- **Effect:** nucleolus 85 → 111 mg/mL (RNA/DNA 6 → 24), nucleoplasm 56 → 60; nucleolus placement 98.6% (was 99.4%); clash rate 2.2e-6.
- **Viewer:** coral pre-ribosomes (in the Ribosomes switch), Show switches for pre-rRNA and RNA polymerases, and an "rRNA genes" guided view that restricts the Show list (new optional `show` field of guided views).
- **Pipeline:** `run_all.sh` now stops when a stage fails (it pipes through tee) and re-encodes the motion modes after s5.
- **Next:** compact nascent RNPs (UTP-A/B on the 5'ETS of short transcripts) instead of free threads; 5S rDNA with Pol III in the spacers; Rabl-connected rDNA on chromosome XII.

**Viewer: motion (Display → Motion: Off / On / Exaggerated).** Illustrative only; see MODEL.md section 7.
- **Conformational motion:** 3,692 AF monomer species move along their AlphaFold 3 pool ensembles, using 4 PCA modes and a conformer path per species. The data is `library/modes.bin`, 29 MB (20 MB gzipped). Built by stage 6 (`scripts/s6_motion.py`).
- **Everything else:** a shader-only rigid wobble, low-pLDDT flexing and free-polymer wiggle. Filaments, pores and nucleosomes are static.
- **Amplitude cap:** none in the shipped data. With a 1.5 nm cap, chain stretching drops (proteins with >5% distorted Cα–Cα pairs at "On": 9% vs 45%) but most motion becomes very small. The cap is `MOTION_CAP` at re-encode (6 s).
- **Possible next step:** superpose the pool models on confident residues only, so floppy termini do not add rigid motion to the core (needs the ~100 min refit).

**Tropomyosin on cables, and the cable cross-section** (`pipeline/cytoskeleton.py`):
- **Spacing fix:** the 8 filaments of a cable were 1 + 7 on a 9 nm ring, which put neighbouring outer filaments only 7.8 nm apart. Packing then left gaps in them (74 short filaments). Cables are now the 8 hexagonal-lattice points nearest the centre at 12 nm spacing (assumed fimbrin bundle spacing): 63 nearly continuous filaments (about 530 of 543 subunits each) in 8 cables.
- **Tropomyosin:** `TPM1_coil` is a coiled-coil dimer, 2 × 200 residues (TPM1 is 199 aa).
  - **Path source:** the path comes from 5JLF (mammalian F-actin–tropomyosin; its coiled coil spans 135 residues, 3.7 nm from the filament axis), transferred onto the 9GO5 actin frame.
  - **Extension:** the filament's own symmetry extends it: the next subunit on the same strand corresponds to 40.0 residues, so 5 actins match one TPM1 dimer.
  - **Smoothing:** helix phase doesn't repeat exactly every 40 residues, so the symmetry copies are blended into one smooth coiled-coil axis. The two chains are wound around it as an ideal left-handed coiled coil (superhelix radius 0.49 nm, pitch 14 nm).
  - **Checks:** Cα spacing is 0.38 nm, the chains are 0.5–1.4 nm apart, the closest approach to actin is 0.44 nm, and consecutive dimers stay 0.52 nm apart.
- **Placement:** one dimer per 10 contiguous filament subunits on each strand. Each is checked residue by residue against the other filaments of its cable, plus obstacles and containment.
  - **Result:** 6,610 dimers (100% of sites), using 13,220 of 35,800 TPM1 chains.
  - **Caveat:** the path is mammalian; no yeast actin–tropomyosin structure exists. TPM2 stays free.
- **Validation:** placed 98.0%, clash rate 1.3e-6 (4 beads sampled), containment 7.0e-5. The cable view renders at 49 ms per frame.

**Full-length mRNA** (`pipeline/mrna.py`):

| Metric | Before | Now |
|---|---|---|
| mRNA nucleotides drawn | 6.8M | 39.7M (target 45M) |
| Polysome mRNA | about 200 nt (ribosome channels and linkers only) | channels and linkers, plus 3′ and 5′ tails, 26.5M nt in all |
| Free mRNAs | 6,864 in-tile threads, about 290 nt each (stuck in crowded tiles) | 6,864 full-length threads, 9.2M nt |

- **Why the old threads were short:** they were walked inside each 200 nm tile after the proteins were packed (target mean 1,500 nt), so they got stuck at about 290 nt.
- **Threading:** all mRNAs are now threaded in the pre-pass, before packing, in one whole-cell grid.
  - **Obstacles:** everything already placed in the cytosol (cables, tropomyosin, patches, microtubules, pores, polysome ribosomes as 100-bead sphere trees), plus all other mRNA.
  - **Lengths:** lognormal, mean 1,500 nt, σ 0.5 (log scale), minimum 300 nt.
  - **Walk:** one bead per nucleotide, 0.6 nm apart, persistence weight 1.0.
  - **Polysome tails:** the 3′ tail leaves the entry of the first ribosome and the 5′ tail leaves the exit of the last one. The 5′ side gets 10–40% of the untranslated remainder. Tails step straight out of the channel until they clear the ribosome envelope.
- **Result:** 84% of threads reach full length, 86% of the target nucleotides overall. Threading takes 70 s.
- **Packing:** the packer avoids the RNA beads, which are bucketed by tile for speed. Placement is unchanged at 98.0%. In-tile mRNA threads are switched off.
- **Encoding:** threads are cut into 64-nt pieces, assigned by midpoint tile. Tiles grew from 412 to 502 MB compressed; the data folder is 536 MB.
- **Validation:** pre-pass mRNA is now part of the clash sample. Clash rate 1.3e-6, containment 4.8e-5.
- **Frame time:** the cytoplasm view takes 10 ms per frame, 15 ms with only mRNA and ribosomes shown in a thicker slab.
- **Caveats:**
  - Ribosomes stay compactly stacked on each mRNA (25–31 nm apart), as in cryo-ET polysomes, so the untranslated length is all in the two tails rather than between ribosomes.
  - The total mRNA count (30,000) is the inventory's assumption.

**ATP synthase dimer rows on cristae rims** (`pipeline/atp_dimers.py`):
- **Structure:** two copies of the yeast monomer 8F39 (the CPX-3281 library entry), superposed on the two monomers of 4B2Q. That's the yeast dimer modelled into a cryo-ET average (37 Å; the superposition deviates 2.3 nm, which is fine for placement).
  - The dimer angle comes out at 83.5°.
  - The two membrane centres lie on a circle of radius 8.8 nm.
  - The monomers stay at least 3.9 nm apart; 8F39 lacks the dimer-specific subunits e and k.
- **Geometry change:** crista half-thickness (the rim tube radius) went from 7 to 9 nm, an 18 nm lumen, to fit that curvature (`cell3d.CRISTA_HALF`).
  - Matrix volume: 0.65 → 0.615 µm³.
  - IMS volume: 0.208 → 0.243 µm³.
  - Inner-membrane area: 28.0 → 28.4 µm².
- **Placement:** dimers centred on the rim tube every 13 nm (assumed row spacing), ridge along the rim. Checks:
  - TM residues within 2.5 nm of the inner-membrane midsurface.
  - Soluble residues in the matrix or IMS, with at most 2% within 0.6 nm of a membrane.
  - No sphere-tree clashes.
- **Result:** 5,881 dimers, using all 11,762 monomers. TM residues sit 0.70 nm (median) from the midsurface, 98.8% within 2 nm.
- **Side effects:**
  - Matrix placement drops 69% → 65%; IMS rises 68% → 86%.
  - Overall placement is 97.9% (99.8% outside mitochondria).
  - Containment 1.3e-4.
  - The mitochondrion view renders at 81 ms per frame.
- **Bug fixed along the way:** `inventory.build()` returned a GO list instead of the geometry, because of a shadowed variable from the GO-fallback change. Only s2 uses that value.

**Respiratory supercomplexes** (`pipeline/supercomplexes.py`):
- **Structure:** III₂IV₂ from 6HU9 (yeast, 3.35 Å), oriented like the other membrane assemblies (score 0.79, 83% of TM residues in the slab).
- **Copies:** min(III₂, IV/2) = 1,561. All complex IV is used; 1,422 III₂ stay free.
- **Validation:** inner-membrane placement 99.3%, clash rate 1.9e-6, containment 1.1e-4. 97.4% of TM residues lie within 2 nm of the midsurface.

**Rendering speed: tried and reverted.**
- **Measurement:** frames are fragment-bound. A quarter of the pixels cuts the mitochondrion view from 80 to 15 ms.
- **Breakdown (mitochondrion view):** background fill 6 ms, lipids about 19 ms (626k lipids × 3 beads, not culled), soluble proteins about 27 ms, membrane proteins about 16 ms, polymers about 11 ms.
- **What was tried:** two-phase Hi-Z occlusion culling. Instances visible last frame are drawn first, a depth pyramid is built, and the rest are tested against it.
- **Why it failed:** it never culled anything. In a molecular scene almost every object's screen footprint includes a background pixel between molecules, so the farthest depth in the footprint stays at infinity.
- **Outcome:** reverted. The code is kept in the session scratchpad.
- **What would actually help:** something that restores early-Z, such as real hemisphere geometry instead of depth-writing impostors, or cheaper lipids.
- **Meanwhile:** adaptive quality already coarsens detail while the user interacts.

**FG-repeat mesh of the nuclear pores** (`pipeline/fg_mesh.py`):
- **What's added:** the disordered domains of the FG-Nups modelled in 7N9F, threaded as residue chains (0.38 nm steps, bead radius 0.3 nm, low persistence).
  - NSP1 residues 1–636 (48 copies per pore), NUP57 1–286 (32), NUP49 1–269 (32), NUP159 388–1102 (16).
  - Anchored at the first modelled residue, pointing at the pore axis.
  - Confined to the pore channel, cytoplasm and nucleoplasm, within 60 nm of the anchor, avoiding the scaffold, chromatin and microtubules.
- **Result:** 14,976 domains on 117 pores; 6.12M of 6.90M residues placed (89%).
- **Packing:** polysomes, mRNA and the packer avoid the FG beads.
- **Viewer:** a "FG repeats" Show checkbox and a picking label.
- **Tile format:** a third polyline channel with a kind per chain (header words 13–14, appended after the membrane data). Free-polymer pick IDs now carry 3 kind bits (15-bit chain index). Kind 4 (pre-rRNA) is reserved.
- **Pipeline:** `mrna.thread_all` is now generic (step, bead radius, persistence, allowed compartments, reach, straight start).
- **Memory fix:** the mRNA threading grid had been sized with 30 cell entries per obstacle sphere (about 7 GB), which pushed the packing workers into swap. It now uses 8; s3 peaks at 9.2 GB.
- **Validation:** clash rate 1.4e-6, containment 4.3e-5. The tiles are 520 MB compressed (data folder 553 MB).
- **Not modelled:** the FG-Nups absent from 7N9F (NUP100, NUP116, NUP145N, NUP42; basket NUP1, NUP60, NUP2), so the nuclear side of the pore has few FG chains.

**Still open:**
- **Cytoskeleton:** fimbrin (SAC6) cross-links between cable filaments, and septins.
- **Nucleolus:** the rDNA array with Pol I and nascent pre-rRNA ("Christmas trees"). The pre-rRNA polyline kind is already in the tile format.
- **Ribosome-associated factors:** NAC, RAC and Ssb, pending a usable structure.
- **Mitochondria:** III₂IV₁ supercomplexes, and supercomplex placement in the flat crista membranes specifically. The ATP synthase dimers lack subunits e and k.
- **Earlier items:** SPB plaques, ER polysomes, the Golgi/endosome split and frame time (patch and cable close-ups take about 45 ms per frame).

**2026-11-08: organisation layer without clustering.**
- **Network neighbourhoods:** STRING (≥ 400, weighted) personalized PageRank profiles (damping 0.5, seed entry dropped), compared by Pearson correlation; complexes take the mean over their subunits (complexes and heterodimers had no subunits before: their UniProt field is empty, so they are now mapped through gene names). Parameters from `scripts/organisation_sweep.py`: co-complex AUROC 0.996, shared localization 0.594 (all variants 0.56–0.61; log profiles slightly higher but not local).
- **Large scale:** the species correlation matrix is the covariance of the log-density fields (128 leading components, 54% of the trace; 128 random fields at 500 nm). Louvain modules are gone. Tile-density correlation follows network similarity (Spearman 0.58 over 86 abundant species; 0.71 vs 0.44 for modules on the 24 most abundant).
- **Small scale:** contact partners from the pool scores × (1 + r⁺)/2; loose partners are the most correlated species (r ≥ 0.4) instead of direct STRING links. Contact partner in contact 8.4% (control 1.3%; modules 6.4%), loose within 15 nm 19.9% (control 14.0%). 24% of molecules placed next to a partner.
- **Validation:** placed 97.8%, clash 2.5e-6 (2 beads across tiles), containment 1.2e-5, crowding CV 0.040.
- **Viewer:** colour "Network neighbourhood (STRING)" from a t-SNE map of the profiles; "Network neighbours" button in the inspector highlights and lists the 30 most correlated species.
- **For later:** the same profiles can drive a denser packing method (L1): they only supply a covariance for where species go and a partner list.

**2026-11-08: literature batch 1 (L2, L3, L5, L6) and an orientation fix.**
- **L2 measured localization (CYCLoPs):** compartment-group shares from the LOC-scores of the three wild-type GFP screens for 2,930 proteins; the annotation keeps deciding within groups. Not used where imaging and annotation share no group (534; GFP mislocalisation), for r-proteins and cell-wall proteins; new groups only from 0.5, never nucleocytoplasm, mitochondria only for the high-confidence proteome. 1.25 M of 42 M copies change compartment (proteasome RPN1 nuclear, COPI coat at the Golgi, actin-patch proteins cortical). Data: `pipeline/data/cyclops/`.
- **L3 mitochondrial proteome:** Vögtle 2017 submitochondrial locations (810 proteins, 252 changed) and Morgenstern 2017 (Table S1: high-confidence proteome removes 28 false mitochondrial assignments, e.g. Glt1, COPI; submitochondrial classes for 20 more). Generic "mitochondrion" copies 0.52 M → 0.02 M. Morgenstern copy numbers agree with PaxDb (Spearman 0.86, same mitochondrial mass share) on a 1.5× larger scale; PaxDb kept.
- **L5 lipids per membrane:** Zinser 1991 Tables 4–5 (phospholipid classes, ergosterol), cardiolipin as its own class (inner membrane 12%), PM ergosterol 40% (assumption) with PS/PE in the cytoplasmic and sphingolipids in the exoplasmic leaflet. The inspector lists the composition of the clicked membrane and leaflet. `web/app/src/lipids.ts`.
- **L6 mitoribosomes on the inner membrane:** 9HM0 74S monosomes, exit tunnel towards the membrane, exit site 5.1 nm above the leaflet; 1,010 of 1,051 placed (96%; free subunits used to fit poorly in the matrix). `pipeline/mitoribosomes.py`.
- **Bug fix, membrane orientation:** single-chain AlphaFold membrane proteins took their side only from "Cytoplasmic" topology domains; 66 proteins with other labels (mostly inner membrane) were oriented by chance and 39 were upside down (TIM50, TIM21, Mgm1 with their IMS domains in the matrix). Rebuilt with the membrane_orient side rule; their motion modes were rotated to match (`work/conformers/pool_modes_before_orientfix.npz` is the old fit).
- **Validation:** placed 98.0% (matrix 64%, was 62%), clash 2.6e-6 (none across tiles), containment 1.1e-5, TM residues 97.3% within 2 nm. Matrix 101 mg/mL (mitoribosomes moved to the membrane), Golgi/endosomes 82 (COPI), IMS 13 (no more generic 10% share).
- **Not done from the batch:** Morgenstern copy numbers (different proteome scale); within-group nucleus/cytosol splits from imaging (classifier unreliable there).

**2026-11-09: literature batch 2 (L7 eisosomes, L8 endocytic invaginations).**
- **Plasma-membrane invaginations:** new `cell3d.T_SLOT` primitive (lumen within r of a rectangle swept radially from a bottom axis), smooth-unioned into the PM field in all four evaluators (cell3d, geometry.ts, WGSL fill, packer via cell3d). Lumen counts as cell wall. Sites: no cortical ER, PM normal within 8° of radial.
- **L7 eisosomes:** 27 furrows (target 30; Karotki 2011), 284 ± 50 nm long, 50 nm deep, U-shaped with a half-cylinder bottom of 10.2 nm membrane radius. Coat from 8QB7/8QB8 (Kefauver 2024): the 7 dimers are screw-related (0.03 nm), giving the lattice (44.8 nm² per dimer); 4,687 dimers on the furrow bottoms (±60°), 12% of Pil1/Lsp1. Inserted residues 99.2% within 2 nm of the furrow mid-surface. Pil1 relocated from "lipid droplet" (UniProt, proteomics) to the PM. `pipeline/eisosomes.py`; "Eisosome" guided view and Show switch.
- **L8 endocytic invaginations:** one tube per actin patch (22), 40 nm wide, tip depth uniform 0–100 nm (Kukulski 2012/2016); patches now grow at these sites (`actin_patches.place(sites=...)`). The network ends up beside the tube more than around it. Coat proteins and Rvs161/167 not placed on the tubes.
- **Pipeline bug found:** numba `cache=True` kept the old inlined `cell3d.mem_sdfs` in pack3d/tiles3d after cell3d changed (furrow lumens packed as cytosol, containment 3e-4). `run_all.sh` now clears the numba cache first.
- **Validation:** placed 98.0%, clash 5.8e-6 (none across tiles), containment 4.2e-5 (7 of 30,000 sampled instances, all grazing nuclear-pore channels; residue beads vs sphere trees), TM residues 97.2% within 2 nm, densities unchanged.
- **Open from L7/L8:** coats hold only 12% of Pil1/Lsp1 (measured copies imply ~2,500 per eisosome); MCC membrane proteins (Sur7, Nce102, Can1) not enriched in furrows; endocytic coat; actin rings around the tubes.
- **L12 bud scar:** one (the cell is taken as a mother that budded once): `cell3d.T_TORUS` ring 1.0 µm across (AFM, Pillet 2014), 50 nm rim (assumed height), plane fitted to the surface under the ring (±6 nm); smooth union with the outer surface only (PM and wall layers unchanged), in cell3d, geometry.ts and the WGSL ray-march (scars in the Geo uniform). Rim filled with procedural chitin microfibrils (free kind 6, 36 GlcNAc per bead, 50% of the wall's 1.5% chitin, assumption); "Bud scar" guided view, "Chitin" Show switch. Birth scar not modelled (no chitin rim). Rebuild: clash 2.0e-6, containment 2.9e-5.
- **L20 peroxisome check (no change):** lumen holds 84k protein copies, close to the ~92k for two glucose peroxisomes in the 2025 absolute peroxisome-proteome preprint (bioRxiv 10.1101/2025.10.20.683451), so the 0.084 µm³ (84 mg/mL) is consistent with the protein load; fewer/smaller peroxisomes would need less cargo. The cargo is dominated by dual-localized GPD1, PNC1, LYS1 (100% peroxisomal: CYCLoPs conflicts keep the annotation) and PET9 (1/3 of the mitochondrial ADP/ATP carrier on peroxisomes from rank weighting of a minor UniProt location). Needs a decision on dual-localization weights.

**2026-11-09 (later): literature batch 3 (L13, L14, L15, L16, L18, part of L10) and two localization rules.**
- **Localization:** (1) the GFP rule keeps locations with experimental evidence from dedicated studies (papers cited as evidence in <40 UniProt entries) alongside the Huh 2003 GFP calls, in UniProt order (GPD1, PNC1, NPL3 regain their cytoplasm; 41 proteins); (2) high-confidence mitochondrial proteins with a mitochondrial primary location get 0.1x on other locations (PET9 5% peroxisomal, not 33%). 57 proteins, 133k copies; peroxisome cargo 130k -> 94k copies (66 mg/mL).
- **L13 Rabl chromosomes:** 32 arms (S288C centromere positions, checked against SGD), each a walk from its kinetochore-MT tip to a telomere site at the NE, with a drift that grows as the arm runs short (`pack3d.arm_walks`, `wc3d.rabl_segments`); 7 telomere foci, angle from the SPB axis by arm length (35-115 deg, assumption); XII-R into and out of the nucleolus. All 32 telomeres within 150 nm of the NE (median 55 nm from target); telomere-SPB distance vs arm length Spearman 0.71. "Chromosome arms" guided view.
- **L14 transcribing Pol II:** Pelechano 2010 per-gene engaged Pol II (691 per cell; Table S1 in `pipeline/data/pelechano2010/`), SGD gene coordinates (`pipeline/data/sgd/`); 552 placed on their genes in the arm walks, 3 nucleosomes evicted around each, DNA through the cleft (1Y1W frame on 1WCM, 0.01 nm), nascent mRNA threaded (52% of the length fits). `pipeline/pol2.py`.
- **L15 RAC on ribosomes:** all 48,900 RAC on ribosomes (27,500 in polysomes, 21,400 free) as `ribosome_80S_RAC`: 7Z3N rRNA on 4V88 (0.27 nm), C. thermophilum RAC as stand-in (yeast 7X3K has no ribosome and does not fit). NAC/Ssb stay free. `pipeline/rac.py`.
- **L16 mtDNA nucleoids:** 15 nucleoids, 35 genomes (3.0 Mb at 2 bp per bead) as 5 kb loop walks in the matrix; 90% of each within 62 nm of its centre. "Mitochondrial nucleoid" view. `pipeline/mtdna.py`.
- **L18 Cvt complex:** one cluster 40 nm from the vacuole: 176 Ape1 dodecamers, 105 Ams1 dimers (half of each, assumption), 195 Atg19; radius 74 nm (= Cvt vesicle size). "Cvt complex" view. `pipeline/cvt.py`.
- **L10 (part): ER ribosome density:** West 2011 measured 550-900 ribosomes/um2 on the cortical ER; free ribosomes added next to the ER/ONM up to 725 with the docked ones: ~570/um2 achieved (73% of the request fits). Fenestrations not built (West 2011 is qualitative about them).
- **Validation:** placed 97.8%, clash 2.3e-6, containment 5.0e-5 (all NPC-channel grazing), TM residues 97.1% within 2 nm.
- **Deferred:** L11 basket (8ZZC is an integrative bead model, needs an IHM reader), L17 polyP (whole-cell 23-250 mM, location disputed), L19 contact sites (no ER near mitochondria in this geometry; NVJ has 409 Nvj1), L9 (Golgi already dispersed cisternae).

**2026-11-09 (evening): L17 polyP, L19 NVJ, open items of L7/L8.**
- **L17 vacuolar polyphosphate:** procedural short chains in the vacuole lumen (viewer worker, free kind 7): 300 mM phosphate units (assumption within "hundreds of millimolar", Bru et al. 2025 eLife), chains lognormal ~24 units (Christ 2020), 4 units per bead, only inside the cut slab. "Polyphosphate" Show switch.
- **L19 nucleus-vacuole junction:** 409 Vac8-Nvj1 (5XJG) on 0.46 um2 of vacuolar membrane within 45 nm of the ONM (`pipeline/nvj.py`); guided view. ERMES etc. impossible here (no cytoplasmic ER).
- **L8 endocytic coat and amphiphysin:** Sla2 40 per site, Sla1 90 and End3 60 declining with depth (Picco 2015): 880/541/324 placed, random orientations; Rvs161-Rvs167 from 50 nm depth to 100 at 100 nm: 411 of 771 (RSA jams before a dense lattice). `pipeline/amphiphysin.py`. Picco counts ~3,000 actin per site at peak (model: 1,500).
- **L7 MCC proteins:** 80% of Sur7 and Nce102 (assumption) onto the furrow walls: 14,874 of 17,144 fit (`pipeline/mcc.py`).
- **L11 basket: not built.** 8ZZC resolves ~1/4 of each Mlp1/2 as coarse beads; AF monomers of Mlp1/2 are compact (39 x 17 x 14 nm), not coiled-coil rods. Documented in MODEL.md limitations.
- **Validation:** placed 97.8%, clash 6.3e-7 (none across tiles), containment 1.5e-5, TM residues 97.3% within 2 nm.

**2026-11-10: L1, proteome total by mass and a squeeze relaxation in the packer.**
- **Total:** proteins scaled to 5.0 pg per cell (measured 4–6 pg: Carroll 2011, Sherman 2002) instead of 4.2e7 molecules (2.7 pg). That is 7.8e7 proteins, 1.85e6 per µm³, and 301,400 ribosomes, matching the ~3e5 growth-law estimate (Kostinski & Reuveni 2021).
- **Squeeze** (`pipeline/squeeze.py`, `wc3d.squeeze_tile`):
  - After random addition, a tile with unplaced copies inserts them in rounds (≤ 15% of the tile's soluble molecules, 4 rounds) and relaxes all soluble molecules as rigid bodies (60 iterations of overlap removal, membranes push back along their distance gradient).
  - It then re-accepts exactly, as random addition does. A round that doesn't gain is undone.
  - Packing takes 83 min, against ~1 min before.
- **Polysomes:** the target mean size is raised so that chains use ≤ 85% of the 30,000 mRNAs. Before, `mrna.place` crashed on more chains than mRNAs.

| Metric | Before L1 | After L1 |
|---|---|---|
| Instances | 21.4 M | 35.8 M |
| Cytosol / nucleoplasm / nucleolus (mg/mL) | 91 / 61 / 111 | 163 / 103 / 137 |
| Cytosol protein + RNA (mg/mL) | 69 + 23 | 127 + 36 |
| Placed overall / outside mitochondria | 97.8% / 99.8% | 94.7% / 98.2% (fails the 99% acceptance) |
| Matrix / IMS / nucleolus / Golgi placed | 63% / 89% / 98.5% / 95% | 40% / 80% / 77% / 81% |
| Clash rate / containment | 6.3e-7 / 1.5e-5 | 9.3e-6 / 1.4e-4 |
| Tile data (gz) | 533 MB | 741 MB (816 MB data in total) |

- **Scaled with the total:** 41 endocytic sites (22 before), 15 actin cables, 758 Vac8–Nvj1, a 180 nm Cvt complex, 8,085 ER-docked ribosomes, 1,948 membrane mitoribosomes and 27,800 pre-ribosomes.
- **Simulation box** (120 nm cytoplasm): hydrated volume fraction 0.198 (0.106 before). GFP-sized D is 40 µm²/s with mean-field hydrodynamics. A stickiness of ~1.5 kT reaches the measured 8.6 µm²/s (2 kT before).
- **New issues:**
  - **Histone surplus:** histones beyond the fixed nucleosome count remain. 80k H2B are dropped (no structure); ~50k free H4 and histone dimers are packed in the nucleoplasm. They should be capped at the nucleosome need plus a small free pool.
  - **Small compartments saturate:** nucleolus 77%, Golgi 81%, matrix 40%.
  - **Polysomes:** only 44% of cytosolic ribosomes end up in chains (70% assumed), because chains run out of space in the denser cytosol.
  - **Validation:** clashes and containment are 10–15× worse than with random addition alone, though still small.

**2026-11-09 (night): simulation branch (`simulation`, worktree `yeast_cell_sim/`).** Simulate a box of the packed cell in the viewer; MODEL.md section 14 has the full description.
- **Model:** rigid-body Brownian dynamics on the GPU (WebGPU compute).
  - Each molecule is a rigid body, using either its sphere tree (schematic, µs) or its solvent-exposed residue beads (residue, ns).
  - D_t and D_r come from a Kirkwood R_h of the residue beads with Stokes–Einstein(–Debye).
  - Mean-field hydrodynamics scales D by the box volume fraction (Mereghetti & Wade 2012; Cichocki *et al.* 1999 factors).
  - Interactions: soft overlap repulsion and an optional cos² stickiness.
  - Membranes are walls, from a distance grid built with the renderer's analytic distance functions.
  - Fixed: membrane proteins, anchored pre-placed structures (encoder flag bit 1), nucleosomes, NPCs and an 8 nm frozen boundary shell.
- **Viewer:**
  - Box placement by click (30–200 nm). In sim mode, whole-cell rendering and streaming stop, and membranes are clipped to the box.
  - Live run, or compute first with a progress bar and then looping playback; scrubbing.
  - A diffusion readout (time-averaged MSD below the confinement scale) compared with yeast measurements.
  - Export of a CG mmCIF, a JSON manifest (sources, poses, R_h, D), the trajectory, `run_openmm.py` and `trajectory_to_dcd.py`. Tested: a 60 nm PM box runs in OpenMM 8.6 on OpenCL.
- **Fixes found by calibration:**
  - The neighbour-grid cap of 24 overflowed. Skipped neighbours gave asymmetric forces that self-propelled clusters, and ribosomes moved faster than in water with stickiness on. The capacity now comes from occupancy.
  - Time-step convergence: dt ≤ 0.2 ns without attraction, ≤ 0.05/ε ns with it. The panel enforces this.
- **Calibration** (120 nm cytoplasm box, pre-L1 data, φ = 0.106), GFP-sized proteins:
  - 77 µm²/s with excluded volume only (0.68 of water);
  - 62 with mean-field hydrodynamics;
  - 13 at ε = 2 kT, close to the measured 8.6 (S. pombe GFP; Delivani 2013).
  - This is the same conclusion as McGuffee & Elcock 2010.
  - **After L1** (merged into the branch; φ = 0.198): 62 µm²/s with excluded volume only, 40 with hydrodynamics, 9.0 at ε = 1.5 kT (measured 8.6).
- **Speed:** 0.49 / 0.32 / 0.08 µs per second for 80 / 120 / 200 nm schematic boxes (pre-L1; the denser 120 nm box runs at 0.18); 0.29 for a 50 nm residue box (Apple laptop GPU).
- **Open:**
  - membrane proteins' lateral diffusion;
  - mRNA tethering of polysomes (they move as free ribosomes);
  - explicit hydrodynamics;
  - electrostatics;
  - re-calibrating on the L1 packing.

**2026-11-10: simulation round 3 (branch `simulation`).** In response to the first hands-on test (molecules looked frozen, box placement clunky, switching models did nothing, residue mode looked better but was slow). MODEL.md section 14 has the details.
- **Speed (engine rewrite):**
  - Changes: two neighbour grids (small and large spheres), cell culling by distance, skipping of cells that hold only the scanning molecule's own spheres, and a parallel per-body force sum (the serial sum was half of residue mode's time).
  - Residue model, 60 nm box: 279 → 2,145 steps/s; with cores, 1,631 (634 with electrostatics).
  - Schematic model, 120 nm box: 643 → 1,546 steps/s.
  - Diffusion unchanged within noise: GFP-sized 41.3 vs 43 µm²/s at ε = 0, and 10.4 vs 10.0 at 1.5 kT.
- **Mid model:** clusters of ~4 residues / 2 nucleotides (k-means); radius = bead radius + 0.8 × rms member distance, which puts the cluster surface 0.13–0.22 nm outside the residue surface. It looks like the residue model, at 0.42 µs/s for a 60 nm box (0.24 with electrostatics). Converged at dt 0.2 ns.
- **Electrostatics:** Debye–Hückel between charged spheres (I = 0.3 M, van Eunen 2010; pH 7.2, Orij 2009; nucleotides −0.25 e, Manning). New pipeline step `pipeline/simdata.py` (called by s5) recovers residue identities from the source structures by geometric matching (98.6% of protein beads) and writes `library/charges.bin`. It costs about 3× in the residue model and 2.5× in the mid model, and changes diffusion little at this salt.
- **Specific binding:** native contacts (C-alpha within 0.8 nm) of 371 confident AlphaFold DB dimers whose partners are both monomer species (`library/contacts.{json,bin}`, median 61 contacts). Gō-like wells; one interface energy (20 kT default) shared over the contact pairs, so it does not depend on the model resolution.
- **Polysome mRNA links:** consecutive polysome ribosomes are recognised from the packing geometry (exit → entry linker, 3–12 nm, cone, entry direction) and joined by a worm-like chain: L_c 62 nm (one ribosome per 156 nt, Arava 2003; 0.49 nm/nt), l_p 1.5 nm (Chen 2012). Drawn as dark brown lines.
- **Calibration revisited (important):**
  - Electrostatics at 0.3 M change diffusion by 0–4%.
  - Without stickiness the three models agree, at 42–54 µm²/s for GFP-sized proteins.
  - A cluster count shows that stickiness slows molecules by aggregating them. In the schematic model at 1.5–2 kT, the values that match GFP, most of an 80 nm box joins one network within 1 µs. The mid model aggregates at 1 kT.
  - The residue model aggregates at 0.5 kT (285 of 446 molecules in one cluster).
  - So the earlier "1.5–2 kT matches GFP" calibration matched the number with a gel. The readout now says so and suggests Viscosity 5–6× as a stand-in for the missing many-body hydrodynamics (Ando & Skolnick 2010: full HI gives GFP D/D0 ≈ 0.08 without attraction).
  - Defaults stay at 0 kT and 1×.
  - In the mid model, stickiness also makes the integration stiff: at 1 kT, D is 17.9 at dt 0.05 ns but 11.8 at 0.01 ns; at 1.5 kT even dt 0.01 ns is unstable (complexes faster than in water). The step limit is now dt ≤ 0.01/ε ns there, and the slider is capped at 2 / 1 / 0.5 kT (schematic / mid / residue).
- **Core spheres:** with surface spheres only, molecules were hollow and strong stickiness pulled them through each other (grid overflow, molecules faster than in water). The fix:
  - each molecule now carries its sphere tree shrunk by 1 nm as core spheres, which meet only other cores;
  - grid capacity doubles on any detected overflow;
  - compute-ahead command buffers are sized to about 0.1 s of GPU work.
- **Tracking:** click molecules (up to 8) to track them. The others fade (screen-door, still clickable) and the tracked centre paths are drawn in colour; the panel lists displacement and per-molecule D. A Playback slider sets 1–60 frames per second.
- **Export:** the manifest carries charges, mRNA links and native contact pairs; `run_openmm.py` adds Debye–Hückel, WLC bonds and Gō wells.
- **Open:**
  - explicit hydrodynamics;
  - flexibility (IDRs);
  - links to ribosomes outside the box;
  - per-pair affinities;
  - electrostatics cost in the residue model;
  - calibration of the native-contact energy against measured affinities.

**2026-11-10: ensemble motion dropped; migration to the HPC cluster.**
- **Dropped:** the pool-ensemble conformational motion (stage 6: `conformers.py`, `pool_modes.py`, `s6_motion.py`, `library/modes.{bin,json}`, the viewer's modes buffer). It needed the foldcomp database of the pool predictions, which is not moving to the cluster. Display → Motion keeps the shader-only rigid wobble, low-pLDDT flexing and polymer wiggle for every molecule. The code is in git history (branch `conformational-dynamics`, last commit before removal b408a2b).
- **Inputs:** `pipeline/config.py` now finds the AFDB tar, the dimers and the pool pair scores in `inputs/` next to the bundle (`YC_INPUTS`), so no Mac paths are needed.
- **Handover:** `HANDOVER.md` in the migration archive describes the setup on the cluster, the state of the work and the history.

**2026-11-11: round H, data-driven validation on the cluster (branches `cluster-setup` → `gpu-bd` → `validation` → `localization-fixes` → `diffusion-model` → `condensates` → `sim-calibration`, linear, not merged).** Plan and pre-registered rules in section H above; MODEL.md 9.1 has the results, 14.2 and 14.5 the simulation side.
- **Cluster setup** (eu-g6-044): 64 cores, 375 GB RAM (not 192), 4 × RTX 4090, inside a Singularity container without Slurm or conda. `YC_NPROC=64`. Packing (s3) takes 18 min (83 min on the Mac). The viewer runs headless on a 4090 (Chromium with Vulkan WebGPU and locally unpacked libraries, `.chromelibs/`, gitignored).
  - Fix: the stages now start their worker processes with one BLAS thread each (`scripts/_common.py`). 64 workers × 64 OpenBLAS threads exceeded the per-user thread limit (4,096), and s5 hung for 6 h with crashed workers.
- **OpenMM on CUDA:** `run_openmm.py` produced NaN on CUDA. Cause: near-duplicate ribosome beads joined by elastic bonds of ~0 length and a WCA core that diverges. Elastic pairs closer than 0.1 nm are dropped and the WCA force is capped below 0.8 σ. A 60 nm cytoplasm box (175 k beads) runs at 700–1,100 ns/day per 4090.
- **`simgpu/`** (new): CuPy rigid-body Brownian dynamics of boxes cut from the packed cell (sphere trees from the library), with mean-field hydrodynamics or the full Ewald-summed polydisperse Rotne–Prager–Yamakawa mobility (Zuk overlap correction; Hasimoto splitting). Tests: the Ewald sum is independent of the splitting parameter to 2e-7, reproduces Hasimoto's periodic self-mobility, and is positive definite.
- **H1 localization benchmark** (`validation/localization.py`): SWAT N-GFP (native and NOP1 promoters) and hyperLOPIT, proteins split by ORF hash. The references themselves agree only 74% (SWAT native vs CYCLoPs), which is the ceiling. Dev error analysis: annotation-only proteins agree 55%, the 625 default-cytosol proteins 40%, peroxisome / lipid droplet / Golgi F1 0.10 / 0.03 / 0.28.
- **H7 nucleus rule (adopted):** the nucleus/cytosol split of nucleocytoplasmic proteins follows CYCLoPs' Nucleus − Cytoplasm score (concentration ratio exp(5 Δ), `inventory.NUC_K`, tuned on dev). Test top-1, assigned shares: SWAT native 64.0 → 66.8%, NOP1 72.6 → 73.8%, hyperLOPIT 84.6 → 84.6%. On the copies placed in the rebuilt cell: 67.0 / 74.1 / 84.5%.
- **H6 membrane anchors (rejected, switched off):** lipid anchors and lipid-binding domains as a membrane-bound fraction. Domain proteins made dev worse (most image as cytosolic); lipid anchors gained 0.1 point on dev and nothing on test. `membrane_anchor.py` stays, `BOUND_F = 0`.
- **H2 diffusion benchmark:** 93 measurements from 25 papers, each with a verbatim quote and DOI (`pipeline/data/validation/diffusion_measurements.tsv`); 35 in-vivo values of soluble molecules in unstressed cells are scored, papers split by DOI hash.
- **H3 far-field hydrodynamics:** full RPY slows GFP-sized molecules in a 120 nm cytoplasm box to 0.56 of water, against 0.62 without hydrodynamics and 0.36 with mean-field factors. So far-field many-body hydrodynamics does not close the gap; the earlier reading of Ando & Skolnick 2010 (D/D₀ 0.08) was at 300 mg/mL and is corrected in MODEL.md 14.5.
- **H4 diffusion model (adopted):** D = D₀(R_h of the construct, IDRs restored, tag added) × crowding(R_h) from BD sweeps of packed boxes (`simgpu/sweep.py`) × a compartment factor fitted on train tracers (cytoplasm 0.365, nucleoplasm 0.26).

  | Model, held-out papers | Ratio RMSE (log10) | Absolute RMSE (log10) |
  |---|---|---|
  | Stokes–Einstein × constant | 0.21 | 0.31 |
  | + crowding(R_h) from BD | **0.11** | **0.23** |
  | + binding (RNA-interactome, pool partners) | 0.18 | 0.28 (rejected) |

  40 nm GEMs: predicted 0.25 µm²/s, measured ~0.3.
- **Viewer:** default Viscosity 2.75× (= 1/0.365); the readout quotes budding-yeast GFP 8.7–21 µm²/s (FCS) and GEM ~0.3. Trimmed IDRs now enlarge R_h in the simulation (`library/rh_scale.json`: 3,102 types, median ×1.14). The inspect panel shows each soluble species' predicted D (`library/diffusion.json`: 3,428 species, median 5.3 µm²/s). Headless test on the 4090: GFP-sized molecules 14.7 µm²/s in an 80 nm schematic box, no console errors.
- **H5 condensates (negative result):** multicomponent second-virial stability from pool PAE likelihood ratios, RNA-interactome capture and IDR stickers; weights set on the nucleus only (nucleolus AUROC 0.68), then frozen. P-body members are enriched in the cytosol (AUROC 0.73), but the pool network adds nothing (RNA binding alone gives the same), the negative controls fail (stress-granule-only 0.72, glycolytic G-bodies 0.91) and the cytosol stays stable. No condensates are placed.
- **Rebuild** (s1–s5 on 64 cores) with the nucleus rule:

  | Metric | Before | After |
  |---|---|---|
  | Species / instances / beads | 5,615 / 35.76 M / 25.04 G | 5,619 / 35.77 M / 25.03 G |
  | Placed overall / outside mitochondria | 94.7% / 98.2% | 94.7% / 98.2% (still fails the 99% acceptance) |
  | Nucleoplasm / nucleolus (mg/mL) | 103 / 137 | 105 / 135 |
  | Clash rate | 9.3e-6 | 8.2e-6 |

- **Open:**
  - membrane proteins in the simulation diffuse 100–10,000× too fast (no lateral friction from the bilayer);
  - peroxisome, lipid-droplet and Golgi localization stay weak; the 625 default-cytosol proteins need another source (perhaps sequence-based predictors, scored on dev);
  - P-body specificity needs valence and network information that generic pair scores lack;
  - the L1 packing shortfalls (matrix 40%, nucleolus 77%, Golgi 81%) remain.

**2026-11-11 (afternoon): round I, imaging, proximity, pair statistics and a topology fix (branch `proximity-microscopy`, from `sim-calibration`; not merged).** Plan and item list in section I; evidence in `yeast_cell_dev_bundle/model_validation_and_tuning.md` (and its PDF in `docs/`).
- **Reviews** (three parallel agents; citations checked): cell biology, microscopy data, proximity data, in `literature_reviews_2026-11-11/`. Unretrievable data in `data_requests.md` (about 35 items, prioritised).
- **Pair-contact atlas** (`pipeline/contacts.py`): for all 35.7 M placed molecules, sphere-tree surface gaps of ≤ 2, 5, 10 nm to every other molecule on 4 GPUs (822 M spheres, about 10 min); per protein pair P(B within g of A) and mean neighbour counts.
- **Pair statistics for kinetic models** (`pipeline/encounters.py`): contact probability, compartment and local concentration of B around A, Smoluchowski encounter rate with in-vivo D; table for 131 kinases × all proteins (`work/contacts/encounters_kinases.parquet`). The well-mixed estimate predicts the snapshot with r = 0.94. Brownian dynamics (`simgpu/encounter_bd.py`): new encounters follow long-time D, re-collisions inside the crowding cage are 2–4× more frequent (short-time D).
- **Proximity tests** (BioGRID 5.0.263; nothing tuned): median per-bait AUROC of the snapshot 0.82 (proximity labelling), 0.82 (PCA), 0.99 (XL-MS), against 0.71/0.78/0.83 for shared compartments; in-vivo kinase–substrate pairs 0.66, in-vitro array pairs 0.52. PPiSeq with 1.47 M tested pairs: at equal abundance compartments 0.64, snapshot 0.63 (the 3D arrangement adds nothing for PCA).
- **Simulated microscopy** (`pipeline/microscopy.py`, `validation/microscopy_gallery.py`): Debye PSF, widefield/confocal, z-stacks; foci per cell vs the Mattiazzi Usaj 2020 screen (idr0078): Sac6 4 vs 5, Vph1 2 vs 1, Snf7 5 vs 2, Sla1 11 vs 6. Errors exposed: Pil1 not punctate (coats hold 7%), Spc42 diffuse (no SPB plaques), Erg6 not on lipid droplets.
- **Ribosomes vs cryo-ET** (`validation/ribosome_cryoet.py`): 17 µM vs 23 µM (Delarue 2018); nearest neighbour 28.6 vs 29.4 nm (CZ cryoET 10007, 249 tomograms); the model's spacing is too regular.
- **Kinase localization test** (`validation/kinase_localization.py`, 15 positives, 6 stress controls from dedicated studies): 10/15 and 6/6 (were 9/15 and 4/6). Failing: TORC1 (Kog1 39% on the vacuole face), Sch9, Bcy1 (39% nuclear), Pho85, Pkh1.
- **Fixes (rebuilt):**
  - soluble proteins of the vacuole, Golgi and peroxisomes now go to new 20 nm shells on the cytosolic face (regions 24–26) unless they can reach the lumen (signal peptide, Cvt cargo, peroxisomal matrix); 680 k copies (V1, EGO, Tco89, Vac8, COPI, clathrin, Arf1/2) had been inside the organelles;
  - no cell-wall share without a signal peptide (132 k Ssa1/2; Fus3, Kss1);
  - the 26S proteasome takes its subunits' nucleus/cytosol split (was 50/50).
- **Also:** OpenBLAS threads per worker capped (s5 had hung at the per-user thread limit); MODEL.md numbers refreshed (chromatin, clashes, structure shares had been from an older build); `scripts/md_to_pdf.py` renders Markdown documents to PDF.

  | Metric | Before | After |
  |---|---|---|
  | Placed overall / outside mitochondria | 94.7% / 98.2% | 94.5% / 98.0% |
  | Golgi lumen / Golgi face placed | 81% / — | 95% / 73% |
  | Vacuole face placed | — | 98% (164 k copies) |
  | Clash rate | 8.2 × 10⁻⁶ | 1.0 × 10⁻⁵ (at the limit) |
  | Cytosol (mg/mL) | 163 | 165 |
  | Localization test top-1 (placed; SWAT native / NOP1 / hyperLOPIT) | 67.0 / 74.1 / 84.5% | 67.0 / 74.1 / 84.5% |
  | PCA / XL-MS per-bait AUROC (snapshot) | 0.80 / 0.97 | 0.82 / 0.99 |

- **Open (section I):** I10 eisosome coats, I11 SPB, I12 kinase locations, I13 polysome compaction and ribosome concentration, I14 membrane diffusion, I15 dilute nucleoplasm, I16 endocytic sites, I17 ER band ribosomes, I18 contact sites; the Golgi coat pool (how much is membrane-bound).

**2026-11-11 (evening): integration of the requested data (branch `data-integration`, from `proximity-microscopy`; not merged).** The user supplied most of `data_requests.md` (`data_requests_validation_2026-11-11.zip`, unpacked into `pipeline/data/validation/`; README there). Two agents extracted numbers from PDFs and Word files into TSVs with verbatim quotes (`membrane_diffusion_measurements.tsv`; `literature/*/` with `literature_numbers_README.md`).
- **Membrane-protein diffusion calibrated (I14):** a factor per membrane on the simulation's Saffman–Delbrück value, from FRAP of transmembrane proteins in unperturbed cells: plasma membrane 5.9 × 10⁻⁴ (Sso1, Snc1 0.0025–0.0029 µm²/s; Valdez-Taubas & Pelham 2003), vacuole 0.016, ER 0.035 (Sec61; Meinema 2013), nuclear envelope 0.079; other membranes take the ER factor (assumption). `pipeline/membrane_diffusion.py` (s5) → `library/membrane_diffusion.json`; viewer option "… calibrated per membrane" (on). A fit, not a test: too few papers to hold out.
- **Diffusion benchmark:** free GFP of Meinema 2013 added (test paper by hash): predicted 12.7 vs 13.6 µm²/s; test absolute RMSE 0.22.
- **Copy numbers:** model (PaxDb) vs Ho 2018 and Das 2026: same ranks (Spearman 0.81–0.93), different dynamic range (10–90%: Ho 18×, model 94×, Das 463×); the input sits closer to the independent proteome than Ho does. Endocytic proteins by FCS (Boeke 2014): Spearman 0.76, model 1.9× higher (FCS sees the free pool).
- **Kinase–substrate (KID, Sharifpoor 2011, without its co-localization evidence):** in-vivo pairs 0.64 / 0.61 (compartments / snapshot), in-vitro-only 0.56 / 0.53, genetic-only 0.54: the PhosphoGRID result replicated.
- **Negative results:** quantitative mitochondrial BioID vs matched control (Singh 2020): no prediction inside mitochondria (AUROC 0.44–0.51); FCCS binding pairs (Boeke 2014): contacts do not predict binding (0.39–0.50).
- **Counts:** translocon-engaged ER ribosomes 7,703 vs 7,715 by cryo-ET (Wu 2026); TORC1 66 vs ~200 (Prouteau 2017).
- **Geometry:** soft X-ray (Chen 2025, G1 cells): vacuole 17%, nucleus 5.9%, lipid droplets 0.37%; relative absorption nucleus 0.90 and mitochondria 1.10 of the cytosol vs model densities 0.64 / 0.65. Mitochondria 4.77% of cell volume in glucose (Di Bartolomeo 2020, CEN.PK) vs the model's 2%: the leading explanation of the 40% matrix fill. Not yet adopted.
- **Not used:** Tarassov's full matrix (binary, only among proteins with ≥ 1 interaction; PPiSeq is the better negative set); Xing 2020 P-body partition coefficients (dcp1Δ strains; for a future P-body build); Slaughter 2007 K_d table (main paper missing).
- **Next:** raise the mitochondrial volume toward the measured fraction and rebuild (I18 → new item); TORC1 copy number and vacuolar placement; nuclear/cytoplasmic intensity test with CYCLoPs intensities.

**2026-11-11 (late): rounds K and M, chromatin and the genome, mRNA (branch `data-integration`; not merged).** Plans and item tables in sections K and M; evidence in `yeast_cell_dev_bundle/model_validation_and_tuning.md` sections 5 and 5a. Three full rebuilds (K, L, M; results of K and L in `work/prev_K`, `work/prev_L`).
- **Chromatin (K1–K7):**
  - every nucleosome carries its genome position (K1);
  - the fibre is built from a measured nucleosome map (K2): Chereji 2018 chemical-cleavage dyads, with Rossi 2021 promoter NDRs removed, giving 64,487 on the arms;
  - per-nucleosome backbone advance and linker reach;
  - Pol II, MCM double hexamers + ORC on 410 OriDB origins, kinetochores (CCAN + 6 MIND) and about 280 ChIP-exo proteins sit at their loci;
  - seven histone marks and H2A.Z are attached as annotation.
- **Fibre: a failure, then a partly tuned choice:**
  - the first fibre (Bystricky 2004 compaction and persistence, 20 nm radius) failed the held-out Hi-C test: cis 0.62 vs 0.37–0.41, P(s) slope −1.6, CEN clustering 13×;
  - nine pre-pass variants were scanned, and the 60 nm radius kept: it is the only one inside the Arbona 2017 R(s) envelope at 20 and 100 kb;
  - Hi-C is labelled "partly tuned" from here on.
- **Telomere steering bug (found in L, fixed in M):**
  - arms ended a median 760–790 nm from their telomere foci: random turns after failed placements erased the drift, and short arms had targets beyond their contour;
  - fix: failure turns keep a pull to the target, the drift cap is 5, and targets lie within 80% of the arm contour;
  - in M the miss is 149 nm, and telomere-to-SPB distance vs arm length has Spearman 0.92 (L 0.52);
  - L's P(s) slope (−1.33, close to log-phase Hi-C) turned out to be an artefact of the bug; M gives −0.92.
- **mRNA (M1–M4, M6):**
  - gene-specific mRNAs with per-gene UTRs (TIF-seq), poly(A) (PAL-seq) and TE-weighted ribosome load;
  - flexible single-strand polysome linkers;
  - polysomes of 170 genes tethered to the mitochondrial outer membrane (Williams 2014);
  - Pab1 and eIF4E on the ends;
  - ribosomes on mRNA rose from 44% to 62% of cytosolic 80S.
- **Clash fix:** in K, mRNA tails ran through the RAC on ribosomes (clash 1.2 × 10⁻⁵, over the limit); the RAC sphere trees are now thread obstacles (M: 7.7 × 10⁻⁶).

  | Held-out test | K | L | M | Data |
  |---|---|---|---|---|
  | Hi-C cis share | 0.62 | 0.53 | 0.51 | 0.37 (G1), 0.41 (log) |
  | P(s) slope 10–200 kb | −1.60 | −1.33 | −0.92 | −1.09, −1.27 |
  | Centromere clustering | 13× | 9.6× | 7.6× | 4.0×, 5.3× |
  | Chromosome-pair Spearman (G1 / log) | 0.12 / 0.18 | 0.06 / 0.19 | 0.27 / 0.04 | — |
  | Nucleosome occupancy AUROC (Oberbeckmann 2019) | — | 0.75 | 0.75 | — |
  | MNase dyads matched within 40 bp (Weiner 2015) | — | 79% | 79% | — |
  | Mito mRNA, Spearman / AUROC (Saint-Georges 2008) | −0.05 / 0.47 | 0.38 / 0.62 | 0.44 / 0.66 | — |
  | Tight polysome fraction (cryo-ET) | 1–4% | 1.4% | 1.4% | 14.5–18.6% |

- **Build M:** 35.7 M instances, 95.9% placed (97.9% outside mitochondria; still below the 99% acceptance), clash 7.7 × 10⁻⁶, cytosol 168 mg/mL. Placed in M:
  - 62,815 nucleosomes, against 57,000–60,000 by absolute occupancy (a population map overcounts a single cell);
  - 358 licensed origins;
  - 22,195 chromatin proteins on peaks;
  - 21,888 Pab1 and 25,124 eIF4E.
  - Proximity benchmarks unchanged.
- **Open:**
  - chromosomes stay too territorial, the P(s) slope is now too flat, and centromeres cluster too tightly (all arms start at the SPB; the kinetochore microtubules may be too short or too aligned);
  - the chromosome-pair pattern is not reproduced (the two maps also disagree with each other);
  - no cohesin, condensin, tRNA genes or Sir/Rap1 at telomeres;
  - marks are annotation only (no readers placed, no viewer colouring);
  - tight polysomes need stacked ribosome orientations;
  - mRNP compaction, eIF4G/closed loop;
  - M5 ER polysomes (docked ribosomes are single);
  - per-gene ribosome density test (Arava 2003 tables not found; in `data_requests.md`);
  - the nuclear Brownian-dynamics run died on a pair-buffer overflow and needs a rerun.

**2026-11-12 (night): round N, ideas from the Jurkat cell model and the open items of K and M (branch `jurkat-ideas`, from `data-integration`; not merged).** Plan in section N; evidence in `yeast_cell_dev_bundle/model_validation_and_tuning.md`. Build N; results of M in `work/prev_M`.
- **Packing (N1):**
  - 2,000 tries per copy and 1,200 failures in a row (was 200 and 25) for all regions; raising the budget only in thin regions let shells take room from the cytosol;
  - half of the tries in thin lumens and shells align the molecule's flattest axis with the membrane normal;
  - placed 96.5% overall and 98.5% outside mitochondria (M 95.9 / 97.9%), still below the 99% target; Golgi face 71 to 76%, envelope periphery 81 to 90%, nucleolus 78 to 81%, matrix 62 to 65%;
  - the pack ran faster (1,214 s vs 1,462 s).
- **ER topology fix (found while checking pools):** soluble ER-annotated proteins without a signal peptide sat in the ER lumen (181,000 copies of 72 proteins: COPII, Sar1, Cdc48, Ypt1, Get3, SRP subunits, Scp160). They now go on the cytosolic face, the rule the other organelles got on 2026-11-11; the lumen keeps Kar2, Pdi1 and the like.
- **Chromatin polymer (N2), tested and not adopted:**
  - OpenMM coarse polymer (2 kb beads, Arbona chain, Rabl tethers, confinement), with the measured nucleosomes back-mapped along it (`chromatin_polymer.py`, `pack3d.guided_walks`, `validation/chromatin_polymer_scan.py`);
  - better contact decay, but more territorial (cis 0.63–0.77) and more centromere clustering (8–20×) in every variant; the fibre walk stays.
- **Hi-C normalisation:** the data measures move with normalisation (log-phase cis share 0.41 to 0.51, telomere clustering 26× to 3.6× after ICE), and the two maps agree on chromosome pairs at only −0.04 / 0.15. hic_test now also scores everything ICE-balanced. Balanced, build N lies inside the two maps' range for P(s) (−1.16), centromere (5.5×) and telomere (2.5×) clustering; only the cis share (0.54 vs 0.39–0.51) is still high.
- **New structures:**
  - SPB plaques (`spb.py`): Spc42-GFP is now one spot (was diffuse);
  - outer kinetochore (16 DASH and 8 Ndc80-complex copies per kinetochore);
  - Pol III and TFIIIC on tRNA genes, Rap1 at chromosome ends, Sir2/3/4 on subtelomeres and HML/HMR (`trna_telomeres.py`; chromatin is crowded, 39% fit);
  - 20 transport receptors per pore channel before the FG mesh (`pore_receptors.py`).
- **mRNA:**
  - ER-translated mRNAs (M5): 2,489 polysomes of 603 Jan 2014 ER-enriched genes tethered to the cortical ER. Held-out Chartron 2016 membrane enrichment: Spearman 0.35, AUROC 0.91 (M: −0.08, 0.49). ER ribosome density about 790 per µm² (West 2011: 550–900; M about 470).
  - Polysome linkers with exponential free lengths, and up to 40 orientations per drawn linker: tight pairs 1.4% to 7.0% of cytosolic ribosomes (cryo-ET 14.5–18.6%); 89% of polysome ribosomes placed (M 88%).
- **Viewer:** seven histone-mark colour modes (marks as per-nucleosome instance flag bits written by s5). Checked headless on the 4090; no console errors.
- **Tools:**
  - `validation/ledger.py`, a per-gene ledger of placed shares by build (snapshots M, N; M to N moved only the intended genes: Nud1, Cnm67, DASH, Sir2);
  - localization share agreement, plain and mass-weighted: test 0.63–0.82 / 0.73–0.90.
- **Nuclear BD rerun:** at 20 nm separation, new encounters are 0.66× the well-mixed long-time prediction in the nucleoplasm and 1.04× in the cytoplasm (chromatin obstruction: D_long/D_short 0.15).
- **Not redone:** diffusion with binding. A partner-binding term was tried before and did not help on held-out rows.
- **Build N:** 36.0 M packed instances, clash 8.4 × 10⁻⁶, cytosol 170 mg/mL. Unchanged:
  - nucleosome tests (AUROC 0.76, 79% MNase match);
  - mitochondrial mRNA test (0.41 / 0.66);
  - proximity (0.84 / 0.81 / 0.99);
  - kinase test (10 of 15, 6 of 6).
- **Open:**
  - tRNA genes: FISH puts them at the nucleolus, the balanced Hi-C maps show no tRNA–tRNA contact excess (1.00–1.05×; model 0.86×); no rule added;
  - stacked orientations of tight polysomes;
  - +TIP and lattice bands on microtubules;
  - condensin;
  - an image-based nucleus/cytoplasm test (needs per-protein images with a nuclear marker);
  - the 99% placement target;
  - the SPB, kinetochore, Rap1 and receptor-per-pore numbers are unchecked (data_requests.md).

**2026-11-12: merged and published (user's request).** All branches of rounds H–N fast-forwarded into `main` (`e6c82fe`), plus `pages-readme` (`e92f435`: the Pages README no longer states a molecule count; it said 21.6 million). Build N exported, smoke-tested headless on the 4090 (3 guided views and a mid-resolution sim run, no console errors) and pushed to `zoidpuff/yeast_cell_viewer` (`a52985d`); site 782 MB, Pages `.git` 2.1 GB (squash soon). The export drops the ensemble-motion files (`modes.bin/json`), as decided on 11-10. Setup notes for pushing from the cluster are in HANDOVER section 6.

**2026-11-12: what an mRNA or a piece of DNA is, in the viewer (branch `polymer-identity`, from `main`; not merged).** The user found that hovering an mRNA or DNA on the live site gave no gene, transcript or position, although the model has them since rounds K and M.
- **Cause:** the tiles hold polymers as anonymous polylines (a picked chain is its index in its tile), and nucleosomes, polysome ribosomes and Pol II as plain instances; the identities stayed in the pre-pass.
- **Fix:** `pipeline/polymer_annot.py` (s5, or `s5_encode.py --annot-only` for an encoded build) writes per tile, in tile order, what each chain and chromatin instance is, plus `library/genome.json` (SGD features, gene layouts, transcripts). Tiles unchanged; 4,842 small files (4.1 MB). The viewer loads them on demand.
- **Checked:** against all 4,842 decoded tiles (chain counts, first points, lengths; 240,448 instances by type and position). Headless on the 4090, hover and inspect panel name genes and positions for mRNA, linker DNA, nucleosomes, ribosomes and Pol II; no console errors.
- **Model error found:** CDS lengths are SGD ORF spans, which include introns. The mRNAs of 208 intron-containing genes, mostly ribosomal-protein genes (18% of mRNA molecules, 17.5% of polysome ribosomes), are too long and carry too many ribosomes; 5.2% of mRNA nucleotides are intron. Documented in MODEL.md (6.5, section 11). The fix is to take the CDS from the summed SGD CDS features in `mrna_genes.gene_table`. It needs a rebuild, and publishing a rebuild rewrites all tiles, so I asked the user first.

**2026-11-12: build O, mRNA layouts without introns (branch `mrna-introns`, merged and published at the user's request).**
- **Fix:** the CDS of each gene is the sum of its SGD coding segments, not the ORF span; TIF-seq 5′ UTRs (genomic distances) lose a 5′ UTR intron where they are longer than it (21 genes). RPL25 843 → 429 nt, RPS8B 5′ UTR 389 → 29 nt.
- **Build O:** 5.8% fewer mRNA nucleotides (27.7 M), 168,017 ribosomes on 25,115 mRNAs (N 171,872 on 25,187; median 5 per mRNA, N 6), ribosome spacing 164 nt of CDS (N 177; Arava 2003 average 156). Tight polysome pairs 7.2% (N 7.0%). Unchanged within noise: placement 96.5% / 98.5% outside mitochondria, densities, overlap rate 7.3 × 10⁻⁶ (N 8.4), Hi-C and nucleosome tests (chromatin identical), mRNA tests (ER 0.35 / 0.91, mitochondria 0.41 / 0.65); the per-gene ledger N → O moves no protein mass between compartments.
- **ER ribosome density** is now `validation/er_ribosome_density.py` (in validate_all): the 20,400 ribosomes within 30 nm of the cortical ER over one face of the ER sheet (26.1 µm²) give 780 per µm² for both N and O; round N's ~790 came from a one-off area estimate.
- **Arava 2003 per-gene ribosome densities** (the natural test of this fix) are still unavailable: Europe PMC says the article is not open access and PMC returns 404.
- **Published:** build O is live (deploy checked). The Pages history was squashed to one commit (`dfe3dbc`) after the user confirmed the force push; the old history is in the backup bundle (HANDOVER section 6).
- **Mistake:** archiving build N, a copy overwrote `work/run_all_N.log` with an older run log; build N's run log is lost (its results are in `work/prev_N/`).

**2026-11-12: round-2 data delivered (`data_requests_round2_2026-11-12.zip`), copied into `pipeline/data/validation/` (README_round2_2026-11-12.md); not yet integrated.** First reading: Joglekar 2006 Table 1 confirms 8 Ndc80 complexes, 16–20 DAM–DASH and 6–7 MIND per kinetochore (the model's numbers); Viswanath 2017 gives the SPB core stoichiometry Spc42:Spc29:Spc110:Cmd1:Cnm67 = 2:1:1:1:1 (Muller 2005) but no copies per SPB; Marcand 1997: Rap1 sites every 18 bp on ~300 bp telomeric tracts (about 17 per end; the model has 15); Swygert 2014: two Sir3 monomers per nucleosome in vitro (the model has one); Berger 2008 published locus maps only as figures; no yeast count of receptors per pore was found. Disome (Meydan & Guydosh 2020, Zhao 2021) and tRNA–nucleolus (Thompson 2003, Haeusler 2008) data are for tests.

**2026-11-12 (evening): round P, build P (branch `round-p`, from `main`; not merged, not published).** Plan and item table in section P.
- **Working mode:** tests that only read the model ran on build O; the model changes (P1-P3, P11, P12) were batched into one rebuild. Two quick checks inside the rebuild caught an SPB bug (the core was first scaled to the whole Spc110 pool, then the proteins drifted from 2:1:1 when the disc filled), so the build was restarted twice within its first minute.
- **Round-2 data used:** Joglekar 2006 (kinetochore counts confirmed), Viswanath 2017 / Muller 2005 (SPB core stoichiometry), Marcand 1997 (Rap1), Swygert 2014 and 2018 (Sir), Zhao 2021 and Meydan & Guydosh 2020 (disomes), Thompson 2003 and Haeusler 2008 (tRNA genes at the nucleolus). Not usable: Berger 2008 (locus maps only as figures); no yeast count of receptors per pore.
- **ETH IBC labs (survey of Michaels, Weis, Neurohr):** Gade 2025 and Kunzi 2026 GEM values as held-out diffusion rows (model 1.2-1.8x too fluid); Kunzi 2026 dry mass consistent, its protein concentration not; Ohsawa 2026 nucleus 5.3-6.3% (model 8%). Open ideas: the Gade perturbations (P10), mRNP diffusion (Joyner 2016), QuPID-MS condensation scores (Narduzzi 2026; table to request), a colloid-osmotic N/C check.
- **Collided disomes:** new pose from 6I7O; 3.9% of cytosolic ribosomes collided vs 5.8% measured; tight pairs 10.4% (7.2% in O) vs 14.5-18.6%.
- **tRNA isoacceptors (user's request):** 42 types, Boltz-2 models (confidence 0.73-0.89; tRNA-Phe 1.2 Å from 1EHZ), split by gene copies; viewer labels, colours, toggle and Mol* model loading work (headless test).
- **Build P:** placement 96.6 / 98.5%, densities unchanged, clash 9.3 × 10⁻⁶ (53 beads; s4 now lists the type pairs: none involve tRNAs or collided pairs), localization, proximity, kinases, nucleosome and mRNA tests unchanged. Hi-C moved with the new chromatin realisation (the SPB is an obstacle for the fibre walk): balanced centromere clustering 7.1× (O 5.5×), slope −1.03 (O −1.16), documented as realisation noise.
- **Open:** P10; tRNA genes at the nucleolus (a rule would consume the test data); SPB absolute copy numbers; Arava 2003 per-gene ribosome loads (not open access).

**2026-11-12 (late): round Q, build Q (branch `round-q`, from `round-p`; not merged, not published).** Plan in section Q. Build Q is archived in `work/prev_Q/` (the working tree moved on to round R).
- **Raw images (Q1, Q5):** YeastRGB C-SWAT archive (9.7 GB, 6,145 strains with a Nup49-mCherry nuclear-envelope marker). `validation/image_nc.py` segments the nucleus per cell from the marker and measures nucleus/cytoplasm intensity; 1,680 strains are above autofluorescence, 1,137 comparable with the model. Held-out half: Spearman 0.81, AUROC 0.98 (model-nuclear against model-cytosolic). Mean intensity against model copies: Spearman 0.74.
- **Silent chromatin (Q2):** Sir proteins on the nucleosomes where Sir3 ChIP-seq is >= 5x input in both replicates, within 30 kb of an end or in HML/HMR (Brothers & Rine 2022; `pipeline/silent_chromatin.py`). Held out: H4K16 acetylation (Weiner 2015) and Sir4 ChIP (Goodnight & Rine 2020).
- **DNase (Q3):** per-base cuts of Zhong 2016 lifted to sacCer3 (`validation/dnase_test.py`).
- **Polysome links (Q4):** entropic springs between neighbouring ribosomes of one mRNA in the GPU engine (`simgpu/engine.py` links; `simgpu/gade_test.py`: 6 boxes x 3 seeds x 10 us, pooled MSDs). Disassembling the polysomes does not speed up 40 nm or 20 nm particles in the simulation (1.02x and 0.98x), against 1.4-1.5x and 1.3x measured (Gade 2025): a failed test. The mRNA is not a simulated body and the simulated lags are microseconds.
- **PDB audit (Q6, Opus subagent; the user's request):** report and table in `pdb_audit_2026-11-12/`. Fixed here: Asc1 and Stm1 booked to the 80S; homo-oligomers Cdc48, Tsa1, Gln1, Mmf1, Hmf1, Pst2, Ynk1 from the PDB. Open: a translating 80S with tRNAs, eEF3, Ssb, eIF5A, 43S/48S, remodelers on nucleosomes, and others.
- **External review of build O (`../yeast_cell_viewer_review.md`; answer in `../yeast_cell_viewer_review_response.md`):** fixed in Q: genotype and condition masks (73 proteins), "spindle pole" only for SPB components, GPI proteins to the wall side, abundance-aware dimer allocation (NAC 122,000 dimers), native histone pairs only, nucleoporin copies per pore from Kim 2018 (the pore had booked one copy of each), stale numbers, files missing from the published source. `validation/anchors.py` reports the stoichiometric anchors.
- **Viewer (the user's request):** mRNA copies per cell and a copy-number histogram on the info card.
- **Build Q against P:** 35.3 M instances; placed 96.7% (98.6% outside mitochondria); clash 3.7 x 10^-6; nuclear-envelope periphery 94% placed (90%); proximity labelling AUROC 0.85 (0.84); mitochondrial mRNA test 0.44 / 0.68 (0.41 / 0.66); chromatin, nucleosome, kinase, localization and image tests unchanged.
- **A crash on the way:** the new dimer allocation left no free actin, and the actin the filaments return had no monomer type; fixed in s3.

**2026-11-12 (night): round R, build R, a consistent G1 cell (branch `round-q`; not merged, not published).** Plan in section R; all of R1 to R7 done. Before/after table for every test: `model_validation_and_tuning.md` section 15a. Build R is the working tree (`work/yc3d`); its validation outputs are also in `work/prev_R/`.
- **What changed (R1-R5):** protein 3.5 pg (5.0), cell 34 µm³ (42), vacuole 17.2% (9%), nucleus 7% (8%), 93 pores by areal density (117), tRNA and mRNA × 0.7; `inventory.CELL_SHARE` carries the factor. The cytosol keeps its density (168 mg/mL), so the diffusion calibration was left alone.
- **Build R:** 5,559 species, 24.9 M instances, 17.9 billion beads; placed 97.5% (98.9% outside mitochondria; matrix 75%, was 66%; Golgi face 84%, was 76%); 211,550 ribosomes, 61% of cytosolic ones on 17,566 mRNAs; site export 599 MB.
- **Better:** anchors (core histones 0.90× need, was 1.29; nucleoporins 1.36×, was 1.5; ribosomes 1.06-1.13× the rRNA counts, was 1.5-1.6); copies against Das 2026 (median 1.23×, was 1.75); matrix density; packing completeness; Cvt complex size (162 nm, was 180).
- **Unchanged:** localization, image N/C test (0.80 / 0.98), proximity, kinases, diffusion benchmark, nucleosome and DNase tests, silent chromatin, ER mRNAs, ribosome concentration (17 µM against 23: cytosol and ribosomes shrank together).
- **Worse:** R(s) at 100-200 kb (376 and 467 nm; Arbona's range starts at 393 and 561); mitochondrial mRNA test 0.35 / 0.61 (0.44 / 0.68; fewer mRNAs per gene); dry-mass sum 240-255 mg/mL (225-240; measured 200-225); translocon-engaged ribosomes 5,507 against 7,715 (7,700 before); TORC1 51 against ~200.
- **Overlap check:** 2.65 × 10⁻⁵ in the first run, above the limit. Four cube seeds showed that the excess is overlap inside the built actin structures (tropomyosin on filaments, filaments in bundles, Arp2/3 at branches; about 0.25% of their beads, in every build). s4 now reports packed molecules (6.4 × 10⁻⁶, within the limit) and the total separately; the total still fails and the actin builders are not fixed.
- **Bugs found on the way:** `FIXED_COPIES` booked 117 pores while the geometry had 93 (now `inventory.N_NPC`); a glucose-repressed protein in LEGACY_OLIGO broke stage 1; cytosol volume 26.62 µm³ was hard-coded in several validation scripts (now `pipeline/volumes.py`).
- **Stale inputs found during the documentation pass (2026-11-12, after the build):** `work/features/protein_features.parquet` had not been regenerated since 2026-11-11, so `copy_numbers.py` and the abundance baselines of the proximity tests read 5.0 pg numbers (my first summary said "copy-number checks unchanged"; that was the stale table). `validation/diffusion.py` had been crashing since round Q on two free-H2B rows whose species (the non-native H4-H2B pair) is no longer built, so "diffusion benchmark unchanged" was also not a re-run. Fixed: `run_all.sh` regenerates the features after stage 1; the two rows use the H2B-H2A.Z dimer; `validate_all.sh` now runs `diffusion.py` and `anchors.py`; `literature_checks.py` reads the organelle fractions from the build. Re-run: proximity and kinase results unchanged, diffusion test error 0.09 / 0.21 (ratio / absolute, log10).
- **Lesson:** a validation script that is not in `validate_all.sh`, or that reads a table the build does not regenerate, silently reports the previous build. Check file dates of every JSON quoted.
- **Documents:** MODEL.md, `model_validation_and_tuning.md` (+ PDF and figures), HANDOVER.md and the review response are at build R.
- **Open:** the cell volume (dry mass asks for 36-38 µm³, the ribosome concentration for about 26); long-range chromatin compaction in the smaller nucleus; tRNA genes at the nucleolus; Gade polysome test; actin-structure overlaps; G1 pore density (Winey 1997 not readable here); the PDB-audit candidates (translating 80S first); review comments 12-17 and 20.

### 2026-11-13: build R merged and published; round-3 data, first look

- **Published (user's request):** `round-q` (rounds P, Q, R) fast-forwarded into `main` (`5138268`); Pages export of build R tested headless on the 4090 (3 views, info cards with histograms, a mid-resolution simulation run; no console errors), pushed as `800c6dd`, deploy succeeded, live `version.json` = local (`20261113100052`). The Pages `.git` now holds two builds (O and R).
- **Round-3 data** (`data_requests_round3_2026-11-13.zip`, filed under `pipeline/data/validation/`; index `README_round3_2026-11-13.md`). First look only; no model change, nothing tuned:
  - **Pores per G1 nucleus (Winey 1997, Table 1, 10 G1 nuclei):** 86 ± 16 (65-113), 11.7 per µm² of envelope. Build R has 93 at 10.8 per µm²: inside the range, so the G1 pore count stands. (Their nuclei are 1.76 µm³ from serial-section EM, against 2.4 µm³ here; shrinkage on fixation not known.) With 86 pores the nucleoporin anchor would be 1.47 instead of 1.36.
  - **Copies per pore (Kim 2018, Supp. Table 2A):** all 31 entries of `ingredients.NPC_COPIES` agree with the table (the table's 32nd row is the second fragment of Nup145; 552 copies per pore with both fragments).
  - **Ribosome concentration (Delarue 2018):** 23 µM is for asynchronous log-phase W303 cells, in cytosol segmented without organelles, and ribosomes fill about 20% of that volume (per the round-3 README; paper not re-read here). So the gap to 17 µM is not explained by cell-cycle stage or by the volume definition.
  - **Translocon-bound ribosomes (Wu 2026):** the 7,715 is for the authors' model of a 35.7 fL haploid G1 cell on galactose, and their total is about 356,000 ribosomes per cell (per the README). So build R's 5,507 in 34 µm³ is a real shortfall (0.71), and two cryo-ET sources now put ribosomes at about 1.4-1.7 times the model's 212,000, against the rRNA-based counts (187-200 thousand) the model sits on. This is the sharpest open conflict in the inventory.
  - **Per-gene ribosome loads (Arava 2003, Supp. Table 3; independent of the ribosome-profiling efficiencies the model uses):** 3,238 genes compared. Ribosomes per mRNA: Spearman 0.44 over all genes, 0.62 for genes with >= 3 and 0.77 for genes with >= 10 translated mRNAs in the model (few mRNAs per gene make the model's per-gene mean noisy). Density per 100 nt: 0.52; given length (residuals of log density on log length) 0.33. Medians 0.58 (model) against 0.47 per 100 nt. Two systematic differences: density falls less steeply with length in the model (log-log slope -0.45 against -0.74), and 97% of the model's mRNAs carry ribosomes while Arava finds a median of 72% of a gene's mRNA on ribosomes. The model has too few ribosome-free mRNAs (570 of 18,136). Not yet a script in `validation/`; computed ad hoc on `prepass.npz` of build R.
  - **Proteome against cell size (for a non-uniform scaling to the G1 cell):** the Lanz 2022 tables are for human cells (HLF, RPE-1), so they give no per-protein slopes for yeast. Neurohr 2019 has yeast proteomes of cdc28-4 cells growing far beyond normal size, which is the wrong regime for a normal G1 cell. A yeast dataset is still needed (added to `data_requests.md`).
- **Review of the technical stack (user's request):** two subagents, one auditing the build code, one surveying other methods; results in the next entry.

### 2026-11-13 (later): round S started; visual check of build R; first findings

- **Visual check of build R against build O (user's question; not done properly before the push):** all 23 guided views,
  whole-cell cuts at 5 planes x 3 axes, close-ups of endosomes, vesicles, cisternae (`stack_review_2026-11-13/visual_check_R/`).
  Vacuoles are worse in R (boxy, carved lobes; plan S9). Endosomes and Golgi unchanged, with an older defect: an empty
  looking layer around their membranes (plan S10). Cristae cross each other in the mitochondrion view of R (to check).
  The rest looks as in O. Rule from the user, saved to memory: every change is inspected visually before it is reported.
- **Not re-run on build R** (answer to the user): stress-condensate negative control (last run 11-11), polysome
  diffusion simulation (build Q boxes), crowding sweep of the diffusion model, copy ledger (plan S7).
- **Audit claims re-run by the main session on build R:** tile seams (volume fraction within 5 nm of a tile face 0.49
  in first-pass tiles, 0.20 in last-pass tiles, interior 0.40), membrane-face copies outside their shell (Golgi 17%,
  envelope 11%, lipid droplets 7%, vacuole 7%, ER 5%, PM 4%), per-species placed share in the matrix (Sod2 37%, Hsp60
  49%, Ssc1 96%), library dates, analytic against built volumes. All confirmed.
- **S1-S5 done** (commits d248f8e to e21f258 on `stack-cleanup`). Build R, new checks: placed 97.5% by copies and
  96.0% by volume; 76 species under 50% placed (335,000 requested copies), 238 under 80%.
- **Found by the first isolated build:** `work/interactions/pair_summary.parquet` (AlphaFold 3 pool pairs, read by
  `features.py`) is an input that lives in `work/` and has no generator in the repository.
- **Found by stage 0** (library rebuilt into a separate folder, compared entry by entry): 5,561 of 5,647 shared entries
  identical; 82 transmembrane proteins are rotated by 180 degrees in the library on disk, which predates the topology
  side rule now in `build_structs.process_af` (61 of them packed in build R, 42,800 copies: Kre2, Ssh1, Ato3, Erv29,
  Ctr1, ...); HTB1 and HTB2 have no entry, so 33,100 free H2B copies were booked as no_structure and never packed.
- **S6, baseline (done):** `main` rebuilt in an isolated run folder (frozen worktree, copied library, 56 instead of 64 workers, `PYTHONHASHSEED=0`) is identical to build R: all 5,912 tiles with the same molecules at the same positions (24,484,185 instances), all 128 pre-pass arrays and the expected totals equal. So the build is deterministic run to run and independent of the worker count and the hash seed (the audit had left this open). `run_all.sh` returned 1 after a complete build (a `[ ] && ...` as last command of its loop); fixed.
- **S9, vacuole (done, commit 5b330d5; sections old against new in `stack_review_2026-11-13/visual_check_R/vacuole_old_vs_new_sections.png`):** one near-spherical vacuole after the soft X-ray tomography of 490 unbudded G1 cells (Chen 2025, file S8: one object in 488 cells, corrected sphericity 0.98 median, ellipsoid elongation 1.10, flatness 1.07). 1 compartment, 16.0 um2 of membrane (was 12 compartments, 40.2 um2, sphericity 0.39), 1.8% carved by nucleus and cortex. Stage 2 runs with it (13 of 14 lipid droplets moved out of its way).
- **S13 (done, commit 8ddaf89):** the concentration rules of the inventory use nucleoplasm 1.65 and cytosol 18.44 um3 (built: 1.56, 18.42; before 1.90, 18.18).
- **Collision shapes, measured (S11):** union volume of the sphere trees against the union of the residue beads: 2.4 times in the matrix, 1.9 in the cytosol, 3.8 in the nucleolus (top 40 species each). Trees of 8 residues per sphere claim 0.74 / 0.87 / 0.59 of the volume of the library's 30-residue trees, at 3.7 times the spheres. Tile test (6 matrix-rich tiles of build R, no neighbours): matrix placed 87.9% with the library trees, 94.9% with 8-residue trees, at 7 times the CPU (5,770 against 828 s); 9 against 13 species under 80%. So finer shapes help but do not remove the biased dropping and are costly: not in the next build; thinning by a pilot (below) first.
- **Requested sphere volume fraction per region (build R):** cytosol 0.38, nucleoplasm 0.39 (99.6% placed); matrix 0.66 (75%), nucleolus 0.66 (80%), Golgi/endosome face 0.70 (84%), envelope face 0.51 (94%).
- **Written, under test (not yet in a build):** tile packing in 6 rounds by size class with the order of the 8 passes reversed every second round (S12; a 216-tile block test reproduces the seam with one round: 0.51 against 0.21 at the faces); shell cap (S10: the Golgi/endosome face keeps 0.55 of every species, 160,700 copies go to the cytosol; envelope face 0.75, 10,100 copies) and the shell test in the attraction path and the squeeze; uniform thinning of over-full compartments by a factor measured in a pilot packing of 10 tiles per region (S11).
- **Data fetched by a subagent** (14 papers; `pipeline/data/validation/literature/g1_cell/README_fetched_2026-11-13.md`); only five paywalled items remain in `data_requests.md`. Two findings for the copy-number round: (1) ribosomal proteins are 15.7% of the model's protein mass, against 29-32% in three proteomics datasets of fast-growing yeast (Metzl-Raz 2017, Björkeroth 2020, Xia 2022; none is BY4741 in YPD): PaxDb's integrated dataset under-weights them for this condition; (2) the 187,000 ribosomes per cell of von der Haar 2008 are not a count but total RNA from six older studies x 83% / rRNA mass; cryo-ET gives 14,000-14,500 per um3 of cytosol (Delarue 2018; Goetz 2025 preprint, which itself converts to 240,000 per cell). Lanz 2024 has per-protein slopes against cell size in yeast (histones about -1, ribosomal proteins about 0).
- **Mitochondria (user's question, 2026-11-13):** entirely procedural, and wrong in look: cristae are discs stacked every 42 nm along the whole tubule, nucleoid DNA is threaded between them. The segmented archive holds no real mitochondria (its networks are synthetic). A subagent is collecting measured ultrastructure for glucose-grown cells and any downloadable segmentation; then a rebuild of the mitochondria builder (new item S16).

### 2026-11-13 (afternoon): round S candidate builds; viewer bug; mitochondria; simulation sub-project started

- **Seed plumbing is neutral:** the plumbing-only commit at seed 0 rebuilt in isolation is identical to build R (5,912 tiles).
- **Candidate build 1** (`runs/S_seed0`; vacuole, library from stage 0, rounds, shell cap, thinning, inventory volumes) against build R, inspected in the viewer (all guided views, cuts, close-ups): sphere volume fraction within 5 nm of tile faces 0.38 (first pass) and 0.37 (last), interior 0.39 (R: 0.49, 0.20, 0.40); membrane-face copies outside their shell 0% everywhere (R: 4-17%); Golgi face 97.7% placed (84%); species under half placed 19 (76); placed 97.5% of the inventory (97.5%); overlaps 2.2e-5 (2.65e-5); containment 9.1e-5 (4.1e-5, to look at). Matrix 73.6% and nucleolus 77.3% of the inventory: the pilot (roomiest tiles) overstated what they hold; pilot tiles are now drawn in proportion to region volume.
- **Viewer bug found by the visual check:** the default region mask had 24 bits, so regions 24-26 (vacuole, Golgi/endosome and peroxisome faces, added later) were packed but not drawn: about 300,000 proteins in build R, live on the site. This, not the packing, was the empty layer around Golgi and endosome membranes. Fixed (`web/app/src/main.ts`), verified in a close-up; goes live with the next push.
- **Found by the build from an empty library:** `RVS161_RVS167_het` (amphiphysin) was a leftover of 2026-11-07 that no stage wrote any more; the builder now makes it (identical beads).
- **Collision shapes for packing (S11):** 15 residues per sphere (`wc3d.PACK_PER`), packing only: 93.1% of the matrix request in the test tiles (library trees 87.9%, 8 residues 94.9%), 3 species under 80% (13), 2.6 times the CPU.
- **Speed (S15):** ATP synthase step 39 s (11 min 49 s; geometry near each crista only; same 6,257 dimers); mRNA-end proteins 2 min (6 min); the serial part of stage 3 is now about 11 min (28).
- **Mitochondria rebuilt (S16; user's question; data collected by a subagent, `literature/mitochondria/`):** one connected curved network (was a dozen straight fragments), outer diameter 333 nm median (Egner 2002: 339), 14 nm between the membranes (Chang 2025: 13-16; was 8), 170 short cristae reaching 50-200 nm in from alternating sides with membrane 0.23 x the outer membrane (measured on the YPD meshes of EMPIAR-12534: 0.2-0.3; was a disc across the tubule every 42 nm, 1.4 x), 14 nucleoid sites in crista-free stretches (1.27 per um, Osman 2015), 15 mtDNA copies (0.45 per fL; was 35). Inner membrane 13.0 um2 (31.5), matrix 0.81 um3 (0.70). Plots in `stack_review_2026-11-13/visual_check_R/`.
- **Candidate build 2** (`runs/S2_seed0`: adds mitochondria, packing shapes, the speed changes): killed in its last round by a session restart, resumed from stage 3. Rule (saved to memory): long jobs are started detached (`setsid nohup`).
- **T1 started** (plan T): scene extraction of the deepest endocytic site of candidate 1 (1,498 actin subunits, 111 Arp2/3, 40 Sla2, 45 Rvs, 3,586 ribosomes); engine; calibration (tracers of 8-40 nm within 0.88-1.09 of the model's diffusion law; implicit part 1.32 x exp(-0.097 R)); first videos.
  - Release, version 1 (placeholders, imposed inward pull) replaced by version 2 after the parameter sheet: the patch as one diffusing body, disassembly from the Picco 2015 curves, nothing imposed; test quantity: the inward drift of the coat and actin centroids (measured about 25-30 nm/s). 12 runs of 6 s under way.
  - Transit: with 5 of the 14 cisternae as TGN acceptors, first contact after a median of 6.3 s; capture at 3 per second within reach gives 78% within 2 min (measured: 66%, Day 2018), so the measurement implies a docking and fusion rate of a few per second of contact, if nothing else holds the vesicle back.
  - Membrane: the solver reproduces the analytic tube (radius 1.01 (kappa/2P)^(1/3); force peak 1.4 and plateau 0.75 of 4 pi P R^2) and gives a tip radius of curvature of 11.4 nm at 1,000 kT and 0.6 MPa (measured 11.6 +- 2.6). With a bare bilayer (20 kT) below a coat of 2,000 nm2 and no Rvs sleeve the membrane pinches at 28 nm depth with 5,200 nm2 above the neck (rvs167 deletion: pinch at about 50 nm; wild-type vesicle 6,380 +- 1,929 nm2): the coat size and the rigidity of the uncoated membrane are being scanned.
- **Structural-model feedback from T1:** the deepest invagination is 92 nm with the coat centre about 70 nm below the membrane; centroid tracking puts the coat at 132-150 nm at scission (tomography: about 100 nm deep, up to 140). The tip of the model's tube is 40 nm wide; measured tip radius 11.6 nm (23 nm across).

### 2026-11-13 (evening): build S (candidate 3) finished, checked and chosen; endocytosis results; next rounds

- **Candidate build 3** (`runs/S3_seed0`, commit 9568da3, from stage 0; the guided views redone with 9408123): placed 99.1% of the inventory by copies and 98.1% by volume (build R 97.5 / 96.0%); 10 of 4,469 species under half placed (76), 60 under 80% (238); sphere volume fraction within 5 nm of tile faces 0.41 (first pass) and 0.40 (last), interior 0.41 (R: 0.49, 0.20, 0.40); membrane-face copies outside their shell 0% for all seven faces (R: 4-17%); localization top-1 on placed copies (test half) SWAT native 0.665, SWAT NOP1 0.742, hyperLOPIT 0.848 (R: 0.674, 0.746, 0.848; inside the CIs). Microtubules: 15 kinetochore, 2 other nuclear, 3 cytoplasmic (candidates 1-2 had lost the cytoplasmic ones: the 80 tries per microtubule could not find a free site among 20 nucleation rings; now 400). 169 cristae, 11 nucleoid sites, 93 pores, none within 75 nm of the nucleus-vacuole contact.
- **Changes after candidate 2, each inspected in the viewer:** nucleolus on the vacuole side of the nucleus (Barbosa 2015, Santana-Sosa 2023, Yang 1989), so the spindle pole body faces the cytoplasm; no pores at the junction; guided views (cell wall clear of the bud scar, cortical ER sheet off the scar, mitochondrion through the densest cristae, Cvt complex framed with the vacuolar membrane, vacuole off the junction, nucleolus on its axis with the envelope in view); lipid detail chosen by the camera's distance to the cut plane, not per tile (the "wedge" on Golgi cisternae was the viewer, packing densities were uniform); captions.
- **Bugs found on the way:** the nucleus and vacuole distance fields are clamped (+-130 and +-90 nm), so a "more than 200 nm from the nucleus" filter never matched; a test server left running from candidate 2 served the first "candidate 3" screenshots (now every server is started with `--directory` and checked); two validation figures failed in isolated builds (they read round-H and calibration files that only the bundle's work folder has; now optional); the condensate test's output files were missing from the run folder after the build although its log was complete (re-run by hand, same numbers; cause not found).
- **User's decisions (2026-11-13, evening):** merge build S and push it; squash the Pages history with a force push; next model round on copy numbers (plan U); the sub-project continues at my discretion, with the modelling strategies made to work together and the membrane-bending proteins (BAR domains, amphipathic helices, ...) included (plan T3).
- **Sub-project results so far (T1; `subprojects/endocytosis_results/`):** release: thermal motion alone moves the vesicle 4-39 nm inward in 6 s (measured 25-104 nm); with the measured actin flow imposed the patch follows it through the crowd (the match is built in; what the run shows is that crowding does not stop it), and diffuses freely after about 5.8 s. Transit: not limited by diffusion (first TGN contact 6.3 s; about 2 captures per second of contact give the measured 66% in 2 min). Budding: the same coat (30 nm preferred radius, 160 kT) makes a full bud with a 5 nm neck on an organelle membrane at 0.02 pN/nm, a cap 37 nm high at 0.8 pN/nm (threshold kappa C0^2/8 = 0.4 pN/nm, Hassinger 2017), and leaves the plasma membrane flat under 0.35 MPa turgor. Scission without the Rvs sleeve: pinch at 28-44 nm (rvs mutants about 50 nm). With the sleeve over the whole tube: the tube closes near its base at 44 nm and the vesicle is 8,830 nm2 (measured 6,380 +- 1,929): not solved, no video presented.

## U. Copy numbers (plan of 2026-11-13, evening)

Asked for by the user ("For the next model round, work on copy-number problems"). Branch `copy-numbers` from `main` after build S.
- U1. Ribosomes: ribosomal proteins are 15.7% of the model's protein mass against 29-32% by proteomics of fast-growing yeast; the 187,000 ribosomes per cell of von der Haar 2008 are derived from total RNA; cryo-ET gives 14,000-14,500 per um3 of cytosol. Decide the number from the measurements that are counts (cryo-ET, rRNA per cell) and make protein mass, ribosomes and free ribosomal proteins consistent with it.
- U2. Scaling from the average cell to the G1 cell per protein (Lanz 2024 slopes: histones about -1, ribosomal proteins about 0), not one factor for all.
- U3. Ribosome-free mRNAs (too few in build S) and the polysome share.
- U4. A copy ledger per species: requested, placed, dropped, thinned, moved by the shell cap; and the inventory's sources per species.
- Every change tested against held-out data that did not set it (the condensate negative control and the localization benchmarks stay as they are).

**Progress, 2026-11-13 (night), branch `copy-numbers` (worktree `dev_copy_numbers/`), test build `runs/U1_seed0`:**
- **What the data say about ribosomes** (`scratchpad` analysis; numbers recomputed from the files on disk): ribosomal proteins are 15.3% of protein mass in PaxDb (the model's input), 15-19% in the mass-spectrometry studies inside Ho 2018 (Lu, Kulak, Lawless), 2-7% in the GFP/TAP studies (which undercount them: Ho's consensus is 7%), 28-32% in Metzl-Raz 2017 and Björkeroth 2020. Counts: 14,000 per um3 of cytoplasm (Delarue 2018, cryo-ET, W303 SCD) and 14,517 +- 3,193 per um3 of cytosol with organelles subtracted (Goetz 2025 preprint); about 220,000 in a G1 cell by serial-section EM (Yamaguchi 2011, abstract only). Build S has 211,550 mature ribosomes, 11,500 per um3 of cytosol (0.8 x cryo-ET). At 3.5 pg of protein, 30% ribosomal protein would be about 460,000 ribosomes, which no count supports, so the proteomic share is not followed; the conflict goes into MODEL.md.
- **U1 (commit 0895e24):** ribosomal proteins scaled together to 14,250 per um3 of the built cytosol: 262,700 r-protein sets (221,500 from PaxDb), other proteins unchanged; protein total 3.60 pg; ribosomal proteins 18.0% of protein mass. The GEM diffusion benchmark (Delarue's GEMs slow with ribosome concentration) is the held-out check.
- **U3 (commits 664b279, f9d431e):** build S had 3% free mRNAs and a mean occupancy of 0.98 (Arava 2003: 0.71; 14% of mRNAs with one ribosome, build S 7%). Now every mRNA is untranslated with one probability (0.27) so that the gene-mean occupancy is 0.71, and the ribosome spacing is the measured mean density (one per 156 nt of ORF, Arava 2003) instead of being fitted to the assumed 70% of ribosomes in polysomes. Consequence, before any build: only about 36% of cytosolic ribosomes are then on mRNAs, with 18,100 mRNAs per cell (Pelechano 2010's 25,900 times the protein share 0.70). Miura 2008 counts 36,100 per cell in YPD, and mRNA concentration is reported to be constant with cell size (which would scale by volume, 0.81, not by protein). Data on the measured polysome share and mRNA per cell are being collected before the mRNA total is changed. New test `validation/mrna_occupancy.py` (per-gene occupancy and density against Arava, ribosomes per mRNA against Riba 2019); build S: occupancy Spearman -0.08 (607 genes), density 0.29, Riba 0.68.
- **U2 (commit d3d77ff):** per-protein size scaling to the G1 cell with the Lanz 2024 slopes (G1 arrest, glucose; 3,853 proteins whose three experiments agree, others the median -0.05), weight (34/42)^slope before the total is set. 98% of proteins move by 0.86-1.14 x; H2B (HTB2) 1.22 x, Pma1 1.10 x.

## Literature candidates (2026-11-08)

A literature survey of what is known about yeast cell biology that the model could add, set against the build with pre-ribosomes and rDNA. "Data" means a usable dataset or structure exists; "consensus" means textbook knowledge without one dataset. Numbers marked (unchecked) were recalled, not confirmed against the paper. Suggested order: L2, L3, L5, L6 first (data-driven, cheap); then L7, L8, L12, L9, L10 (visible structures); then L13, L14 (nucleus); L1 last, as its own project. Already planned: IDRs, log-phase P-bodies.

**Density**
- **L1. Proteome total and ribosome count (high impact, high effort). Done 2026-11-10 (5 pg, squeeze relaxation; see Progress).** Now 4.2e7 proteins = 2.7 pg; 162,600 ribosomes in 26.6 µm³ of cytosol = 10 µM. BioNumbers: 4–6 pg protein per haploid cell (BNID 110550, 105082). In-cell cryo-ET (Delarue *et al.* 2018 Cell): ribosomes about 23 µM, about 20% of cytosolic volume (the 23 µM is from a commentary). So the model is about 2× too dilute, which is most of the gap to the ~200 mg/mL dry mass. Needs a packer that goes beyond random sequential addition (compression or relaxation), not just new numbers.

**Data that replaces assumptions**
- **L2. Measured localization fractions.** Now rank weights 1, ½, ¼ from annotations; 625 proteins default to cytosol. CYCLoPs (Koh *et al.* 2015 G3; 16 compartments, >20 M cells, cyclops.ccbr.utoronto.ca); yeast hyperLOPIT (Nightingale *et al.* 2019, unchecked) would also split Golgi from endosomes. Low to medium effort.
- **L3. Mitochondrial proteome.** Now fixed 75/25 inner/outer membrane and 90/10 matrix/IMS splits; matrix 62% filled. Morgenstern *et al.* 2017 Cell Rep. 19:2836: 901 high-confidence proteins with copy numbers and submitochondrial localization. Also tests whether the matrix shortfall is volume or over-assignment. Low effort.
- **L4. mRNA identities and ribosome loading.** Now lognormal lengths, polysomes 1 + Poisson(4), 70% in polysomes assumed. Per-gene mRNA abundances and lengths, and ribosome density per mRNA (Arava *et al.* 2003 PNAS: well below maximal packing, falls with ORF length, most mRNAs ribosome-bound). Each polysome becomes a real gene. Medium effort.
- **L5. Lipid composition per organelle.** Now one whole-cell composition (Ejsing 2009). Zinser *et al.* 1991 (PM rich in ergosterol and PS), cardiolipin in the inner mitochondrial membrane, PS/PE in the PM inner leaflet (consensus), lipid droplets as a monolayer around a neutral-lipid core. Low effort (procedural).

**Missing structures with data**
- **L6. Mitoribosomes on the inner membrane.** Pfeffer *et al.* 2015 Nat Commun 6:6019: the vast majority of yeast mitoribosomes are tethered to the inner membrane, in polysome-like clusters, through Mba1 and the 21S rRNA expansion segment 96-ES1. Now loose in the matrix, few fit. Better accuracy and relieves the matrix. Medium effort; best value.
- **L7. Eisosomes (MCC).** 30–50 furrows per cell, about 50 nm deep, 200–300 nm long (Strádalová 2009 and reviews). Near-native Pil1/Lsp1 lattice with bound PI(4,5)P2, PS and sterol: Kefauver *et al.* 2024 Nature 632:664, PDB 8QB7–8QBG. Now the PM is flat. Medium effort.
- **L8. Membrane invaginations at endocytic sites.** Now 22 actin patches on flat membrane. Kukulski *et al.* 2012 Cell: 211 intermediates by correlative ET, timed; clathrin arrives on flat membrane, invagination starts with actin, amphiphysin on the tube. Picco *et al.* 2015 eLife: molecules per site for coat and actin proteins (unchecked). Medium effort.
- **L9. Golgi as scattered cisternae.** Now Golgi/endosomes as 1% primitives. Not stacked; about 6–10 early and as many late/TGN cisternae (Papanikou & Glick 2009); mostly flat disks in tomograms. Then endosomes and vesicles. Medium effort.
- **L10. Cortical ER as fenestrated sheets and tubules.** Now continuous sheets over 48% of the PM, no cytoplasmic ER. West *et al.* 2011 JCB 193:333: tubules and fenestrated cisternae at the PM, plus cytoplasmic cisternae linking to the nuclear envelope. Medium effort.
- **L11. Nuclear pore basket and the missing FG-Nups.** Akey *et al.* 2022 Cell, 2023 Mol Cell 83:3283; Mlp1:Mlp2:Pml39 = 4:2:1 per spoke (J Cell Sci). Medium effort.
- **L12. Bud and birth scars.** Every unbudded cell has a birth scar; a mother also has one chitin bud-scar ring per division (consensus). Low effort; on the cell exterior.

**Nucleus**
- **L13. Chromosome layout from Hi-C.** Now walks start at the kinetochore microtubules but telomeres are free. Rabl configuration with telomeres at the envelope in about 6–8 clusters (Gotta 1996, unchecked); Duan *et al.* 2010 Nature 465:363 3D genome model and later Hi-C. Includes Rabl-connected rDNA on chromosome XII. Medium effort.
- **L14. Transcribing Pol II and nucleosome-free promoters.** Now all ~30,000 Pol II free (Borggrefe 2001, BNID 104516) and nucleosomes evenly spaced. Only a few % to ~15% of Pol II is engaged (run-on, PMC2982843; Struhl 2007); promoters are nucleosome-depleted (genome maps: Brogaard 2012). Engaged Pol II with nascent mRNA at per-gene rates. Medium to high effort.

**Smaller additions**
- **L15. RAC, Ssb and NAC on ribosomes**, with nascent chains at the tunnel exit. Native *S. cerevisiae* 80S–RAC (Nat Commun 2022, 10.1038/s41467-022-31127-4); 7Z3N/7Z3O (*C. thermophilum*). Medium.
- **L16. mtDNA nucleoids.** 10–20 Abf2-packed nucleoids per haploid cell, 25–50+ copies of the ~86 kb genome. Medium.
- **L17. Vacuole contents.** Polyphosphate (VTC polymerase 7YTJ, 8I6V; in-cell ³¹P-NMR quantification, PMC9911596); nucleus–vacuole junction (Nvj1–Vac8). Low to medium; polyP could be procedural.
- **L18. Cvt complex.** prApe1 dodecamers (4R8F) with Atg19 form a large cargo complex in growing cells. Low.
- **L19. Contact sites and seipin.** ERMES, NVJ, ER–PM tethers (tricalbins, Ist2); seipin Sei1/Ldb16 (7OXP, Klug *et al.* 2021 Nat Commun 12:5892) at ER–lipid droplet necks. Medium.
- **L20. Peroxisome count check.** 1–2 (to several) very small peroxisomes in glucose (PMC3245845), against 0.2% of volume here. Low.

## Where the model stands (before v3)

| Area | Current state | Gap |
|---|---|---|
| Assemblies | 42 PDB types. Everything else (5,441 types, 49% of mass) is a single AlphaFold chain. | Most complexes and oligomers are placed as separate monomers. |
| Homo-oligomers | 23 PDB homo-oligomers (from a top-40 abundance cut) | UniProt says 360 more placed types are oligomers (9.1% of mass). |
| Complex Portal | 18 of 623 multi-subunit complexes are modeled. Selection only tried the top 40 and skipped every membrane complex. | Of the 91 largest unmodeled complexes, 41 have an RCSB entry covering ≥ 80% of subunit mass, 19 of them membrane complexes. |
| Membrane proteins | 1,229 TM genes. Only Pma1, V-ATPase and NPC are assemblies. | 1,204 of 1,229 are placed as monomers. 61.6 k TM copies annotated as cytosolic are dropped. |
| Cytoskeleton | ACT1 (111 k), tubulin, septins, Arp2/3, cofilin, profilin and tropomyosin are scattered as free monomers. | No filaments, patches, cables, spindle or SPB. |
| RNA | rRNA (ribosome), tRNA, schematic mRNA, nucleosomal DNA. | No polysomes, and RNA is stripped from all automatically added complexes (snoRNPs, UTP-A). |
| Disorder | Residues with pLDDT < 50 are trimmed. | 3.1% of mass removed overall, but 16% in the nucleoplasm and 15% in the nucleolus. |

## 0. Prerequisites (before any model change)

- **0.1 Make the pipeline re-runnable.** Replace the hard-coded `/workspace/yc*` paths with a `config.py`, and port the notebook cells to scripts (handover §2, §3.2). Every item below triggers a re-run of about 10 min.
- **0.2 Make the viewer tolerate new type indices.** Type indices change on every repack. `SPECIAL` in `data.ts` already looks types up by id string; check that nothing else hard-codes them.
- **0.3 Use local inputs.** These are already on this machine:
  - AFDB proteome: `~/Downloads/UP000002311_559292_YEAST_v6.tar`
  - IntAct: `~/Downloads/IntActYeast.txt`
  - Alliance paralog pairs: `~/Downloads/alliance_yeast_paralog_pairs.csv`
  - AF dimers: `~/Projects/bioemu_docker/afdb_yeast_dimers`

## A. Assemblies: use the structures that exist (highest impact per effort)

### A1. Homo-oligomers from the AF dimer set
- **What it covers.** 604 of the 648 homodimer models match placed AF monomer types: 4.7 M copies, 11.7% of mass. Examples are FBA1 (809 k copies), PGI1, TAL1, URA1 and ALD6.
- **Annotation conflicts.** For 425 of the 604, UniProt says monomer or unknown. For 29, UniProt says the protein forms a larger oligomer than a dimer.
- **Work:**
  - Filter by confidence (ipTM, ipSAE, pDockQ2). The set runs from ipTM 0.37 to 0.97, so pick thresholds by checking a sample against PDB.
  - Prefer a PDB homo-assembly where one exists. Extend `rcsb_homo_assembly` from the top 40 to all 360 annotated oligomers.
  - For tetramers and higher, use the PDB assembly or predict the full oligomer. The dimer alone is not enough.
- **Pipeline changes.** Add the dimer as a library type. Halve the copies in `expected_totals` via `oligo_map`, as the PDB oligomers already do.

### A2. Large complexes from RCSB (answers "are we using cryo-EM fully?" — no)

`extend_complexes.py` only tried the 40 largest soluble Complex Portal complexes. It required ≥ 80% coverage in a single assembly and dropped every nucleic-acid chain. The current RCSB check (`analysis/rcsb_coverage.csv`) finds these ≥ 80%-coverage entries for complexes that are not yet modeled:

| Group | Complexes (best entry) |
|---|---|
| Soluble, abundant | 6-phosphofructokinase 3O8O (850 kDa, 12.7 k copies), eIF3 8CAH, eIF2 6FYY, RAC 7X3K, eRF1–eRF3 4CRN, CDC48–SHP1 8UB4, CDC48–NPL4–UFD1 8DAR, FACT 7NKY, PeBoW 8V87, RQC 9OFV, Box H/ACA snoRNP 9G25, IDH 3BLX, RNR 2CVY |
| Mitochondrial ribosome | 54S 8UXA, 37S 9HM0 |
| Membrane | ATP synthase 8F39, complex III 8YIO, complex IV 9ETZ / 6T0B, complex II 9QDM, TOM 9J99, SEC complex 7AFT, Sec61 7KAJ, Ssh1 2WW9, OST 8AGE / 7OCI, EMC 6WB9, β-1,3-glucan synthase FKS1 9PE1 / FKS2 9LU8, PMT1–PMT2 6P2R, VTC 7YTJ, SPOTS 8C81, SNARE 3B5N |

- Run the selection over **all** 623 complexes, not the top 40.
- **Partial structures (50–80% coverage)**, such as SAGA 6T9I, INO80 9C9Z, TFIID 6T9K and CCR4–NOT: fill the missing subunits from AF monomers or AF dimers docked at known interfaces, or place the core and leave the rest as monomers.
- **Keep RNA chains** that belong to the complex: snoRNA in Box C/D and H/ACA, rRNA in mito ribosomes and pre-ribosomes, RNase P/MRP, SRP 7SL. Right now `extend_complexes` drops all nucleic acids.
- **Ribosome.** 4V88 is a 2011 crystal structure. Consider a cryo-EM translating 80S that includes tRNAs and eEF1A/eEF2, and a separate ER-bound form (see C1).
- **Overlapping complexes.** Complex Portal lists variants (eEF1, RSC, CURI) that share subunits. Draw subunits from the pool by stoichiometry so total copies are conserved. The current code already subtracts consumed subunits; it just needs to handle many more complexes.

### A3. Heterodimers from the AF set, with filtering
- **What it covers.** For 1,050 of the 1,250 models, both partners are placed AF monomers, up to about 11.6% of mass in total (pairs overlap).
- **Filters needed.** The top pairs by abundance show why:
  - **Different compartments:** KAR2 (ER lumen) with SSE1 (cytosol) is impossible. Require that both partners share a localization.
  - **Paralog pairs:** SSA1–SSA2 and HXK1–HXK2 are paralogs with high ipTM, but they are probably not real heterocomplexes. Some paralog pairs are real (PFK1–PFK2, FRS1–FRS2, BMH1–BMH2), so flag paralog pairs and require supporting evidence rather than dropping them all. Sources: IntAct, Complex Portal, PDB.
  - **Hubs:** SKP1 (17 partners), CDC28 and PHO85 (15 each) and the LSm proteins (12) bind partners in mutually exclusive ways. Allocate partner pools by copy number, for example greedily by confidence or as a small linear program, instead of pairing every model.
  - **Overlap with modeled complexes:** 214 pairs touch a protein already in a PDB complex. Skip them, or merge them where the dimer extends a complex.
- **Interactions that are real and visible:** ACT1–PFY1 (ipTM 0.90, profilin–actin), ACT1–COF1, BMH1–BMH2, SAM1–SAM2, ASN1–ASN2, ZUO1–SSZ1 (part of RAC).

## B. Membranes

### B1. Automatic membrane orientation for multi-chain assemblies
This is the prerequisite for the 19 membrane complexes in A2 and for membrane heterodimers.
- Today only NPC, Pma1 and V-ATPase are oriented, each by hand-written code in `build_structs.process_complex`.
- Generalise it: map UniProt TRANSMEM ranges onto chains, fit the slab normal, then set the sign from the TOPO_DOM "Cytoplasmic" annotation. This is the method `process_af` already uses for monomers. Alternatively, take pre-oriented coordinates from OPM or PDBTM.

### B2. Sub-membrane localization
- **Mitochondria.** TM proteins go to the inner and outer membranes by a fixed 75/25 split, and soluble ones to matrix and IMS by 90/10 (`wc3d.expected_totals`). UniProt text distinguishes "inner membrane", "outer membrane", "intermembrane space" and "matrix", but the `RULES` regex collapses them all into "mitochondri". Parse them separately, and use a mitochondrial sub-compartment proteome as a fallback.
- **Other organelles.** Make the same split for ER vs NE, Golgi vs endosome, and the lipid-droplet monolayer (hairpin proteins). Today LD TM proteins are sent to the ER.
- **Dropped proteins.** 61.6 k TM copies are annotated "cytosol" and dropped. Re-localize them from a GFP atlas (see E2).

### B3. Membrane-specific lipids
- One composition (whole-cell, Ejsing 2009) is used for every membrane (`worker.ts:217`).
- Use organelle compositions: cardiolipin in the mitochondrial inner membrane; ergosterol and sphingolipid enrichment in the PM; leaflet asymmetry (PS and PE on the inner leaflet of the PM).

### B4. Membrane supra-structures
- **Eisosomes.** Pil1 (40 k) and Lsp1 (36 k) are BAR-domain proteins that form lattice-coated furrows in the PM. Right now they are scattered monomers. Build furrows with Sur7 and Nce102 in the MCC domains.
- **Mitochondrial cristae.** ATP synthase dimer rows at the cristae rims; respiratory supercomplexes (III₂IV₂) in the flat cristae membrane.
- **ER-bound ribosomes** on Sec61 or the SEC complex (see C1).

## C. Nucleic acids, translation, chromatin

### C1. Polysomes and ribosome state
- **Polysomes.** mRNA threads are packed after the proteins, through leftover cytosol (`pack3d.py:303`), so ribosomes and mRNA are never connected. Place ribosomes along mRNA paths at a polysome spacing, with the mRNA entering and leaving each ribosome's channel.
- **Ribosome count.** 167 k placed against an estimated ~3 × 10⁵ for YPD (handover §9). Reconcile the count. The only ER fraction is a hard-coded 3% (`ingredients.COMPLEXES`); find a literature value and dock those ribosomes on Sec61.
- **Ribosome-associated factors.** NAC (66 k copies), RAC/Ssb, eEF1A (TEF1/TEF2, among the most abundant proteins) and eIF3 should sit on ribosomes, not be scattered.

### C2. Chromatin organisation
- **Rabl configuration:** centromeres clustered at the SPB, telomeres at the NE, arms reaching across the nucleus. Add 16 chromosome territories instead of anonymous walks.
- ~~**rDNA array** in the nucleolus, with Pol I (7.7 k copies) on active repeats.~~ Done 2026-11-07 (active genes and nucleosomal repeats; not connected into one array).
- **Gene-level features:** nucleosome-depleted regions at promoters, and Pol II on transcribed genes.

### C3. Other RNPs and the nucleolus
- ~~**Pre-ribosomes.** 90S, pre-60S and pre-40S particles make up much of the nucleolar mass.~~ Done 2026-11-07: five stages with their RNA (see Progress).
- **Lower-copy RNPs:** snRNPs and spliceosome, telomerase, RNase P/MRP, SRP.
- **Mitochondrial DNA nucleoids:** about 20–50 copies of the 86 kb mtDNA, packaged by Abf2, with mito ribosomes at the inner membrane.

## D. Cytoskeleton (currently all free monomers)

| Structure | Make it from |
|---|---|
| Actin patches (endocytic sites) | Branched networks: F-actin plus Arp2/3, capping protein CAP1/2, Abp1, Crn1. Place at the PM, polarized if the cell is budded. |
| Actin cables | Formin-nucleated (BNI1/BNR1) bundles with tropomyosin TPM1/2 and fimbrin SAC6, along the cortex |
| G-actin pool | Profilin–actin and cofilin–actin pairs (both in the AF heterodimer set). Split F/G actin using a literature ratio. |
| Microtubules | SPB embedded in the NE (Spc42 layer, Spc110/γ-TuRC), a few cytoplasmic MTs, and short intranuclear MTs or a spindle depending on cell-cycle stage. Build 13-protofilament lattices from the tubulin dimer (5W3J is already in the library). |
| Septins | CDC3/10/11/12/SHS1 filaments at the bud neck. Requires a bud (E3). |

Filaments need a "polymer instance" type in both the packer and the viewer. Placed as rigid long instances they would break the sphere-tree assumptions. The procedural-polymer path already used for mRNA and glucan is the natural template.

## E. Inventory and geometry

- **E1. Copy numbers.** PaxDb whole-organism integrated (about 2-fold uncertainty), rescaled to the Ho et al. 2018 total. Use the Ho 2018 per-protein unified abundances directly, which are closer to log phase in rich medium. Re-check the total protein density (handover §9: 1.0 × 10⁶/µm³ vs a published 1–3 × 10⁶).
- **E2. Localization source.** 625 proteins default to cytosol. Add a genome-wide GFP localization table (CYCLoPs, LoQAtE or the full YeastGFP set) as the primary source, with UniProt and GO as fallbacks.
- **E3. Geometry.**
  - **Nucleus:** 7% of cell volume vs 10.5% by EM (handover §9).
  - **Cell-cycle state:** decide whether this is an unbudded G1 cell or a budded one. That decision drives septins, actin polarity and spindle state.
  - **Mitochondria:** consider a network geometry taken from tomography or soft X-ray data.
- **E4. Packing saturation.** The mitochondrial matrix is 55% filled and the IMS 42%. Options: small-first backfill, multi-pass RSA, or a denser packer. Assemblies (A1–A3) reduce the instance count, which may help.

## F. Disordered regions and condensates

- **F1. Represent trimmed IDRs** instead of deleting them (about 3.1% of mass; 16% in the nucleoplasm). Two options:
  - excluded-volume random-coil beads tethered to the folded domain, generated procedurally like mRNA;
  - sampled conformations, for example with BioEmu.
- **F2. NPC FG-repeat mesh.** Nsp1, Nup100/116/145 and others fill the pore channel with FG repeats. This is the IDR case with the highest visual and functional impact.
- **F3. Nucleolar and P-body condensates.** Higher local density where the data support it.

## G. Validation to add alongside the changes

- Copy conservation per gene after each assembly step (subunits consumed ≤ pool).
- Fraction of mass in assemblies, per compartment, before and after.
- Interface sanity: no clash inside a predicted dimer. Spot-check confidence thresholds against PDB complexes.
- Compare with cryo-ET: ribosome density and nearest-neighbour distances, ER-bound ribosome fraction, and ATP synthase spacing on cristae.

## Suggested order

| Step | Items | Why now | Effort |
|---|---|---|---|
| 1 | 0.1, 0.2 | Unblocks every re-run | 0.5–1 day |
| 2 | A1, A2 (soluble), keep RNA | Biggest gain in accuracy for the effort, using data already in hand | 1–2 days |
| 3 | B1, then A2 (membrane), B2 | Membrane proteome as assemblies, in the right membranes | 2–3 days |
| 4 | A3 | Interaction layer with filtering and pool allocation | 2 days |
| 5 | C1 | Polysomes and ER-bound ribosomes; a visible change and a reconciled ribosome count | 2–3 days |
| 6 | E1, E2, E3 | Inventory and geometry rescale, then a full re-run | 1–2 days |
| 7 | D, C2, B4 | Cytoskeleton, chromatin and membrane supra-structures; needs new polymer machinery | 1–2 weeks |
| 8 | F, C3, B3 | Disorder, RNPs, lipids | open-ended |

Effort figures are rough guesses for one person working with an agent.
