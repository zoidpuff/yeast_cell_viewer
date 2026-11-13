"""Stage 0: the base structure library (pipeline/lib, or $YC_LIB): the entries every later stage assumes to exist.

Until 2026-11-13 these were made only by cells 3, 4 and 8 of the historical notebook, so the library on disk was state
that no script could rebuild and a change to the sphere-tree recipe or the pLDDT cut-off never reached it.

  1. curated complexes (ingredients.COMPLEXES) from their RCSB assembly, membrane complexes oriented (build_structs)
  2. one AlphaFold DB model per protein of the inventory (AFDB_TAR), residues with pLDDT < 50 left out,
     transmembrane proteins oriented from the UniProt topology
  3. the curated PDB homo-oligomers of inventory v2 (LEGACY_OLIGO_PDB; paralogs reuse the entry of their paralog)
  4. tRNA-Phe (1EHZ), the reference for the tRNA species of trna_species.py
Everything else in the library is written by stages 1 and 3 (Complex Portal assemblies, heterodimers, homo-oligomers,
barrel monomers, the model() of each builder, tRNA species).

  python scripts/s0_library.py [--only complexes,af,oligo,trna]
Needs the PDB assemblies in $YC_PDB (downloaded from RCSB when missing) and the AlphaFold DB proteome tar.
"""
from _common import *
import gemmi
import inventory as inv, ingredients as ingm, build_structs as bs, structlib as sl

# gene -> RCSB entry of its homo-oligomer (chosen in inventory v2 by structlib.rcsb_homo_assembly for the most abundant
# soluble oligomeric enzymes; TDH1/2, ENO2, ADH2 and PDC5 use the structure of their paralog). Read back from the
# library entries of 2026-11-03, which record their entry
LEGACY_OLIGO_PDB = dict(TDH3="3PYM", PDC1="1QPB", CDC19="1A3W", ADH1="4W6Z", TPI1="1NEY", GPM1="1QHF", ENO1="2AL1", AHP1="4H86",
                        SOD1="1F1G", HSP82="2WER", IPP1="1E9G", GND1="2P4Q", TRR1="3ITJ", HOM6="1EBF", YDJ1="1XAO", CYS3="1N8P",
                        TKL1="1GPU", KRS1="7F6W", TDH1="3PYM", TDH2="3PYM", ENO2="2AL1", ADH2="4W6Z", PDC5="1QPB")


def fetch(pdb):
    path = f"{C.PDB}/{pdb}-assembly1.cif.gz"
    if not os.path.exists(path):
        import requests
        r = requests.get(f"https://files.rcsb.org/download/{pdb}-assembly1.cif.gz", timeout=300); r.raise_for_status()
        open(path, "wb").write(r.content)
    return path


def main():
    only = set(sys.argv[sys.argv.index("--only") + 1].split(",")) if "--only" in sys.argv else {"complexes", "af", "oligo", "trna"}
    os.makedirs(C.LIB, exist_ok=True)
    df, geo = inv.build()
    ing, _ = ingm.build_ingredients(df)
    topo = pd.read_csv(f"{C.DATA}/uniprot_topology.tsv", sep="\t"); topo.columns = ["acc", "topo", "tm", "intramem"]
    topo = topo.set_index("acc")
    if "complexes" in only:
        pma1_tm = bs.parse_ranges(topo.loc[df.set_index("gene").loc["PMA1", "uniprot"], "tm"], "TRANSMEM")
        for r in ing[ing.kind == "complex"].itertuples():
            fetch(r.pdb)
            n, ns, rmax = bs.process_complex(r.ing_id, r.pdb, pma1_tm)
            log(f"complex {r.ing_id:16s} {r.pdb}: {n} beads, {ns} spheres, radius {rmax:.1f} nm")
    if "af" in only:
        res, missing = bs.process_all_af(ing, topo)
        err = {k: v[1] for k, v in res.items() if v[0] == "ERR"}
        log(f"AlphaFold models: {len(res) - len(err)} written, {len(err)} failed, {len(missing)} proteins without a model in the tar")
        for k, v in list(err.items())[:20]:
            log("  failed", k, v)
        json.dump(dict(written=len(res) - len(err), failed=err, missing=sorted(missing)), open(f"{C.RESULTS}/library_af.json", "w"), indent=0)
    if "oligo" in only:
        i1 = ing.set_index("ing_id")
        for g, pdb in LEGACY_OLIGO_PDB.items():
            if g not in i1.index:
                continue
            fetch(pdb)
            n, ns, rmax = bs.process_homo_pdb(g + "_oligo", pdb, i1.loc[g, "uniprot"])
            log(f"homo-oligomer {g + '_oligo':12s} {pdb}: {n} beads, {ns} spheres")
    if "trna" in only:
        st = gemmi.read_structure(fetch("1EHZ")); st.remove_ligands_and_waters()
        xyz, knd, _ = sl.beads_from_structure(st); xyz = xyz - xyz.mean(0)
        n, ns, rmax = bs.save("tRNA", xyz, knd, per=15)
        log(f"tRNA (1EHZ): {n} beads, {ns} spheres")
    log("library:", len([f for f in os.listdir(C.LIB) if f.endswith(".npz")]), "entries in", C.LIB)


if __name__ == "__main__":
    main()
