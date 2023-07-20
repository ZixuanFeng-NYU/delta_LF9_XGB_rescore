import os
import pathlib
import sys
import pandas as pd
import fileinput
from featureSASA_zf_edit import sasa

def cal_SASA(tup):#out,fn,lig,pro,datadir):
    # fn: pdbid
    # lig: ligand part should be path of pdbid_posesrank.pdb
    # pro: protein part
    pro = tup[0]
    lig = pathlib.Path(tup[1])
    fn = lig.name.partition('_')[0]
    #pro = os.path.join(datadir,pro)
    #lig = os.path.join(datadir,lig)
    outpath = pathlib.Path(f"sasa/{lig.name.replace('.pdb','')}.csv")
    datadir = "."
    sasa_features = sasa(datadir,str(pro),str(lig))
    sasa_com = sasa_features.sasa
    sasa_pro = sasa_features.sasa_pro
    sasa_lig = sasa_features.sasa_lig
    with open(outpath,"w") as out:
        out.write(fn + "," +  ",".join([str(round(i,2) )for i in sasa_com]) + "," + ",".join([str(round(i,2) )for i in sasa_lig]) + "," + ",".join([str(round(i,2) )for i in sasa_pro]))


if __name__ == "__main__":
    # Pairs of protein.ligand structures at Path
    protein_path = pathlib.Path("/scratch/zf2012/zf_research/Research_part5_Spring2023/05-12-2023_PCBA_DUDE/pubchem_HS90A/diffdock")
    lig_path = pathlib.Path("/scratch/zf2012/zf_research/Research_part5_Spring2023/05-12-2023_PCBA_DUDE/pubchem_HS90A/diffdock/poses/")
    prot_files = os.listdir(protein_path)
    collection_of_pairs = []
    for lig_file in os.listdir(lig_path):
        print(lig_file)
        lig_file = pathlib.Path(lig_file)
        pdb_id, sep, tail = lig_file.name.partition("_")
        prot_name = "AF-P07900-F1-model_v4.pdb"
        collection_of_pairs.append((pathlib.Path(protein_path/prot_name), pathlib.Path(os.path.join(lig_path,lig_file))))
    print("entering multi")

    for collection in collection_of_pairs:
        cal_SASA(collection)
