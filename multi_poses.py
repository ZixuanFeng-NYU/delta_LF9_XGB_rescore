import multiprocessing
import os

os.mkdir('poses')
def process_files(file_name):
    for j in range(1, 6):
        output_file = f"folder_combined/poses/{file_name}_{j}.pdb"
        if not os.path.isfile(output_file):
            input_file = f"folder_combined/docked_ligand_pdbqt/{file_name}/{file_name}_docked_ligand_rank{j}.pdbqt"
            command = f"/scratch/zf2012/zf_research/Lin_F9_test-master/smina.static -r protein_pdbqt/ -l {input_file} --local_only --scoring Lin_F9 -o {output_file}"
            os.system(command)

if __name__ == '__main__':
    file_names = os.listdir("folder_combined/docked_ligand_pdbqt")
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    pool.map(process_files, file_names)
    pool.close()
    pool.join()
