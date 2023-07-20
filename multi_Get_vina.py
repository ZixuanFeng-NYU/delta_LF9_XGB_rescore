import multiprocessing
import os
from functools import partial
import subprocess

protein_path = 'protein_pdbqt'

def process_file(i, rank, protein_path):
    input_file = f'docked_ligand_pdbqt/{i}/{i}_docked_ligand_rank{rank}.pdbqt'
    output_file = f'vina_feature/{i}_{rank}.csv'

    # Assuming 'cal_vina_features.py' is a standalone script, you can call it like this:
    command = f"python cal_vina_features.py {protein_path} {input_file} {output_file} {i}"
    subprocess.run(command, shell=True)

def main():
    pool = multiprocessing.Pool()
    func = partial(process_file, protein_path)

    ligand_dirs = os.listdir('docked_ligand_pdbqt')

    for i in ligand_dirs:
        for rank in range(1, 6):
            pool.apply_async(func, (i, rank))

    pool.close()
    pool.join()

if __name__ == '__main__':
    main()
