import multiprocessing
from functools import partial
from cal_vina_features import cal_vina_features

protein_path = 'protein_pdbqt'

def process_file(protein_path, i, rank):
    input_file = f''
    output_file = f''
    return cal_vina_features(protein_path, input_file, output_file, i)

def main():
    pool = multiprocessing.Pool()
    func = partial(process_file, protein_path)

    for i in os.listdir('docked_ligand_pdbqt'):
        for rank in range(1, 6):
            pool.apply_async(func, (i, rank))

    pool.close()
    pool.join()

if __name__ == '__main__':
    main()

