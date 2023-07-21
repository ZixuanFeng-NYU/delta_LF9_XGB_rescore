import os
import multiprocessing

def process_files(i):
    folder_path = f"docked_ligand_pdbqt/{i}"
    os.makedirs(folder_path, exist_ok=True)

    for j in range(1, 6):
        sdf_file = f"folder_combined/index{i}_*/rank{j}_confidence*.sdf"
        if not os.path.isfile(sdf_file):
            output_file = f"{folder_path}/{i}_docked_ligand_rank{j}_added_h.sdf"
            command = f"obabel {sdf_file} -h -osdf -O {output_file}"
            os.system(command)

    for j in range(1, 6):
        pdbqt_file = f"docked_ligand_pdbqt/{i}/{i}_docked_ligand_rank{j}.pdbqt"
        if not os.path.isfile(pdbqt_file):
            sdf_file = f"docked_ligand_pdbqt/{i}/{i}_docked_ligand_rank{j}_added_h.sdf"
            output_file = f"docked_ligand_pdbqt/{i}/{i}_docked_ligand_rank{j}.pdbqt"
            command = f"/home/zf2012/Meeko/scripts/mk_prepare_ligand.py -i {sdf_file} -o {output_file}"
            os.system(command)

if __name__ == '__main__':
    num_processes = multiprocessing.cpu_count()  # Number of CPU cores
    pool = multiprocessing.Pool(processes=num_processes)

    ligand_id_list=os.listdir('folder_combined')
    id_list=[ligand_id.split('_')[0].split('x')[1] for ligand_id in ligand_id_list]
    # Generate the range of values for 'i'
    i_values = id_list

    # Map the function to the list of 'i' values
    pool.map(process_files, i_values)

    pool.close()
    pool.join()
