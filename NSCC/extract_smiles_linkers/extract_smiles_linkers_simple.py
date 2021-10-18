from icecream import ic
from openbabel import openbabel
from openbabel import pybel
import os
import pandas as pd
from pathlib import Path
import sys
import tempfile

# Defining global parameters
root_repo_folder = Path(os.path.realpath(__file__)).parents[1]
dataset = 'pretest'
dataset_name = 'mofid_{}_combined'.format(dataset)
dataset_dataframe = pd.read_csv(os.path.join(root_repo_folder, 'output_cleaned', '{}.csv'.format(dataset_name)), index_col=0)
dataset_dataframe.info()
smiles_linkers_full_dict = []
total_smiles_index = 0


# Load Linkers' SMILES into a DataFrame
for index, mof in dataset_dataframe.iterrows():
    # ic(mof['cifname'])
    smiles_linkers_list = [smiles_linker.strip(" \'") for smiles_linker in mof['smiles_linkers'].strip('[]').split(',')]
    cif_smiles_index = 0
    for smiles_linker in smiles_linkers_list:
        smiles_linkers_full_dict.append({'cifname': mof['cifname'],
                                         'cif_smiles_index': cif_smiles_index, 'smiles': smiles_linker})
        cif_smiles_index += 1
        total_smiles_index += 1

linkers_dataframe = pd.DataFrame(smiles_linkers_full_dict)
linkers_dataframe.info()
ic(linkers_dataframe.head())


smiles_xyz_counter = 0
smiles_xyz_total = linkers_dataframe['smiles'].size

# Now write every linker SMILES into xyz files
for index, smiles in linkers_dataframe.iterrows():
    mymol_smiles = pybel.readstring("smi", smiles['smiles'])
    mymol_smiles.addh()
    mymol_smiles.make3D(forcefield='gaff', steps=1000)
    if mymol_smiles.spin != 1:
        ic(index)
        ic(mymol_smiles.charge)
        ic(mymol_smiles.spin)
    ic('smiles_linker_{}_{}_{}.xyz'.format(dataset, smiles_xyz_counter, smiles_xyz_total))
    mymol_smiles.write(format='orcainp', filename='smiles_linker_{}_{}_{}.xyz'.format(dataset, smiles_xyz_counter,
                                                                                  smiles_xyz_total), overwrite=True)
    smiles_xyz_counter += 1


# # Code for one structure, know to work
# smiles_linker_example = smiles_linkers_full_list[0][0]
# mymol = pybel.readstring("smi", smiles_linker_example)
# ic(smiles_linker_example)
# mymol.make3D(forcefield='gaff', steps=1000)
# mymol.write(format='xyz', filename='smiles_linker_example_befH.xyz', overwrite=True)
# mymol = pybel.readstring("smi", smiles_linker_example)
# mymol.addh()
# mymol.make3D(forcefield='gaff', steps=1000)
# mymol.write(format='xyz', filename='smiles_linker_example_aftH.xyz', overwrite=True)


# # Loop the xyz generation on every linkers' smiles in each mof cif
# for smiles_linkers in smiles_linkers_full_list:
#     for smiles_linker in smiles_linkers:
#         obConversion.ReadString(mol, smiles_linker)
#         ic(mol.NumAtoms())
#         mol.AddHydrogens()
#         ic(mol.NumAtoms())
#
#     for smiles_linker in smiles_linkers_list:
#         with tempfile.NamedTemporaryFile('w', encoding='utf8') as tmp:
#             ic(tmp.name)
#             ic(smiles_linker)
#             tmp.write(smiles_linker)
