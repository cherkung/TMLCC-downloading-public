from porE.hea.HEA import HEA
from pore import porosity as p
from pore import psd
from porE.io.ase2pore import *
from icecream import ic
import os
import pandas as pd
from pathlib import Path
import sys
import tempfile

# Defining global parameters
root_repo_folder = Path(os.path.realpath(__file__)).parents[2]
dataset_name = 'test'
dataset = pd.read_csv(os.path.join(root_repo_folder, 'data', '{}.csv'.format(dataset_name)))
pores_table_presort = []


def collect_result(pores):
    # Callback function (for parallelization)
    pores_table_presort.append(pores)


def one_mof(MOFname):
    mof_name = MOFname + '.cif'
    ic(mof_name)
    mof_cif = os.path.join(root_repo_folder, 'data', 'mof_cif_{}'.format(dataset_name), mof_name)
    result = {'MOFname': MOFname}
    # Execute HEA, using the cif file
    Phi_HEA = HEA(mof_cif)[0]
    result['Phi_HEA'] = Phi_HEA
    # convertes cif to porE xyz (i.e., pypore.xyz)
    ase2pore(mof_cif)
    structure = 'pypore.xyz'
    #
    # Execute porosity evaluation, using the overlapping sphere approach (OSA)
    #
    print('-----------')
    print('\nRun OSA\n')
    print('-----------')
    Phi, density, poreV, V_total, V_vdwSum, V_overlap = p.osa(structure)
    result['poreV_OSE'] = poreV
    #
    # Execute an analyis of the pore size distribution (PSD)
    # First number   -> Number of starting points
    # Second number  -> Number of MC steps
    #
    print('-----------')
    print('\nRun PSD\n')
    print('-----------')
    no_pores, tmp1, tmp2, tmp3, tmp4 = psd.get_psd(structure, 200, 1000)
    pores = tmp1[0:no_pores]
    distr = tmp2[0:no_pores]
    pore_pos_cart = tmp3[0:no_pores]
    pore_pos_frac = tmp4[0:no_pores]
    #
    # Execute porosity evaluation, using the grid point apporach (GPA)
    # a probe radius for the accessible porosity needs to be provided
    #
    probe_R = 1.20
    #
    # Here, explicitly provide the full grid (grid points along each cell vector)
    # FullGrid
    grid_a = grid_b = grid_c = 30
    print('-----------------------------------')
    print('\nRun GPA: grid_a, grid_b, grid_c\n')
    print('-----------------------------------')
    Phi_void, Phi_acc, density, poreV_void, poreV_acc = p.gpa_fullgrid(structure, probe_R, grid_a, grid_b, grid_c)
    #
    # Here, a grid point density per A is provided instead
    # GridPerA
    grid_density = 2.0
    print('-------------------------')
    print('\nRun GPA: grid_density\n')
    print('-------------------------')
    Phi_void, Phi_acc, density, poreV_void, poreV_acc = p.gpa_gridpera(structure, probe_R, grid_density)
    result['Phi_void'] = Phi_void
    result['Phi_acc'] = Phi_acc
    ic(result)
    return result


def main():
    # Load pores_{}_missing.csv
    missing_MOFs = pd.read_csv(os.path.join(root_repo_folder, 'output_cleaned',
                                            'pores_{}_missing.csv'.format(dataset_name)), index_col=0, header=0,
                               names=['MOFname'])
    ic(missing_MOFs)
    missing_MOFs.info()

    # Implement the extraction function
    for MOFname in missing_MOFs['MOFname']:
        ic(MOFname)
        collect_result(one_mof(MOFname))

    pores_table_final = pd.DataFrame(pores_table_presort).sort_values(by='MOFname', ignore_index=True)
    pores_table_final.to_csv(os.path.join(root_repo_folder, 'data', 'pores_{}_correction.csv'.format(dataset_name)))


if __name__ == "__main__":
    main()
