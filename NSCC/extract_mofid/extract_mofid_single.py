from icecream import ic
from mofid.run_mofid import cif2mofid
import os
import pandas as pd
from pathlib import Path
import sys


# Defining global parameters
root_repo_folder = Path(os.path.realpath(__file__)).parents[2]
dataset_name = 'train'
dataset = pd.read_csv(os.path.join(root_repo_folder, 'data', '{}.csv'.format(dataset_name)))
mofid_table_presort = []



def collect_result(mofid):
    # Callback function (for parallelization later)
    mofid_table_presort.append(mofid)


def one_mof(MOFname):
    # Adjust the mof information to correspond to the round_idx given, and direct round_idx to the CIF file
    mof_name = "{}.cif".format(MOFname)
    ic(mof_name)
    mof_cif = os.path.join(root_repo_folder, 'data', 'mof_cif_{}'.format(dataset_name), mof_name)
    # Call the cif2mofid function as described in github.com/snurr-group/mofid
    mofid = cif2mofid(mof_cif, output_path=os.path.join(root_repo_folder, 'data', 'mofid'))
    return mofid


def main():
    # Load mofid_train_missing.csv
    missing_MOFs = pd.read_csv(os.path.join(root_repo_folder, 'output_cleaned', 'mofid_train_missing.csv'),
                               index_col=0, header=0, names=['MOFname'])
    ic(missing_MOFs)
    missing_MOFs.info()

    # Implement the extraction function
    for MOFname in missing_MOFs['MOFname']:
        collect_result(one_mof(MOFname))

    mofid_table_final = pd.DataFrame(mofid_table_presort).sort_values(by='cifname', ignore_index=True)
    mofid_table_final.to_csv(os.path.join(root_repo_folder, 'data', 'mofid_train_correction.csv'))


if __name__ == '__main__':
    main()
