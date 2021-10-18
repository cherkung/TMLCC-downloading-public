from icecream import ic
from mofid.run_mofid import cif2mofid
import os
import pandas as pd
from pathlib import Path
import sys


# Defining global parameters
root_repo_folder = Path(os.path.realpath(__file__)).parents[2]
dataset_name = 'pretest'
dataset = pd.read_csv(os.path.join(root_repo_folder, 'data', '{}.csv'.format(dataset_name)))
mofid_table_presort = []



def collect_result(mofid):
    # Callback function (for parallelization later)
    mofid_table_presort.append(mofid)


def one_mof(round_idx):
    # Adjust the mof information to correspond to the round_idx given, and direct round_idx to the CIF file
    mof_idx = round_idx + 1
    if dataset_name == 'pretest':
        mof_name = 'mof_unit_pretest_{}.cif'.format(mof_idx)
    else:
        mof_name = 'mof_unit_{}.cif'.format(mof_idx)
    ic(mof_name)
    mof_cif = os.path.join(root_repo_folder, 'data', 'mof_cif_{}'.format(dataset_name), mof_name)
    # Call the cif2mofid function as described in github.com/snurr-group/mofid
    mofid = cif2mofid(mof_cif, output_path=os.path.join(root_repo_folder, 'data', 'mofid'))
    return mofid


def main():
    # Parsing arguments from the job submission script
    try:
        if int(sys.argv[2]) == 1:
            raise IndexError
        round_number = int(sys.argv[1])
        max_round = int(sys.argv[2])
        if round_number != max_round:
            round_start_idx = len(dataset.index) // (max_round - 1) * (round_number - 1)
            round_end_idx = round_start_idx + len(dataset.index) // (max_round - 1)
        else:
            round_start_idx = len(dataset.index) // (max_round - 1) * (round_number - 1)
            round_end_idx = len(dataset.index)
        pass
    except IndexError:
        round_number = 0
        round_start_idx = 0
        round_end_idx = len(dataset.index)
        max_round = 1

    # Implement the extraction function, preparing for parallelization later
    for round_idx in range(round_start_idx, round_end_idx):
        collect_result(one_mof(round_idx))
        if round_idx % 10 == 0:
            mofid_table_temp = pd.DataFrame(mofid_table_presort).sort_values(by='cifname', ignore_index=True)
            mofid_table_temp.to_csv(
                os.path.join(root_repo_folder, 'data', 'mofid_temp_{}_{}_{}.csv'.format(dataset_name, round_number,
                                                                                        max_round)))

    mofid_table_final = pd.DataFrame(mofid_table_presort).sort_values(by='cifname', ignore_index=True)
    mofid_table_final.to_csv(os.path.join(root_repo_folder, 'data', 'mofid_{}_{}_{}.csv'.format(dataset_name,
                                                                                                round_number,
                                                                                                max_round)))


if __name__ == '__main__':
    main()
