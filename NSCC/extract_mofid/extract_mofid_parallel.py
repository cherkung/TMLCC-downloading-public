from icecream import ic
from mofid.run_mofid import cif2mofid
import multiprocessing as mp
import os
import pandas as pd
from pathlib import Path
import sys


# Defining global parameters
root_repo_folder = Path(os.path.realpath(__file__)).parents[2]
dataset_name = 'test'
dataset = pd.read_csv(os.path.join(root_repo_folder, 'data', '{}.csv'.format(dataset_name)))
mofid_table_presort = []


def collect_result(mofid):
    # Callback function
    mofid_table_presort.append(mofid)


def one_mof(round_idx):
    mof_name = dataset.iloc[round_idx]['MOFname'] + '.cif'
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
        round_number = 1
        round_start_idx = 0
        round_end_idx = len(dataset.index)
        max_round = 1

    # Implement the extraction function with Parallelize
    pool = mp.Pool(mp.cpu_count())
    for round_idx in range(round_start_idx, round_end_idx):
        pool.apply_async(one_mof, args=(round_idx, ), callback=collect_result)
    pool.close()
    pool.join()

    # sort the mofid table and close those that we don't need
    mofid_table = pd.DataFrame(mofid_table_presort).sort_values(by='cifname', ignore_index=True)
    mofid_table.to_csv(os.path.join(root_repo_folder, 'data', 'mofid_{}_{}_{}.csv'.format(dataset_name, round_number,
                                                                                          max_round)))


if __name__ == '__main__':
    main()
