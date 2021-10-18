import pandas as pd
from icecream import ic
import os
from pathlib import Path
import numpy as np

# Define global variable
root_folder = Path(os.path.realpath(__file__)).parents[1]
ic(root_folder)
content_path = 'output_cleaned/data'


# It is established that csv and table loaded here is the same
# Will try to find if those existed in the csv files

# # Load the csv output
#
# mofid_pretest_csv = pd.read_csv('clean_train_combined.csv')
# missing_csv = []
# for i in range(2000):
#     if "mof_unit_{}.cif".format(i+1) not in np.array(mofid_pretest_csv["MOFname"]):
#         missing_csv.append("mof_unit_{}.cif".format(i+1))
#
# # Create an array of adjacent mof position
# adjacent_MOFs = []
# for missing_mof in missing_csv:
#     adjacent_left = 'mof_unit_{}.cif'.format(int(missing_mof.split('_')[2].split('.')[0])-1)
#     adjacent_MOFs.append(adjacent_left)
#     adjacent_right = 'mof_unit_{}.cif'.format(int(missing_mof.split('_')[2].split('.')[0])+1)
#     adjacent_MOFs.append(adjacent_right)
# ic(adjacent_MOFs)
#
# # Finally, find_surface in pretest set
# combined_dataframe = []
# for file in os.listdir('data/find_surface'):
#     if file.startswith("clean_") and file.endswith(".csv") and 'pretest' not in file:
#         ic(file)
#         df = pd.read_csv(os.path.join(root_folder, content_path, 'find_surface', file), index_col=None, header=0)
#         for adjacent_MOF in adjacent_MOFs:
#             if adjacent_MOF in np.array(df['MOFname']):
#                 ic(adjacent_MOF)
#                 ic(file)
#         combined_dataframe.append(df)
#
# frame = pd.concat(combined_dataframe)
# frame = frame.drop(columns='Unnamed: 0', axis=1).sort_values(by='MOFname', ignore_index=True)
# frame.to_csv('clean_train_combined.csv')
#
# missing = []
# for i in range(2000):
#     if "mof_unit_{}.cif".format(i+1) not in np.array(frame["MOFname"]):
#         missing.append("mof_unit_{}.cif".format(i+1))

# Finally, find_surface in test set
combined_dataframe = []
for file in os.listdir(os.path.join(root_folder, 'output_cleaned', 'data', 'find_surface')):
    if file.startswith("clean_test_") and file.endswith("20.csv"):
        ic(file)
        df = pd.read_csv(os.path.join(root_folder, content_path, 'find_surface', file), index_col=None, header=0)
        combined_dataframe.append(df)

frame = pd.concat(combined_dataframe)
frame = frame.drop(columns='Unnamed: 0', axis=1).sort_values(by='MOFname', ignore_index=True)
frame.to_csv('clean_test_combined.csv')

clean_test_csv = pd.read_csv('clean_test_combined.csv', index_col=0)
missing_csv = []
for index, mof in clean_test_csv.iterrows():
    if mof['MOFname'] not in np.array(frame['MOFname']):
        missing_csv.append(mof['MOFname'])
