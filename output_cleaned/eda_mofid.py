import pandas as pd
from icecream import ic
import os
from pathlib import Path
import numpy as np

# Define global variable
root_folder = Path(os.path.realpath(__file__)).parents[1]
ic(root_folder)
content_path = 'output_cleaned/data'
dataset_name = 'train'


# # It is established that csv and table loaded here is the same
# # Will try to find if those existed in the csv files
#
# # Load the train.csv given
# origin_train = pd.read_csv(os.path.join(root_folder, 'data', 'train.csv'))
#
# # Then, mofid in train set
# combined_dataframe = []
# for file in os.listdir('data/extract_mofid'):
#     if file.startswith('mofid_train_'):
#         ic(file)
#         df = pd.read_csv(os.path.join(root_folder, content_path, 'extract_mofid', file), index_col=None, header=0)
#         combined_dataframe.append(df)
#
# frame = pd.concat(combined_dataframe)
# frame = frame.drop(columns='Unnamed: 0', axis=1).sort_values(by='cifname', ignore_index=True)
# frame_dd = frame.drop_duplicates(subset='cifname', ignore_index=True)
# frame_dd.to_csv('mofid_train_combined.csv')
#
# # Perform the comparison
# missing = []
# for i in range(len(origin_train)):
#     if origin_train.iloc[i]['MOFname'] not in np.array(frame['cifname']):
#         missing.append(origin_train.iloc[i]['MOFname'])
# ic(missing)
#
# # Write it down
# missing_frame = pd.DataFrame(missing)
# missing_frame.to_csv('mofid_train_missing.csv')
#
#
# origin_pretest = pd.read_csv(os.path.join(root_folder, 'data', 'pretest.csv'))
#
# # Then, mofid in pretest set
# combined_dataframe = []
# for file in os.listdir('data/extract_mofid'):
#     if file.startswith('mofid_pretest_') and file.endswith('_10.csv'):
#         ic(file)
#         df = pd.read_csv(os.path.join(root_folder, content_path, 'extract_mofid', file), index_col=None, header=0)
#         combined_dataframe.append(df)
#
# file='mofid_pretest_correction.csv'
# df = pd.read_csv(os.path.join(root_folder, content_path, 'extract_mofid', file), index_col=None, header=0)
# combined_dataframe.append(df)
#
# frame = pd.concat(combined_dataframe)
# frame = frame.drop(columns='Unnamed: 0', axis=1).sort_values(by='cifname', ignore_index=True)
# frame_dd = frame.drop_duplicates(subset='cifname', ignore_index=True)
# frame.to_csv('mofid_pretest_combined.csv')
#
# # Perform the comparison
# missing = []
# for i in range(len(origin_pretest)):
#     if origin_pretest.iloc[i]['MOFname'] not in np.array(frame['cifname']):
#         missing.append(origin_pretest.iloc[i]['MOFname'])
# # Write it down
# missing_frame = pd.DataFrame(missing)
# missing_frame.to_csv('mofid_pretest_missing.csv')

# Then, mofid in test set
combined_dataframe = []
for file in os.listdir('output_cleaned/data/extract_mofid'):
    if file.startswith('mofid_test') and file.endswith('_20.csv'):
        ic(file)
        df = pd.read_csv(os.path.join(root_folder, content_path, 'extract_mofid', file), index_col=None, header=0)
        combined_dataframe.append(df)

frame = pd.concat(combined_dataframe)
frame = frame.drop(columns='Unnamed: 0', axis=1).sort_values(by='cifname', ignore_index=True)
frame_dd = frame.drop_duplicates(subset='cifname', ignore_index=True)
frame.to_csv('mofid_test_combined.csv')

# Perform the comparison

mofid_test_csv = pd.read_csv('mofid_test_combined.csv', index_col=0)
missing = []
for index, mof in mofid_test_csv.iterrows():
    if mof['cifname'] not in np.array(frame['cifname']):
        missing.append(mof['cifname'])
