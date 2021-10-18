import pandas as pd
from icecream import ic
import os
from pathlib import Path

# Define global variable first
root_folder = Path(os.path.realpath(__file__)).parents[1]
ic(root_folder)
content_path = 'output_cleaned/data'


combined_dataframe = []
# First, mofid in pretest set
for file in os.listdir('data/extract_mofid'):
    if file.endswith("_80.csv") and "mofid_pretest" in file:
        ic(file)
        df = pd.read_csv(os.path.join(root_folder, content_path, 'extract_mofid', file), index_col=None, header=0)
        combined_dataframe.append(df)

frame = pd.concat(combined_dataframe)
frame = frame.drop(columns='Unnamed: 0', axis=1).sort_values(by='cifname', ignore_index=True)
frame.to_csv('mofid_pretest_combined.csv')

# Then, mofid in train set
combined_dataframe = []
for file in os.listdir('data/extract_mofid'):
    if file.endswith(".csv") and "mofid_train" in file:
        ic(file)
        df = pd.read_csv(os.path.join(root_folder, content_path, 'extract_mofid', file), index_col=None, header=0)
        combined_dataframe.append(df)

frame = pd.concat(combined_dataframe)
frame = frame.drop(columns='Unnamed: 0', axis=1).sort_values(by='cifname', ignore_index=True)
frame.to_csv('mofid_train_combined.csv')

# Now, find_surface in train set
combined_dataframe = []
for file in os.listdir('data/find_surface'):
    if file.startswith("clean") and file.endswith(".csv") and 'pretest' not in file:
        ic(file)
        df = pd.read_csv(os.path.join(root_folder, content_path, 'find_surface', file), index_col=None, header=0)
        combined_dataframe.append(df)

frame = pd.concat(combined_dataframe)
frame = frame.drop(columns='Unnamed: 0', axis=1).sort_values(by='MOFname', ignore_index=True)
frame.to_csv('clean_train_combined.csv')

# Finally, find_surface in pretest set
combined_dataframe = []
for file in os.listdir('data/find_surface'):
    if file.startswith("clean_pretest_") and file.endswith(".csv"):
        ic(file)
        df = pd.read_csv(os.path.join(root_folder, content_path, 'find_surface', file), index_col=None, header=0)
        combined_dataframe.append(df)

frame = pd.concat(combined_dataframe)
frame = frame.drop(columns='Unnamed: 0', axis=1).sort_values(by='MOFname', ignore_index=True)
frame.to_csv('clean_pretest_combined.csv')
