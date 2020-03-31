import os
import pandas as pd


def merge_summary_data(folder):
    """
    This function will merge all files in a folder ending with the suffix "summary.csv".

    Arguments:
        folder: folder containing the files to be merged.

    Returns:
        A single csv file containing the merged data.
    """
    file_list = [f for f in os.listdir(folder) if f.endswith('summary.csv')]
    combined_csv = pd.concat([pd.read_csv(folder + f, header=0) for f in file_list])
    combined_csv.to_csv(folder + "agg_data.csv", index=False)
