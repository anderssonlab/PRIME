#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: zmk214
"""

import os
import numpy as np
import pandas as pd
import re
import logging
from scipy.sparse import load_npz
import glob
from collections import defaultdict


def log_also_console(message, level="info"):
    print(message)  # Console
    if level == "info":
        logging.info(message)
    elif level == "warning":
        logging.warning(message)
    elif level == "error":
        logging.error(message)
    elif level == "debug":
        logging.debug(message)
    else:
        logging.info(message)


def extract_filenames(dir_path, keep_word=None, drop_word=None):

    filename_ls = os.listdir(dir_path)

    if keep_word:
        filename_ls = [name for name in filename_ls
                       if all(substring in name
                              for substring in keep_word)]
    elif keep_word is not None:
        raise ValueError(
            "Invalid parameter type for keep_word. Expected list of strings.")

    if drop_word:
        filename_ls = [name for name in filename_ls
                       if all(substring not in name
                              for substring in drop_word)]
    elif drop_word is not None:
        raise ValueError(
            "Invalid parameter type for drop_word. Expected list of strings.")

    if '.DS_Store' in filename_ls:
        filename_ls.remove('.DS_Store')

    filename_ls.sort()

    return filename_ls


def strip_suffix(filename, suffix):
    """
    Removes a given suffix and extension from a filename.
    Example: strip_suffix("sample_chr1_profiles_subtnorm.parquet",
                          "_profiles_subtnorm")
             → "sample_chr1"
    """
    pattern = re.escape(suffix) + r"\.[^.]+$"
    return re.sub(pattern, "", filename)


def check_matching_filenames(profile_filelist, metadata_filelist,
                             profile_suffix="_profiles_subtnorm",
                             metadata_suffix="_metadata"):
    profile_bases = [strip_suffix(f, profile_suffix)
                     for f in profile_filelist]
    metadata_bases = [strip_suffix(f, metadata_suffix)
                      for f in metadata_filelist]

    if profile_bases != metadata_bases:
        mismatches = [(p, m) for p, m in zip(profile_bases, metadata_bases)
                      if p != m]
        logging.error(
            f"❌ Filename mismatch detected! Mismatched pairs: {mismatches}")
        raise ValueError(
            f"Filename mismatch detected! Mismatched pairs: {mismatches}")

    logging.info("✅ Profile and metadata filenames match.")


def set_index_if_exists(df, index_name='rownames'):
    """
    Sets a specified column as the index of the DataFrame if it exists,
    and removes the index name to prevent it from appearing as an additional
    column in the DataFrame.

    Parameters:
    df (pd.DataFrame): The DataFrame in which the index needs to be set.
    index_name (str): The name of the column to set as the index.

    Returns:
    pd.DataFrame: The DataFrame with the specified column set as the index.
    """
    if index_name in df.columns:
        df = df.set_index(index_name)
        df.index.name = None
    return df


def load_metadata(metadata_path, ext):
    if ext == ".parquet":
        logging.info(
            f"📥 Reading metadata file: {metadata_path} (format: {ext})")
        df = pd.read_parquet(metadata_path)
    elif ext == ".csv":
        logging.info(
            f"📥 Reading metadata file: {metadata_path} (format: {ext})")
        df = pd.read_csv(metadata_path, header=0, index_col=None)
    else:
        raise ValueError(f"❌ Unsupported metadata file format: {ext}")

    return set_index_if_exists(df)


def extract_profiles(df):
    """
    Extract profile vectors from a DataFrame.

    Parameters:
    df (pandas.DataFrame): DataFrame containing profile data.

    Returns:
    numpy.ndarray: NumPy array representing profile vectors.
    """
    # Convert DataFrame to NumPy array
    numpy_array = df.values

    # Convert NaN values to 0
    numpy_array[np.isnan(numpy_array)] = 0

    return numpy_array


def adjust_genomic_positions_for_bed(ranges_df):
    """
    Adjusts genomic positions for BED file format.

    Parameters:
    ranges_df (pd.DataFrame):
    DataFrame with columns 'chrom', 'chromStart', 'chromEnd', 'strand'.

    Returns:
    pd.DataFrame: Adjusted DataFrame ready for BED file output.
    """
    ranges_df['chromStart'] = ranges_df['chromStart'] - 1
    # Convert to 0-based start
    return ranges_df


def move_metadata_columns_to_ranges(ranges_df, metadata_df, columns_to_keep):
    """
    Moves columns from the metadata DataFrame to the ranges DataFrame,
    excluding the columns specified in `columns_to_keep`.

    Parameters:
    ranges_df (pd.DataFrame): The DataFrame containing genomic ranges.
    metadata_df (pd.DataFrame): The DataFrame containing additional metadata.
    columns_to_keep (list): List of columns that should not be moved
    from metadata_df.

    Returns:
    pd.DataFrame: The updated ranges DataFrame
    with additional metadata columns.
    """
    columns_to_move = [col for col in metadata_df.columns
                       if col not in columns_to_keep]
    for col in columns_to_move:
        ranges_df[col] = metadata_df[col]
    return ranges_df


def load_profile_matrix(profile_path, ext, model):
    if ext == ".parquet":
        logging.info(f"   Reading profile file: {profile_path} (format: {ext})")
        df = pd.read_parquet(profile_path)
        df = set_index_if_exists(df)
        return extract_profiles(df)

    elif ext == ".csv":
        logging.info(f"   Reading profile file: {profile_path} (format: {ext})")
        df = pd.read_csv(profile_path, header=0, index_col=None)
        df = set_index_if_exists(df)
        return extract_profiles(df)

    elif ext == ".npz":
        logging.info(f"   Reading profile file: {profile_path} (format: {ext})")
        sparse_mat = load_npz(profile_path)
        try:
            model.predict_proba(sparse_mat[0])
            logging.info("✅ Model accepts sparse input — using sparse matrix.")
            return sparse_mat
        except Exception as e:
            logging.warning(
                "⚠️ Model does not support sparse input — converting to dense."
                f"Reason: {e}"
            )
            return sparse_mat.toarray()

    else:
        raise ValueError(f"❌ Unsupported profile file format: {ext}")


def prepare_ranges_df(metadata_df, profile_base_name, y_proba):
    ranges_df = metadata_df.copy()

    ranges_df.rename(columns={
        'seqnames': 'chrom',
        'start': 'chromStart',
        'end': 'chromEnd'
    }, inplace=True)

    ranges_df['strand'] = ranges_df['strand'].astype(str)
    ranges_df['score'] = y_proba[:, 1]
    ranges_df['name'] = profile_base_name

    ranges_df = ranges_df.reindex(columns=['chrom',
                                           'chromStart',
                                           'chromEnd',
                                           'name',
                                           'score',
                                           'strand'])
    ranges_df['chromStart'] = ranges_df['chromStart'].astype(int)
    ranges_df['chromEnd'] = ranges_df['chromEnd'].astype(int)

    ranges_df = ranges_df.sort_values(by=['chrom', 'chromStart'])
    ranges_df = adjust_genomic_positions_for_bed(ranges_df)

    return ranges_df


def combine_bed_files(input_dir, output_dir):
    log_also_console(f"🔍 Combining BED files from: {input_dir}")
    os.makedirs(output_dir, exist_ok=True)

    bed_files = glob.glob(os.path.join(input_dir, "*.bed"))
    if not bed_files:
        log_also_console("⚠️ No .bed files found to combine.")
        return

    prefix_groups = defaultdict(list)

    # Group files by prefix (everything before _chr*)
    for bed_file in bed_files:
        filename = os.path.basename(bed_file)
        match = re.match(r"(.+)(_chr[^_]+)_profiles_subtnorm\.bed", filename)

        if match:
            prefix = match.group(1)
            prefix_groups[prefix].append(bed_file)
        else:
            log_also_console(
                f"❗ Skipping file (unexpected format): {filename}")

    if not prefix_groups:
        log_also_console("⚠️ No groups matched the expected pattern (_chr*).")
        return

    for prefix, files in prefix_groups.items():
        log_also_console(f"🔧 Combining files for prefix: {prefix}")
        combined_path = os.path.join(output_dir, f"{prefix}_combined.bed")
        combined_df = None

        for i, file in enumerate(sorted(files)):
            logging.info(f"   Reading {os.path.basename(file)}")
            df = pd.read_csv(file, sep="\t", header=0)
            if combined_df is None:
                combined_df = df
            else:
                if df.columns.equals(combined_df.columns):
                    combined_df = pd.concat([combined_df, df.iloc[1:]],
                                            ignore_index=True)
                else:
                    log_also_console(
                        f"⚠️ Column mismatch in {file}, including full file."
                    )
                    combined_df = pd.concat([combined_df, df],
                                            ignore_index=True)

        combined_df.to_csv(combined_path, sep="\t", header=True, index=False)
        log_also_console(f"✅ Combined files into {combined_path}")
