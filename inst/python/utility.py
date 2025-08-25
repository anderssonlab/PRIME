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
import sys


def setup_logger(log_file=None):
    log_level = logging.DEBUG
    logger = logging.getLogger()
    logger.setLevel(log_level)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    # Always log to stderr
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(formatter)
    logger.addHandler(sh)

    return logger


def extract_filenames(dir_path, keep_word=None, drop_word=None):

    allowed_extensions = (".npz", ".csv", ".parquet")

    dir_path = os.path.normpath(dir_path)
    filename_ls = [f for f in os.listdir(dir_path)
                   if f.endswith(allowed_extensions)]

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


def build_matched_file_dict(profile_dir, metadata_dir,
                            profile_suffix="_profiles_subtnorm",
                            metadata_suffix="_metadata"):

    profile_files = extract_filenames(profile_dir)
    metadata_files = extract_filenames(metadata_dir)

    check_matching_filenames(profile_files, metadata_files,
                             profile_suffix=profile_suffix,
                             metadata_suffix=metadata_suffix)

    matched = {}

    for pf in profile_files:
        base, ext = os.path.splitext(pf)
        base = base.replace(profile_suffix, "")
        matched.setdefault(base, {})["profile_path"] = os.path.join(profile_dir, pf)
        matched[base]["profile_ext"] = ext

    for mf in metadata_files:
        base, ext = os.path.splitext(mf)
        base = base.replace(metadata_suffix, "")
        if base not in matched:
            continue  # or log missing profile
        matched[base]["metadata_path"] = os.path.join(metadata_dir, mf)
        matched[base]["metadata_ext"] = ext

    # Optionally validate
    unmatched = [k for k, v in matched.items() if "profile_path" not in v or "metadata_path" not in v]
    if unmatched:
        raise ValueError(f"❌ Unmatched samples found: {unmatched}")

    return matched


def define_num_cores(requested_core, total_jobs):

    cpu_count = os.cpu_count() or 1
    auto_core = min(25, cpu_count // 2)
    num_core = requested_core or max(1, auto_core)
    num_core = min(num_core, total_jobs)

    print(f"⚙️ Using {num_core} core{'s' if num_core > 1 else ''} for {total_jobs} job{'s' if total_jobs != 1 else ''}")
    return num_core


def extract_profiles(df):
    # Convert DataFrame to NumPy array
    numpy_array = df.values
    # Convert NaN values to 0
    numpy_array[np.isnan(numpy_array)] = 0
    return numpy_array


def set_index_if_exists(df, index_name='rownames'):
    if index_name in df.columns:
        df = df.set_index(index_name)
        df.index.name = None
    return df


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
                "⚠️ Model does not support sparse input — converting input profile to dense."
                f"Reason: {e}"
            )
            return sparse_mat.toarray()

    else:
        raise ValueError(f"❌ Unsupported profile file format: {ext}")


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


def adjust_genomic_positions_for_bed(ranges_df):
    ranges_df['chromStart'] = ranges_df['chromStart'] - 1
    # Convert to 0-based start
    return ranges_df


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


def move_metadata_columns_to_ranges(ranges_df, metadata_df, columns_to_keep):
    columns_to_move = [col for col in metadata_df.columns
                       if col not in columns_to_keep]
    for col in columns_to_move:
        ranges_df[col] = metadata_df[col]
    return ranges_df


def combine_bed_files(input_dir, output_dir):
    print(f"🔍 Combining BED files from: {input_dir}")
    os.makedirs(output_dir, exist_ok=True)

    bed_files = [f for f in glob.glob(os.path.join(input_dir, "*.bed"))
                 if not f.endswith("_combined.bed")]

    #bed_files = glob.glob(os.path.join(input_dir, "*.bed"))

    if not bed_files:
        print("⚠️ No .bed files found to combine.")
        return

    prefix_groups = defaultdict(list)

    # Group files by prefix (everything before _chr*)
    for bed_file in bed_files:
        filename = os.path.basename(bed_file)
        match = re.match(r"^(.*)_chr[^_\.\-]+(?=[_\.\-]|\.bed$)", filename)

        if match:
            prefix = match.group(1)
            prefix_groups[prefix].append(bed_file)
        else:
            print(
                f"❗ Skipping file (unexpected format): {filename}")

    if not prefix_groups:
        print("⚠️ No groups matched the expected pattern (_chr*).")
        return

    for prefix, files in prefix_groups.items():
        print(f"🔧 Combining files for prefix: {prefix}")
        combined_path = os.path.join(output_dir, f"{prefix}_combined.bed")
        combined_df = None

        for i, file in enumerate(sorted(files)):
            logging.info(f"   Reading {os.path.basename(file)}")
            df = pd.read_csv(file, sep="\t", header=0, comment="#")
            if combined_df is None:
                combined_df = df
            else:
                combined_df = pd.concat([combined_df, df],
                                        ignore_index=True)

        combined_df.to_csv(combined_path, sep="\t", header=True, index=False)

        print(f"✅ Combined files into {combined_path}")
