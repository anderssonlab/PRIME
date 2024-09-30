import os
import numpy as np
import pickle
import pandas as pd
import pyarrow.parquet as pq
from lightgbm import LGBMClassifier

import argparse





def extract_filenames(dir_path, keep_word=None, drop_word=None):

    filename_ls = os.listdir(dir_path)

    if keep_word:
        filename_ls = [name for name in filename_ls if all(substring in name for substring in keep_word)]
    elif keep_word is not None:
        raise ValueError("Invalid parameter type for keep_word. Expected list of strings.")

    if drop_word:
        filename_ls = [name for name in filename_ls if all(substring not in name for substring in drop_word)]
    elif drop_word is not None:
        raise ValueError("Invalid parameter type for drop_word. Expected list of strings.")

    if '.DS_Store' in filename_ls:
        filename_ls.remove('.DS_Store')

    filename_ls.sort()

    return(filename_ls)

def extract_ranges(df):

    """
    Extracts chromosome, start, end, and strand from DataFrame index.

    Args:
    df (pandas.DataFrame): DataFrame with index containing strings in the format "chrom:start-end;strand".

    Returns:
    pandas.DataFrame: DataFrame with 'chrom', 'chromStart', 'chromEnd', and 'strand' columns.
    """
    # Split index and expand into separate columns
    split_index = df.index.str.split(':|-|;', expand=True)

    # Untuple the split_index
    chrom, chromStart, chromEnd, strand = zip(*split_index.values)

    # Create a DataFrame with column names
    result_df = pd.DataFrame({
        'chrom': chrom,
        'chromStart': chromStart,
        'chromEnd': chromEnd,
        'strand': strand
    })

    return result_df

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
        df.index.name = None  # Remove the name of the index to avoid an extra column named 'rownames'
    return df


def adjust_genomic_positions_for_bed(ranges_df):
    """
    Adjusts genomic positions for BED file format.

    Parameters:
    ranges_df (pd.DataFrame): DataFrame with columns 'chrom', 'chromStart', 'chromEnd', 'strand'.

    Returns:
    pd.DataFrame: Adjusted DataFrame ready for BED file output.
    """
    ranges_df['chromStart'] = ranges_df['chromStart'] - 1  # Convert to 0-based start
    return ranges_df


def move_metadata_columns_to_ranges(ranges_df, metadata_df, columns_to_keep):
    """
    Moves columns from the metadata DataFrame to the ranges DataFrame,
    excluding the columns specified in `columns_to_keep`.

    Parameters:
    ranges_df (pd.DataFrame): The DataFrame containing genomic ranges.
    metadata_df (pd.DataFrame): The DataFrame containing additional metadata.
    columns_to_keep (list): List of columns that should not be moved from metadata_df.

    Returns:
    pd.DataFrame: The updated ranges DataFrame with additional metadata columns.
    """
    columns_to_move = [col for col in metadata_df.columns if col not in columns_to_keep]
    for col in columns_to_move:
        ranges_df[col] = metadata_df[col]
    return ranges_df

# Function to load environment variables from .env file
def load_env_vars(script_dir, profile_main_dir, subdir_name, model_path):

    # Set the working directory
    # os.chdir(script_dir)

    # Import model
    model = pickle.load(open(model_path, 'rb'))

    # Ensure the output directory exists
    if not os.path.exists(profile_main_dir+"/predictions"):
        os.makedirs(profile_main_dir+"/predictions")
        print(f"Directory '{profile_main_dir}/predictions' created successfully.")
    else:
        print(f"Directory '{profile_main_dir}/predictions' already exists.")

    # Ensure the output sub-directory exists
    if not os.path.exists(profile_main_dir+"/predictions/"+subdir_name):
        os.makedirs(profile_main_dir+"/predictions/"+subdir_name)
        print(f"Directory '{profile_main_dir}/predictions/{subdir_name}' created successfully.")
    else:
        print(f"Directory '{profile_main_dir}/predictions/{subdir_name}' already exists.") 

    return script_dir, profile_main_dir, subdir_name, model


def wrapup_model_prediction(script_dir, profile_dir, profile_filename, metadata_dir, metadata_filename, model, output_dir, name_prefix, threshold, file_format='parquet'):

    # Example usage of extract_filenames
    filenames_without_extensions = os.path.splitext(profile_filename)[0]
    print(filenames_without_extensions)

    # Define the input paths
    input_profiles_subtnorm_path = os.path.join(profile_dir, profile_filename)
    input_metadata_path = os.path.join(metadata_dir, metadata_filename)

    # Read files based on the specified format
    if file_format == 'parquet':
        subtnorm_df = pd.read_parquet(input_profiles_subtnorm_path)
        metadata_df = pd.read_parquet(input_metadata_path)
    elif file_format == 'csv':
        subtnorm_df = pd.read_csv(input_profiles_subtnorm_path, header=0, index_col=None)
        metadata_df = pd.read_csv(input_metadata_path, header=0, index_col=None)
    else:
        raise ValueError("Unsupported file format. Please use 'parquet' or 'csv'.")

    # Set "rownames" as index if it exists
    subtnorm_df = set_index_if_exists(subtnorm_df)
    metadata_df = set_index_if_exists(metadata_df)

    # Ensure lengths match
    if len(subtnorm_df) != len(metadata_df):
        raise ValueError("The metadata and profile subtnorm dataframes do not have the same length")
        sys.exit(1)

    # Extract ranges_df and profiles
    ranges_df = extract_ranges(subtnorm_df)
    subtnorm_np = extract_profiles(subtnorm_df)

    # Predict probabilities
    y_proba = model.predict_proba(subtnorm_np)
    ranges_df['score'] = y_proba[:, 1]

    ranges_df['name'] = filenames_without_extensions

    # Reorder columns and ensure chromStart and chromEnd are integers
    ranges_df = ranges_df.reindex(columns=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'])
    ranges_df['chromStart'] = ranges_df['chromStart'].astype(int)
    ranges_df['chromEnd'] = ranges_df['chromEnd'].astype(int)

    # Sort the DataFrame by 'chrom' and 'chromStart'
    ranges_df = ranges_df.sort_values(by=['chrom', 'chromStart'])

    # Adjust genomic positions for BED format
    ranges_df = adjust_genomic_positions_for_bed(ranges_df)

    # Move extra metadata columns to ranges_df
    columns_to_keep = ["seqnames", "start", "end", "width", "strand", 'chrom', 'chromStart', 'chromEnd', 'name', 'score']
    ranges_df = move_metadata_columns_to_ranges(ranges_df, metadata_df, columns_to_keep)

    # Save results to .bed files
    output_all_results = os.path.join(output_dir, f'{name_prefix}_pred_all_{filenames_without_extensions}.bed') 
    ranges_df.to_csv(output_all_results, sep='\t', header=True, index=False)

    # Save selected results if score >= threshold
    output_slt_results = os.path.join(output_dir, f'{name_prefix}_pred_slt{threshold}_{filenames_without_extensions}.bed')
    selected_ranges = ranges_df[ranges_df['score'] >= threshold]
    selected_ranges.to_csv(output_slt_results, sep='\t', header=True, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script to predict genomewide TC normalized data.')

    parser.add_argument('-w', '--script_dir', default=".", type=str, 
                        help='Path to the script directory')
    parser.add_argument('-p', '--profile_main_dir', type=str, required=True, 
                        help='Path to the profile directory')
    parser.add_argument('-r', '--profile_sub_dir', type=str, default="tcs", 
                        help='Sub-directory name in the main profile directory')
    parser.add_argument('-m', '--model_path', type=str, required=True, 
                        help='Path to the input model')
    parser.add_argument('-n', '--name_prefix', type=str, required=True, 
                        help='Name added to the output files, indicate model name and library (celltype) name') 
    parser.add_argument('-t', '--threshold', type=float, default=0.5,
                        help='Threshold value for prediction')
    parser.add_argument('-f', '--file_format', type=str, default='parquet', choices=['parquet', 'csv'],
                        help='File format for input files (default: parquet)')

    args = parser.parse_args()

    # Load environment variables
    script_dir = args.script_dir
    profile_main_dir = args.profile_main_dir
    profile_sub_dir = args.profile_sub_dir.split(',')
    model_path = args.model_path

    name_prefix = args.name_prefix
    file_format = args.file_format

    threshold = args.threshold


    for subdir_name in profile_sub_dir:

        # Load environment variables and model
        script_dir, profile_main_dir, subdir_name, model = load_env_vars(script_dir, profile_main_dir, subdir_name, model_path)

        profile_dir = profile_main_dir+"/profiles_subtnorm/"+subdir_name
        metadata_dir = profile_main_dir+"/metadata/"+subdir_name
        output_dir = profile_main_dir+"/predictions/"+subdir_name

        profile_file_ls = extract_filenames(profile_dir)
        metadata_file_ls = extract_filenames(metadata_dir)

        # Execute data processing  
        for i in range(len(profile_file_ls)):
            wrapup_model_prediction(script_dir,
                                    profile_dir,
                                    profile_file_ls[i],
                                    metadata_dir,
                                    metadata_file_ls[i],
                                    model,
                                    output_dir,
                                    name_prefix,
                                    threshold,
                                    file_format=file_format)