import os
import pickle
import logging
import argparse
import numpy as np


from utility import (
    extract_filenames,
    check_matching_filenames,
    log_also_console,
    load_metadata,
    load_profile_matrix,
    prepare_ranges_df,
    move_metadata_columns_to_ranges
)


def load_env_vars(script_dir, profile_main_dir, model_path):
    os.chdir(script_dir)

    with open(model_path, 'rb') as f:
        model = pickle.load(f)

    prediction_dir = os.path.join(profile_main_dir, "predictions")
    os.makedirs(prediction_dir, exist_ok=True)

    print(f"📁 Prediction output directory ready: {prediction_dir}")
    return script_dir, profile_main_dir, model


def wrapup_model_prediction(script_dir,
                            profile_dir,
                            profile_filename,
                            metadata_dir,
                            metadata_filename,
                            model,
                            output_dir,
                            name_prefix):
    input_profiles_subtnorm_path = os.path.join(profile_dir, profile_filename)
    input_metadata_path = os.path.join(metadata_dir, metadata_filename)

    split_profile = os.path.splitext(profile_filename)
    split_meta = os.path.splitext(metadata_filename)

    profile_base_name = split_profile[0]
    meta_base_name = split_meta[0]

    profile_ext = split_profile[1].lower()
    meta_ext = split_meta[1].lower()

    log_also_console(
        f"🔄 Predicting for: {profile_base_name+profile_ext} "
        f"(with metadata: {meta_base_name+meta_ext})"
    )

    metadata_df = load_metadata(input_metadata_path, meta_ext)
    subtnorm_np = load_profile_matrix(input_profiles_subtnorm_path,
                                      profile_ext,
                                      model)

    if subtnorm_np.shape[0] != len(metadata_df):
        raise ValueError("Mismatch between profile and metadata lengths.")

    y_proba = model.predict_proba(subtnorm_np)
    ranges_df = prepare_ranges_df(metadata_df, profile_base_name, y_proba)

    columns_to_keep = ["seqnames", "start", "end", "width", "strand",
                       'chrom', 'chromStart', 'chromEnd', 'name', 'score']
    ranges_df = move_metadata_columns_to_ranges(ranges_df,
                                                metadata_df,
                                                columns_to_keep)

    output_bed_file = os.path.join(output_dir,
                                   f"{name_prefix}_pred_all_"
                                   f"{profile_base_name}.bed")

    # Convert all numeric columns to string without scientific notation
    for col in ranges_df.columns:
        if np.issubdtype(ranges_df[col].dtype, np.number):
            ranges_df[col] = ranges_df[col].apply(lambda x: str(x))

    ranges_df.to_csv(output_bed_file, sep='\t', header=True, index=False)

    log_also_console(f"📄 Prediction saved to: {output_bed_file}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--script_dir", required=True)
    parser.add_argument("--profile_main_dir", required=True)
    parser.add_argument("--model_path", required=True)
    parser.add_argument("--log_file", default="prediction.log")
    parser.add_argument("--name_prefix", default="run")
    parser.add_argument("--single_prediction", action="store_true")
    parser.add_argument("--profile_filename",
                        help="Only used in single prediction mode")
    parser.add_argument("--metadata_filename",
                        help="Only used in single prediction mode")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="🔍 [%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[logging.FileHandler(args.log_file)]
    )

    script_dir, profile_main_dir, model = load_env_vars(
        args.script_dir,
        args.profile_main_dir,
        args.model_path
    )

    profile_dir = os.path.join(profile_main_dir, "profiles_subtnorm")
    metadata_dir = os.path.join(profile_main_dir, "metadata")
    output_dir = os.path.join(profile_main_dir, "predictions")

    if args.single_prediction:
        if not args.profile_filename or not args.metadata_filename:
            raise ValueError(
                "Both --profile_filename and --metadata_filename are required "
                "in single prediction mode."
            )
        wrapup_model_prediction(
            script_dir=script_dir,
            profile_dir=profile_dir,
            profile_filename=args.profile_filename,
            metadata_dir=metadata_dir,
            metadata_filename=args.metadata_filename,
            model=model,
            output_dir=output_dir,
            name_prefix=args.name_prefix
        )
    else:
        profile_files = extract_filenames(profile_dir)
        metadata_files = extract_filenames(metadata_dir)
        check_matching_filenames(profile_files, metadata_files)

        for p_file, m_file in zip(profile_files, metadata_files):
            wrapup_model_prediction(
                script_dir=script_dir,
                profile_dir=profile_dir,
                profile_filename=p_file,
                metadata_dir=metadata_dir,
                metadata_filename=m_file,
                model=model,
                output_dir=output_dir,
                name_prefix=args.name_prefix
            )


if __name__ == "__main__":
    main()
