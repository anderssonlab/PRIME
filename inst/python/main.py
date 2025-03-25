import os
import pickle
import logging
import argparse
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed


from utility import (
    extract_filenames,
    check_matching_filenames,
    log_also_console,
    load_metadata,
    load_profile_matrix,
    prepare_ranges_df,
    move_metadata_columns_to_ranges,
    combine_bed_files
)


def load_env_vars(script_dir, profile_main_dir, model_path):
    os.chdir(script_dir)

    with open(model_path, 'rb') as f:
        model = pickle.load(f)

    prediction_dir = os.path.join(profile_main_dir, "predictions")
    os.makedirs(prediction_dir, exist_ok=True)

    log_also_console(f"📁 Prediction output directory ready: {prediction_dir}")
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

    logging.info(
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
        if pd.api.types.is_numeric_dtype(ranges_df[col]):
            ranges_df[col] = ranges_df[col].apply(lambda x: str(x))

    ranges_df.to_csv(output_bed_file, sep='\t', header=True, index=False)

    logging.info(f"📄 Prediction saved to: {output_bed_file}")


def worker_wrap(arg):
    wrapup_model_prediction(*arg)


def parallel_prediction(task_args, num_workers):
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(worker_wrap, arg) for arg in task_args]
        for future in as_completed(futures):
            future.result()  # Raise any exceptions here


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--script_dir", required=True)
    parser.add_argument("--profile_main_dir", required=True)
    parser.add_argument("--combined_outdir",
                        help=(
                            "Output dir to save combined BED files "
                            "(batch mode only)"
                        ))
    parser.add_argument("--model_path", required=True)
    parser.add_argument("--log_file", default="prediction.log")
    parser.add_argument("--name_prefix", default="run")
    parser.add_argument("--single_prediction", action="store_true")
    parser.add_argument("--profile_filename",
                        help="Only used in single prediction mode")
    parser.add_argument("--metadata_filename",
                        help="Only used in single prediction mode")
    parser.add_argument("--num_core", type=int, default=None,
                        help="Number of parallel processes to use (optional)")

    args = parser.parse_args()

    # Require --combined_outdir if not using --single_prediction
    if not args.single_prediction and not args.combined_outdir:
        log_also_console("❌ Error: --combined_outdir is required "
                         "when not in --single_prediction mode.")
        raise ValueError("--combined_outdir is required"
                         "when not in --single_prediction mode.")

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

        task_args = [
            (script_dir, profile_dir, p_file,
             metadata_dir, m_file, model,
             output_dir, args.name_prefix)
            for p_file, m_file in zip(profile_files, metadata_files)
        ]

        cpu_count = os.cpu_count()
        num_workers = args.num_core or max(
            1,
            min(25, (cpu_count // 2) if cpu_count else 1)
        )
        log_also_console(f"⚙️  Running with {num_workers} parallel workers.")

        parallel_prediction(task_args, num_workers)

        log_also_console("📦 Combining BED files...")
        combine_bed_files(
            input_dir=output_dir,
            output_dir=args.combined_outdir
        )

        log_also_console("🔚 All predictions completed.")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log_also_console(f"❌ Uncaught error: {e}")
        raise
