import argparse
import os
import sys
import pickle
import pandas as pd
from multiprocessing import Pool


from utility import (
    setup_logger,
    build_matched_file_dict,
    define_num_cores,
    load_profile_matrix,
    load_metadata,
    prepare_ranges_df,
    move_metadata_columns_to_ranges,
    combine_bed_files
)


def parse_args():
    parser = argparse.ArgumentParser(description="Run PRIMEloci prediction.")
    parser.add_argument("--script_dir", required=True)
    parser.add_argument("--profile_main_dir", required=True)
    parser.add_argument("--combined_outdir", required=True)
    parser.add_argument("--model_path", required=True)
    parser.add_argument("--log_file", required=True)
    parser.add_argument("--name_prefix", default="PRIMEloci")
    parser.add_argument("--num_core", type=int, default=None,
                        help="Number of cores to use. Defaults to "
                        "min(25, half of available CPU cores).")
    return parser.parse_args()


def predict_worker(matched_entry,
                   model_path,
                   output_dir,
                   name_prefix,
                   log_file=None):
    try:
        with open(model_path, "rb") as f:
            model = pickle.load(f)

        profile_np = load_profile_matrix(matched_entry["profile_path"],
                                         matched_entry["profile_ext"], model)
        metadata_df = load_metadata(matched_entry["metadata_path"],
                                    matched_entry["metadata_ext"])

        if profile_np.shape[0] != len(metadata_df):
            raise ValueError("Mismatch between profile and metadata lengths.")

        y_proba = model.predict_proba(profile_np)

        ranges_df = prepare_ranges_df(metadata_df,
                                      matched_entry['index'],
                                      y_proba)
        columns_to_keep = ["seqnames", "start", "end", "width", "strand",
                           'chrom', 'chromStart', 'chromEnd', 'name', 'score']
        ranges_df = move_metadata_columns_to_ranges(ranges_df,
                                                    metadata_df,
                                                    columns_to_keep)

        output_bed_file = os.path.join(output_dir,
                                       f"{name_prefix}_pred_all_"
                                       f"{matched_entry['index']}.bed")

        for col in ranges_df.columns:
            if pd.api.types.is_numeric_dtype(ranges_df[col]):
                ranges_df[col] = ranges_df[col].apply(lambda x: str(x))
        ranges_df.to_csv(output_bed_file, sep='\t', header=True, index=False)

        print(f"📄 Prediction saved to: {output_bed_file}")
        return {"status": "success", "sample": matched_entry["index"]}

    except Exception as e:
        error_msg = f"[ERROR] Prediction failed for {matched_entry['index']}: {e}"
        print(error_msg)
        if log_file is not None:
            with open(log_file, "a") as logf:
                logf.write(error_msg + "\n")
        return {"status": "error",
                "sample": matched_entry["index"],
                "error": str(e)}


def run_parallel_prediction(args):
    profile_dir = os.path.join(args.profile_main_dir,
                               "profiles_subtnorm")
    metadata_dir = os.path.join(args.profile_main_dir,
                                "metadata")
    output_dir = os.path.join(args.profile_main_dir,
                              "predictions")
    matched_profile = build_matched_file_dict(profile_dir, metadata_dir,
                                              profile_suffix="_profiles_subtnorm",
                                              metadata_suffix="_metadata")

    if not matched_profile:
        print("❌ No profile files found.", file=sys.stderr)
        sys.exit(1)

    print(f"🔍 Found {len(matched_profile)} profile files")
    print(f"🧠 Model: {args.model_path}")

    total_jobs = len(matched_profile)
    def_num_core = define_num_cores(args.num_core, total_jobs)

    print(matched_profile)

    task_list = [
        ({**entry, "index": index_name},
         args.model_path, output_dir,
         args.name_prefix, args.log_file)
        for index_name, entry in matched_profile.items()
    ]

    with Pool(processes=def_num_core) as pool:
        results = pool.starmap(predict_worker, task_list)

    errors = [res for res in results if res["status"] == "error"]

    if errors:
        print("\n❌ Some predictions failed:\n")
        for err in errors:
            print(f"Sample: {err['sample']}\nError: {err['error']}\n")

        if args.log_file:
            with open(args.log_file, "a") as f:
                f.write("\n❌ Some predictions failed:\n")
                for err in errors:
                    f.write(f"Sample: {err['sample']}\nError: {err['error']}\n")

        sys.exit(1)
    else:
        print("✅ All predictions completed successfully.")

    for line in results:
        print(line)


def main():
    args = parse_args()
    logger = setup_logger(args.log_file)
    run_parallel_prediction(args)
    combine_bed_files(input_dir=os.path.join(args.profile_main_dir,
                                             "predictions"),
                      output_dir=args.combined_outdir)


if __name__ == "__main__":
    main()
