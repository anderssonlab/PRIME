import os
import subprocess

# 🔧 CONFIGURATION
script_dir = "/Users/natsudanav/Desktop/PRIME/python/"
profile_main_dir_pq = (
    "/Users/natsudanav/Desktop/PRIME/tests/testoutput/"
    "PRIMEloci_tmp/PRIMEloci_profiles"
)
profile_main_dir_npz = "/Users/natsudanav/Desktop/testoutput23Mar/PRIMEloci_tmp/PRIMEloci_profiles"
combine_outdir = "/Users/natsudanav/Desktop/testoutput23Mar/"
model_path = "/Users/natsudanav/Documents/data_PRIMEloci_dev/model_PRIMEloci/PRIMEloci_GM12878_model_1.0.sav"
log_file = combine_outdir+"PRIMEloci_log/PRIMEloci-5.log"


# 🟨 Test #1: Single prediction (.npz)
single_npz_profile = "K562_C1test__chr1_profiles_subtnorm.npz"
single_metadata = "K562_C1test__chr1_metadata.parquet"

single_cmd_npz = [
    "python3", os.path.join(script_dir, "main_prediction.py"),
    "--script_dir", script_dir,
    "--profile_main_dir", profile_main_dir_npz,
    "--model_path", model_path,
    "--log_file", log_file,
    "--name_prefix", "single_npz",
    "--single_prediction",
    "--profile_filename", single_npz_profile,
    "--metadata_filename", single_metadata
]


# 🟨 Test #2: Single prediction (.parquet)
single_parquet_profile = "K562_C1_test_chr1_profiles_subtnorm.parquet"
single_metadata_2 = "K562_C1_test_chr1_metadata.parquet"

single_cmd_parquet = [
    "python3", os.path.join(script_dir, "main_prediction.py"),
    "--script_dir", script_dir,
    "--profile_main_dir", profile_main_dir_pq,
    "--model_path", model_path,
    "--log_file", log_file,
    "--name_prefix", "single_parquet",
    "--single_prediction",
    "--profile_filename", single_parquet_profile,
    "--metadata_filename", single_metadata_2
]


# 🟨 Test #3: Batch mode (auto match)
batch_cmd = [
    "python3", os.path.join(script_dir, "main_prediction.py"),
    "--script_dir", script_dir,
    "--profile_main_dir", profile_main_dir_npz,
    "--combined_outdir", combine_outdir,
    "--model_path", model_path,
    "--log_file", log_file,
    "--name_prefix", "batch_run",
    "--num_core", "2"
]


def run_case(name, cmd):
    print(f"\n🔹 Running: {name}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    print(result.stdout)

    if result.stderr:
        print("⚠️ STDERR:", result.stderr)

    if result.returncode != 0:
        print(f"❌ {name} failed with exit code {result.returncode}")
    else:
        print(f"✅ {name} completed successfully.")

    print("🔸 Done.\n" + "-" * 50)


# run_case("Single Prediction (.npz)", single_cmd_npz)
# run_case("Single Prediction (.parquet)", single_cmd_parquet)


run_case("Batch Prediction (all matched)", batch_cmd)
