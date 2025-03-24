import os
import subprocess

# 🔧 CONFIGURATION
script_dir = "/Users/natsudanav/Desktop/PRIME/python/main_prediction.py"
profile_main_dir = "/Users/natsudanav/Desktop/PRIME/tests/testoutput/PRIMEloci_tmp/PRIMEloci_profiles"
model_path = "/Users/natsudanav/Documents/data_PRIMEloci_dev/model_PRIMEloci/PRIMEloci_GM12878_model_1.0.sav"
log_file = "test_prediction.log"

# 🟨 Test #1: Single prediction (.npz)
single_npz_profile = "test_sample_01_profiles_subtnorm.npz"
single_metadata = "test_sample_01_metadata.parquet"

single_cmd_npz = [
    "python", os.path.join(script_dir, "main_prediction.py"),
    "--script_dir", script_dir,
    "--profile_main_dir", profile_main_dir,
    "--model_path", model_path,
    "--log_file", log_file,
    "--name_prefix", "single_npz",
    "--single_prediction",
    "--profile_filename", single_npz_profile,
    "--metadata_filename", single_metadata
]

# 🟨 Test #2: Single prediction (.parquet)
single_parquet_profile = "test_sample_02_profiles_subtnorm.parquet"
single_metadata_2 = "test_sample_02_metadata.parquet"

single_cmd_parquet = [
    "python", os.path.join(script_dir, "main_prediction.py"),
    "--script_dir", script_dir,
    "--profile_main_dir", profile_main_dir,
    "--model_path", model_path,
    "--log_file", log_file,
    "--name_prefix", "single_parquet",
    "--single_prediction",
    "--profile_filename", single_parquet_profile,
    "--metadata_filename", single_metadata_2
]

# 🟨 Test #3: Batch mode (auto match)
batch_cmd = [
    "python", os.path.join(script_dir, "main_prediction.py"),
    "--script_dir", script_dir,
    "--profile_main_dir", profile_main_dir,
    "--model_path", model_path,
    "--log_file", log_file,
    "--name_prefix", "batch_run"
]

# 🚀 RUN TESTS
def run_case(name, cmd):
    print(f"🔹 Running: {name}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.stderr:
        print("⚠️ STDERR:", result.stderr)

run_case("Single Prediction (.npz)", single_cmd_npz)
run_case("Single Prediction (.parquet)", single_cmd_parquet)
run_case("Batch Prediction (all matched)", batch_cmd)
