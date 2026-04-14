"""
Manually-run tests for verifying there are no (unexpected) accuracy regressions
on some of the example models. These are a little bit slower (e.g., a couple of
minutes instead of a couple of seconds), so they are not currently including
in the automatically run tests.

To run:
pytest test/manual_accuracy.py
"""

from mrpast.result import load_json_pandas
from typing import List, Optional
import glob
import os
import pandas
import shutil
import subprocess

MRPAST = "mrpast"

THISDIR = os.path.dirname(os.path.realpath(__file__))
MODEL_DIR = os.path.join(THISDIR, "..", "examples")
JOBS = 6

# mrpast should be completely DETERMINISTIC, given these parameters. The only thing that
# is not completely deterministic is the timeout handling, but that should be very close
# (e.g., the error should not decrease by much if we get one more solver iteration in before
#  the timeout applies)
DEFAULT_SOLVER_TO = 30  # 30 seconds to solver, so we keep this fast
DEFAULT_SOLVER_REPS = 15  # 15 replicates is a little low, but should be enough
DEFAULT_TIMESLICES = 200  # Necessary for accuracy on complex models!


def scrape_for(output: List[str], prefix: str) -> Optional[str]:
    for line in output:
        if line.startswith(prefix):
            return line.split(prefix)[1]
    return None


def simulate(model: str, out_dir: str, prefix: str, reps: int = 10) -> str:
    if os.path.exists(out_dir):
        print(f"WARNING: {out_dir} already exists; removing before running tests")
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)
    out_args = os.path.join(out_dir, prefix)
    cmd = [
        MRPAST,
        "simulate",
        "-j",
        str(JOBS),
        "-r",
        str(reps),
        "--seq-len",
        str(100000000),
        "-e",
        str(1e-8),
        "-n",
        str(10),
        os.path.join(MODEL_DIR, model),
        out_args,
    ]
    print(f"Running {cmd}")
    subprocess.check_call(cmd)
    return out_args


def process(
    model: str,
    arg_prefix: str,
    reps: int = 20,
    time_slices: int = DEFAULT_TIMESLICES,
    solver_timeout: int = DEFAULT_SOLVER_TO,
) -> str:
    outdir = os.path.dirname(arg_prefix)
    cmd = [
        MRPAST,
        "process",
        "-j",
        str(JOBS),
        "--num-times",
        str(time_slices),
        "-r",
        str(reps),
        "--suffix",
        "test",
        "--out-dir",
        outdir,
        os.path.join(MODEL_DIR, model),
        arg_prefix,
    ]
    print(f"Running command: {cmd}")
    subprocess.check_call(cmd)

    # Now run the solver with a timeout of 1 minute
    cmd = [
        MRPAST,
        "solve",
        "-j",
        str(JOBS),
        "--timeout",
        str(solver_timeout),
    ] + list(glob.glob(os.path.join(outdir, "*.json")))
    print(f"Running command: {cmd}")
    output = subprocess.check_output(cmd).decode("utf-8")
    print(output)
    best_file = scrape_for(
        output.split("\n"), "The output with the highest likelihood is "
    )
    return best_file


def sum_err_by_prefix(result_df: pandas.DataFrame, prefix: str, is_rel: bool) -> float:
    total_err = 0
    for _, row in result_df.iterrows():
        if row["description"].startswith(prefix):
            err = abs(row["final"] - row["ground_truth"])
            if is_rel:
                err = err / row["ground_truth"]
            total_err += err
    return total_err


def check_abs_errors(df, max_epoch, max_migrate, max_coalrate, max_grate, max_admix=-1):
    abs_epoch_err = sum_err_by_prefix(df, "Epoch ", False)
    assert abs_epoch_err > 0
    assert abs_epoch_err <= max_epoch, df
    abs_migrate_err = sum_err_by_prefix(df, "Migration rate ", False)
    assert abs_migrate_err > 0
    assert abs_migrate_err <= max_migrate, df
    abs_coalrate_err = sum_err_by_prefix(df, "Coalescence rate ", False)
    assert abs_coalrate_err > 0
    assert abs_coalrate_err <= max_coalrate, df
    abs_grate_err = sum_err_by_prefix(df, "Growth rate ", False)
    assert abs_grate_err > 0
    assert abs_grate_err <= max_grate, df
    if max_admix >= 0:
        abs_admix_err = sum_err_by_prefix(df, "Admixture proportion ", False)
        assert abs_admix_err > 0
        assert abs_admix_err <= max_admix, df


def test_ooa3_sim():
    MAX_ABS_EPOCH_ERR = 40
    MAX_ABS_MIGRATE_ERR = 1e-4
    MAX_ABS_COALRATE_ERR = 8e-6
    MAX_ABS_GRATE_ERR = 5e-5

    arg_prefix = simulate("ooa_3g09.yaml", "test.ooa3", "test_ooa3")
    best_result = process("ooa_3g09.yaml", arg_prefix)
    assert best_result is not None, "No best result found in solver output!"
    df = load_json_pandas(best_result)

    # Check absolute errors
    check_abs_errors(
        df,
        MAX_ABS_EPOCH_ERR,
        MAX_ABS_MIGRATE_ERR,
        MAX_ABS_COALRATE_ERR,
        MAX_ABS_GRATE_ERR,
    )


def test_aa5_sim():
    MAX_ABS_EPOCH_ERR = 300
    MAX_ABS_MIGRATE_ERR = 5e-5
    MAX_ABS_COALRATE_ERR = 5e-5
    MAX_ABS_GRATE_ERR = 0.0005
    MAX_ABS_ADMIX_ERR = 0.02

    arg_prefix = simulate("aa5.yaml", "test.aa5", "test_aa5")
    best_result = process("aa5.yaml", arg_prefix)
    assert best_result is not None, "No best result found in solver output!"
    df = load_json_pandas(best_result)
    df.to_csv("tmp.csv", index=False)

    # Check absolute errors
    check_abs_errors(
        df,
        MAX_ABS_EPOCH_ERR,
        MAX_ABS_MIGRATE_ERR,
        MAX_ABS_COALRATE_ERR,
        MAX_ABS_GRATE_ERR,
        MAX_ABS_ADMIX_ERR,
    )
