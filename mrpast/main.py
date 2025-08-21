# Migration Rate and Population Size Across Space and Time (mrpast)
# Copyright (C) 2025 April Wei
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# with this program.  If not, see <https://www.gnu.org/licenses/>.
from collections import defaultdict
from copy import deepcopy
from enum import Enum
from multiprocessing import Pool
from typing import Dict, List, Tuple, Optional, Set
from numpy.typing import NDArray

import argparse
import glob
import json
import os
import pandas as pd
import re
import subprocess
import sys
import time
import uuid

try:
    import demes
except ImportError:
    demes = None  # type: ignore

from mrpast.polarize import polarize_vcf
from mrpast.from_demes import convert_from_demes
from mrpast.helpers import (
    dump_model_yaml,
    get_best_output,
    make_zarr,
    remove_ext,
    which,
    load_old_mrpast,
    one_file_or_one_per_chrom,
)
from mrpast.simulate import (
    run_simulation,
    build_demography,
)
from mrpast.arginfer import (
    infer_arg,
    ArgTool,
    DEFAULT_SAMPLE_SUFFIX,
    DEFAULT_NUM_SAMPLES as DEFAULT_ARG_SAMPLES,
)
from mrpast.model import (
    UserModel,
    ModelSolverInput,
    print_model_warnings,
)
from mrpast.simprocess import (
    BootstrapCoalSampler,
    JackknifeCoalSampler,
    StandardCoalSampler,
    get_coalescence_list,
    get_time_slices,
    merge_coals,
    sample_coal_matrices,
    simulate_muts_and_export,
)
from mrpast.result import (
    load_json_pandas,
    tab_show,
)

DEFAULT_INDIVIDUALS = 10  # per population
DEFAULT_MAX_GENERATION = 1_000_000.0
DEFAULT_MAX_SUBSETS = 1
DEFAULT_MIN_TIME_UNIT = 1.0
DEFAULT_MUT_RATE = 1.2e-8
DEFAULT_RANDOM_SEED = 42
DEFAULT_RECOMB_RATE = 1e-8
DEFAULT_SEQ_LENGTH = 100_000_000
DEFAULT_SIM_REPLICATES = 20
DEFAULT_SOLVE_REPS_PER_EPOCH = 10
DEFAULT_TIME_SLICES = 20
DEFAULT_TREE_SAMPLE_RATE = 125_000
# Everything up to the sample identifier is what we group by. I.e., this groups "across" all
# samples for the same prefix (such as chromosome identifier). For example, if you run on
# chromosomes 1 and 2 and get 100 ARG samples, then this will produce two groups (one for each
# chromosome) of 100 samples each.
CHROM_GROUP_BY = f"(.+)\.vcf_.*{DEFAULT_SAMPLE_SUFFIX}"
# This groups by sample, for when the coal matrices should summarize all chromosomes across
# the ARG samples.
SAMPLE_GROUP_BY = f"{DEFAULT_SAMPLE_SUFFIX}([0-9]+)(_v[0-9]+)?.trees"

CMD_SIMULATE = "simulate"
CMD_PROCESS = "process"
CMD_SOLVE = "solve"
CMD_SIM2VCF = "sim2vcf"
CMD_ARGINFER = "arginfer"
CMD_MODEL = "model"
CMD_INIT = "init"
CMD_CONFIDENCE = "confidence"
CMD_POLARIZE = "polarize"
CMD_SHOW = "show"
CMD_SELECT = "select"


class BootstrapOpt(Enum):
    none = ""
    coalcounts = "coalcounts"
    jackknife = "jackknife"

    def __str__(self):
        return self.value


def time_call(command: List[str], **kwargs) -> float:
    start_time = time.time()
    subprocess.check_call(command, **kwargs)
    return time.time() - start_time


solver_exe = which("mrp-solver", required=True)
eval_exe = which("mrp-eval", required=True)


def convert_simulation(
    arg_file: str,
    mut_rate: float,
    leave_out: Set[int] = set(),
    seed: int = DEFAULT_RANDOM_SEED,
    emit_zarr: bool = False,
) -> str:
    """
    Convert simulation results (a tree-sequence) into data that can be used for
    ARG inference.

    :param arg_file: The tree-sequence filename.
    :param mut_rate: The mutation rate. If the simulation has already simulated mutations
        then this should be set to 0.0.
    :param leave_out: The set of population indexes to leave out of the conversion. If there
        are D populations then this number should be between 0 and (D-1).
    :param seed: The random seed to use for simulating mutations.
    :param emit_zarr: Set to True to emit ZARR/VCF for use with tsinfer.
    :return: The VCF (or ZARR/VCF) filename.
    """
    pop_map = simulate_muts_and_export(
        arg_file,
        arg_file + ".vcf",
        seed=seed,
        mut_rate=mut_rate,
        leave_out_pops=leave_out,
    )
    pop_map_file = f"{arg_file}.popmap.json"
    with open(pop_map_file, "w") as f:
        f.write(pop_map.to_json(indent=2))
    vcf_file = f"{arg_file}.vcf"
    if emit_zarr:
        out_file = make_zarr(vcf_file, delete_orig=True)
    else:
        out_file = vcf_file
    print(f"Wrote {out_file} and {pop_map_file}")
    return out_file


def solve_single(
    input_file: str, timeout: Optional[float] = None, matrix: Optional[int] = None
) -> Tuple[str, float]:
    if matrix is not None:
        suffix = f".mat{matrix}"
    else:
        suffix = ""
    output_file = ".".join(input_file.split(".")[:-1]) + suffix + ".out.json"
    cmd = [solver_exe, input_file, output_file]
    if timeout is not None:
        cmd.extend(["--timeout", str(timeout)])
    if matrix is not None:
        cmd.extend(["--select-matrix", str(matrix)])
    elapsed = time_call(cmd)
    return output_file, elapsed


def solve(
    solver_inputs: List[str],
    jobs: int = 1,
    timeout: Optional[float] = None,
    specific_matrices: List[int] = [],
) -> List[Tuple[str, float]]:
    assert (
        solver_exe is not None
    ), "No mrp-solver executable found. Invalid installation or bad $PATH."
    arguments = []
    for matrix in specific_matrices or [None]:  # type: ignore
        for input_file in solver_inputs:
            arguments.append((input_file, timeout, matrix))
    outputs = []
    if len(solver_inputs) == 1 or jobs == 1:
        for arg_list in arguments:
            outputs.append(solve_single(*arg_list))
    else:
        with Pool(jobs) as pool:
            outputs = pool.starmap(solve_single, arguments)
    return outputs


def get_coaldist_from_arg(
    arg_prefix,
    jobs=1,
    leave_out_pops=[],
    min_time_unit=1.0,
    tree_sample_rate: int = DEFAULT_TREE_SAMPLE_RATE,
    rate_maps: Optional[str] = None,
    rate_map_threshold: float = 1.0,
) -> List[str]:
    # Load the ARG from the tree-sequence(s) and bin the coalescence times.
    glb = f"{arg_prefix}*.trees"
    tree_files = list(sorted(glob.glob(glb)))
    if not tree_files:
        print("No tree files matched glob {glb}!", file=sys.stderr)
        exit(1)
    print(f"Found {len(tree_files)} tree files")

    expanded_rate_maps: List[Optional[str]] = []
    if rate_maps is not None:
        # Group by "chromosome" (for simulation, this is just numbered replicated, but you can think of them
        # as being independently evolving chromosomes under the same demographic constraints). Once we have
        # the group we can match the rate_maps to the group.
        group_regex = re.compile(CHROM_GROUP_BY)
        grouped = defaultdict(list)
        for fn in tree_files:
            m = group_regex.search(fn)
            if m is None:
                grouped[fn].append(fn)
            else:
                group_id = m.group(1)
                grouped[group_id].append(fn)

        group_ids = list(sorted(grouped.keys()))
        rate_map_list = one_file_or_one_per_chrom(
            rate_maps, group_ids, ".txt", desc="recombination map"
        )
        assert len(rate_map_list) == len(grouped)

        tree_files = []
        expanded_rate_maps = []
        for group_key, rate_map in zip(group_ids, rate_map_list):
            for tree_file in grouped[group_key]:
                tree_files.append(tree_file)
                expanded_rate_maps.append(rate_map)
    else:
        expanded_rate_maps = [None for _ in tree_files]
    params = [
        [
            fn,
            f"{fn}-coal.txt",
            tree_sample_rate,
            min_time_unit,
            leave_out_pops,
            rm,
            rate_map_threshold,
        ]
        for fn, rm in zip(tree_files, expanded_rate_maps)
    ]
    if jobs == 1:
        coal_filenames = [get_coalescence_list(*p) for p in params]
    else:
        with Pool(jobs) as pool:
            coal_filenames = pool.starmap(get_coalescence_list, params)
    return coal_filenames


# Turn fine-grained coalescence counts into a summary coal matrix (time x state).
# A solver input can have many coal matrices, each one representing a sample. There
# are different ways to get a sampled coal matrix:
# * From a single ARG per chromosome (either simulated or one sample from inference):
#   * Bootstrap: sample tree coal counts with replacement to produce I matrices.
#   * Jackknife: sample blocks of trees to produce I matrices (where each has 1 block left out)
#   * "none": just produce a single coal matrix
# * From multiple ARGs per chromosome (N samples from inference):
#   * "none": produces N coal matrices -- these are tagged as "ARG samples"
#   * Bootstrap and Jackknife: sum the coal counts across N ARGs first, and then generate I
#     samples as per normal.
#
# The above is implemented by a combination of the "bootstrap" and "group_by" flags. "group_by"
# produces N>=1 groups that correspond to the ARG samples.
def get_coal_counts(
    model_file: str,
    grouped_filenames: Dict[str, List[str]],
    time_slices: List[float],
    bootstrap: BootstrapOpt,
    bootstrap_iters: int,
    jobs: int = 1,
    seed: int = DEFAULT_RANDOM_SEED,
) -> Tuple[List[NDArray], str, List[str]]:
    """
    Produce a list of coal matrices that will be used as part of the solver input.

    :return: A tuple (matrices, description) where matrices is a list of coal matrices
        and description is a string describing what the samples mean (ARG samples, boostrapped, etc.)
    """
    # Load the configuration and create the symbolic model from it.
    model = UserModel.from_file(model_file)

    description = None
    if bootstrap == BootstrapOpt.none:
        sampler = StandardCoalSampler()
        description = "arg_samples"
    elif bootstrap == BootstrapOpt.coalcounts:
        sampler = BootstrapCoalSampler(bootstrap_iters)
        description = "bootstrap"
    elif bootstrap == BootstrapOpt.jackknife:
        sampler = JackknifeCoalSampler(bootstrap_iters)
        description = "jackknife"
    else:
        assert False, f"Unrecognized bootstrap option {bootstrap}"

    number_of_groups = len(grouped_filenames)
    deme_pair_index0 = model.get_pair_ordering()
    coal_matrices = []
    if number_of_groups > 1:
        print(f"Using {number_of_groups} ARG sample with sampling method {description}")
        if bootstrap == BootstrapOpt.none:
            for _, group_files in sorted(grouped_filenames.items(), key=lambda t: t[0]):
                group_sampler = deepcopy(sampler)
                sample_coal_matrices(
                    group_files,
                    time_slices,
                    deme_pair_index0,
                    group_sampler,
                    jobs=jobs,
                    seed=seed,
                )
                assert len(group_sampler.coal_matrices) == 1
                coal_matrices.append(group_sampler.coal_matrices[0])
        else:
            coal_files = []
            for group, grouped_files in grouped_filenames.items():
                merged_file = f"{group}-avg-coal.txt"
                print(f"Merging {grouped_files} -> {merged_file}")
                merge_coals(grouped_files, merged_file)
                coal_files.append(merged_file)
            sample_coal_matrices(
                coal_files,
                time_slices,
                deme_pair_index0,
                sampler,
                jobs=jobs,
                seed=seed,
            )
            coal_matrices = list(sampler.coal_matrices)
    else:
        print(f"Using a single ARG sample with sampling method {description}")
        coal_files = list(grouped_filenames.values())[0]
        sample_coal_matrices(
            coal_files,
            time_slices,
            deme_pair_index0,
            sampler,
            jobs=jobs,
            seed=seed,
        )
        coal_matrices = list(sampler.coal_matrices)
    return coal_matrices, description, sampler.sample_hashes()


def generate_solver_input(
    model_file,
    coal_count_matrices: List[NDArray],
    time_slices,
    replicates: Optional[int] = None,
    sampling_description: Optional[str] = None,
    sampling_hashes: Optional[List[str]] = None,
    generate_ground_truth: bool = False,
) -> Tuple[List[str], List[str]]:
    # Load the configuration and create the symbolic model from it.
    user_model = UserModel.from_file(model_file)
    solver_input = user_model.to_solver_model(
        generate_ground_truth=generate_ground_truth
    )
    # Add the coalescence and time discretization information to the solver input.
    solver_input.coal_count_matrices = coal_count_matrices
    solver_input.time_slices_gen = time_slices
    solver_input.sampling_description = sampling_description
    solver_input.sampling_hashes = sampling_hashes

    solver_in_ids = []
    solver_in_data = []
    # If requested, we create an extra solver input that has all parameters initialized with
    # their ground truth -- this typically lets the solver find the true minima for comparison,
    # whereas the random init replicants may get stuck in other local minima.
    if generate_ground_truth:
        solver_in_data.append(solver_input.to_json(indent=2))
        solver_in_ids.append("truth")

    if replicates is None:
        replicates = DEFAULT_SOLVE_REPS_PER_EPOCH * user_model.num_epochs
    assert replicates is not None
    for i in range(replicates):
        repl_input = solver_input.randomize()
        solver_in_data.append(repl_input.to_json(indent=2))
        solver_in_ids.append(str(i))
    return solver_in_ids, solver_in_data


def process_solver_outputs(
    output_data: List[Tuple[str, float]],
    verbose: bool = False,
) -> str:
    """
    Find and return the best output from the solver.
    """
    WORST_NEGLL = 2**64
    best_negLL = WORST_NEGLL
    worst_negLL = 0
    avg_elapsed = 0.0
    avg_negLL = 0.0
    best_output = None
    assert len(output_data) > 0, "Must pass at least one output"
    for filename, elapsed_time in output_data:
        with open(filename) as f:
            result = json.load(f)
        avg_elapsed += elapsed_time
        negLL = result.get("negLL") or WORST_NEGLL
        avg_negLL += negLL
        if negLL < best_negLL:
            best_negLL = negLL
            best_output = filename
        if negLL > worst_negLL:
            worst_negLL = negLL
    avg_negLL /= len(output_data)
    avg_elapsed /= len(output_data)
    if verbose:
        print(f"Worst (negative log) likelihood: {worst_negLL}")
        print(f"Best (negative log) likelihood: {best_negLL}")
        print(f"Average (negative log) likelihood: {avg_negLL}")
        print(f"Average solve time: {avg_elapsed} seconds")
    assert best_output is not None
    return best_output


def add_common(parser):
    parser.add_argument(
        "--jobs",
        "-j",
        type=int,
        default=1,
        help="Number of jobs (threads) to use. Defaults to 1.",
    )
    parser.add_argument(
        "--seed", type=int, default=DEFAULT_RANDOM_SEED, help="Set the random seed."
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Verbose output, including timing information.",
    )


def parse_intlist(argument: str) -> List[int]:
    try:
        return [int(a) for a in argument.split(",") if a]
    except ValueError:
        raise RuntimeError(f"Invalid integer list: {argument}")


def suffixed_int(arg_value):
    assert len(arg_value) > 0
    try:
        return (int(arg_value), None)
    except ValueError:
        suffix = arg_value[-1]
        return (int(arg_value[:-1]), suffix)


def float_or_filename(arg_value):
    try:
        return float(arg_value)
    except ValueError:
        rv = str(arg_value)
        assert os.path.isfile(rv), f"Invalid filename: {rv}"
        return rv


def float_or_str(arg_value):
    try:
        return float(arg_value)
    except ValueError:
        return str(arg_value)


def time_slice_list(time_slice_str: Optional[str]) -> Tuple[List[float], bool]:
    if time_slice_str is None:
        return ([], False)
    is_extra = False
    if time_slice_str.startswith("+"):
        time_slice_str = time_slice_str[1:]
        is_extra = True
    time_list = list(map(float, time_slice_str.split(",")))
    return (time_list, is_extra)


def process_ARGs(
    model: str,
    arg_prefix: str,
    suffix: str,
    jobs: int,
    do_solve: bool = False,
    out_dir: str = ".",
    add_ground_truth: bool = False,
    num_times: int = DEFAULT_TIME_SLICES,
    replicates: Optional[int] = None,
    bootstrap: BootstrapOpt = BootstrapOpt.none,
    bootstrap_iter: int = 100,
    max_generation: float = DEFAULT_MAX_GENERATION,
    tree_sample_rate: int = DEFAULT_TREE_SAMPLE_RATE,
    min_time_unit: float = DEFAULT_MIN_TIME_UNIT,
    leave_out: str = "",
    group_by: Optional[str] = None,
    time_slice_str: Optional[str] = None,
    verbose: bool = False,
    rate_maps: Optional[str] = None,
    rate_map_threshold: float = 1.0,
    left_skew_times: bool = False,
    seed: int = DEFAULT_RANDOM_SEED,
) -> List[str]:
    """
    Given a set of ARGs, extract the pair-wise coalescence information and turn it into a concrete
    model input. Optionally, run the maximum likelihood solver to get estimated parameter values.
    """
    model_base = ".".join(os.path.basename(model).split(".")[:-1])
    suffix = str(uuid.uuid4())[:8] if suffix is None else suffix
    generic_prefix = f"{model_base}.{suffix}"
    solve_in_prefix = f"{generic_prefix}.solve_in"

    # Step 1: collect coalescence distributions from the ARG, or load it from a previously
    # generated JSON file.
    coal_filenames = get_coaldist_from_arg(
        arg_prefix,
        jobs=jobs,
        leave_out_pops=parse_intlist(leave_out),
        min_time_unit=min_time_unit,
        tree_sample_rate=tree_sample_rate,
        rate_maps=rate_maps,
        rate_map_threshold=rate_map_threshold,
    )
    if verbose:
        print(f"Wrote coalescences to {coal_filenames}")

    unmatched = []
    groups = defaultdict(list)
    if group_by is None:
        if bootstrap == BootstrapOpt.none:
            group_by = SAMPLE_GROUP_BY
        else:
            group_by = CHROM_GROUP_BY
    group_regex = re.compile(group_by)
    for fn in coal_filenames:
        m = group_regex.search(fn)
        if m is None:
            unmatched.append(fn)
        else:
            group_id = m.group(1)
            groups[group_id].append(fn)
    if unmatched and len(groups) > 0:
        print(f"Some coal files matched --group-by regex, others did not.")
        exit(1)
    if unmatched:
        print(f"WARNING: No matches for --group-by {group_by}", file=sys.stderr)
        groups["all"] = unmatched

    # Get time slices and coal matrices - this may produce many matrices if we are using
    # a sampling strategy like bootstrap or jackknife.
    ts_list, ts_is_extra = time_slice_list(time_slice_str)
    if not ts_list or ts_is_extra:
        time_slices = get_time_slices(
            coal_filenames, num_times, max_generation, left_skewed=left_skew_times
        )
        time_slices = sorted(time_slices + ts_list)
    else:
        time_slices = ts_list

    coal_matrices, sampling_description, sampling_hashes = get_coal_counts(
        model,
        groups,
        time_slices,
        bootstrap,
        bootstrap_iter,
        jobs=jobs,
        seed=seed,
    )
    # Step 2: generate the solver input file.
    solver_input_ids, solver_input_data = generate_solver_input(
        model,
        coal_matrices,
        time_slices,
        replicates=replicates,
        sampling_description=sampling_description,
        sampling_hashes=sampling_hashes,
        generate_ground_truth=add_ground_truth,
    )

    solver_inputs = []
    for ident, data in zip(solver_input_ids, solver_input_data):
        filename = os.path.join(
            out_dir, f"{solve_in_prefix}.{sampling_description}.{ident}.json"
        )
        with open(filename, "w") as f:
            f.write(data)
        solver_inputs.append(filename)
    print(f"Created solver inputs: {solver_inputs}")

    if do_solve:
        solver_outputs = solve(solver_inputs, jobs)
        print(f"Created outputs: {[fn for (fn, _) in solver_outputs]}")
        best_output = process_solver_outputs(solver_outputs, verbose=True)
        print(f"The output with the highest likelihood is {best_output}")
    return solver_inputs


def main():
    parser = argparse.ArgumentParser(
        description="Infer demographic history using ARGs."
    )
    subparsers = parser.add_subparsers(dest="command")

    simulate_parser = subparsers.add_parser(
        CMD_SIMULATE,
        help="Simulate data from your model in order to validate the model performance.",
    )
    add_common(simulate_parser)
    simulate_parser.add_argument(
        "model", help="The input YAML file specifying the model"
    )
    simulate_parser.add_argument(
        "--replicates",
        "-r",
        type=int,
        default=DEFAULT_SIM_REPLICATES,
        help=f"Number of simulation replications to perform. Defaults to {DEFAULT_SIM_REPLICATES}.",
    )
    simulate_parser.add_argument(
        "arg_prefix", help="The prefix for the output tree-sequence files"
    )
    simulate_parser.add_argument(
        "--seq-len",
        "-s",
        type=int,
        default=DEFAULT_SEQ_LENGTH,
        help=f"Length of sequences in base-pairs. Default to {DEFAULT_SEQ_LENGTH}.",
    )
    simulate_parser.add_argument(
        "--recomb-rate",
        "-e",
        type=float_or_str,
        default=DEFAULT_RECOMB_RATE,
        help=f"Rate of recombination, or filename/prefix for recombination map. A prefix will match '<prefix>*.txt'. Defaults to {DEFAULT_RECOMB_RATE}.",
    )
    simulate_parser.add_argument(
        "--individuals",
        "-n",
        type=int,
        default=DEFAULT_INDIVIDUALS,
        help=f"Number of individuals per population. Defaults to {DEFAULT_INDIVIDUALS}.",
    )
    simulate_parser.add_argument(
        "--debug-demo",
        "-d",
        action="store_true",
        help="Output results from msprime demography debugger.",
    )

    process_parser = subparsers.add_parser(
        CMD_PROCESS, help="Process an ARG to generate solver input."
    )
    add_common(process_parser)
    process_parser.add_argument(
        "model", help="The input YAML file specifying the model"
    )
    process_parser.add_argument(
        "arg_prefix",
        help="The prefix of the input tree-seq file(s) specifying the ARG. Assumes .trees file extension.",
    )
    process_parser.add_argument(
        "--replicates",
        "-r",
        type=int,
        default=None,
        help=f"Number of solver replications to perform. Defaults to {DEFAULT_SOLVE_REPS_PER_EPOCH} * num_epochs.",
    )
    process_parser.add_argument(
        "--num-times",
        "-t",
        type=suffixed_int,
        default=(DEFAULT_TIME_SLICES, None),
        help=f"Number of time slices to use. Defaults to {DEFAULT_TIME_SLICES}. Use the suffix 'l' or 'L' to use left-skewed time slices.",
    )
    process_parser.add_argument(
        "--solve",
        "-s",
        action="store_true",
        help="Solve the model after generating the solver inputs.",
    )
    process_parser.add_argument(
        "--add-ground-truth",
        "-g",
        action="store_true",
        help="Generate an additional solver input(s) using the ground-truth parameter values.",
    )
    process_parser.add_argument(
        "--suffix",
        help="Filenames will use the provided suffix instead of a random one.",
    )
    process_parser.add_argument(
        "--out-dir", "-o", help="Output directory.", default="."
    )
    process_parser.add_argument(
        "--min-time-unit",
        "-u",
        default=DEFAULT_MIN_TIME_UNIT,
        type=float,
        help=f"The minimum time unit for distinguishing between coalescence events (default: {DEFAULT_MIN_TIME_UNIT} generation)",
    )
    process_parser.add_argument(
        "--max-generation",
        "-m",
        default=DEFAULT_MAX_GENERATION,
        type=float,
        help=f"Ignore all coalescence events occuring after the given generation (default: {DEFAULT_MAX_GENERATION})",
    )
    process_parser.add_argument(
        "--tree-sample-rate",
        "-b",
        default=DEFAULT_TREE_SAMPLE_RATE,
        type=int,
        help=f"Sample a tree from the ARG every tree-sample-rate base pairs (default: {DEFAULT_TREE_SAMPLE_RATE} bp)",
    )
    process_parser.add_argument(
        "--leave-out",
        default="",
        type=str,
        help=f"Comma-separated list of population IDs to leave out when counting coalescence",
    )
    process_parser.add_argument(
        "--bootstrap",
        type=BootstrapOpt,
        default=BootstrapOpt.none,
        choices=list(BootstrapOpt),
        help=f"Bootstrap the sampled trees to create more than once coalescent matrix.\n"
        "  coalcounts: standard bootstrap of over marginal trees."
        "  jackknife: leave-one-out jacktree over blocks of marginal trees.",
    )
    process_parser.add_argument(
        "--bootstrap-iter",
        "-i",
        default=100,
        type=int,
        help=f"How many blocks to split the trees in for jackknifing, number of reps for standard bootstrap. Default: 100",
    )
    process_parser.add_argument(
        "--group-by",
        default=None,
        help=f"Regex to group ARGs or coal files by. By default group by chromosome for bootstrapping and "
        "by sample otherwise.",
    )
    process_parser.add_argument(
        "--time-slices",
        default=None,
        help=f"The comma-separated list of time slice values instead of computing them from coalescence counts. "
        "Or, if prefixed with '+', the list of time slices to append to the auto-generated time slices.",
    )
    process_parser.add_argument(
        "--rate-maps",
        default=None,
        help=f"A filename prefix for tskit-style RateMap files, whose lexicographic sort order matches the input ARGs "
        "lexicographic sort order. Generates a glob '<prefix>*.txt'. Used for determining tree sampling (see --rate-map-threshold)",
    )
    process_parser.add_argument(
        "--rate-map-threshold",
        default=1e-9,
        type=float,
        help=f"Only sample trees from regions with a recombination rate <= to this. Requires --rate-maps",
    )

    solve_parser = subparsers.add_parser(
        CMD_SOLVE, help="Infer the demographic parameters of the model."
    )
    solve_parser.add_argument(
        "solver_inputs",
        nargs="+",
        help="The solver input JSON files. The output filenames will be derived from the input filenames.",
    )
    solve_parser.add_argument(
        "--timeout",
        default=None,
        type=float,
        help="Timeout in seconds. Solver returns the current best result upon timeout.",
    )
    add_common(solve_parser)

    sim2vcf_parser = subparsers.add_parser(
        CMD_SIM2VCF, help="Convert a simulation result to VCF and a population map."
    )
    sim2vcf_parser.add_argument("arg_file", help="The ARG (.trees) file to process.")
    sim2vcf_parser.add_argument(
        "--prefix",
        "-p",
        action="store_true",
        help="Treat arg_file as a prefix, and search for all <arg_prefix>*.trees files",
    )
    sim2vcf_parser.add_argument(
        "--leave-out",
        default="",
        type=str,
        help=f"Comma-separated list of population IDs to leave out when converting to VCF",
    )
    sim2vcf_parser.add_argument(
        "--mut-rate",
        default=DEFAULT_MUT_RATE,
        type=float,
        help=f"The mutation rate, for simulating mutations on existing trees.",
    )
    sim2vcf_parser.add_argument(
        "--zarr",
        "-z",
        action="store_true",
        help=f"Output VCF/ZARR files, required for tsinfer usage.",
    )
    add_common(sim2vcf_parser)

    arginfer_parser = subparsers.add_parser(
        CMD_ARGINFER, help="Infer an ARG from a VCF file."
    )
    arginfer_parser.add_argument(
        "vcf_prefix",
        help='The prefix of VCF file(s) to process. Generates a glob "<vcf_prefix>*.vcf"',
    )
    arginfer_parser.add_argument(
        "arg_prefix",
        help="The prefix to use when writing the resulting ARGs to disk (.trees files)",
    )
    arginfer_parser.add_argument(
        "pop_map", help="The file containing the population map (*.popmap.json)"
    )
    arginfer_parser.add_argument(
        "--ne-override",
        "-N",
        default="auto",
        type=float_or_str,
        help="Provide an override for the auto-calculated (diploid) effective population size.",
    )
    arginfer_parser.add_argument(
        "--mut-rate",
        "-m",
        default=DEFAULT_MUT_RATE,
        type=float,
        help=f"Expected mutation rate. Default {DEFAULT_MUT_RATE}.",
    )
    arginfer_parser.add_argument(
        "--recomb-rate",
        "-r",
        default=DEFAULT_RECOMB_RATE,
        type=float_or_str,
        help=f"Expected recombination rate, or recombination map filename. Default {DEFAULT_RECOMB_RATE}.",
    )
    arginfer_parser.add_argument(
        "--samples",
        "-s",
        type=int,
        default=DEFAULT_ARG_SAMPLES,
        help=f"How many ARGS to sample. Default {DEFAULT_ARG_SAMPLES}.",
    )
    arginfer_parser.add_argument(
        "--thin",
        "-t",
        type=int,
        default=None,
        help=f"How many MC/MC iterations between samples. Default depends on the inference tool.",
    )
    arginfer_parser.add_argument(
        "--dry-run",
        "-d",
        action="store_true",
        help=f"Just emit the arguments that would be used when running SINGER.",
    )
    arginfer_parser.add_argument(
        "--tool",
        type=ArgTool,
        choices=list(ArgTool),
        default=ArgTool.RELATE,
        help='Which ARG inference tool to run: "tsinfer" (default), "relate", or "singer"',
    )
    arginfer_parser.add_argument(
        "--ancestral",
        "-a",
        default=None,
        help="The ancestral FASTA file (input). Assumes the positions start counting at 1.",
    )

    add_common(arginfer_parser)

    model_parser = subparsers.add_parser(
        CMD_MODEL, help="Process an input model and produce information about it"
    )
    model_parser.add_argument("model", help="The model YAML file")
    model_parser.add_argument(
        "--to-demes",
        "-d",
        type=str,
        help="Write a Demes YAML file representing the model.",
    )
    model_parser.add_argument(
        "--debug",
        action="store_true",
        help="Emit msprime demography debugger output for the given model",
    )

    create_parser = subparsers.add_parser(
        CMD_INIT, help="Create an initial model, to be edited by the user."
    )
    init_choices = create_parser.add_mutually_exclusive_group()
    init_choices.add_argument(
        "--from-demes",
        "-d",
        type=str,
        help="Convert a Demes YAML file into a mrpast model.",
    )
    init_choices.add_argument(
        "--from-old-mrpast",
        type=str,
        help="Convert an old matrix-based mrpast YAML file into a current mrpast model.",
    )

    confidence_parser = subparsers.add_parser(
        CMD_CONFIDENCE, help="Calculate parameter confidence intervals."
    )
    confidence_parser.add_argument(
        "solved_result", help="A JSON file output by the solver."
    )
    confidence_parser.add_argument(
        "--bootstrap",
        "-b",
        action="store_true",
        help="Solve for all bootstrapped samples instead of using GIM (theoretical).",
    )
    confidence_parser.add_argument(
        "--replicates",
        "-r",
        type=int,
        default=None,
        help=(
            f"Number of solver replications to perform per bootstrap sample. "
            f"Defaults to {DEFAULT_SOLVE_REPS_PER_EPOCH} * num_epochs."
        ),
    )
    add_common(confidence_parser)

    polarize_parser = subparsers.add_parser(CMD_POLARIZE, help="Polarize a VCF file.")
    polarize_parser.add_argument("vcf_file", help="The input VCF file.")
    polarize_parser.add_argument(
        "ancestral",
        help="The ancestral FASTA file (input). Assumes the positions start counting at 1.",
    )
    polarize_parser.add_argument(
        "out_prefix", help="The output prefix for the haps/sample results."
    )

    show_parser = subparsers.add_parser(CMD_SHOW, help="Show solver results.")
    show_parser.add_argument("solved_result", help="A JSON file output by the solver.")
    show_parser.add_argument(
        "--sort-by", "-s", default="Index", help="Sort parameters by the column name."
    )

    select_parser = subparsers.add_parser(CMD_SELECT, help="AIC-based model selection.")
    select_parser.add_argument(
        "solved_results",
        nargs="+",
        help="Two or more JSON file output by the solver.",
    )
    select_parser.add_argument(
        "--bootstrap",
        "-b",
        action="store_true",
        help="Emit the distribution of AIC values for all bootstrapped samples. Requires "
        "that you have previously run 'mrpast confidence --bootstrap' to produce a .csv "
        "for each of the solved_results.",
    )

    args = parser.parse_args()

    if args.command == CMD_SIMULATE:
        total_trees = run_simulation(
            args.model,
            args.arg_prefix,
            args.seq_len,
            args.replicates,
            recomb_rate=args.recomb_rate,
            samples_per_pop=args.individuals,
            debug_demo=args.debug_demo,
            seed=args.seed,
            jobs=args.jobs,
        )
        print(f"Wrote {total_trees} total marginal trees")
    elif args.command == CMD_PROCESS:
        num_times, left_skew = args.num_times
        if left_skew is not None:
            assert left_skew.lower() == "l"
            left_skew = True
        else:
            left_skew = False
        process_ARGs(
            args.model,
            args.arg_prefix,
            args.suffix,
            args.jobs,
            do_solve=args.solve,
            out_dir=args.out_dir,
            add_ground_truth=args.add_ground_truth,
            num_times=num_times,
            replicates=args.replicates,
            bootstrap=args.bootstrap,
            bootstrap_iter=args.bootstrap_iter,
            max_generation=args.max_generation,
            tree_sample_rate=args.tree_sample_rate,
            min_time_unit=args.min_time_unit,
            leave_out=args.leave_out,
            group_by=args.group_by,
            time_slice_str=args.time_slices,
            verbose=args.verbose,
            rate_maps=args.rate_maps,
            rate_map_threshold=args.rate_map_threshold,
            left_skew_times=left_skew,
            seed=args.seed,
        )
    elif args.command == CMD_SOLVE:
        solver_outputs = solve(args.solver_inputs, args.jobs, args.timeout)
        print(f"Created outputs: {[fn for (fn, elapsed) in solver_outputs]}")
        best_output = process_solver_outputs(solver_outputs, verbose=True)
        print(f"The output with the highest likelihood is {best_output}")
    elif args.command == CMD_SIM2VCF:
        if args.prefix:
            arg_files = list(glob.glob(f"{args.arg_file}*.trees"))
            print(f"Found {len(arg_files)} ARGs")
        else:
            arg_files = [args.arg_file]
        arguments = []
        for i, arg_file in enumerate(arg_files):
            arguments.append(
                [
                    arg_file,
                    args.mut_rate,
                    set(parse_intlist(args.leave_out)),
                    args.seed + i,
                    args.zarr,
                ]
            )
        with Pool(args.jobs) as p:
            p.starmap(convert_simulation, arguments)
    elif args.command == CMD_ARGINFER:
        infer_arg(
            args.tool,
            args.vcf_prefix,
            args.arg_prefix,
            args.pop_map,
            Ne=args.ne_override,
            mu=args.mut_rate,
            recomb=args.recomb_rate,
            fasta=args.ancestral,
            samples=args.samples,
            jobs=args.jobs,
            seed=args.seed,
            thin=args.thin,
            dry_run=args.dry_run,
        )
    elif args.command == CMD_MODEL:
        print(f"Model {args.model} is valid")
        print_model_warnings(args.model)
        if args.to_demes or args.debug:
            assert (
                demes is not None
            ), 'Could not find demes module; try "pip install demes"'
            model = UserModel.from_file(args.model)
            demography, _ = build_demography(model)
            print(demography.debug())
            if args.to_demes:
                demes_graph = demography.to_demes()
                outfile = args.to_demes
                demes.dump(demes_graph, outfile)
                print(f"Wrote {outfile}")
    elif args.command == CMD_INIT:
        if args.from_demes is not None:
            print(
                "WARNING: Conversion from Demes model is still pretty rough, closely inspect "
                "the resulting model.",
                file=sys.stderr,
            )
            dump_model_yaml(convert_from_demes(args.from_demes), sys.stdout)
        elif args.from_old_mrpast is not None:
            model = load_old_mrpast(args.from_old_mrpast)
            dump_model_yaml(model, sys.stdout)
        else:
            print(
                'You must pass in additional options to the "init" command. See --help',
                file=sys.stderr,
            )
            exit(1)
    elif args.command == CMD_CONFIDENCE:
        if not args.bootstrap:
            result = subprocess.check_output(
                [eval_exe, "intervals", args.solved_result]
            )
            resulting_obj = json.loads(result.decode("utf-8"))
            print(json.dumps(resulting_obj, indent=2))
        else:
            base_name = remove_ext(os.path.basename(args.solved_result))
            out_dir = f"{base_name}.bootstrap.out"
            if os.path.exists(out_dir):
                print(
                    f"{out_dir} already exists; remove it if you want to rerun the analysis",
                    file=sys.stderr,
                )
                exit(2)
            os.mkdir(out_dir)
            with open(args.solved_result) as f:
                base_input = ModelSolverInput.from_json(f.read())
            if args.replicates is not None:
                reps = args.replicates
            else:
                reps = DEFAULT_SOLVE_REPS_PER_EPOCH * base_input.num_epochs
            inputs = []
            for i in range(reps):
                rep_input = base_input.randomize()
                rep_fn = os.path.join(out_dir, f"input.rep{i}.json")
                with open(rep_fn, "w") as fout:
                    fout.write(rep_input.to_json(indent=2))
                inputs.append(rep_fn)

            result_df = pd.DataFrame()
            samples = len(base_input.coal_count_matrices)
            outputs_with_times = solve(inputs, args.jobs, None, list(range(samples)))
            for i in range(samples):
                group = []
                for out_file, _ in outputs_with_times:
                    if f".mat{i}." in out_file:
                        group.append(out_file)
                best_output, best_negLL = get_best_output(group)
                print(f"Best result for matrix {i}: {best_output}")
                if best_output is not None:
                    best_df = load_json_pandas(best_output)
                    sample_column = [i for _ in range(len(best_df))]
                    best_df["sample"] = sample_column
                    negLL_column = [best_negLL for _ in range(len(best_df))]
                    best_df["negLL"] = negLL_column
                    result_df = pd.concat([result_df, best_df])
            csv_file = f"{base_name}.bootstrap.csv"
            result_df.to_csv(csv_file)
            print(f"Wrote DataFrame to {csv_file}")
    elif args.command == CMD_POLARIZE:
        out_file = args.out_prefix + ".vcf"
        if os.path.exists(out_file):
            raise RuntimeError(f"Output file {out_file} already exists")
        with open(args.vcf_file) as fin, open(out_file, "w") as fout:
            stats = polarize_vcf(fin, fout, args.ancestral)
        print(f"Finished polarizing. Wrote {out_file}")
        stats.print(sys.stdout, prefix="  ")
    elif args.command == CMD_SHOW:
        tab_show(args.solved_result, args.sort_by)
    elif args.command == CMD_SELECT:
        cmd = [eval_exe, "select"]
        if args.bootstrap:
            cmd.append("--bootstrap")
        if len(args.solved_results) < 2:
            print(
                "You must pass at least two solver output JSON files as input",
                file=sys.stderr,
            )
            exit(1)
        result = subprocess.check_output(cmd + args.solved_results)
        print(result)
    else:
        parser.print_help()
        exit(1)


if __name__ == "__main__":
    main()
