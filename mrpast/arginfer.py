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
from enum import Enum
from multiprocessing import Pool
from typing import Optional, List, Tuple, Union
import glob
import itertools
import os
import shutil
import sys
import tempfile
import tskit
import uuid

from mrpast.helpers import which, run, remove_ext, load_ratemap, MAX_CHROM_SIZE
from mrpast.model import PopMap


class ArgTool(Enum):
    SINGER = "singer"
    RELATE = "relate"
    TSINFER = "tsinfer"

    def __str__(self):
        return self.value


# We've found that larger thinning improves our results with SINGER
DEFAULT_SINGER_THIN = 400
DEFAULT_NUM_SAMPLES = 10
DEFAULT_SAMPLE_SUFFIX = "sample"
DEFAULT_SINGER_POLAR = 0.99  # Assume polarized


def _check_version_at_least(
    version_string: str,
    error_msg: str,
    min_major: int = 0,
    min_minor: int = 0,
    min_patch: int = 0,
):
    parts = version_string.split(".")
    if len(parts) == 2:
        major, minor = parts
        patch = "0"
    else:
        assert len(parts) == 3, f"Invalid version string: {version_string}"
        major, minor, patch = parts
    try:
        act_tuple = (int(major), int(minor), int(patch))
    except ValueError:
        raise RuntimeError(f"Invalid version string: {version_string}")
    min_tuple = (min_major, min_minor, min_patch)
    assert (
        min_tuple <= act_tuple
    ), f"{error_msg}, got version {act_tuple} but needed at least {min_tuple}"


def tsinfer_run(
    vcz_path: str,
    arg_prefix: str,
    mut_rate: float,
    recomb_rate_or_map: Union[float, str],
    pop_map_file: str,
    fasta_file: Optional[str] = None,
    dry_run: bool = False,
    ploidy: int = 2,
    jobs: int = 1,
):
    # We delay import of tsinfer, because it tries to precompile some Python code (I think, via numba)
    # and it can make import extremely slow.
    try:
        import tsinfer
        import tsdate
        import zarr
    except ImportError:
        assert (
            False
        ), "tsinfer/tsdate/zarr not found; try 'pip install tsinfer tsdate zarr'"

    # tsinfer changed the way they handle input; make sure we have the ZARR/VCF version.
    _check_version_at_least(
        tsinfer.__version__, "Invalid tsinfer version", min_minor=4, min_patch=1
    )

    with open(pop_map_file) as f:
        pop_map = PopMap.from_json(f.read())
    if isinstance(recomb_rate_or_map, str):
        ratemap = load_ratemap(recomb_rate_or_map)
    else:
        ratemap = tskit.RateMap(position=[0, MAX_CHROM_SIZE], rate=[recomb_rate_or_map])

    if dry_run:
        return

    vcf_zarr = zarr.open(vcz_path)
    if fasta_file is None:
        ancestral_state = vcf_zarr["variant_allele"][:, 0]
    else:
        try:
            import pyfaidx
        except ImportError as e:
            print(
                "Polarizing with tsinfer requires pyfaidx; try 'pip install pyfaidx'",
                file=sys.stderr,
            )
            raise e
        fasta_reader = pyfaidx.Fasta(fasta_file)
        assert (
            len(fasta_reader.values()) == 1
        ), "Your FASTA file has more than one contig; we only support a single contig."
        ancestral_str = "X" + str(list(fasta_reader.values())[0])
        if "ancestral_state" in vcf_zarr:
            del vcf_zarr["ancestral_state"]
        tsinfer.add_ancestral_state_array(vcf_zarr, ancestral_str)
        ancestral_state = "ancestral_state"

    vdata = tsinfer.VariantData(vcz_path, ancestral_state)
    print("Running tsinfer", file=sys.stderr)
    inferred_ts = tsinfer.infer(vdata, recombination_rate=ratemap, num_threads=jobs)

    print("Running tsdate", file=sys.stderr)
    simplified_ts = tsdate.preprocess_ts(inferred_ts)

    redated_ts = tsdate.date(simplified_ts, mutation_rate=mut_rate)

    with_pops = attach_populations_ts(redated_ts, pop_map, ploidy)

    out_file = f"{arg_prefix}.tsdate.trees"
    with_pops.dump(out_file)


def split_anc(filename: str) -> List[str]:
    """
    When RELATE samples multiple ARGs (just the branch lengths, no topologies) it produces
    a .anc file with multiple samples in it. Unfortunately this cannot be converted to .trees,
    so we have to split that into multiple .anc files first.

    :param filename: The filename of the .anc file.
    """
    assert filename.endswith(".anc")
    root = filename[:-4]
    result_filenames = []
    sample_count = 1
    with open(filename) as f:
        header_lines = []
        for i, line in enumerate(f):
            if not line.startswith("NUM_"):
                break
            if line.startswith("NUM_SAMPLES_PER_TREE"):
                sample_count = int(line.split()[-1])
                break
            header_lines.append(line.strip())
        try:
            output_files = []
            for i in range(sample_count):
                outfn = f"{root}_sample{i}.anc"
                result_filenames.append(outfn)
                output_files.append(open(outfn, "w"))
                output_files[i].write("\n".join(header_lines) + "\n")
            for line in f:
                pos, remainder = line.split(": ")
                for i in range(sample_count):
                    output_files[i].write(f"{pos}:")
                for item in remainder.strip().split(") "):
                    idx, data = item.split(":(")
                    for i in range(sample_count):
                        output_files[i].write(f" {idx}:(")
                    parts = data.split(" ")
                    assert len(parts) > sample_count
                    non_lengths = " ".join(parts[sample_count:])
                    for i in range(sample_count):
                        output_files[i].write(f"{parts[i]} {non_lengths})")
                for i in range(sample_count):
                    output_files[i].write("\n")
        finally:
            for o in output_files:
                o.close()
    return result_filenames


def get_samples(sample_file: str) -> List[str]:
    """
    Get the sample identifiers from the .sample file (Relate).
    """
    samples = []
    with open(sample_file) as f:
        for i, line in enumerate(f):
            # Skip the first two rows: first is column names, second is unused.
            # (see https://myersgroup.github.io/relate/input_data.html#FileFormat)
            if i > 1:
                line = line.strip()
                samples.append(line.split()[0].strip())
    return samples


def get_vcf_stats(vcf_file: str):
    """
    :return: Tuple (range, number of sites, number of individuals). The range is a pair
        of the first and last base-pair position observed.
    """
    first_pos = None
    individuals = 0
    sites = 0
    last_pos = None
    with open(vcf_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            data = line.split("\t")
            if individuals <= 0:
                individuals = len(data[9:])
            samples = itertools.chain.from_iterable(
                map(lambda x: x.split("|"), data[9:])
            )
            unique = len(set(samples))
            pos = int(data[1])
            if first_pos is None:
                first_pos = pos
            if pos != last_pos and unique > 1:
                sites += 1
            last_pos = pos
    return ((first_pos, last_pos), sites, individuals)


def watterson_ne(mu: float, sites: int, samples: int, bp_len: int) -> float:
    """
    Estimate the (panmictic) effective population size.
    """
    denom = 0.0
    for i in range(1, samples):
        denom += 1 / i
    theta = int(sites / denom) / bp_len
    return theta / (4 * mu)


def attach_populations_ts(
    tree_seq: tskit.TreeSequence, pop_map: PopMap, ploidy: int = 2
) -> tskit.TreeSequence:
    """
    Update the given tree-sequence file with the population information. There are two "types"
    of tree-sequences:

    * Ones made via Relate or SINGER. These have the individual order such that the
      sample nodes (0, 1) are for individual 0, (2, 3) are for individual 1, etc. There is
      no "individual table" in the tree file.
    * Ones made via tsinfer. These have an individual table, and assume the order of individuals
      matches the original order in the .vcz data.

    :param tree_seq: The tree-sequence.
    :param pop_map: A PopMap object describing individuals and their populations.
    :param ploidy: Default value of 2. The ploidy of the individuals.
    """
    tables = tree_seq.dump_tables()
    tables.populations.clear()
    # ARGs from tsinfer will already have a good schema, but ones from SINGER and Relate
    # will not, so we add the basic JSON schema in those cases.
    tables.populations.metadata_schema = tskit.MetadataSchema.permissive_json()

    indiv2pop = {}
    for pop_id in range(pop_map.num_pops):
        tables.populations.add_row(metadata={"name": pop_map.names[pop_id]})
        for indiv in pop_map.mapping[pop_id]:
            indiv2pop[indiv] = pop_id

    mapped = {}
    try:
        for id, row in enumerate(tables.nodes):
            if not (bool)(row.flags & tskit.NODE_IS_SAMPLE):
                continue
            if row.individual == tskit.NULL:
                individual = id // ploidy
            else:
                individual = row.individual
            mapped[id] = indiv2pop[individual]
    except:
        print(f"Mapped: {mapped}")
        print(f"indiv2pop: {indiv2pop}")

    for id, value in mapped.items():
        tables.nodes[id] = tables.nodes[id].replace(population=value)

    tables.sort()
    return tables.tree_sequence()


def attach_populations(tree_seq_file: str, pop_map: PopMap, ploidy: int = 2):
    """
    Update the given tree-sequence file with the population information. This function
    only works with tree-sequences created where the individual order is such that the
    sample nodes (0, 1) are for individual 0, (2, 3) are for individual 1, etc.

    :param tree_seq_file: The filename of the tree-sequence.
    :param pop_map: A PopMap object describing individuals and their populations.
    """
    orig_ts = tskit.load(tree_seq_file)
    ts = attach_populations_ts(orig_ts, pop_map, ploidy)
    with tempfile.TemporaryDirectory() as tmpdirname:
        filename = f"tmp.{str(uuid.uuid4())[:8]}.trees"
        ts.dump(filename)
        run(["mv", filename, tree_seq_file])


def singer_run(
    vcf_file: str,
    arg_prefix: str,
    Ne: float,
    mut_rate: Union[float, str],
    recomb: Union[float, str],
    jobs: int,
    samples: int,
    burnin_iter: int,
    seed: int,
    pop_map_file: str,
    thin: int = DEFAULT_SINGER_THIN,
    bp_range: Tuple[int, int] = (0, 500_000_000),
    sample_suffix: str = DEFAULT_SAMPLE_SUFFIX,
    dry_run: bool = False,
):
    """
    Run SINGER on a single VCF file.
    """
    parallel_singer = which("parallel_singer", required=True)
    singer_master = which("singer_master", required=True)
    convert_to_tskit = which("convert_to_tskit", required=True)

    if dry_run:

        def _run(args, **kwargs):
            print(args)

    else:
        _run = run

    burnin_samples = (burnin_iter + (thin - 1)) // thin

    if isinstance(recomb, float):
        assert isinstance(mut_rate, float)
        ratio = recomb / mut_rate
        recomb_map = None
    else:
        recomb_map = str(recomb)
        ratio = 1

        # Have to use a mutation map when using recombination map.
        if not isinstance(mut_rate, str):
            with open("mutmap.txt", "w") as f:
                f.write(f"0\t{bp_range[1]+1}\t{mut_rate}")
            mut_rate = os.path.abspath("mutmap.txt")

    vcf_file = os.path.abspath(vcf_file)
    with tempfile.TemporaryDirectory() as tmpdirname:
        orig_dir = os.getcwd()
        os.chdir(tmpdirname)
        print(f"Running in directory {tmpdirname}")
        try:
            vcf_prefix = remove_ext(os.path.basename(vcf_file), ext="vcf")
            with open(pop_map_file) as f:
                pop_map = PopMap.from_json(f.read())
            _run(f'cp "{vcf_file}" {vcf_prefix}.vcf', shell=True)
            vcf_index = 0

            total_samples = samples + burnin_samples
            common_args = [
                "-vcf",
                vcf_prefix,
                "-n",
                total_samples,
                "-thin",
                thin,
                "-polar",
                str(DEFAULT_SINGER_POLAR),
            ]
            if isinstance(mut_rate, str):
                common_args.extend(["-mut_map", mut_rate])
            else:
                common_args.extend(["-ratio", ratio, "-m", mut_rate])
            if Ne > 0.0:
                common_args.extend(["-Ne", Ne])
            if recomb_map is not None:
                common_args.extend(["-recomb_map", recomb_map])
            # XXX parallel_singer doesn't really work with jobs=1
            if jobs == 1:
                assert bp_range is not None
                _run(
                    [
                        singer_master,
                        "-start",
                        bp_range[0],
                        "-end",
                        bp_range[1],
                        "-output",
                        "TEMP_ARG_OUTPUT",
                    ]
                    + common_args
                )
                _run(
                    [
                        convert_to_tskit,
                        "-input",
                        "TEMP_ARG_OUTPUT",
                        "-output",
                        arg_prefix,
                        "-start",
                        1,
                        "-end",
                        total_samples,
                        "-step",
                        1,
                    ]
                )
            else:
                _run(
                    [
                        parallel_singer,
                        "-output",
                        arg_prefix,
                        "-num_cores",
                        jobs,
                    ]
                    + common_args
                )
            for i, sampled_tree in enumerate(glob.glob(f"{arg_prefix}*.trees")):
                index = int(sampled_tree.split("_")[-1].replace(".trees", ""))
                base = ".".join(os.path.basename(sampled_tree).split(".")[:-1])
                if index >= burnin_samples:
                    print(f"Attaching populations to inferred ARG tree-sequence")
                    attach_populations(sampled_tree, pop_map)
                    _run(
                        [
                            "mv",
                            sampled_tree,
                            f"{orig_dir}/{base}.{sample_suffix}{index:03d}_v{vcf_index}.trees",
                        ]
                    )
            assert dry_run or i > 0, "No sampled ARGs found"
        finally:
            os.chdir(orig_dir)


def relate_run(
    relate_root: str,
    vcf_file: str,
    arg_prefix: str,
    haploid_ne: float,
    mut_rate: float,
    recomb_rate_or_map: Union[float, str],
    num_samples: int,
    seed: int,
    pop_map_file: str,
    thin: Optional[int] = None,
    ancestral_fasta: Optional[str] = None,
    dry_run: bool = False,
):
    """
    Run Relate on a single VCF file.
    """
    assert vcf_file.endswith(".vcf"), f"Invalid VCF file (bad extension): {vcf_file}"
    assert haploid_ne > 0.0

    if dry_run:

        def _run(args, **kwargs):
            print(args)

    else:
        _run = run

    RELATE_FILE_FORMATS = os.path.join(relate_root, "bin", "RelateFileFormats")
    assert os.path.isfile(RELATE_FILE_FORMATS)
    RELATE_PREP_INPUTS = os.path.join(
        relate_root, "scripts", "PrepareInputFiles", "PrepareInputFiles.sh"
    )
    RELATE_MAIN = os.path.join(relate_root, "bin", "Relate")
    RELATE_POPSIZE = os.path.join(
        relate_root, "scripts", "EstimatePopulationSize", "EstimatePopulationSize.sh"
    )
    RELATE_BRANCH_LENS = os.path.join(
        relate_root, "scripts", "SampleBranchLengths", "SampleBranchLengths.sh"
    )

    prefix = vcf_file[:-4]
    base_prefix = os.path.basename(prefix)
    pop_prefix = f"{base_prefix}.popsize"

    out_dir = os.path.abspath(os.path.dirname(arg_prefix))
    assert os.path.isdir(out_dir), f"Output directory {out_dir} does not exist"
    arg_prefix = os.path.basename(arg_prefix)

    with open(pop_map_file) as f:
        pop_map = PopMap.from_json(f.read())

    orig_dir = os.getcwd()
    try:
        with tempfile.TemporaryDirectory() as tmpdirname:
            print(f"Running in tempdir: {tmpdirname}")
            os.chdir(tmpdirname)

            if isinstance(recomb_rate_or_map, str):
                rec_map = recomb_rate_or_map
            else:
                # Relate always wants a recombination map, so we provide a generic constant rate one
                # when we have a constant rec rate.
                r = float(recomb_rate_or_map)
                rec_map = "recmap.txt"
                with open(rec_map, "w") as f:
                    f.write("pos COMBINED_rate Genetic_Map\n")
                    f.write(f"0 {r*1e8} 0.0\n")
                    f.write(f"{MAX_CHROM_SIZE} {r*1e8} {r * 1e2 * MAX_CHROM_SIZE}\n")

            _run(
                [
                    RELATE_FILE_FORMATS,
                    "--mode",
                    "ConvertFromVcf",
                    "-i",
                    prefix,
                    "--haps",
                    f"{base_prefix}.haps",
                    "--sample",
                    f"{base_prefix}.sample",
                ]
            )

            if ancestral_fasta is not None:
                _run(
                    [
                        RELATE_PREP_INPUTS,
                        "--haps",
                        f"{base_prefix}.haps",
                        "--sample",
                        f"{base_prefix}.sample",
                        "-o",
                        "prepared",
                        "--ancestor",
                        ancestral_fasta,
                    ]
                )
                os.remove(f"{base_prefix}.haps")
                os.remove(f"{base_prefix}.sample")
                shutil.copyfile("prepared.haps.gz", f"{base_prefix}.haps.gz")
                shutil.copyfile("prepared.sample.gz", f"{base_prefix}.sample.gz")
                _run(["gunzip", f"{base_prefix}.haps.gz"])
                _run(["gunzip", f"{base_prefix}.sample.gz"])
                assert not os.path.isfile("prepared.dist")

            # We don't run the "annotate" script, because we don't need RELATE to know which samples
            # are associated with which populations. We attach the populations to the .trees file at
            # the end, and RELATE does everything agnostic to the structured demography.

            sample_ids = get_samples(f"{base_prefix}.sample")
            pop_file = f"{base_prefix}.poplabels"
            with open(pop_file, "w") as f:
                f.write("sample population group sex\n")
                for smpl in sample_ids:
                    f.write(f"{smpl} PAN PAN NA\n")

            intermediate_arg = f"{arg_prefix}.intermediate"
            _run(
                [
                    RELATE_MAIN,
                    "--mode",
                    "All",
                    "-m",
                    mut_rate,
                    "--haps",
                    f"{base_prefix}.haps",
                    "--sample",
                    f"{base_prefix}.sample",
                    "--map",
                    rec_map,
                    "--seed",
                    seed,
                    "-o",
                    intermediate_arg,
                    "--effectiveN",
                    haploid_ne,
                ]
            )

            _run(
                [
                    RELATE_POPSIZE,
                    "-i",
                    intermediate_arg,
                    "-m",
                    mut_rate,
                    "--poplabels",
                    pop_file,
                    "--seed",
                    seed,
                    "-o",
                    pop_prefix,
                ]
            )

            sample_cmd = [
                RELATE_BRANCH_LENS,
                "--mode",
                "All",
                "-m",
                mut_rate,
                "-i",
                intermediate_arg,
                "--map",
                rec_map,
                "--num_samples",
                num_samples,
                "--format",
                "a",
                "--coal",
                f"{pop_prefix}.coal",
                "-o",
                arg_prefix,
            ]
            if thin is not None:
                sample_cmd.extend(["--num_proposals", thin])
            _run(sample_cmd)

            full_anc = f"{arg_prefix}.anc"
            assert os.path.isfile(full_anc)
            anc_files = split_anc(full_anc)
            for anc in anc_files:
                assert anc.endswith(".anc")
                sample_prefix = anc[:-4]
                shutil.copyfile(f"{intermediate_arg}.mut", f"{sample_prefix}.mut")
                _run(
                    [
                        RELATE_FILE_FORMATS,
                        "--mode",
                        "ConvertToTreeSequence",
                        "-i",
                        sample_prefix,
                        "-o",
                        sample_prefix,
                    ]
                )
                tree_file = f"{sample_prefix}.trees"
                assert os.path.isfile(tree_file)
                attach_populations(tree_file, pop_map)
                shutil.copyfile(tree_file, os.path.join(out_dir, tree_file))
    finally:
        os.chdir(orig_dir)


def run_one_vcf_file(
    tool: ArgTool,
    vcf_file: str,
    arg_prefix: str,
    Ne: Union[str, float],
    mut_rate: float,
    recomb_rate_or_map: Union[float, str],
    num_samples: int,
    seed: int,
    pop_map_file: str,
    fasta_file: Optional[str] = None,
    burnin_iter: int = 1000,
    ploidy: int = 2,
    jobs: int = 1,
    thin: Optional[int] = None,
    dry_run: bool = False,
):
    def do_vcf_calculations():
        (start_bp, end_bp), sites, individuals = get_vcf_stats(vcf_file)
        if Ne == "auto":
            haploid_Ne = watterson_ne(
                mut_rate, sites, individuals * ploidy, (end_bp - start_bp)
            )
            print(f"Computed haploid Ne = {haploid_Ne}")
            diploid_ne = haploid_Ne / ploidy
        else:
            assert isinstance(Ne, float)
            diploid_ne = Ne
        return diploid_ne, start_bp, end_bp

    sys.stdout.flush()

    if tool == ArgTool.TSINFER:
        tsinfer_run(
            vcf_file,
            arg_prefix,
            mut_rate,
            recomb_rate_or_map,
            pop_map_file,
            fasta_file=fasta_file,
            dry_run=dry_run,
            jobs=jobs,
        )
    elif tool == ArgTool.RELATE:
        relate_root = os.environ.get("RELATE_ROOT")
        assert (
            relate_root is not None and len(relate_root) > 0
        ), f"RELATE_ROOT not set; required for running RELATE"
        diploid_ne, _, _ = do_vcf_calculations()
        relate_run(
            relate_root,
            vcf_file,
            arg_prefix,
            diploid_ne * 2,
            mut_rate,
            recomb_rate_or_map,
            num_samples,
            seed,
            pop_map_file,
            thin,
            fasta_file,
            dry_run,
        )
    else:
        assert tool == ArgTool.SINGER
        which("parallel", required=True)
        which("parallel_singer", required=True)
        which("convert_to_tskit", required=True)
        assert (
            fasta_file is None
        ), 'SINGER does not accept an ancestral FASTA file; try "mrpast polarize" instead.'
        if thin is None:
            thin = DEFAULT_SINGER_THIN
        diploid_ne, start_bp, end_bp = do_vcf_calculations()
        singer_run(
            vcf_file,
            arg_prefix,
            diploid_ne,  # SINGER uses diploid Ne
            mut_rate,
            recomb_rate_or_map,
            jobs,
            num_samples,
            burnin_iter,
            seed,
            pop_map_file,
            bp_range=(start_bp, end_bp),
            thin=thin,
            dry_run=dry_run,
        )


def infer_arg(
    tool: ArgTool,
    vcf_glob_prefix: str,
    arg_prefix: str,
    pop_map_file: str,
    Ne: Union[float, str] = "auto",
    mu: float = 1.2e-8,
    recomb: Union[float, str] = 1e-8,
    fasta: Optional[str] = None,
    samples: int = DEFAULT_NUM_SAMPLES,
    burnin_iter: int = 1000,
    ploidy: int = 2,
    jobs: int = 1,
    seed: int = 42,
    thin: Optional[int] = None,
    dry_run: bool = False,
):
    """
    Run SINGER on one or more VCF files matching a prefix pattern.
    """
    if tool == ArgTool.TSINFER:
        suffix = "vcz"
        Ne = 0  # Unused for tsinfer/tsdate
    else:
        suffix = "vcf"
    # The VCF files and the recombination map files must match in sorted order!
    vcf_files = list(
        sorted(map(os.path.abspath, glob.glob(f"{vcf_glob_prefix}*.{suffix}")))
    )

    recomb_list: List[Union[float, str]] = []
    if isinstance(recomb, str):
        if recomb.endswith(".txt") and os.path.isfile(recomb):
            recomb_list = [os.path.abspath(recomb)]
        else:
            recomb_list = list(
                map(os.path.abspath, sorted(glob.glob(f"{recomb}*.txt")))
            )
    else:
        recomb_list = [recomb]
    if len(recomb_list) == 1:
        recomb_list = recomb_list * len(vcf_files)
    assert len(vcf_files) == len(
        recomb_list
    ), "You must provide either a single recombination rate/map, or exactly one for each VCF/VCZ"
    fasta_list: List[Optional[str]] = []
    if isinstance(fasta, str):
        # Matched by lexicographic ordering. I.e., chr1.fa matches to chr1.vcf via order of the name.
        fasta_list = list(map(os.path.abspath, sorted(glob.glob(f"{fasta}*.fa"))))
        assert len(fasta_list) == len(
            vcf_files
        ), "You must provide an ancestral FASTA file pattern than matches one file for each VCF/VCZ"
    else:
        fasta_list = [fasta] * len(vcf_files)

    pop_map_file = os.path.abspath(pop_map_file)
    print(f"Running {str(tool)} on {len(vcf_files)} VCF/VCZ files")
    orig_dir = os.getcwd()
    out_dir = os.path.dirname(arg_prefix)
    if out_dir:
        os.chdir(out_dir)
        arg_prefix = os.path.basename(arg_prefix)
    outer_jobs = jobs if tool == ArgTool.RELATE else 1
    try:
        arguments = []
        for vcf_file, recomb_rate_or_map, fasta_file in zip(
            vcf_files, recomb_list, fasta_list
        ):
            suffix = os.path.basename(vcf_file[len(vcf_glob_prefix) :])
            arguments.append(
                (
                    tool,
                    vcf_file,
                    arg_prefix + suffix,
                    Ne,
                    mu,
                    recomb_rate_or_map,
                    samples,
                    seed,
                    pop_map_file,
                    fasta_file,
                    burnin_iter,
                    ploidy,
                    jobs,
                    thin,
                    dry_run,
                )
            )
        if outer_jobs == 1:
            for args in arguments:
                run_one_vcf_file(*args)
        else:
            with Pool(outer_jobs) as p:
                p.starmap(run_one_vcf_file, arguments)
    finally:
        os.chdir(orig_dir)
