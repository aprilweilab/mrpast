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
from typing import Optional, List, Union, Dict, Any, Tuple
from yaml import dump
import json
import os
import subprocess
import sys
import tempfile
import msprime

try:
    from yaml import CDumper as Dumper  # type: ignore
except ImportError:
    from yaml import Dumper  # type: ignore


def which(exe: str, required=False) -> Optional[str]:
    """
    Find the named executable, first via system PATH and then via the Python PATH.

        :param exe: The executable name.
    :param required: If True, throw an exception when not found instead of returning None.
    :return: None if the executable is not found.
    """
    try:
        result = (
            subprocess.check_output(["which", exe], stderr=subprocess.STDOUT)
            .decode("utf-8")
            .strip()
        )
    except subprocess.CalledProcessError:
        result = None
    if result is None:
        sys_paths = [os.path.join(p, "mrpast") for p in sys.path]
        for p in sys_paths + [os.path.realpath(os.path.dirname(__file__))]:
            p = os.path.join(p, exe)
            if os.path.isfile(p):
                result = p
                break
    if required and result is None:
        raise RuntimeError(f"Could not find executable {exe}")
    return result


def run(cmd: Union[str, List[str]], shell: bool = False, verbose: bool = False):
    if verbose:
        print(f"Running: {cmd}")
    if shell:
        subprocess.check_call(cmd, shell=True)
    else:
        subprocess.check_call([str(c) for c in cmd])


def remove_ext(filename: str, ext: Optional[str] = None) -> str:
    file_ext = filename.split(".")[-1]
    removed = ".".join(filename.split(".")[:-1])
    assert len(file_ext) < len(filename), "Filename has no extension"
    if ext is not None:
        assert ext == file_ext, f"Unexpected file extension on {filename}"
    return removed


def count_lines(filename: str) -> int:
    BUF_SIZE = 1024 * 100
    new_lines = 0
    with open(filename, "rb") as f:
        while True:
            data = f.read(BUF_SIZE)
            new_lines += data.count(b"\n")
            if len(data) < BUF_SIZE:
                break
    return new_lines


def dump_model_yaml(model: Dict[str, Any], out):
    """
    Generic YAML dumpers tend to produce JSON or not very readable YAML.

    This is specific to how we want the models layed out for readability, but
    falls back to a YAML dumper otherwise.
    """

    def dump_matrices(matrices: List[List[List[Union[int, float]]]], indent=2):
        prefix = " " * indent
        print(f"{prefix}matrices:", file=out)
        for matrix in matrices:
            for i, row in enumerate(matrix):
                if i == 0:
                    row_prefix = f"{prefix}- - "
                else:
                    row_prefix = f"{prefix}  - "
                print(f"{row_prefix}{json.dumps(row)}", file=out)

    def dump_vectors(vectors: List[List[Union[int, float]]], indent=2):
        prefix = " " * indent
        print(f"{prefix}vectors:", file=out)
        for vector in vectors:
            print(f"{prefix}- {json.dumps(vector)}", file=out)

    def dump_parameters(parameters: List[Dict[str, float]], indent=2, no_label=False):
        prefix = " " * indent
        if not no_label:
            print(f"{prefix}parameters:", file=out)
        for param in parameters:
            for i, key in enumerate(sorted(param.keys())):
                if i == 0:
                    line_prefix = f"{prefix}- "
                else:
                    line_prefix = f"{prefix}  "
                print(f"{line_prefix}{key}: {json.dumps(param[key])}", file=out)

    def dump_vectors_raw(vectors: List[List[Union[int, float]]], indent=2):
        prefix = " " * indent
        for vector in vectors:
            line_prefix = f"{prefix}- "
            print(f"{line_prefix}{json.dumps(vector)}", file=out)

    keys = set(model.keys())
    for key in (
        "ploidy",
        "pop_names",
    ):
        if key in keys:
            keys.remove(key)
            out.write(dump({key: model[key]}, Dumper=Dumper))
    for key in ("coalescence", "growth"):
        if key in keys:
            keys.remove(key)
            print(f"{key}:", file=out)
            dump_vectors(model[key]["vectors"])
            dump_parameters(model[key]["parameters"])
    for key in ("migration",):
        if key in keys:
            keys.remove(key)
            print(f"{key}:", file=out)
            dump_matrices(model[key]["matrices"])
            dump_parameters(model[key]["parameters"])
    for key in ("epochTimeSplit",):
        if key in keys:
            keys.remove(key)
            print(f"{key}:", file=out)
            dump_parameters(model[key], indent=0, no_label=True)
    for key in ("populationConversion",):
        if key in keys:
            keys.remove(key)
            print(f"{key}:", file=out)
            dump_vectors_raw(model[key], indent=0)
    remaining = {}
    for key in keys:
        remaining[key] = model[key]
    if remaining:
        out.write(dump(remaining, Dumper=Dumper))


def haps2vcf(input_prefix, output_prefix, ploidy=2):
    haps_file = f"{input_prefix}.haps"
    sample_file = f"{input_prefix}.sample"

    indivs = []
    with open(sample_file) as f:
        for i, line in enumerate(f):
            if i <= 1:
                continue
            id1, id2, missing = line.strip().split()
            if id1 == id2:
                ident = id1
            else:
                ident = f"{id1}_{id2}"
            indivs.append(ident)

    with open(haps_file) as f, open(f"{output_prefix}.vcf", "w") as fout:
        for i, line in enumerate(f):
            line = line.strip().split()
            chrom, varid, pos, ancestral, alt = line[:5]
            alleles = line[5:]
            assert len(indivs) == len(alleles) / ploidy
            if i == 0:
                print("##fileformat=VCFv4.2", file=fout)
                print("##source=mrpast", file=fout)
                print(f"##contig=<ID={chrom}>", file=fout)
                print(
                    f'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                    file=fout,
                )
                indivs_str = "\t".join(indivs)
                print(
                    f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{indivs_str}",
                    file=fout,
                )
            if ploidy == 1:
                paired_alleles = "\t".join(alleles)
            else:
                paired_alleles = []
                for i in range(0, len(alleles), ploidy):
                    paired_alleles.append("|".join(alleles[i : i + ploidy]))
                paired_alleles = "\t".join(paired_alleles)
            print(
                f"{chrom}\t{pos}\t{varid}\t{ancestral}\t{alt}\t.\t.\t.\tGT\t{paired_alleles}",
                file=fout,
            )


def relate_polarize(
    relate_root: str,
    vcf_file: str,
    ancestral_fasta: str,
    out_prefix: str,
):
    """
    Polarize the given VCF file using Relate's scripts, and then re-export it back to a VCF file.
    """
    assert vcf_file.endswith(".vcf"), f"Invalid VCF file (bad extension): {vcf_file}"
    vcf_file = os.path.abspath(vcf_file)
    ancestral_fasta = os.path.abspath(ancestral_fasta)

    RELATE_FILE_FORMATS = os.path.join(relate_root, "bin", "RelateFileFormats")
    assert os.path.isfile(RELATE_FILE_FORMATS)
    RELATE_PREP_INPUTS = os.path.join(
        relate_root, "scripts", "PrepareInputFiles", "PrepareInputFiles.sh"
    )

    prefix = vcf_file[:-4]
    base_prefix = os.path.basename(prefix)
    out_prefix = os.path.abspath(out_prefix)

    orig_dir = os.getcwd()
    try:
        with tempfile.TemporaryDirectory() as tmpdirname:
            print(f"Running in tempdir: {tmpdirname}")
            os.chdir(tmpdirname)

            run(
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

            run(
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
            run(["gunzip", f"prepared.haps.gz"])
            run(["gunzip", f"prepared.sample.gz"])
            assert not os.path.isfile("prepared.dist")

            haps2vcf("prepared", out_prefix)
    finally:
        os.chdir(orig_dir)


def make_zarr(vcf_file: str, delete_orig: bool = False) -> str:
    """
    This method uses command-line tools to convert from VCF (uncompressed) to the ZARR/VCF that is required
    for input to tsinfer. This process can be unfortunately quite slow. Requires tools:
    * bgzip
    * bcftools or tabix (just one)
    * vcf2zarr
    Result will be the same name/directory as the original file, but with a .vcz extension.
    """
    dir = os.path.dirname(vcf_file)
    base = remove_ext(os.path.basename(vcf_file), ext="vcf")
    vcz_file = os.path.join(dir, f"{base}.vcz")
    if os.path.exists(vcz_file):
        raise FileExistsError(
            f"Output {vcz_file} already exists; remove and try again."
        )
    bgzip = which("bgzip", required=True)
    bcftools = which("bcftools")
    tabix = which("tabix")
    assert (
        tabix is not None or bcftools is not None
    ), "bcftools or tabix is required for conversion to ZARR/VCF. Please install one of the two."
    vcf2zarr = which("vcf2zarr", required=True)
    assert vcf2zarr is not None

    with tempfile.TemporaryDirectory() as tmpdirname:
        intermediate = os.path.join(tmpdirname, f"{base}.vcf.gz")
        run(
            f'{bgzip} "{vcf_file}" --stdout > "{intermediate}"',
            shell=True,
            verbose=True,
        )
        if tabix:
            run([tabix, intermediate], verbose=True)
        else:
            assert bcftools
            run([bcftools, "index", intermediate], verbose=True)
        run([vcf2zarr, "convert", intermediate, vcz_file], verbose=True)
    if delete_orig:
        os.remove(vcf_file)
    return vcz_file


MIN_REC_RATE = 1e-12
MAX_CHROM_SIZE = 500_000_000


def load_ratemap(ratemap_file: str) -> msprime.RateMap:
    """
    Load a text file as a tskit.RateMap. Don't allow recombination rates that
    are too low (below 1e-12), as these tend to cause numerical issues later.
    """
    positions = []
    rates = []
    with open(ratemap_file) as f:
        for i, line in enumerate(f):
            line = line.strip()
            if line:
                start, end, rate_str = line.split()
                rate = float(rate_str)
                if rate < MIN_REC_RATE:
                    rate = MIN_REC_RATE
                if i == 0:
                    positions.append(0.0)
                    rates.append(float(rate))
                positions.append(float(start))
                rates.append(float(rate))
    positions.append(max(float(end), MAX_CHROM_SIZE))
    return msprime.RateMap(position=positions, rate=rates)


def get_best_output(filenames: List[str]) -> Tuple[Optional[str], float]:
    """
    Given a list of solver output filenames, pick the one with the lowest negative log-likelihood
    and return that filename plus the negative log-likelihood value.

    :return: Tuple (best_filename, best_negLL)
    """
    bestLL = 2**64
    best = None
    for fn in filenames:
        with open(fn) as f:
            data = json.load(f)
        negLL = data["negLL"]
        if negLL is None:
            negLL = 2**64
        if negLL < bestLL:
            bestLL = negLL
            best = fn
    return (best, bestLL)
