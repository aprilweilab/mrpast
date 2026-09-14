"""
Use grapp and scikit-allel to extract site-frequency-spectrums from datasets
in GRG format, and emit the results in fastsimcoal2 format.
"""
from grapp.popgen import allele_counts, pop_allele_counts, population_pairs
import pygrgl
import allel
import argparse
import os
import sys


def create_fsc2_spectrums(
    grg: pygrgl.GRG,
    directory: str,
    model_name: str,
    folded: bool = False,
    verbose: bool = False,
):
    af_prefix = "DAF" if not folded else "MAF"
    assert not os.path.exists(directory), f"{directory} already exists; remove and try again"
    os.makedirs(directory)

    def _write_sfs(pop, ac):
        if folded:
            sfs = allel.sfs_folded(ac)
        else:
            sfs = allel.sfs(ac)
        popstr = f"pop{pop}"
        fname = f"{model_name}_{af_prefix}{popstr}.obs"
        with open(os.path.join(directory, fname), "w") as fout:
            print("1 observations", file=fout)
            print("\t".join(map(lambda i: f"d{pop}_{i}", range(len(sfs)))), file=fout)
            print(" ".join(map(str, sfs)), file=fout)

    # We have no populations - treat it as a single population, and exit early.
    pop_names = grg.get_populations()
    if not pop_names:
        if verbose:
            print(f"No population information in GRG; emitting a single SFS", file=sys.stderr)
        ac = allele_counts(grg, impute_missing=True, return_ref=folded)
        _write_sfs(0, ac)
        return
    if verbose:
        print(f"Emitting SFS for these populations (in order):", file=sys.stderr)
        for i, p in enumerate(pop_names):
            print(f"  {i}: {p}", file=sys.stderr)

    pop_ac = pop_allele_counts(grg, impute_missing=True, return_ref=folded)
    for pop in range(pop_ac.shape[0]):
        assert pop < len(pop_names)
        _write_sfs(pop, pop_ac[pop])

    for popA, popB in population_pairs(grg):
        assert popA > popB
        if folded:
            jsfs = allel.joint_sfs_folded(pop_ac[popA], pop_ac[popB])
        else:
            jsfs = allel.joint_sfs(pop_ac[popA], pop_ac[popB])
        fname = f"{model_name}_joint{af_prefix}pop{popA}_{popB}.obs"
        with open(os.path.join(directory, fname), "w") as fout:
            print("1 observations", file=fout)
            pidA = lambda i: f"d{popA}_{i}"
            pidB = lambda i: f"d{popB}_{i}"            
            print("\t".join(map(pidB, range(jsfs.shape[1]))), file=fout)
            for i in range(jsfs.shape[0]):
                print(pidA(i) + " " + " ".join(map(str, jsfs[i, :])), file=fout)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create fsc2-style SFS files."
    )
    parser.add_argument("grg", help="The input dataset in GRG format.")
    parser.add_argument("prefix", help="The output file prefix: usually the model prefix.")
    parser.add_argument("outdir", help="The directory to place outputs in; must not yet exist.")
    parser.add_argument(
        "--folded",
        action="store_true",
        help="Emit the folded SFS, instead of the full (derived allele known) SFS.",
    )
    args = parser.parse_args()
 
    grg = pygrgl.load_immutable_grg(args.grg)
    create_fsc2_spectrums(grg, args.outdir, args.prefix, folded=args.folded, verbose=True)