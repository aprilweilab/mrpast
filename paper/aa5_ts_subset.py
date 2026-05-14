# Helper script to downsample ARGs by removing a population that we want to treat
# as unsampled.
import tskit
import sys

# Haplotypes by population. We know this mapping because this matches the order in the
# aa5.yaml model, which is what we simulated the data from.
# Also, when running "mrpast process" it will emit a mapping of populations from the
# model to the ARG, which can be used to validate.
SAMPLES_50 = (
    list(range(100))          # AFR
    + list(range(100, 200))   # EUR
    + list(range(200, 300))   # ASIA
                              # N_AMER: skip 300-400
    + list(range(400, 500))   # ADMIX
)


for chrom in range(0, 20):
    samples = SAMPLES_50
    suffix = "_50"

    infile = f"aa5{suffix}_tsi_aa5{suffix}_sim__{chrom}-0.trees.vcz.tsdate.trees"
    outfile = f"aa5{suffix}_tsi_unsamp_{chrom}-0.downsample.trees"

    ts = tskit.load(infile)
    ts = ts.simplify(samples=samples, filter_populations=False)
    ts.dump(outfile)
    print(f"Saved {outfile}")
