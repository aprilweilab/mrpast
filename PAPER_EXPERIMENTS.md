# Paper Experiments

Walk-through for running the main experiments from the _mrpast_ paper. Some of these commands, especially the ARG inference steps, can take hours to run, so we have documented them here instead of providing a more interactive tutorial (e.g., Jupyter notebooks). We marked each step with the order of time it takes (seconds, minutes, hours). See [the documentation](https://mrpast.readthedocs.io/) for more details and examples.

Throughout we use a lot of threads (20 or more), so the assumption is that these commands are being run on a large-ish server. To run locally (e.g., on a laptop), just reduce the `-j` parameter and wait longer.

## Recombination rate maps, ancestral sequences, etc.

### Regular rate maps

Almost all _mrpast_ commands that use a recombination rate map use the [msprime.RateMap](https://tskit.dev/msprime/docs/stable/rate_maps.html) format, which is just a three column text file with an entry per line and no header.  The columns are **(1)** base position (start), **(2)** base position (end), **(3)** rate.
The one exception is `mrpast arginfer` _only when using Relate_, which uses [the HapMap-style format](https://myersgroup.github.io/relate/input_data.html), though Relate is not used in any of the experiments described in this document.

The rate maps for the simulated data in these experiments were created by:
1. Downloading the HapMap recombination maps from [stdpopsim](https://popsim-consortium.github.io/stdpopsim-docs/stable/index.html).
2. Running [make_rate_map.py](scripts/make_rate_map.py) for chromosome 1, e.g.: `python scripts/make_rate_map.py genetic_map_Hg38_chr1.txt > ratemap1.chr1.txt`
3. Using the same ratemap for all 20 "simulated chromosomes", by duplicating `ratemap1.chr1.txt` to `ratemap1.chr2.txt` ... `ratemap1.chr19.txt`

### Background-selection modified rate maps

The rate maps for real data (1,000 Genomes Project) were created from the same HapMap maps, by updating them to take into account background selection.

1. Running [make_rate_map.py](scripts/make_rate_map.py) for all 22 autosomes, e.g.: `python scripts/make_rate_map.py genetic_map_Hg38_chr1.txt > ratemap.chr1.txt`
2. Modifying the rate maps so that any regions with elevated B-scores have extremely high recombination rates, so that the `--rate-map-threshold` option can avoid selecting trees from those regions.
    * Download the background selection BED file from [here](https://github.com/gmcvicker/bkgd/tree/master/data/hg38) and unzip it
    * Use the [bscore_modify.py](scripts/bscore_modify.py) script: `bscore_modify.py bkgd_hg38.bed ratemap.chr1.txt 1 > ratemap.withb.chr1.txt` (repeat for all autosomes)

In our initial testing, background selection did not appear to substantially change _mrpast_'s inferrence results, so it is possible that our results can be replicated using just the original rate maps.

### GRCh38 ancestral sequence

Simulated data is already "polarized", so does not need an ancestral sequence. For real data we used the [Ensembl Release 112 GRCh38 ancestral sequences](https://ftp.ensembl.org/pub/release-112/fasta/ancestral_alleles/) for [homo sapiens](https://ftp.ensembl.org/pub/release-112/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh38.tar.gz).
After untarring that package, _mrpast_ just needs a directory and file prefix for all FASTA files: one for each autosome used in the ARG
inference. See [here](https://mrpast.readthedocs.io/en/latest/concepts.html) for information about order and numbering of chromosome-based information.

### Population maps

We manually constructed the [mrpast population maps](https://mrpast.readthedocs.io/en/latest/concepts.html#population-maps) for the 1,000 Genomes Project analyses, after ensuring that the input (VCF) datasets ordered the samples correctly.

* Out-of-Africa 3-population model: [ooa3.i35.popmap.json](paper/ooa3.i35.popmap.json)
  * Using FIN instead of CEU: [ooa3fin.i35.popmap.json](paper/ooa3fin.i35.popmap.json)
* American Admixture 5-population (one unsampled) models:
  * [admix_clm.i35.popmap.json](paper/admix_clm.i35.popmap.json)
  * [admix_mxl.i35.popmap.json](paper/admix_mxl.i35.popmap.json)
  * [admix_pur.i35.popmap.json](paper/admix_pur.i35.popmap.json)

## Stepping stone model

We ran two variations of the stepping stone model. The first uses the [20deme1epoch.yaml](examples/20deme1epoch.yaml), which has parameters for migrations rates _and_ coalescent rates.
The second uses the [20d1e.fixed_ne.yaml](examples/20d1e.fixed_ne.yaml) model, where all the coalescent rates have been fixed to constants.

### With coalescent rate parameters

```
mkdir -p 20d1e.simdata/
mkdir -p 20d1e_10.output/

# Simulate the model (minutes to hours)
mrpast simulate -j 20 --individuals 10 --recomb-rate ratemap1.chr1.txt examples/20deme1epoch.yaml 20d1e.simdata/20d1e_10_sim_

# Process the model (extract coalescent count matrix) (minutes)
mrpast process -j 20 --rate-maps ratemap1.chr --tree-sample-rate 125000 --bootstrap coalcounts --suffix sim --out-dir 20d1e_10.output examples/20deme1epoch.yaml 20d1e.simdata/20d1e_10_sim_

# The 20d1e.output directory now contains our solver inputs and we can solve them all. We'll
# give a solver timeout of 48 hours (the default is 24 hours).
mrpast solve -j 10 --timeout 172800 20d1e_10.output/20deme1epoch.sim.solve_in.bootstrap.*.json
```

### Fixed coalescent rates

Same model as above, but the coalescent rate parameters are fixed instead of free. We already have the data simulated (see above), so we just need to process and solve.

```
mkdir -p 20d1e_10_fixed_ne.output/

# Process the model (minutes)
mrpast process -j 20 --rate-maps ratemap1.chr --tree-sample-rate 125000 --bootstrap coalcounts --suffix sim --out-dir 20d1e_10_fixed_ne.output 20d1e.fixed_ne.yaml 20d1e.simdata/20d1e_10_sim__

# The 20d1e_10_fixed_ne.output directory now contains our solver inputs and we can solve them all. We'll
# give a solver timeout of 48 hours (the default is 24 hours).
mrpast solve -j 10 --timeout 172800 20d1e_10_fixed_ne.output/20d1e.fixed_ne.sim.solve_in.bootstrap.*.json
```

## Single population growth model

The follow commands were run 50 times each, to get 50 separate "best results" on different simulated replicates. The `--seed` option ensures that we get different replicates, and is left as a variable `$REP` in the commands below, since it should be set to the replicate number (`0`, `1`, ...).

### Simulated ARGs, extra time slices

In the main figure, we show results with adding some extra time slices to cover the growth epoch.

```
mkdir -p 1t12.simdata/${REP}/
mkdir -p 1t12.output/${REP}/

# Simulate the model (seconds to minutes)
mrpast simulate -j 20 --individuals 50 --seed ${REP} examples/1t12.gr.yaml 1t12.simdata/${REP}/1t12_sim_

# Process AND SOLVE the model (minutes). Only on big models do we need to solve
# separately (because the solve is so slow and we need to set a timeout)
mrpast process -j 20 --solve --time-slices +100,150,200 --num-times 200 --tree-sample-rate 125000 --bootstrap coalcounts --suffix sim --out-dir 1t12.output/${REP}/ examples/1t12.gr.yaml 1t12.simdata/${REP}/1t12_sim_
```

### Simulated ARGs, standard time slices

In the supplement, we show results WITHOUT adding some extra time slices to cover the growth epoch.

```
mkdir -p 1t12_nots.output/${REP}/

mrpast process -j 20 --solve --num-times 200 --tree-sample-rate 125000 --bootstrap coalcounts --suffix sim --out-dir 1t12_nots.output/${REP}/ examples/1t12.gr.yaml 1t12.simdata/${REP}/1t12_sim_
```

### Inferred ARGs

We only inferred ARGs for the first (0th) replicate of simulated data.

```
mkdir -p 1t12.tsi.output/

# Convert simulated data to ZARR/VCF (minutes)
mrpast sim2vcf -j 20 -p --zarr 1t12.simdata/${REP}/1t12_sim_

# Infer the ARGs (hours)
mkdir -p 1t12.0.tsinfer/
mrpast arginfer -j 40 --tool tsinfer 1t12.simdata/0/1t12_sim_  1t12.0.tsinfer/1t12_tsi_ 1t12.simdata/0/1t12_sim__0-0.trees.popmap.json

# Process and solve the model (minutes)
mrpast process -j 20 --solve --time-slices +100,150,200 --num-times 200 --tree-sample-rate 125000 --bootstrap coalcounts --suffix tsi --out-dir 1t12.tsi.output/ examples/1t12.gr.yaml 1t12.0.tsinfer/1t12_tsi_
```

## Out-of-Africa 3-population model

### Simulated ARGs

```
mkdir -p ooa3.simdata/
mkdir -p ooa3.output/

mrpast simulate -j 20 --individuals 50 examples/ooa_3g09.yaml --recomb-rate ratemap1.chr ooa3.simdata/ooa3_50_sim_

# Process and solve the model (minutes).
mrpast process -j 20 --num-times 200 --solve --suffix simarg_ts200 --out-dir ooa3.output --tree-sample-rate 125000 --rate-maps ratemap1.chr --rate-map-threshold "1e-9" --bootstrap coalcounts examples/ooa3_3g09.yaml ooa3.simdata/ooa3_50_sim_

# Generate confidence intervals via bootstrapping local trees (minutes to hours). Here we assume that 
# ooa3.output/ooa_3g09.simarg_ts200.solve_in.bootstrap.31.out.json was the best result, but it may differ for your run
# We restrict to 20 replicates for speed (and it doesn't really affect results), but this can be omitted.
# Each replicate is given a solver timeout of 60 seconds to prevent really slow (almost always suboptimal) solutions
# from becoming the bottleneck. This command is running 100 bootstraps @ 20 replicates each which is 2000 solver runs.
mrpast confidence -j 20 --replicates 20 --timeout 60 ooa3.output/ooa_3g09.simarg_ts200.solve_in.bootstrap.31.out.json
```

### Inferred ARGs

```
mkdir -p ooa3.tsi/

# Convert to zarr (minutes to hours)
mrpast sim2vcf -j 20 --zarr --mut-rate 2.35e-8 -p ooa3.simdata/ooa3_50_sim_

# Infer ARGs (hours)
mrpast arginfer -j 40 --mut-rate 2.35e-8 --recomb-rate ratemap1.chr --tool tsinfer ooa3.simdata/ooa3_50_sim_ ooa3.tsi/ooa3_50_tsi_ ooa3.simdata/ooa3_50_sim__0-0.trees.popmap.json

# Process and solve (minutes)
mrpast process -j 20 --num-times 50L --solve --suffix tsi_ts50L --out-dir ooa3.output --tree-sample-rate 125000 --rate-maps ratemap1.chr --rate-map-threshold "1e-9" --bootstrap coalcounts examples/ooa_3g09.yaml ooa3.tsi/ooa3_50_tsi_

# Generate confidence intervals via bootstrapping local trees (minutes to hours).
mrpast confidence -j 20 --replicates 20 --timeout 60 ooa3.output/ooa_3g09.tsi_ts50L.solve_in.bootstrap.31.out.json
```

### 1,000 Genomes Project ARGs

**Data preparation:**

We used the high coverage 1,000 Genomes Project data from [here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/)

* Use the [kgp_sample_sets.sh](paper/kgp_sample_sets.sh) script to build the different mixtures of samples
  * See [ped_no_families.py](paper/ped_no_families.py) for an example of how we identified unrelated samples.
* Use the [kgp_vcfs.sh](paper/kgp_vcfs.sh) script to create the .vcf.gz files for all our datasets, for all 22 autosomes
* Run `vcf2zarr` from [bio2zarr](https://sgkit-dev.github.io/bio2zarr/vcf2zarr/overview.html) (`pip install bio2zarr[vcf]`)

**ARG inference and mrpast results:**

```
mkdir -p kgp_ooa3_i35.tsi/
mkdir -p kgp.ooa3.output/

# Infer the ARGs (takes 2-4 hours)
mrpast arginfer -j 30 --mut-rate 2.35e-8 --recomb-rate ratemap.withb.chr --tool tsinfer --ancestral homo_sapiens_ancestor_ kgp_ooa3_i35/1kGP_high_coverage_Illumina.chr kgp.ooa3.output/kgp_ooa3_i35_ ooa3.i35.popmap.json

# Process and solve the model (takes a few minutes at most)
mrpast process -j 20 --num-times 50L --suffix tsi_50L --out-dir kgp.ooa3.output --rate-maps ratemap.withb.chr --rate-map-threshold 1e-9 --solve --bootstrap coalcounts ooa_3g09.yaml kgp_ooa3_i35.tsi/kgp_ooa3_i35_1kGP_high_coverage_Illumina.chr
```

The comparison between the symmetric and asymmetric model was done on the above 1,000 genomes ARGs, but fixing the epoch times in both models.
Replace the models in the `mrpast process` commands with [ooa_3g09.fixedepochs.yaml](examples/ooa_3g09.fixedepochs.yaml) and [ooa_3g09.fe.asym.yaml](examples/ooa_3g09.fe.asym.yaml) to perform this comparison.


## American admixture 5-population model

### Simulated ARGs

```
mkdir -p aa5.simdata/
mkdir -p aa5.output/

# Simulate (minutes)
mrpast simulate -j 20 --individuals 50 examples/aa5.yaml --recomb-rate ratemap1.chr aa5.simdata/aa5_50_sim_

# Process and solve (minutes):
# - We add a few extra time slices to cover the very recent epoch
# - We leave population 3 (NAT) out of the coalescence count matrix, because we are treating that population
#   as unsampled.
#mrpast process -j 20 --leave-out 3 --num-times 150 --time-slices +5,10,15,20 --solve --suffix simarg_unsamp_ts150p --out-dir aa5.output --tree-sample-rate 125000 --rate-maps ratemap1.chr --rate-map-threshold "1e-9" --bootstrap coalcounts examples/aa5.yaml aa5.simdata/aa5_50_sim_

# Generate confidence intervals via bootstrapping local trees (minutes to hours).
mrpast confidence -j 20 --replicates 20 --timeout 60 aa5.output/aa5.simarg_unsamp_ts150p.solve_in.bootstrap.28.out.json
```

### Inferred ARGs

We want to treat NAT as unsampled, so we'll use [aa5_ts_subset.py](paper/aa5_ts_subset.py) to remove them
from the ARG.

```
mkdir -p aa5.tsi/

# Export to ZARR/VCF (minutes to hours)
mrpast sim2vcf -j 20 --zarr --mut-rate 2.35e-8 -p aa5.simdata/aa5_50_sim_

# Infer ARGs (hours)
mrpast arginfer -j 40 --mut-rate 2.35e-8 --recomb-rate ratemap1.chr --tool tsinfer aa5.simdata/aa5_50_sim_ aa5.tsi/aa5_50_tsi_ aa5.simdata/aa5_50_sim__0-0.trees.popmap.json

# Remove NAT population from ARGs
python aa5_ts_subset.py

# Process and solve (minutes)
mrpast process -j 20 --num-times 150 --time-slices +5,10,15,20 --solve --suffix tsi_unsamp_ts150p --out-dir aa5.output --tree-sample-rate 125000 --rate-maps ratemap1.chr --rate-map-threshold "1e-9" --bootstrap coalcounts examples/aa5.yaml aa5.tsi/aa5_50_tsi_unsamp_

# Generate confidence intervals via bootstrapping local trees (minutes to hours).
mrpast confidence -j 20 --replicates 20 --timeout 60 aa5.output/aa5.tsi_unsamp_ts150p.solve_in.bootstrap.56.out.json
```

### 1,000 Genomes Project ARGs

We ran the same model (`aa5.yaml`) with three different populations for the ADMIX population: CLM, MXL, and PUR.
Below we show the commands for running against the CLM data; just replace "_clm" with the other population identifiers as desired.

**Data preparation:**

Data preparation steps are covered in the section above (OOA3 model). Those scripts emit the data for this model as well.

**ARG inference and mrpast results:**

```
mkdir -p kgp_admix_clm_i35.tsi
mkdir -p aa5_clm.output

# Infer ARGs (hours)
mrpast arginfer -j 40 --mut-rate 2.35e-8 --recomb-rate ratemap.withb.chr --tool tsinfer --ancestral homo_sapiens_ancestor_ kgp_admix_clm_i35/1kGP_high_coverage_Illumina.chr kgp_admix_clm_i35.tsi/kgp_admix_clm_i35_ admix_clm.i35.popmap.json

# Process and solve (minutes)
# --map-pops is used to map population 3 (0-based) from the ARG to population 4 (0-based) in the model,
# which is the ADMIX population. This is because the ARG only has 4 populations, but the model have 5,
# and this leaves population 3 (NAT) in the model unsampled.
mrpast process -j 40 --num-times 150 --time-slices +5,10,15,20 --suffix full_clm_150p --out-dir aa5_clm.output --rate-maps ratemap.withb.chr --rate-map-threshold 1e-9 --map-pops 3:4 --solve --bootstrap coalcounts examples/aa5.yaml kgp_admix_clm_i35.tsi/kgp_admix_clm_i35_1kGP_high_coverage_Illumina.chr

# Generate confidence intervals via bootstrapping local trees (minutes to hours).
mrpast confidence -j 20 --replicates 20 --timeout 60 aa5_clm.output/aa5.full_clm_150p.solve_in.bootstrap.52.out.json
```