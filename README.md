![](https://github.com/aprilweilab/mrpast/actions/workflows/python-package.yml/badge.svg)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/mrpast/README.html)

# mrpast

Infer demographic parameters from Ancestral Recombination Graphs (ARGs), by extracting pairwise coalescence counts for use in a maximum likelihood model. For detail on the method and results, see the [preprint](https://www.biorxiv.org/content/10.1101/2025.10.07.680347v1):

> DeHaas, Drew, Zhibai Jia, Leo Speidel, and Xinzhu Wei. "Inference of complex demographic history using composite likelihood based on whole-genome genealogies." bioRxiv (2025): 2025-10.

See [the documentation](https://mrpast.readthedocs.io/en/latest/) for commands, examples, and concepts.

See [PAPER_EXPERIMENTS.md](PAPER_EXPERIMENTS.md) for the configurations that were used in the paper.

![](docs/diagram.png)

## Install

Python 3.8 or newer is supported.

Install from PyPi:
```
pip install mrpast
```

On Linux, this will use prebuilt binaries. On MacOS, this will trigger a source code build, which requires [CMake](https://cmake.org/download/) and `gcc` or `clang` (C++17 support required).

You can also install the [conda package](https://bioconda.github.io/recipes/mrpast/README.html) via the [bioconda](https://bioconda.github.io/) channel: `conda install mrpast`.

## Build/Install from repository

Recommend using a virtual environment, the below creates and activates one:
```
python3 -m venv MyEnv
source MyEnv/bin/activate
```

Clone repo, then build and install:
```
git clone --recursive https://github.com/aprilweilab/mrpast.git
pip install mrpast/
```

## Usage

There are three primary subcommands to `mrpast`, and they are usually run in this order:
1. `mrpast simulate`
2. `mrpast process`
3. `mrpast solve`

These steps describe the "Simulated ARG" workflow, where no ARG inference is performed. See
[the documentation](https://mrpast.readthedocs.io/en/latest/) for workflows making use of inferred ARGs.

### Simulation

In order to test out a demographic model, it is recommended that you start out
by simulating that model and verifying that `mrpast` can recover the model
parameters with the necessary accuracy. The simulation is done via
[msprime](https://tskit.dev/msprime/docs/stable/intro.html) and produces
ancestral recombination graphs (ARGs) in the form of a [tskit](https://tskit.dev/tskit/docs/stable/introduction.html)
tree-sequence file (`.trees`).

Example:
```
# Simulate the model 10 times, using a DNA sequence length of 100Kbp and the default recombination rate (1e-8)
mrpast simulate --replicates 10 --seq-len 100000 --debug-demo examples/5deme1epoch.yaml 5de1
```

This creates 10 tree-sequence files (ARGs) that are named like `5de1*.trees`, using the given model.

### Processing

Given an ARG in tree-sequence format, either from simulation (see above) or from
ARG inference, we then extract coalescence information.

Example:
```
# Use 10 CPU threads to process the data and produce 10 replicates (expanded models) to be solved (later).
# `--bootstrap` creates 100 bootstrap samples by default, the average of which is used for input the maximum
# likelihood function
mrpast process --jobs 10 --replicates 10 --suffix trial1 --bootstrap coalcounts examples/5deme1epoch.yaml 5de1
```

See `mrpast process --help` for more options that control time discretization, distance between sampled trees, etc.

If we want, we could use `--solve` to run the solver as soon as processing completed. Otherwise, see the next section.

### Solving

If you didn't pass `--solve` to `mrpast process` then you can run the solver via:
```
mrpast solve --jobs 10 5deme1epoch.*.solve_in.*.json
```

The resulting output files will be listed, and the best output (maximum likelihood)
will be listed as well. The JSON files for the output contains the parameter
values, their bounds, their initialized values, and (if present) their ground
truth values. Assuming the best result was `5deme1epoch.trial1.solve_in.0.out.json`,
we can quickly view the results via

```
mrpast show -n 5deme1epoch.trial1.solve_in.0.out.json
```

### Other workflows

#### Simulated Data, Inferred ARG

The simulated data, inferred ARG workflow is:
1. `mrpast simulate`: Simulate your model with some ground-truth parameter values.
2. `mrpast sim2vcf -p`: Convert all .trees files with the given prefix to VCF files, and emit the corresponding .popmap.json files (which maps each sample to a population).
3. `mrpast arginfer`: Infer ARG from the VCF files, and then attach the population IDs to the ARG (.trees files) using the .popmap.json
4. `mrpast process`: Process and solve the inferred ARGs

#### Real Data, Inferred ARG

The real data workflow is:
1. Manually create a .popmap.json file for your VCF dataset. See the documentation for more details.
2. `mrpast arginfer`: Infer ARG from the VCF files, and then attach the population IDs to the ARG (.trees files) using the .popmap.json
3. `mrpast process`: Process and solve the inferred ARGs

## Modeling

The demographic model is specified via [YAML](https://yaml.org/). See the [examples directory](https://github.com/aprilweilab/mrpast/tree/main/examples) for example models. See [the documentation](https://mrpast.readthedocs.io/en/latest/) for details on model syntax and behavior.

## Alternative installation options

1. Compile for the native CPU; this can speed up the numerical solver, but makes the resulting package less portable.
```
MRPAST_ENABLE_NATIVE=1 pip install mrpast/
```

2. Build the solver in debug mode, so GDB can be attached.
```
MRPAST_DEBUG=1 pip install mrpast/
```

