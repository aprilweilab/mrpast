![](https://github.com/aprilweilab/mrpast/actions/workflows/python-package.yml/badge.svg)

# mrpast

Infer demographic parameters from Ancestral Recombination Graphs (ARGs).

## Install

_Not yet published to PyPi_.

## Build/Install from repository

Requires Python 3.9 or newer, CMake, and a version of gcc or clang that supports C++17.

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

### Alternative installation options

1. Compile for the native CPU; this can drastically speed up the numerical solver, but makes the resulting package less portable.
```
MRPAST_ENABLE_NATIVE=1 pip install mrpast/ 
```

2. Build the solver in debug mode, so GDB can be attached.
```
MRPAST_DEBUG=1 pip install mrpast/
```

## IMPORTANT: NEW MODEL FORMAT

If you used MrPast prior to August 14, 2025, you may have models in the "old" format. The new format attempts
to be more user friendly. See [the examples](https://github.com/aprilweilab/mrpast/tree/main/examples) for the
new format.

To convert an "old style" model, `foo.yaml` to the new style, just do:
```
mrpast init --from-old-mrpast foo.yaml > foo.new.yaml
```

The documentation is still not updated w/r/t the new model format, but hopefully the examples are sufficient to
explain the changes.

## Usage

There are three primary subcommands to `mrpast`, and they are usually run in this order:
1. `mrpast simulate`
2. `mrpast process`
3. `mrpast solve`

These steps describe the "Simulated ARG" workflow, where no ARG inference is performed. See
the documentation for workflows making use of inferred ARGs.

### Simulation

In order to test out a demographic model, it is recommended that you start out
by simulating that model and verifying that `mrpast` can recover the model
parameters with the necessary accuracy. The simulation is done via
[msprime](https://tskit.dev/msprime/docs/stable/intro.html) and produces an
ancestral recombination graph (ARG) in the form of a tree-sequence file
(`.trees`).

Example:
```
# Simulate the model 10 times, using a DNA sequence length of 100Kbp and the default recombination rate
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

The resulting output files will be listed, and the best output (best likelihood)
will be listed as well. The JSON files for the output contains the parameter
values, their bounds, their initialized values, and (if present) their ground
truth values.

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

The demographic model is specified via [YAML](https://yaml.org/). See the [examples directory](https://github.com/aprilweilab/mrpast/tree/main/examples) for example models. See the documentation for details on model syntax and behavior.