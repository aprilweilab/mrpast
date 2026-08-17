Command Reference
=================

Command options can change; use ``mrpast <command> -h`` to see the most up-to-date documentation. This page
just gives you an overview and may be slightly out of date.

``mrpast simulate``
-------------------

Given a model *with ground truth*, simulate it via `msprime <https://tskit.dev/msprime/docs/stable/intro.html>`_.

::

    usage: mrpast simulate [-h] [--jobs JOBS] [--seed SEED] [--verbose] [--replicates REPLICATES] [--seq-len SEQ_LEN] [--recomb-rate RECOMB_RATE] [--individuals INDIVIDUALS] [--debug-demo] model arg_prefix

    positional arguments:
    model                 The input YAML file specifying the model
    arg_prefix            The prefix for the output tree-sequence files

    options:
    -h, --help            show this help message and exit
    --jobs JOBS, -j JOBS  Number of jobs (threads) to use. Defaults to 1.
    --seed SEED           Set the random seed.
    --verbose, -v         Verbose output, including timing information.
    --replicates REPLICATES, -r REPLICATES
                            Number of simulation replications to perform. Defaults to 20.
    --seq-len SEQ_LEN, -s SEQ_LEN
                            Length of sequences in base-pairs. Default to 100000000.
    --recomb-rate RECOMB_RATE, -e RECOMB_RATE
                            Rate of recombination, or filename/prefix for recombination map. A prefix will match '<prefix>*.txt'. Defaults to 1e-08.
    --individuals INDIVIDUALS, -n INDIVIDUALS
                            Number of individuals per population. Defaults to 10.
    --debug-demo, -d      Output results from msprime demography debugger.


``mrpast process``
------------------

Process a set of ARGs into coalescence matrices. Optionally, run the solver and infer model parameters as well.

::

    usage: mrpast process [-h] [--jobs JOBS] [--seed SEED] [--verbose] [--replicates REPLICATES] [--num-times NUM_TIMES] [--solve] [--add-ground-truth] [--suffix SUFFIX] [--out-dir OUT_DIR]
                        [--min-time-unit MIN_TIME_UNIT] [--max-generation MAX_GENERATION] [--tree-sample-rate TREE_SAMPLE_RATE] [--leave-out LEAVE_OUT] [--bootstrap {none,coalcounts,jackknife}]
                        [--bootstrap-iter BOOTSTRAP_ITER] [--group-by GROUP_BY] [--time-slices TIME_SLICES] [--rate-maps RATE_MAPS] [--rate-map-threshold RATE_MAP_THRESHOLD] [--map-pops MAP_POPS]
                        [--timeout TIMEOUT]
                        model arg_prefix

    positional arguments:
    model                 The input YAML file specifying the model
    arg_prefix            The prefix of the input tree-seq file(s) specifying the ARG. Assumes .trees file extension.

    options:
    -h, --help            show this help message and exit
    --jobs JOBS, -j JOBS  Number of jobs (threads) to use. Defaults to 1.
    --seed SEED           Set the random seed.
    --verbose, -v         Verbose output, including timing information.
    --replicates REPLICATES, -r REPLICATES
                            Number of solver replications to perform. Defaults to 10 * num_epochs.
    --num-times NUM_TIMES, -t NUM_TIMES
                            Number of time slices to use. Defaults to 20. Use the suffix 'l' or 'L' to use left-skewed time slices.
    --solve, -s           Solve the model after generating the solver inputs.
    --add-ground-truth, -g
                            Generate an additional solver input(s) using the ground-truth parameter values.
    --suffix SUFFIX       Filenames will use the provided suffix instead of a random one.
    --out-dir OUT_DIR, -o OUT_DIR
                            Output directory.
    --min-time-unit MIN_TIME_UNIT, -u MIN_TIME_UNIT
                            The minimum time unit for distinguishing between coalescence events (default: 1.0 generation)
    --max-generation MAX_GENERATION, -m MAX_GENERATION
                            Ignore all coalescence events occuring after the given generation (default: 1000000.0)
    --tree-sample-rate TREE_SAMPLE_RATE, -b TREE_SAMPLE_RATE
                            Sample a tree from the ARG every tree-sample-rate base pairs (default: 125000 bp)
    --leave-out LEAVE_OUT
                            Comma-separated list of population IDs to leave out when counting coalescence
    --bootstrap {none,coalcounts,jackknife}
                            Bootstrap the sampled trees to create more than once coalescent matrix. coalcounts: standard bootstrap of over marginal trees. jackknife: leave-one-out jacktree over blocks of
                            marginal trees.
    --bootstrap-iter BOOTSTRAP_ITER, -i BOOTSTRAP_ITER
                            How many blocks to split the trees in for jackknifing, number of reps for standard bootstrap. Default: 100
    --group-by GROUP_BY   Regex to group ARGs or coal files by. By default group by chromosome for bootstrapping and by sample otherwise.
    --time-slices TIME_SLICES
                            The comma-separated list of time slice values instead of computing them from coalescence counts. Or, if prefixed with '+', the list of time slices to append to the auto-generated time
                            slices.
    --rate-maps RATE_MAPS
                            A filename prefix for tskit-style RateMap files, whose lexicographic sort order matches the input ARGs lexicographic sort order. Generates a glob '<prefix>*.txt'. Used for determining
                            tree sampling (see --rate-map-threshold)
    --rate-map-threshold RATE_MAP_THRESHOLD
                            Only sample trees from regions with a recombination rate <= to this. Requires --rate-maps
    --map-pops MAP_POPS   A list of <idx1>:<idx2>, comma-separated, which maps a particular population to another population, based on their 0-based indices. Useful for when the ARG populations are in a
                            different order (or not sampled) compared to the model.
    --timeout TIMEOUT     Solver timeout in seconds. Solver returns the current best result upon timeout.


``mrpast solve``
-----------------

Run after ``mrpast process``, on the ``.json`` files it created. Or just use the ``--solve`` option when running ``mrpast process`` and then you don't have to run this command.

::

    usage: mrpast solve [-h] [--timeout TIMEOUT] [--jobs JOBS] [--seed SEED] [--verbose] solver_inputs [solver_inputs ...]

    positional arguments:
    solver_inputs         The solver input JSON files. The output filenames will be derived from the input filenames.

    options:
    -h, --help            show this help message and exit
    --timeout TIMEOUT     Timeout in seconds. Solver returns the current best result upon timeout.
    --jobs JOBS, -j JOBS  Number of jobs (threads) to use. Defaults to 1.
    --seed SEED           Set the random seed.
    --verbose, -v         Verbose output, including timing information.


``mrpast sim2vcf``
------------------

Create ``.vcf`` files (or VCF/ZARR directories, for *tsinfer*) from ARGs that were simulated by ``mrpast``. Helpful because it also emits the population map ``.json`` file, which ``mrpast arginfer`` will need.
Usually used with the ``--prefix`` flag.

::

    usage: mrpast sim2vcf [-h] [--prefix] [--leave-out LEAVE_OUT] [--mut-rate MUT_RATE] [--zarr] [--jobs JOBS] [--seed SEED] [--verbose] arg_file

    positional arguments:
    arg_file              The ARG (.trees) file to process.

    options:
    -h, --help            show this help message and exit
    --prefix, -p          Treat arg_file as a prefix, and search for all <arg_prefix>*.trees files
    --leave-out LEAVE_OUT
                            Comma-separated list of population IDs to leave out when converting to VCF
    --mut-rate MUT_RATE   The mutation rate, for simulating mutations on existing trees.
    --zarr, -z            Output VCF/ZARR files, required for tsinfer usage.
    --jobs JOBS, -j JOBS  Number of jobs (threads) to use. Defaults to 1.
    --seed SEED           Set the random seed.
    --verbose, -v         Verbose output, including timing information.


``mrpast arginfer``
-------------------

Infer ARGs from genotype data. Does a lot of helpful things like determining `Ne` for Relate and SINGER, handling ancestral state for *tsinfer*, etc.

::

    usage: mrpast arginfer [-h] [--ne-override NE_OVERRIDE] [--mut-rate MUT_RATE] [--recomb-rate RECOMB_RATE] [--samples SAMPLES] [--thin THIN] [--dry-run] [--tool {singer,relate,tsinfer}]
                       [--ancestral ANCESTRAL] [--jobs JOBS] [--seed SEED] [--verbose]
                       vcf_prefix arg_prefix pop_map

    positional arguments:
    vcf_prefix            The prefix of VCF file(s) to process. Generates a glob "<vcf_prefix>*.vcf"
    arg_prefix            The prefix to use when writing the resulting ARGs to disk (.trees files)
    pop_map               The file containing the population map (*.popmap.json)

    options:
    -h, --help            show this help message and exit
    --ne-override NE_OVERRIDE, -N NE_OVERRIDE
                            Provide an override for the auto-calculated (diploid) effective population size.
    --mut-rate MUT_RATE, -m MUT_RATE
                            Expected mutation rate. Default 1.2e-08.
    --recomb-rate RECOMB_RATE, -r RECOMB_RATE
                            Expected recombination rate, or recombination map filename. Default 1e-08.
    --samples SAMPLES, -s SAMPLES
                            How many ARGS to sample. Default 10.
    --thin THIN, -t THIN  How many MC/MC iterations between samples. Default depends on the inference tool.
    --dry-run, -d         Just emit the arguments that would be used when running SINGER.
    --tool {singer,relate,tsinfer}
                            Which ARG inference tool to run: "tsinfer" (default), "relate", or "singer"
    --ancestral ANCESTRAL, -a ANCESTRAL
                            The ancestral FASTA file (input). Assumes the positions start counting at 1.
    --jobs JOBS, -j JOBS  Number of jobs (threads) to use. Defaults to 1.
    --seed SEED           Set the random seed.
    --verbose, -v         Verbose output, including timing information.


``mrpast model``
----------------

Validate, view, and/or export a ``mrpast`` model.

::

    usage: mrpast model [-h] [--to-demes TO_DEMES] [--debug] model

    positional arguments:
    model                 The model YAML file

    options:
    -h, --help            show this help message and exit
    --to-demes TO_DEMES, -d TO_DEMES
                            Write a Demes YAML file representing the model.
    --debug               Emit msprime demography debugger output for the given model


``mrpast show``
---------------

Show solver result. Just a helpful utility - often the `Python API <python_api.html>`_ is more useful for deeply analyzing results.

::

    usage: mrpast show [-h] [--sort-by SORT_BY] [--show-ne] solved_result

    positional arguments:
    solved_result         A JSON file output by the solver.

    options:
    -h, --help            show this help message and exit
    --sort-by SORT_BY, -s SORT_BY
                            Sort parameters by the column name.
    --show-ne, -n         Convert coalescence rates to Ne (effective population sizes).

``mrpast confidence``
---------------------

Performs two operations: generates parameter confidence intervals information, and generates bootstrap sample results that are needed for ``mrpast select``.

.. note::

    The default mode runs bootstrapping, but ``--gim`` does not. Bootstrapping can be very slow! Use ``--timeout``, ``--jobs``, and ``--replicates`` to make
    bootstrapping take less time.

::

    usage: mrpast confidence [-h] [--simple-expect] [--bootstrap] [--gim] [--replicates REPLICATES] [--timeout TIMEOUT] [--jobs JOBS] [--seed SEED] [--verbose] solved_result

    positional arguments:
    solved_result         A JSON file output by the solver.

    options:
    -h, --help            show this help message and exit
    --simple-expect       NOT RECOMMENDED. Calculate the Jacobian from the MLE, instead of averaging the gradients over many samples (bootstraps or ARG samples). This can help with some numerical issues, but
                            the resulting parameter intervals are likely to be over-confident.
    --bootstrap, -b       Solve MLE for all bootstrapped samples; can be very slow on large models!
    --gim, -g             Use the GIM method of computing confidence intervals. Faster, but possibly less accurate than bootstrapping.
    --replicates REPLICATES, -r REPLICATES
                            Number of solver replications to perform per bootstrap sample. Defaults to 10 * num_epochs.
    --timeout TIMEOUT     Solver timeout in seconds. Solver returns the current best result upon timeout.
    --jobs JOBS, -j JOBS  Number of jobs (threads) to use. Defaults to 1.
    --seed SEED           Set the random seed.
    --verbose, -v         Verbose output, including timing information.

``mrpast select``
-----------------

Generate JSON output that contains three model comparison statistics: ``AIC``, ``AIC_cl`` (preferred statistic: composite-likelihood adjusted AIC), and ``cl`` (composite likelihood).
Typical usage is ``mrpast select ... > selection_results.json``, and then loading the JSON file in an editor or with Python to examine it.

.. warning::

    Unless you use ``mrpast select --bootstrap``, you are only comparing a single solution from each model that is being compared. It is recommended to use
    ``--bootstrap`` to get a distribution of ``AIC_cl`` values to compare.

::

    usage: mrpast select [-h] [--bootstrap] solved_results [solved_results ...]

    positional arguments:
    solved_results   Two or more JSON file output by the solver.

    options:
    -h, --help       show this help message and exit
    --bootstrap, -b  Emit the distribution of AIC values for all bootstrapped samples. Requires that you have previously run 'mrpast confidence --bootstrap' to produce a .csv for each of the solved_results.

``mrpast pops``
---------------

Show the population information that is attached to an existing ARG(s):

::

    usage: mrpast pops show [-h] arg_prefix

    positional arguments:
    arg_prefix  The filename prefix for finding the input ARGs (.trees files)

    options:
    -h, --help  show this help message and exit

Attach a `population map <../overview/concepts.html#population-maps>`_ to an existing ARG(s):

::

    usage: mrpast pops attach [-h] [--ploidy PLOIDY] arg_prefix out_prefix pop_map

    positional arguments:
    arg_prefix       The filename prefix for finding the input ARGs (.trees files)
    out_prefix       The output prefix for writing the ARGs (now containing population info).
    pop_map          The file containing the population map (*.popmap.json)

    options:
    -h, --help       show this help message and exit
    --ploidy PLOIDY  The ploidy of individuals. Default: 2

``mrpast coalplot``
-------------------

Plot the coalescence distributions of closely related models/datasets. For example, plot simulated data against inferred
data. The model underlying the comparison does not have to be the same for all data, but the demes represented by those
models must be the same.

::

    usage: mrpast coalplot [-h] [--pmf] [--labels LABELS [LABELS ...]] [--keep-df KEEP_DF] output_file model result_jsons [result_jsons ...]

    positional arguments:
    output_file           Filename for the output image. Passed directly to matplotlib.pyplot.savefig().
    model                 A mrpast model YAML file that can be used to get deme names.
    result_jsons          One or more JSON files output by 'mrpast process' or 'mrpast solve'.

    options:
    -h, --help            show this help message and exit
    --pmf                 Instead of plotting the CDF (default), plot the proportion of coalescences that occur within each time slice as a line plot (PMF).
    --labels LABELS [LABELS ...]
                            Label the input JSON files
    --keep-df KEEP_DF     Save the underlying pandas.DataFrame in the given filename.

