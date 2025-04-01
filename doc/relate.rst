
.. _relate:

Relate
======

This page walks you through an example using mrpast's direct integration with
`Relate <https://myersgroup.github.io/relate/index.html>`_.

.. note::
  You don't have to use mrpast's ARG inference integrations: you can use
  anything that produces a `tskit tree-sequence <https://tskit.dev/learn/>`_. However we assume that ARGs we
  process have the populations table setup properly to map individuals to
  populations. When using your own ARGs you may need to add this information
  yourself using the tskit APIs. See `attach_populations_ts() <https://github.com/aprilweilab/mrpast/blob/main/mrpast/arginfer.py>`_
  for an example of how to do this.


Installing dependencies
~~~~~~~~~~~~~~~~~~~~~~~

- Download and unpack one of the Relate installations from the `relate webpage <https://myersgroup.github.io/relate/index.html>`_.
- Create an environment variable ``RELATE_ROOT`` that points at your Relate installation's top-level directory. For Linux, this can be done via ``export RELATE_ROOT=/path/to/your/relate/`` (except with the correct path).
- Now look for ``EstimatePopulationSize.sh`` script in your relate directory under ``relate/scripts/EstimatePopulationSize/EstimationPopulationSize.sh`` and comment out or remove the last line (it should start with ``Rscript``). The alternative is to ensure that you have `R installed <https://www.r-project.org/>`_ with the relevant modules.
- Similarly, remove/comment out the last line of ``scripts/PrepareInputFiles/PrepareInputFiles.sh``

See a Relate `Pull Request for fixing these issues <https://github.com/MyersGroup/relate/pull/4>`_ for more details.

An Example
~~~~~~~~~~

.. note::
  Throughout this example we use a small number of simulation replicates ("chromosomes") to make this example
  take less time to run through. For accuracy, you would want to run with the default number of chromosomes (20)
  when testing out a new model.

Step 1: Simulate some data
--------------------------

We want to do a more realistic simulation, so we are going to use a recombination rate map.
Download the HapMap recombination map from stdpopsim (https://stdpopsim.s3.us-west-2.amazonaws.com/genetic_maps/HomSap/HapMapII_GRCh38.tar.gz), and untar the file.

::

  wget https://stdpopsim.s3.us-west-2.amazonaws.com/genetic_maps/HomSap/HapMapII_GRCh38.tar.gz
  tar -xf HapMapII_GRCh38.tar.gz


Now convert chromosome1's map to a file based on
`msprime.RateMap <https://tskit.dev/msprime/docs/stable/api.html#msprime.RateMap>`_, using the
`make_rate_map.py <https://github.com/aprilweilab/mrpast/blob/main/scripts/make_rate_map.py>`_
script in the mrpast repository.

::
  
  wget https://raw.githubusercontent.com/aprilweilab/mrpast/refs/heads/main/scripts/make_rate_map.py
  python make_rate_map.py genetic_map_Hg38_chr1.txt > ratemap_Hg38_chr1.txt

We'll simulate through mrpast directly, so we run:

::

  # Grab the example 5-deme synthetic model
  wget https://raw.githubusercontent.com/aprilweilab/mrpast/refs/heads/main/examples/5deme1epoch.yaml
  mkdir -p 5d1e.simdata/
  mrpast simulate -j 2 --replicates 2 --recomb-rate ratemap_Hg38_chr1.txt 5deme1epoch.yaml 5d1e.simdata/5d1e_msprime_

I've used ``-j 2`` here to use 2 threads, you can reduce or increase as computing
resources allow. The above will simulate 2 replicates ("chromosomes"), each of
which is 100Mbp in length, and with 10 diploid individuals per deme that the
model specifies. In this case, our model has 5 demes, so we'll have 50
individuals, which means 100 haplotypes. See ``mrpast simulate --help`` for more
options.

If you look in the ``5d1e.simdata/`` directory, you'll see a bunch of ``*.trees``
files that contain our simulated chromosomes.

.. note::
  For recombination maps, mrpast can take either a single file (with a ``.txt`` extension) or a file prefix.
  When using a single file, it will use the same map for all chromosomes being processed. When using a file
  prefix, the number of files matching that prefix must match the number of chromosomes. And the recombination
  maps are paired up according to lexicographic order: so it is expected that each file has the same prefix and
  then a suffix like ``chr1``, ``chr2``, etc., prior to the file extension.

Step 2: Export the data
-----------------------

Since we want to use inferred ARGs, we need to export our simulated data from
tree-sequences to input that can be used by ARG inference. Relate uses .haps/.sample
input files, but the mrpast front-end for Relate takes VCF. The ``sim2vcf`` subcommand
performs this conversion for us:

::

  mrpast sim2vcf --prefix --mut-rate 1.2e-8 --jobs 2 5d1e.simdata/5d1e_msprime_

Here ``--prefix`` means to convert any tree-sequence starting with the given
prefix (``5d1e.simdata/5d1e_msprime_``) We specify the mutation rate, because
our initial simulation was just coalescent trees in the ARG and not the actual
mutations. When we export to a genotype matrix (like VCF) we need the actual
mutations.

After running this command, the ``5d1e.simdata/`` directory should contain a
bunch of ``.vcf`` files, named similarly to the tree-sequence files.  There will
also be some .json files: these are the population maps, that specify how each
individual in the genotype data maps to each population. Here is an
example of one of these files:

::

  {
    "mapping": [
      [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
      [10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
      [20, 21, 22, 23, 24, 25, 26, 27, 28, 29],
      [30, 31, 32, 33, 34, 35, 36, 37, 38, 39],
      [40, 41, 42, 43, 44, 45, 46, 47, 48, 49]
    ],
    "names": [
      "pop_0",
      "pop_1",
      "pop_2",
      "pop_3",
      "pop_4"
    ]
  }

Each row in "mapping" corresponds to a population. So the 1st row (row 0) is the
first population (population 0), as defined in the original model file
(``mrpast/examples/5deme1epoch.yaml``). Each of the numbers in the row is an
individual ID (or "index") that is in the given population. So individuals 0-9
are in population 0, individuals 10-19 are in population 1, etc. The names of
the populations, in the same order, are "pop_0", "pop_1", etc. Since this is just
a synthetic example, the populations don't have meaningful names.


Step 3: Infer ARGs with Relate
------------------------------

Now we can run Relate on our data.

::

  mkdir -p 5d1e.relate/
  mrpast arginfer -j 2 --mut-rate 1.2e-8 --recomb-rate genetic_map_Hg38_chr1.txt --tool relate 5d1e.simdata/5d1e_msprime_ 5d1e.relate/5d1e_rel_ 5d1e.simdata/5d1e_msprime__0-0.trees.popmap.json


Even with only two chromosomes, the above can take some time (more than 30 minutes on my laptop).

.. note::
  For all recombination and mutation rate maps, mrpast uses files based on ``tskit.RateMap`` with no header row, EXCEPT that
  Relate takes a standard 3-column style map input. The Relate recombination rate file assumes columns 
  like ``Chromosome      Position(bp)    Rate(cM/Mb)``, with a header row. Hence, we use
  ``genetic_map_Hg38_chr1.txt`` above instead of ``ratemap_Hg38_chr1.txt``.

At the end of this step, we now have two sets of ARGs:

1. Simulated ARGs (``.trees`` files) in 5d1e.simdata/
2. Inferred ARGs (``.trees`` files) in 5d1e.relate/

Step 4: Process the ARGs
------------------------

We can now process the ARGs and solve for our model parameters. Lets first solve using the simulated ARGs:

::

  mkdir -p 5d1e.simarg.output/
  mrpast process -j 6 --solve --out-dir 5d1e.simarg.output/ --bootstrap coalcounts 5deme1epoch.yaml 5d1e.simdata/5d1e_msprime_


When processing completes, it will print something like "The output with the
highest likelihood is 5d1e.simarg.output/5deme1epoch.478104b8.solve_in.bootstrap.8.out.json". We can
then examine the result via:

::

  mrpast show 5d1e.simarg.output/5deme1epoch.478104b8.solve_in.bootstrap.8.out.json

Output will be something like:

::

    Index  Description                    Relative Error    Absolute Error        Truth        Final  Epochs
  -------  ---------------------------  ----------------  ----------------  -----------  -----------  --------
        0  Migration rate from 0->1           0.937115         0.000149021  0.000159021  1e-05        [0]
        1  Migration rate from 0->3           0.11433          0.000273001  0.00238783   0.00266083   [0]
        2  Migration rate from 0->4           0.0283911        5.07069e-05  0.00178602   0.00173531   [0]
        3  Migration rate from 1->0           1.32452          0.000527682  0.000398394  0.000926076  [0]
        4  Migration rate from 1->3           0.261859         0.000145552  0.00055584   0.000410289  [0]
        5  Migration rate from 1->4           0.329364         0.000155512  0.00047216   0.000316648  [0]
        6  Migration rate from 2->3           0.451319         3.71508e-05  8.23161e-05  0.000119467  [0]
        7  Migration rate from 2->4           0.0146486        1.85477e-06  0.000126617  0.000124762  [0]
        8  Migration rate from 3->0           0.0532306        1.47234e-05  0.000276597  0.00029132   [0]
        9  Migration rate from 3->1           0.47047          4.3004e-05   9.14063e-05  0.00013441   [0]
       10  Migration rate from 3->2           0.0437323        6.75181e-05  0.00154389   0.00161141   [0]
       11  Migration rate from 4->0           0.718232         5.56921e-05  7.75405e-05  2.18484e-05  [0]
       12  Migration rate from 4->1           0.0833728        7.46222e-06  8.95042e-05  9.69665e-05  [0]
       13  Migration rate from 4->2           0.516272         3.5834e-05   6.94091e-05  0.000105243  [0]
       14  Coalescence rate for deme 0        0.082526         4.19061e-05  0.000507792  0.000465886  [0]
       15  Coalescence rate for deme 1        0.233479         0.00019529   0.000836434  0.00103172   [0]
       16  Coalescence rate for deme 2        0.0132266        1.16129e-06  8.77998e-05  8.89611e-05  [0]
       17  Coalescence rate for deme 3        0.0497271        2.7004e-06   5.43044e-05  5.1604e-05   [0]
       18  Coalescence rate for deme 4        0.00781619       4.24454e-06  0.000543044  0.000547289  [0]

Now lets process the inferred ARGs:

::

  mkdir -p 5d1e.relarg.output/
  mrpast process -j 6 --solve --out-dir 5d1e.relarg.output/ --bootstrap coalcounts 5deme1epoch.yaml 5d1e.relate/5d1e_rel_


And again examine the result:

::

  mrpast show 5d1e.relarg.output/5deme1epoch.0dc1d8ef.solve_in.bootstrap.0.out.json

Output will be something like:

::

    Index  Description                    Relative Error    Absolute Error        Truth        Final  Epochs
  -------  ---------------------------  ----------------  ----------------  -----------  -----------  --------
        0  Migration rate from 0->1            61.8848         0.00984098   0.000159021  0.01         [0]
        1  Migration rate from 0->3             0.920164       0.00219719   0.00238783   0.000190633  [0]
        2  Migration rate from 0->4             0.994401       0.00177602   0.00178602   1e-05        [0]
        3  Migration rate from 1->0            24.1008         0.00960161   0.000398394  0.01         [0]
        4  Migration rate from 1->3             0.468887       0.000260626  0.00055584   0.000295214  [0]
        5  Migration rate from 1->4             0.267925       0.000126503  0.00047216   0.000345657  [0]
        6  Migration rate from 2->3           120.483          0.00991768   8.23161e-05  0.01         [0]
        7  Migration rate from 2->4             9.81543        0.0012428    0.000126617  0.00136942   [0]
        8  Migration rate from 3->0             0.963846       0.000266596  0.000276597  1e-05        [0]
        9  Migration rate from 3->1             1.26084        0.000115248  9.14063e-05  0.000206655  [0]
       10  Migration rate from 3->2             5.47713        0.00845611   0.00154389   0.01         [0]
       11  Migration rate from 4->0             0.871035       6.75405e-05  7.75405e-05  1e-05        [0]
       12  Migration rate from 4->1            18.3003         0.00163796   8.95042e-05  0.00172746   [0]
       13  Migration rate from 4->2             0.463705       3.21853e-05  6.94091e-05  0.000101594  [0]
       14  Coalescence rate for deme 0          0.514556       0.000261288  0.000507792  0.000246505  [0]
       15  Coalescence rate for deme 1          0.657592       0.000550033  0.000836434  0.000286401  [0]
       16  Coalescence rate for deme 2         14.1415         0.00124162   8.77998e-05  0.00132942   [0]
       17  Coalescence rate for deme 3         15.9959         0.000868646  5.43044e-05  0.000922951  [0]
       18  Coalescence rate for deme 4          0.784535       0.000426037  0.000543044  0.000117007  [0]

You can see that the overall relative error is significantly higher with the inferred ARGs than the simulated ARGs.
Both of these methods (simulated ARGs and inferred ARGs) have higher relative error than they would if we had used
more data (such as 20 chromosomes in our simulation), but the amount of data seems to especially affect inferred ARG
results. Another thing that can improve inferred ARG results is using the ``----rate-maps`` and
``--rate-map-threshold`` which lets you specify a recombination map and then only sample trees in regions with
recombination rate below the given threshold (``1e-9`` is usually a good threshold).

Take a look at the `analyzing real data <real_data.html>`_ section for more hints about improving results on larger datasets.