
.. _tsinfer:

tsinfer and tsdate
==================

This page walks you through an example using mrpast's direct integration with
`tsinfer <https://tskit.dev/tsinfer/docs/stable/introduction.html>`_ and
`tsdate <https://tskit.dev/tsdate/docs/latest/>`_.

.. note::
  You don't have to use mrpast's ARG inference integrations: you can use
  anything that produces a `tskit tree-sequence <https://tskit.dev/learn/>`_. However we assume that ARGs we
  process have the populations table setup properly to map individuals to
  populations. When using your own ARGs you may need to add this information
  yourself using the tskit APIs. See `attach_populations_ts() <https://github.com/aprilweilab/mrpast/blob/main/mrpast/arginfer.py>`_
  for an example of how to do this.


Installing dependencies
~~~~~~~~~~~~~~~~~~~~~~~

Install tsinfer and tsdate: ``pip install tsinfer tsdate bio2zarr[vcf] --upgrade``

The ``--upgrade`` is to ensure you have at least version ``0.4.1`` of tsinfer, which uses a new
file format compared to older versions.

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

We're going to use the `OutOfAfrica_3G09 model <https://popsim-consortium.github.io/stdpopsim-docs/stable/catalog.html#sec_catalog_homsap_models_outofafrica_3g09>`_
from stdpopsim. We have already converted this model from the Demes model
(see `the modeling section <modeling.html>`_), so we'll simulate through mrpast directly:

::

  # Grab the example OOA3 model
  wget https://raw.githubusercontent.com/aprilweilab/mrpast/refs/heads/main/examples/ooa_3g09.yaml
  mkdir -p ooa3.simdata/
  mrpast simulate -j 2 --replicates 2 --recomb-rate ratemap_Hg38_chr1.txt ooa_3g09.yaml ooa3.simdata/ooa3_msprime_

I've used ``-j 2`` here to use 2 threads, you can reduce or increase as computing
resources allow. The above will simulate 2 replicates ("chromosomes"), each of
which is 100Mbp in length, and with 10 diploid individuals per deme that the
model specifies. In this case, our model has 3 demes, so we'll have 30
individuals, which means 60 haplotypes. See ``mrpast simulate --help`` for more
options.

If you look in the ``ooa3.simdata/`` directory, you'll see a bunch of ``*.trees``
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
tree-sequences to input that can be used by ARG inference. tsinfer uses VCF/ZARR
input files and the ``sim2vcf`` subcommand performs this conversion for us:

::

  mrpast sim2vcf --zarr --prefix --mut-rate 2.35e-8 --jobs 2 ooa3.simdata/ooa3_msprime_

Here ``--prefix`` means to convert any tree-sequence starting with the given
prefix (``ooa3.simdata/ooa3_msprime_``) We specify the mutation rate, because
our initial simulation was just coalescent trees in the ARG and not the actual
mutations. When we export to a genotype matrix (like VCF/ZARR) we need the actual
mutations.

After running this command, the ``ooa3.simdata/`` directory should contain a
bunch of ``.vcz`` directories, named similarly to the tree-sequence files.  There will
also be some .json files: these are the population maps, that specify how each
individual in the genotype data maps to each population. Here is an
example of one of these files:

::

  {
    "mapping": [
      [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
      [10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
      [20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
    ],
    "names": [
      "YRI",
      "CEU",
      "CHB"
    ]
  }

Each row in "mapping" corresponds to a population. So the 1st row (row 0) is the
first population (population 0), as defined in the original model file
(`mrpast/examples/ooa_3g09.yaml`). Each of the numbers in the row is an
individual ID (or "index") that is in the given population. So individuals 0-9
are in population 0, individuals 10-19 are in population 1, etc. The names of
the populations, in the same order, are "YRI", "CEU", etc.


Step 3: Infer ARGs with tsinfer and tsdate
------------------------------------------

Now we can run tsinfer and tsdate on our data.

::

  mkdir -p ooa3.tsinfer/
  mrpast arginfer -j 4 --mut-rate 2.35e-8 --recomb-rate ratemap_Hg38_chr1.txt --tool tsinfer ooa3.simdata/ooa3_msprime_ ooa3.tsinfer/ooa3_ts_ ooa3.simdata/ooa3_msprime__0-0.trees.popmap.json


Even with only two chromosomes, the above can take some time (on the order of 20 minutes).
At the end of this step, we now have two sets of ARGs:

1. Simulated ARGs (``.trees`` files) in ooa3.simdata/
2. Inferred ARGs (``.trees`` files) in ooa3.tsinfer/

Step 4: Process the ARGs
------------------------

We can now process the ARGs and solve for our model parameters. Lets first solve using the simulated ARGs:

::

  mkdir -p ooa3.simarg.output/
  mrpast process -j 4 --num-times 50L --solve --out-dir ooa3.simarg.output/ --bootstrap coalcounts ooa_3g09.yaml ooa3.simdata/ooa3_msprime_


When processing completes, it will print something like "The output with the highest likelihood is ooa3.simarg.output/ooa_3g09.b0f8fc9b.solve_in.bootstrap.10.out.json".
We can then examine the result via:

::

  mrpast show ooa3.simarg.output/ooa_3g09.b0f8fc9b.solve_in.bootstrap.10.out.json

Which gives output something like:

::

    Index  Description                    Relative Error    Absolute Error           Truth           Final  Epochs
  -------  ---------------------------  ----------------  ----------------  --------------  --------------  ---------
        0  Epoch 0->1                         0.00783693       6.64571       848             841.354        []
        1  Epoch 1->2                         0.0733242      410.615        5600            6010.62         []
        2  Epoch 2->3                         0.0705646      620.968        8800            9420.97         []
        3  Migration rate from 0->1           0.0268311        6.70776e-06     0.00025         0.000256708  [1]
        4  Migration rate from 0->1           0.0970003        2.91001e-06     3e-05           3.291e-05    [0]
        5  Migration rate from 0->2           0.463635         8.80906e-06     1.9e-05         2.78091e-05  [0]
        6  Migration rate from 1->2           0.623214         5.98285e-05     9.6e-05         0.000155829  [0]
        7  Coalescence rate for deme 0        0.119665         8.1962e-06      6.84932e-05     6.02969e-05  [3]
        8  Coalescence rate for deme 0        0.00297262       1.20838e-07     4.06504e-05     4.07712e-05  [0, 1, 2]
        9  Coalescence rate for deme 1        0.00117098       2.78804e-07     0.000238095     0.000238374  [1]
       10  Coalescence rate for deme 1        0.00362068       6.09022e-08     1.68207e-05     1.67598e-05  [0]
       11  Coalescence rate for deme 2        0.00221614       2.04856e-08     9.2438e-06      9.22331e-06  [0]
       12  Growth rate for deme 1             0.0079936        3.19744e-05     0.004           0.00396803   [0]
       13  Growth rate for deme 2             0.0158041        8.69227e-05     0.0055          0.00558692   [0]

Now lets process the inferred ARGs:

::

  mkdir -p ooa3.tsarg.output/
  mrpast process -j 4 --num-times 50L --solve --out-dir ooa3.tsarg.output/ --bootstrap coalcounts ooa_3g09.yaml ooa3.tsinfer/ooa3_ts_


And again examine the result:

::

  mrpast show ooa3.tsarg.output/ooa_3g09.b499d52d.solve_in.bootstrap.24.out.json

Which gives output something like:

::

    Index  Description                    Relative Error    Absolute Error           Truth            Final  Epochs
  -------  ---------------------------  ----------------  ----------------  --------------  ---------------  ---------
        0  Epoch 0->1                          0.246024      208.628         848             1056.63         []
        1  Epoch 1->2                          0.285714     1600            5600             7200            []
        2  Epoch 2->3                          0.16793      1477.78         8800            10277.8          []
        3  Migration rate from 0->1            0.278297        6.95743e-05     0.00025          0.000180426  [1]
        4  Migration rate from 0->1            1.13211         3.39634e-05     3e-05            6.39634e-05  [0]
        5  Migration rate from 0->2            1.2886          2.44834e-05     1.9e-05          4.34834e-05  [0]
        6  Migration rate from 1->2            2.14897         0.000206301     9.6e-05          0.000302301  [0]
        7  Coalescence rate for deme 0         0.105603        7.23311e-06     6.84932e-05      6.126e-05    [3]
        8  Coalescence rate for deme 0         0.0417301       1.69634e-06     4.06504e-05      4.23468e-05  [0, 1, 2]
        9  Coalescence rate for deme 1         0.322068        7.66829e-05     0.000238095      0.000161412  [1]
       10  Coalescence rate for deme 1         0.0971438       1.63402e-06     1.68207e-05      1.84547e-05  [0]
       11  Coalescence rate for deme 2         0.625862        5.78534e-06     9.2438e-06       1.50291e-05  [0]
       12  Growth rate for deme 1              0.242275        0.0009691       0.004            0.0030309    [0]
       13  Growth rate for deme 2              0.32707         0.00179889      0.0055           0.00370111   [0]

You can see that the overall relative error is higher with the inferred ARGs than the simulated ARGs.
Both of these methods (simulated ARGs and inferred ARGs) have higher relative error than they would if we had used
more data (such as 20 chromosomes in our simulation). Another thing that can improve inferred ARG results is using the ``----rate-maps`` and
``--rate-map-threshold`` which lets you specify a recombination map and then only sample trees in regions with
recombination rate below the given threshold (``1e-9`` is usually a good threshold).

Take a look at the `analyzing real data <real_data.html>`_ section for more hints about improving results on larger datasets.