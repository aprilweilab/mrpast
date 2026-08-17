
Externally-created ARGs
=======================

This page walks you through an example using mrpast with `tskit ARGs <https://tskit.dev/learn/>`_ that
were inferred outside of ``mrpast`` (e.g., without using ``mrpast arginfer``). 

.. warning::

  On the models we have tested, using ``mrpast`` with ``SINGER`` or ``Relate`` was not as accurate as ``tsinfer+tsdate``. It may be worth comparing multiple ARG inference methods on your model(s).

``mrpast process`` is the command that takes your ARGs as input, and it requires the following:

* One or more ``.trees`` files with the same prefix. The expectation is that the suffix encodes two pieces of information: (1) the chromosome number, and the (2) MC/MC sample number (if applicable). The chromosome numbering should follow the same ordering as any `Rate Maps <../overview/concepts.html#rate-maps>`_ that will be used.
* Proper sample-to-population mapping in the ``.trees`` files.

ARG filenames and sampling
~~~~~~~~~~~~~~~~~~~~~~~~~~

``mrpast`` works best with as many local trees as possible, which means usually you should use an ARG for each of
the autosomes of the population(s) you are studying. The primary ordering of the ARG filenames is expected to
be based on the *chromosome number*. For example, for humans you might have ARG filenames that look like:

::

  myarg_1.trees
  myarg_2.trees
  ...
  myarg_22.trees

With such ARGs, you can do the following use-cases:

1. Multiple chromosomes, single MC/MC replicate, no bootstrapping: ``mrpast process --bootstrap none <model> myarg_``. This will produce a coalescence matrix from all sampled trees, across all chromosomes. Subsequent model solving will produce a single point estimate.
2. Multiple chromosomes, single MC/MC replicate, with bootstrapping (the default): ``mrpast process <model> myarg_``. This will produce ``100`` coalescence matrices from all sampled trees, across all chromosomes, by bootstrap sampling said trees. Subsequent model solving will produce a single point estimate, based on the average of the bootstrapped matrices. Future calls to ``mrpast confidence`` and ``mrpast select`` can use the bootstrapped matrices in more detail.

You will see output messages like ``Using a single ARG sample with sampling method arg_samples`` and
``Using a single ARG sample with sampling method bootstrap``, respectively, for the above use cases.

See `bootstrapping <../overview/bootstrapping.html>`_ for more details on how bootstrapping works.

If you have MC/MC samples from an ARG (e.g., branch length MC/MC samples from ``Relate`` or full ARG MC/MC samples from ``SINGER``)
you can use those *instead of bootstrapping*. This has two effects: (a) the point estimates are derived from the coalescence
matrix that is the average of the coalescence matrices for each MC/MC sample, and (b) downstream confidence intervals from ``mrpast confidence``
will use the MC/MC samples instead of bootstrap samples as input for estimating parameter variance.

.. warning::

  In our limited testing, MC/MC samples from ``Relate`` and ``SINGER`` produced *less variation in mrpast results*
  than tree-based bootstrapping did. That is, you are likely to get wider confidence intervals (and better coverage)
  with bootstrapping.


Extending our human example above, if we had 2 MC/MC samples per chromosome (i.e., 44 ARGs), the naming should look like:

::

  myarg_1.sample0.trees
  myarg_1.sample1.trees
  myarg_2.sample0.trees
  myarg_2.sample1.trees
  ...
  myarg_22.sample0.trees
  myarg_22.sample1.trees


You can make use of these MC/MC samples via ``mrpast process --bootstrap none <model> myarg_``. You should
see an output message like ``Using 2 ARG samples with sampling method arg_samples``. Theoretically, you can
combine bootstrapping and MC/MC samples, but this is not well-tested and is somewhat complicated - it is
recommended to use ``--bootstrap none`` if MC/MC samples are used. If you want bootstrapping with an MC/MC-based
ARG tool, it is recommended to just use the last MC/MC sample for each chromosome, instead of all MC/MC samples.

.. note::

  The default regular expression for finding MC/MC samples assumes that everything before ``sample`` is the
  chromosome ordering, and everything after is sample ordering (regex is ``sample([0-9]+)(_v[0-9]+)?.trees``).
  Advanced users can change this regex with the ``--group-by`` option.


Attaching and checking populations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can attach population information to your ARGs in a few ways.

1. You can write your own script(s). See `attach_populations_ts() <https://github.com/aprilweilab/mrpast/blob/main/mrpast/arginfer.py>`_ for an example of how to do this. The main requirements are that (a) every sample should be assigned to a population, and (b) the `population metadata <https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Population.metadata>`_ should have a ``name`` field.
2. You can create a mrpast-style ``.popmap.json`` file (see `here <../overview/concepts.html#population-maps>`_) and use the ``mrpast pops attach`` command.

``mrpast pops show <args prefix>`` can be used to show the population information in a set of ARGs.
