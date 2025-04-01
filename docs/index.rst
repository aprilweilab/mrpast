mrpast Overview
===============

mrpast is a tool for inferring demographic parameters from phased genomic data.
It uses Ancestral Recombination Graphs (for an overview of ARGs, see `here <https://pmc.ncbi.nlm.nih.gov/articles/PMC10796009/>`_).
You can either infer ARGs from your data yourself, or use one of mrpast's integrations with
`Relate <https://myersgroup.github.io/relate/index.html>`_,
`tsinfer <https://tskit.dev/tsinfer/docs/stable/introduction.html>`_, or 
`SINGER <https://github.com/popgenmethods/SINGER>`_.


.. contents::
   :depth: 2

Usage Modes
~~~~~~~~~~~

mrpast takes a user-specified model as input. We have a bunch of
`example models <https://github.com/aprilweilab/mrpast/tree/main/examples>`_
that are a good starting place, and you can look at the `modeling <modeling.html>`_ section of the docs as well.
Models specify a ``ground_truth`` value for each parameter, which can be set to a random value initially (if you
don't know an expected value). The upside of providing a ``ground_truth`` value is that you can then simulate the
input model.

Simulating input models lets you explore the identifiability of the model with respect to mrpast. You can simulate
ARGs from a model, and if mrpast cannot accurately identify the model parameters then it is unlikely to be a useful
model on real data. Once you have verified that the simulated ARGs can be useful for mrpast, you can convert
the simulated ARGs to raw data (``mrpast sim2vcf``) and then re-infer ARGs (``mrpast arginfer``) to see if the
results are still accurate on inferred ARGs.

If a model behaves well with both simulated ARGs and inferred ARGs on simulated data, then you can move on to
`real data <real_data.html>`_.

Simplest Usage Example
~~~~~~~~~~~~~~~~~~~~~~

Below is the simplest usage for mrpast, which generates and processes simulated ARGs with a constant recombination rate.
This example only takes a couple minutes to run end-to-end on a typical laptop.

::

  # Install
  pip install mrpast

  # Download one of the example models
  wget https://github.com/aprilweilab/mrpast/raw/refs/heads/main/examples/5deme1epoch.yaml

  # Simulate some ARGs
  mkdir -p 5d1e.simdata && mrpast simulate -j 6 5deme1epoch.yaml 5d1e.simdata/5d1e_

  # Process the coalescence distribution from the ARGs and solve the maximum likelihood parameters
  mkdir -p 5d1e.output && mrpast process -j 6 --solve --out-dir 5d1e.output/ --bootstrap coalcounts 5deme1epoch.yaml 5d1e.simdata/5d1e_


The last command above will emit a message like
"The output with the highest likelihood is 5d1e.output/5deme1epoch.bfc276c6.solve_in.bootstrap.3.out.json"
at the end. We can then use the "mrpast show" command to see the results from that output:

::

  # Show the parameter values for the best result
  mrpast show 5d1e.output/5deme1epoch.bfc276c6.solve_in.bootstrap.3.out.json

Further Information
~~~~~~~~~~~~~~~~~~~

Use the links below or on the left menu-bar to learn more.

.. toctree::
  :maxdepth: 2

  Installation <installation>
  Concepts <concepts>
  Relate workflows <relate>
  tsinfer workflows <tsinfer>
  SINGER workflows <singer>
  Examining results <results>
  Using real data <real_data>
  Creating/modifying models <modeling>
  python_api

