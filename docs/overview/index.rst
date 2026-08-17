Overview
========

``mrpast`` is a tool for inferring demographic parameters from Ancestral Recombination Graphs (created from phased genomic data).

For an end-to-end example, see the `tsinfer workflow <../workflows/tsinfer.html>`_.

Recommendations for use
-----------------------

``mrpast`` can be used with many preexisting models, such as those available in `stdpopsim <https://popsim-consortium.github.io/stdpopsim-docs/stable/catalog.html>`_.
However, users are often developing their own model to explore the demographic history of a (potentially new) set of samples. Here we briefly describe
some of the ``mrpast`` features that aid in this process, and how to use them.

Process
~~~~~~~

The steps for creating a new model might look something like this:

* Hypothesize a model ``M1``, perhaps based on existing knowledge of the history of the populations/organisms involved.

  * This involves `creating a model <modeling.html>`_ with demes, migration rates, coalescence rates, and other parameters.

* Generate some reasonable "ground truth" for the model parameters.

  * Often this comes from running ``mrpast`` on real data, and letting it infer parameter values.

  * You may also have external knowledge about certain parameters; this can inform the ``ground_truth`` values directly, or the parameter lower/upper bounds in the model.

* Given a model with ``ground_truth`` you can now simulate the model with ``mrpast simulate``

* Given simulated ARGs, is ``mrpast process`` able to accurately re-infer the ground truth values? How much uncertainty is there in the estimates?

  * Poor results here may indicate that the model has identifiability problems, or ``mrpast`` has difficulty producing stable/robust results for the model.

* Now use ``mrpast sim2vcf`` to export the simulated data, ``mrpast arginfer`` to infer ARGs, and ``mrpast process`` to infer model parameters from simulated data ARGs.

  * Is ``mrpast`` still able to accurately recover the ground truth parameters?

  * Poor results here may indicate loss of information (or bias) in the ARG inference process. The accuracy here gives a sense of the best possible accuracy you can achieve on the real data, *if the model is in fact "true"*

* Alternate models ``M2``, ``M3``, etc., may be created (e.g., by modifying ``M1``). How do we choose between them?

  * Repeat the steps above with alternate models, to see if they are all equally identifiable. Models are approximating the true history, and the choice of what is modeled thoroughly (i.e., how detailed the model is) may require trade-offs with what is actually identifiable/inferrable.

  * Once you have a few reasonably inferrable models, you can use ``mrpast select`` to perform model comparison. This ``AIC_cl`` heuristic can give you some idea of which model better fits the data. Just like the inference process, model comparison should be done on both simulated and real data -- if the correct model cannot be determined from simulated data, then it is unlikely that model comparison on real data will be useful.

* Finally, on the chosen model, parameter confidence intervals can be computed with ``mrpast confidence``. Like everything else, this should be done on both simulated and real data - some parameters may have confidence intervals that are too tightly bounded, even in simulation.

* The number of time slices (``mrpast process --num-times <T>``) should be consistent for all steps above. The rule of thumb is that models with one or two epochs can use a few time slices (``20-40``), but models with more epochs should use more (``100-200``).

* Understanding `bootstrapping <bootstrapping.html>`_ is helpful for both confidence intervals and model comparison.

Usage Modes
-----------

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
`real data <../workflows/real_data.html>`_.

Simplest Usage Example
----------------------

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


Use the links below or on the menu-bar to learn more.

.. toctree::
  :maxdepth: 2

  Installation <installation>
  Concepts <concepts>
  Creating/modifying models <modeling>
  Bootstrapping <bootstrapping>


