.. _examples:

Examples
========

This page shows some very brief examples for how to use the ``mrpast`` command line in certain scenarios.

ARG has too many populations
----------------------------

In this example, we'll simulate some ARGs with 3 populations and use them in a model with 2 populations.

::

  mkdir ooa3.simdata/
  mkdir ooa3_output/
  mkdir ooa2_output/

  # Simulate our ARGs first (2 ARGs, 10MB length each, 10 individuals per population, 3 population model)
  mrpast simulate -r 2 --seq-len 10000000 -n 10 examples/ooa_3g09.yaml ooa3.simdata/ooa3_

  # Use them in the same 3-population model
  mrpast process --out-dir ooa3_output examples/ooa_3g09.yaml ooa3.simdata/ooa3__

  # Now use them in the 2-population model by removing the 3rd population from the ARG
  mrpast process --leave-out 2 --out-dir ooa2_output examples/ooa2_2t12.yaml ooa3.simdata/ooa3__

  # Same thing, but lets also rearrange the order of the populations that we map. So we'll let the 2nd
  # and 3rd population from the ARG be used in our model, and drop the 1st population.
  mrpast process --leave-out 0 --map-pops 1:0,2:1 --out-dir ooa2_output examples/ooa2_2t12.yaml ooa3.simdata/ooa3__

ARG has too few populations
---------------------------

In the opposite scenario from above, we'll use the same simulated ARGs (3 populations) and use them in a model
with 4 populations, which means we need to treat one of them as unsampled (we won't have any coalescences from
the ARG to use).

::

  mkdir ooa4_output/

  # The model has 4 populations (0, 1, 2, 3) and we choose to map the 3 ARG populations to the first
  # 3 populations of that model. For the last mapping, we just need to choose a number that is NOT in
  # the ARG's populations and map it to "3" (the 4th population). Any number will do - here we use 100.
  mrpast process --map-pops 0:0,1:1,2:2,100:3 --out-dir ooa4_output/ ../examples/ooa_4j17.yaml ooa3.simdata/ooa3__