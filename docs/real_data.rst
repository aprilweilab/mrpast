.. _real_data:

Using real data
===============

The workflow for analyzing real (non-simulated data) is typically:

1. ``mrpast polarize`` the data; this requires Relate to be installed, since we currently make use of their polarization scripts.
2. ``mrpast arginfer`` on the polarized data to produce ARGs.
3. ``mrpast process --solve`` to process and solve the maximum likelihood problem.
4. ``mrpast confidence`` to generate confidence intervals on the parameters.

There are some additional considerations:

- It is often best to pass these options to ``mrpast process``: ``--rate-maps ratemap.chr --rate-map-threshold 1e-9``. This requires the ARG to only be sampled from regions with a recombination rate less than ``1e-9``. We have found that ARG inference tends to be more accurate in such regions.
- If you are concerned about particular regions of the genome (either the quality of sequencing, or things like selection influencing results), you can modify your rate maps to set the recombination rate really high in those regions. Then the above recombination rate threshold will prevent sampling from those regions.
