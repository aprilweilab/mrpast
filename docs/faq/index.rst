.. _faq_howto:

FAQ/HOWTO
=========

Retaining coalescence counts
----------------------------

By default, the coalescence counts are temporarily kept in files in the system temp directory, and then removed
prior to ``mrpast process`` completing. If you want to keep these counts for some reason, then set
``MRP_COAL_DIR=<path>`` in your environment when running ``mrpast``.

System temp dir problems
------------------------

If you are having problems with your system temp directory (e.g., on a shared compute cluster), you
can tell ``mrpast`` to use another location by setting ``MRP_TMP_DIR=<path>`` in your environment.

GIM vs. Bootstrapping
---------------------

There are two scenarios where the Godambe information matrix (GIM) is used:

* Optionally, for parameter confidence intervals. ``mrpast confidence --gim`` estimates the confidence interval for each parameter from the GIM (based on first and second derivatives near the maximum likelihood estimate). This is *instead of* bootstrapping, and is much faster than bootstrapping, and for really large models (e.g., 20 demes) is the only practical option.
* For model comparison with `Akaike Information Criterion <https://en.wikipedia.org/wiki/Akaike_information_criterion>`_. This is independent from parameter confidence intervals. Model comparison can be performed on a *single instance* or *N bootstrap samples* - we recommend the latter. See `bootstrapping <../overview/bootstrapping.html>`_ for more details.