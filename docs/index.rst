mrpast Documentation
====================

``mrpast`` is a tool for inferring demographic parameters from phased genomic data.
It uses Ancestral Recombination Graphs (for an overview of ARGs, see `here <https://pmc.ncbi.nlm.nih.gov/articles/PMC10796009/>`_).
You can either infer ARGs from your data yourself, or use one of mrpast's integrations with
`tsinfer <https://tskit.dev/tsinfer/docs/stable/introduction.html>`_ (preferred),
`Relate <https://myersgroup.github.io/relate/index.html>`_, or
`SINGER <https://github.com/popgenmethods/SINGER>`_.

In our experience, ``mrpast`` obtains the best results with ``tsinfer+tsdate`` (the latter using the ``variational_gamma`` method).
For some models, the results can be substantially better than when using other ARG inference tools.

Use the links below or on the menu-bar to learn more.

.. toctree::
  :maxdepth: 2

  overview/index
  reference/index
  workflows/index
  faq/index


