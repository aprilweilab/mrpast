
.. _singer:

SINGER
======

.. warning::
  SINGER support is still experimental, and our integration is not as well tested or supported as tsinfer or Relate.

This walkthrough uses `this fork of SINGER <https://github.com/dcdehaas/SINGER>`_ that has a few minor changes from
the `main fork <https://github.com/popgenmethods/SINGER>`_.

Installation:

::

  git clone https://github.com/dcdehaas/SINGER.git
  cd SINGER
  mkdir build && cd build
  cmake .. -DCMAKE_BUILD_TYPE=Release
  make -j 4
  cd ..
  ./make_release.sh
  export PATH=$PATH:$PWD/releases/singer-0.1.8-beta-linux-x86_64/

After running the above build and setting ``$PATH`` to include the SINGER executables you should have everything
needed to run SINGER. You only need to do the build once, but you'll need to set ``$PATH`` everytime you want
to use SINGER (or add it to your login script).

Steps 1-2: Same as Relate
-------------------------

See the `Relate walkthrough <relate.html>`_ and do steps 1 and 2, which should result in some VCF files and
population map JSON files.

Step 3: Infer ARGs with SINGER
------------------------------

Now we can run SINGER on our data.

::

  mkdir -p 5d1e.singer/
  mrpast arginfer -j 2 --mut-rate 1.2e-8 --recomb-rate ratemap_Hg38_chr1.txt --tool singer 5d1e.simdata/5d1e_msprime_ 5d1e.relate/5d1e_rel_ 5d1e.simdata/5d1e_msprime__0-0.trees.popmap.json

.. note::
  If you get error ``RuntimeError: Could not find executable parallel_singer`` then you don't have SINGER setup
  properly on your ``$PATH``. See installation instructions above.

Step 4-5: Same as Relate
------------------------

See the `Relate walkthrough <relate.html>`_ and do steps 4 and 5, except instead of ``5d1e.relate`` use ``5d1e.singer``.

.. note::
  SINGER performs MC/MC on both the topology of the ARGs and the branch lengths. You may need to increase the number
  of samples for ``mrpast arginfer`` (see ``--samples``) to get an adequate exploration of the ARG space.
