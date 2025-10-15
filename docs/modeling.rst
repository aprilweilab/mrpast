.. _modeling:

Creating/Modifying Models
=========================

mrpast models are written in `YAML <https://yaml.org/>`_. They can be written fairly easily by hand, or you
can convert a `Demes <https://popsim-consortium.github.io/demes-spec-docs/main/introduction.html>`_ models
into a mrpast model (and then edit it as needed).

Converting from Demes
~~~~~~~~~~~~~~~~~~~~~

Example of converting from a Demes model:

::

  mrpast init --from-demes my_demes.yaml > my_mrpast.yaml

The resulting model may need to be modified to change lower/upper bounds on parameters. The conversion also
cannot tell when you want two parameters to be "the same" -- e.g., if the Demes model has symmetric migration
between demes ``A`` and ``B`` it will show up as two parameters in the mrpast model. You can remove one of them
and update the migration rate matrix to ensure there is only a single parameter for the migration.

Converting to Demes
~~~~~~~~~~~~~~~~~~~

You can convert a mrpast model to a Demes model via:

::

  mrpast model my_mrpast.yaml --to-demes my_demes.yaml

Components of a mrpast model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lets look at the `ooa_2t12.yaml <https://github.com/aprilweilab/mrpast/blob/main/examples/ooa_2t12.yaml>`_ mrpast model.
It starts with some comments and the ploidy of the individuals taking part in the model.

::

  # stdpopsim model OutOfAfrica_2T12, according to which:
  # "Model parameters are taken from Fig. S5 in Fu et al. (2013)"

  # Diploid individuals.
  ploidy: 2

Next we have a list of (optional) population names. Through-out, all lists and matrices that refer to populations
(or "demes") must have the same cardinality. The number of populations is never explicitly given, it is computed
from the lengths of these lists and matrices. We use :math:`D` to refer to the nubmer of populations.

::

  # The names of population0 and population1, respectively.
  pop_names:
  - AFR
  - EUR

For each population, we need an effective population size. In mrpast this is modeled as a coalescence rate, which
is just :math:`\frac{1}{ploidy \times N_e}` where `N_e` is the effective population size. At the heart of mrpast
models is a "parameter", which just a 1-based index associated with a bounded unknown value. Parameters can have
ground truth associated with them, for the purposes of simulation and quantifying error rates. The coalescence
rates are made up of two components: a list of coalesence rate definitions, and a list of parameters. Each parameter
has an ``index:`` field that specifies the unique, non-zero index it uses.

::

  # The coalescence rate parameterization, which is REQUIRED.
  coalescence:
    # Demes can be referenced by number (0, 1) or name (AFR, EUR). Rate values can be either
    # constants (just a float value) or a parameter reference of the form {"param": <index>}
    entries:
    - {"epoch": 0, "deme": "AFR", "rate": {"param": 3}}
    - {"epoch": 0, "deme": "EUR", "rate": {"param": 6}}
    - {"epoch": 1, "deme": "AFR", "rate": {"param": 2}}
    - {"epoch": 1, "deme": "EUR", "rate": {"param": 5}}
    - {"epoch": 2, "deme": "AFR", "rate": {"param": 2}}
    - {"epoch": 2, "deme": "EUR", "rate": {"param": 4}}
    - {"epoch": 3, "deme": "AFR", "rate": {"param": 2}}
    - {"epoch": 4, "deme": "AFR", "rate": {"param": 1}}

The corresponding parameters are next, and are nested within the ``coalescence`` object. Each object (``coalescence``,
``migration``, ``growth``) has its own parameter space within which indexes must be unique.

::

    # These indexes specified by these parameters matches up with the parameter references above.
    parameters:
    - {"ground_truth": 6.839945280437756e-05, "lb": 1e-07, "ub": 0.01, "index": 1}
    - {"ground_truth": 3.45447008428907e-05, "lb": 1e-07, "ub": 0.01, "index": 2}
    - {"ground_truth": 1.1570737191765288e-06, "lb": 1e-07, "ub": 0.01, "index": 3}
    - {"ground_truth": 0.00026867275658248256, "lb": 1e-07, "ub": 0.01, "index": 4}
    - {"ground_truth": 5.388388380070718e-05, "lb": 1e-07, "ub": 0.01, "index": 5}
    - {"ground_truth": 9.971355417745619e-07, "lb": 1e-07, "ub": 0.01, "index": 6}

So above, e.g., the ``{"param": 3}`` matches up with the ``{..., "index": 3}``, etc.

Next, an optional ``growth`` section, defining the rate at which each population grows. The growth rate is applied
backwards in time: the coalescence rate :math:`C(i)` given in the first section defines an effective population size :math:`N_e(i)`,
which is the population size for population :math:`i` at the *start* (nearest to current time) of the particular epoch.
The population size through-out the epoch (and at it's end, furthest away from current time) is defined as
:math:`\frac{1}{2 C(i) e^{\alpha(i) t}}` where :math:`\alpha(i)` is the growth rate and :math:`t` is the
number of generations backwards in time.

::

  # Growth rates are optional.
  growth:
    entries:
    - {"epoch": 0, "deme": "AFR", "rate": {"param": 1}}
    - {"epoch": 0, "deme": "EUR", "rate": {"param": 3}}
    - {"epoch": 1, "deme": "EUR", "rate": {"param": 2}}
    parameters:
    - {"ground_truth": 0.0166, "lb": 0.001, "ub": 0.05, "index": 1}
    - {"ground_truth": 0.0030700000000000002, "lb": 0.001, "ub": 0.05, "index": 2}
    - {"ground_truth": 0.0195, "lb": 0.001, "ub": 0.05, "index": 3}


Next we have the (optional) ``migration`` section which follows a similar parameterization, except that there are two demes
(source and destination) in the definition.

::

  # Migration rates are optional.
  migration:
    # The source/dest are with respect to BACKWARDS IN TIME. Everything in MrPast models is
    # backwards in time. So "source" is where lineages are migrating _from_ backwards in time
    # and where lineages are migrating _to_ forwards in time.
    entries:
    - {"epoch": 0, "source": "AFR", "dest": "EUR", "rate": {"param": 2}}
    - {"epoch": 0, "source": "EUR", "dest": "AFR", "rate": {"param": 2}}
    - {"epoch": 1, "source": "AFR", "dest": "EUR", "rate": {"param": 2}}
    - {"epoch": 1, "source": "EUR", "dest": "AFR", "rate": {"param": 2}}
    - {"epoch": 2, "source": "AFR", "dest": "EUR", "rate": {"param": 1}}
    - {"epoch": 2, "source": "EUR", "dest": "AFR", "rate": {"param": 1}}
    parameters:
    - {"ground_truth": 0.00015, "lb": 1e-05, "ub": 0.01, "index": 1}
    - {"ground_truth": 2.5e-05, "lb": 1e-05, "ub": 0.01, "index": 2}

Next have the times between each epoch. These can be fixed to a specific value by setting the lower
bound, upper bound, and ground truth all to the same value. 

::

  # The epoch splits are required, unless there is only a single epoch. If you have 4 splits
  # then there are 5 epochs.
  epochTimeSplit:
  - {"ground_truth": 204.6, "lb": 100.0, "ub": 300.0, "index": 1}
  - {"ground_truth": 920.0, "lb": 800.0, "ub": 1100.0, "index": 2}
  - {"ground_truth": 2040.0, "lb": 1900.0, "ub": 2200.0, "index": 3}
  - {"ground_truth": 5920.0, "lb": 5700.0, "ub": 6100.0, "index": 4}


Lastly, we define the conversion relationship between populations. Forwards in time you can think of this as populations
starting by splitting off from another population. Backwards in time this can be thought of as populations merging into
one. This is called "admixture", though you can define a split (just copy from one population to another)
by setting the proportion to be ``1.0``, as we do below. See the ``aa5.yaml`` example model
for demonstration on how to parameterize admixture when the proportions are not ``1.0``.

::

  # This admixture is just a population split, since the proportion is 1.0. After a split
  # occurs (backward in time) the derived population ceases to exist and should not be
  # referenced in any of the epochs where it is "dead"
  admixture:
    entries:
    - {"epoch": 3, "ancestral": "AFR", "derived": "EUR", "proportion": 1.0}
