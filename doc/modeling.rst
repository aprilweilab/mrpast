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

  # Informational only: the names of population0 and population1, respectively.
  pop_names:
  - AFR
  - EUR

For each population, we need an effective population size. In mrpast this is modeled as a coalescence rate, which
is just :math:`\frac{1}{ploidy \times N_e}` where `N_e` is the effective population size. At the heart of mrpast
models is a "parameter", which just a 1-based index associated with a bounded unknown value. Parameters can have
ground truth associated with them, for the purposes of simulation and quantifying error rates. The coalescence
rates are made up of two components: a list of vectors of parameter indexes, and a list of parameters. Each parameter
has an ``index:`` field that specifies the unique, non-zero index it uses. A zero value indicates no parameter. For
coalescence rates, no parameter is only valid when that population is inactive (e.g., it hasn't been created yet).

::

  # The coalescence rate parameterization.
  coalescence:
    # Each vector corresponds to an epoch, and we must have exactly one vector per epoch. Here
    # we have 5 epochs and 2 demes (populations).
    # Parameters are indexed 1-based.
    vectors:
    - [3, 6]  # Parameter 3 defines the coalescent rate of population0, and param 6 similarly is for population1, in epoch0
    - [2, 5]
    - [2, 4]
    - [2, 0]  # Here, in epoch3, there is no coalescent rate parameter for population1. This is because that population is inactive (not yet created)
    - [1, 0]

The corresponding parameters are next, and are nested within the ``coalescence`` object. Each object (``coalescence``,
``migration``, ``growth``) has its own parameter space within which indexes must be unique.

::

  parameters:
    # The ground-truth is used for simulation and for comparing results.
    - ground_truth: 6.839945280437756e-05
      # The lower and upper bound are to restrict the search space when solving maximum likelihood.
      lb: 1.0e-07
      ub: 0.01
      # The index here is how we can refer to this parameter in the coalescence vectors. If we put a "1" in the coalescence
      # vectors above, we are referring to this parameter.
      index: 1
    - ground_truth: 3.45447008428907e-05
      lb: 1.0e-07
      ub: 0.01
      index: 2
    - ground_truth: 1.1570737191765288e-06
      lb: 1.0e-07
      ub: 0.01
      index: 3
    - ground_truth: 0.00026867275658248256
      lb: 1.0e-07
      ub: 0.01
      index: 4
    - ground_truth: 5.388388380070718e-05
      lb: 1.0e-07
      ub: 0.01
      index: 5
    - ground_truth: 9.971355417745619e-07
      lb: 1.0e-07
      ub: 0.01
      index: 6

Next we have the ``migration`` section which follows a similar parameterization, except that each epoch is described
by a :math:`D \times D` matrix instead of a :math:`D`-length vector. Here, a 0 (no parameter) means there is no possible
migration between the two populations in the given epoch.

::

  migration:
    # Since migration rate is between a pair of populations, we have a single matrix per epoch instead of just a vector.
    matrices:
    # This first matrix is [ [0, 2], [2, 0] ]. The "2" values refer to the parameter with "index: 2" below. Parameters are
    # scoped only within their particular parameterization, so "2" under "migration:" is different from "2" under "coalescence:".
    - - [0, 2]
      - [2, 0] # The 0s on the diagonal are because migration between a population and itself does not make sense.
    - - [0, 2]
      - [2, 0]
    - - [0, 1]
      - [1, 0]
    - - [0, 0] # These last two epochs have no migration at all.
      - [0, 0]
    - - [0, 0]
      - [0, 0]
    parameters:
    - ground_truth: 0.00015
      lb: 1.0e-05
      ub: 0.01
      index: 1
    - ground_truth: 2.5e-05
      lb: 1.0e-05
      ub: 0.01
      index: 2

Next, an optional ``growth`` section, defining the rate at which each population grows. The growth rate is applies
backwards in time: the coalescence rate :math:`C(i)` given in the first section defines an effective population size :math:`N_e(i)`,
which is the population size for population :math:`i` at the *start* (nearest to current time) of the particular epoch.
The population size through-out the epoch (and at it's end, furthest away from current time) is defined as
:math:`\frac{1}{2 C(i) e^{\alpha(i) t}}` where :math:`\alpha(i)` is the growth rate and :math:`t` is the
number of generations backwards in time.

::

  # Growth rate parameterization for populations.
  growth:
    vectors:
    - [1, 3]  # Epoch0: Both population0 and population1 are growing
    - [0, 2]  # Epoch1: Only population1 is growing
    - [0, 0]  # Epoch2: No growth
    - [0, 0]  # Epoch3: No growth
    - [0, 0]  # Epoch4: No growth
    parameters:
    - ground_truth: 0.0166
      lb: 0.001
      ub: 0.05
      index: 1
    - ground_truth: 0.0030700000000000002
      lb: 0.001
      ub: 0.05
      index: 2
    - ground_truth: 0.0195
      lb: 0.001
      ub: 0.05
      index: 3


Our final set of parameters is for the times between each epoch. mrpast does not support fixed epoch times per-se, they
are always bounded parameters. If you need "fixed" epoch times, just set the bounds to be very tight.

::

  # The parameters for the transition time between each epoch. There are 5 epochs, so we need 4 of
  # these parameters (for each "in between")
  epochTimeSplit:
  # Split between epoch0 -> epoch1
  - ground_truth: 204.6
    lb: 100.0
    ub: 562.0
    index: 1
  # Split between epoch1 -> epoch2
  - ground_truth: 920.0
    lb: 562.0
    ub: 1480.0
    index: 2
  # Split between epoch2 -> epoch3
  - ground_truth: 2040.0
    lb: 1480.0
    ub: 3980.0
    index: 3
  # Split between epoch3 -> epoch4
  - ground_truth: 5920.0
    lb: 3980.0
    ub: 7910.0
    index: 4


Lastly, we define the conversion relationship between populations. Forwards in time you can think of this as populations
starting by splitting off from another population. Backwards in time this can be thought of as populations merging into
one. The values are 0-based, so :math:`0` is the first population and :math:`D-1` is the last population. If you want no
splits or merges in your model, each vector should just be :math:`0, 1, ..., D-1`.

::

  # How the populations map to each other in between epochs. When the position and the number in it match, there
  # is no splitting of the population. This is viewed backwards in time as a merge of the populations, so if we
  # have v[i] = j, then backwards in time population i merges into population j in between the relevant epochs.
  populationConversion:
  - [0, 1]  # epoch0 -> epoch1: v[0] == 0 and v[1] == 1, so there are no splits
  - [0, 1]  # epoch1 -> epoch2: v[0] == 0 and v[1] == 1, so there are no splits
  - [0, 0]  # epoch2 -> epoch3: v[0] == 0, but v[1] == 0, so here we have population1 splitting off from population0 at the start of epoch2 / end of epoch3
  - [0, 0]  # epoch3 -> epoch4: no change from the previous scenario, so there are no new splits

