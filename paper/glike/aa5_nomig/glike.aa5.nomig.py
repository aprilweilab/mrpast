import glike
import itertools
import math
import numpy
import os
import random
import sys
import tskit

# Use the same population map for both glike and mrpast
from mrpast.model import PopMap

ADMIX_LOWER_BOUND = 0.05

## Parameters: ground truth and bounds values.
PARAMETERS = {
 "c_afr_e5": {"truth": 6.84931506849315e-05, "bounds": (5.0e-06, 0.0025)},
 "c_afr_e04": {"truth": 4.065040650406504e-05, "bounds": (5.0e-06, 0.0025)},
 "c_eur_e3": {"truth": 0.0002380952380952381, "bounds": (5.0e-06, 0.0025)},
 "c_eur_e2": {"truth": 0.000215857756459556, "bounds": (1.0e-06, 0.01)},
 "c_asia_e2": {"truth": 0.00030888182875505216, "bounds": (1.0e-06, 0.01)},
 "c_eur_e01": {"truth": 1.682066345910231e-05, "bounds": (5.0e-06, 0.0025)},
 "c_asia_e01": {"truth": 9.243796257772543e-06, "bounds": (5.0e-06, 0.0025)},
 "c_nat_e01": {"truth": 0.00047, "bounds": (5.0e-06, 0.0025)},
 "c_admix_e0": {"truth": 5.0e-05, "bounds": (5.0e-06, 0.0025)},
 "g_eur_e1": {"truth": 0.004, "bounds": (0.001, 0.05)},
 "g_asia_e1": {"truth": 0.0055, "bounds": (0.001, 0.05)},
 "t_e0": {"truth": 12.0, "bounds": (5.0, 50.0)},
 "t_e1": {"truth": 650.0, "bounds": (500.0, 800.0)},
 "t_e2": {"truth": 848.0, "bounds": (800.0, 1000.0)},
 "t_e3": {"truth": 5600.0, "bounds": (5400.0, 5800.0)},
 "t_e4": {"truth": 8800.0, "bounds": (8600.0, 9000.0)},
 "a_afr_admix": {"truth": 0.104687, "bounds": (ADMIX_LOWER_BOUND, 0.80)},
 "a_eur_admix": {"truth": 0.776327, "bounds": (ADMIX_LOWER_BOUND, 0.80)},
}

# When we run with the tsinfer ARG, we skip trees with these indices because they cause
# a bug in gLike that throws an exception.
SKIP_TSINFER_TREES = (12373, 17347, 34123, 42254, 48294)

# For all the gLike demos, they use ploidy=1, and presumably they assume their Ne values are also haploid
# effective pop sizes. In mrpast, when dealing with humans, we always use ploidy=2. So the Ne value would
# be 1/2*coal_rate. However, all our parameters are already specified at coal rate, so we don't need to
# transform it here.
def CR(coal_rate):
    return coal_rate

def aa5_demography(c_afr_e5, c_afr_e04, c_eur_e3, c_eur_e2, c_asia_e2, c_eur_e01, c_asia_e01, c_nat_e01, c_admix_e0,
                   g_eur_e1, g_asia_e1,
                   t_e0, t_e1, t_e2, t_e3, t_e4,
                   a_afr_admix, a_eur_admix):
    demo = glike.Demo()
    # Epoch 0
    demo.add_phase(glike.Phase(
       0, t_e0, # Simulated Truth: 0 - 12
       [CR(c_afr_e04), CR(c_eur_e01), CR(c_asia_e01), CR(c_nat_e01), CR(c_admix_e0)],
       grs = [0, 0, 0, 0, 0],
       populations = ["AFR", "EUR", "ASIA", "NAT", "ADMIX"]))

    # Epoch 1 (admixture event)
    a_nat_admix = max(0, 1 - (a_afr_admix + a_eur_admix))  # Use max() here to avoid numeric errors that may result in small negative
    P_ADMIXture = numpy.array([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [a_afr_admix, a_eur_admix, 0, a_nat_admix]])
    demo.add_phase(glike.Phase(
       t_e0, t_e1, # Simulated Truth: 12 - 650
       [CR(c_afr_e04), CR(c_eur_e01), CR(c_asia_e01) , CR(c_nat_e01)],
       grs = [0, g_eur_e1, g_asia_e1, 0],
       P = P_ADMIXture,
       populations = ["AFR", "EUR", "ASIA", "NAT"]),
     discretize = 250)  # Balancing time (this model is really slow) with accuracy, since I think this is the tricky epoch
                        # (because of length and growth), we give it more discretization

    # Epoch 2 (NAT split from ASIA)
    P_NAT_split = numpy.array([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [0, 0, 1]])
    demo.add_phase(glike.Phase(
        t_e1, t_e2,  # Simulated Truth: 650 - 848
        [CR(c_afr_e04), CR(c_eur_e2), CR(c_asia_e2)],
        [0, 0, 0],
        P = P_NAT_split,
        populations = ["AFR", "EUR", "ASIA"]))

    # Epoch 3 (ASIA split from EUR)
    P_ASIA_split = numpy.array([
        [1, 0],
        [0, 1],
        [0, 1]])
    demo.add_phase(glike.Phase(
        t_e2, t_e3,
        [CR(c_afr_e04), CR(c_eur_e3)],
        P = P_ASIA_split,
        populations = ["AFR", "EUR"]))

    # Epoch 4 (EUR split from AFR)
    P_ooa_split = numpy.array([
      [1],
      [1]])
    demo.add_phase(glike.Phase(
      t_e3, t_e4,
      [CR(c_afr_e04)],
      P = P_ooa_split,
      populations = ["AFR"]))

    # Epoch 5
    demo.add_phase(glike.Phase(
      t_e4, math.inf,
      [CR(c_afr_e5)],
      populations = ["AFR"]))

    return demo


# Set this to True to use gLike-example-style bounds on the admixture proportions
# so that they cannot exceed each other. Otherwise, we use mrpast-style bounds
# which do allow overlap. This is just due to a limitation of the mrpast solver
# which currently can't do symbolic bounds.
USE_SYMBOLIC_BOUNDS = True


if __name__ == "__main__":
    if sys.argv[1] == "DEMO_DUMP":
        truth = {k: v["truth"] for k, v in PARAMETERS.items()}
        demo = aa5_demography(**truth)
        #demography = glike.demo_to_demography(demo)
        #print("##### msprime Demography (not accurate - the raw demography is)")
        #print(demography.debug())
        #print()
        print("##### Raw Demography")
        demo.print()
        exit(0)
    tree_file = sys.argv[1]
    ts = tskit.load(tree_file)
    if "TREES10" in os.environ:
        SPACING = 3000000  # Roughly 10 trees
    elif "TREES25" in os.environ:
        SPACING = 1200000  # Roughly 25 trees
    elif "TREES50" in os.environ:
        SPACING = 600000  # Roughly 50 trees
    else:
        assert "TREES100" in os.environ
        SPACING = 300000  # Roughly 100 trees
    trees = [ts.at(pos).copy() for pos in range(SPACING, 30000000+SPACING, SPACING)]
    if "TSINFER" in os.environ:
        new_trees = list(filter(lambda t: t.index not in SKIP_TSINFER_TREES, trees))
        if len(new_trees) != len(trees):
            print(f"Skipped {len(trees) - len(new_trees)} trees due to gLike bug")
        trees = new_trees
    print(f"Using {len(trees)} trees", file=sys.stderr)

    # This assumes that the samples are the first 2N nodes in the tree sequence.
    pop_map_file = sys.argv[2]
    with open(pop_map_file) as f:
        pop_map = PopMap.from_json(f.read())
    print("POPULATION MAP", file=sys.stderr)
    glike_samples = {}
    for pop_id in range(pop_map.num_pops):
        name = pop_map.names[pop_id]
        sample_nodes = list(sorted(itertools.chain.from_iterable([(indiv*2, indiv*2+1) for indiv in pop_map.mapping[pop_id]])))
        print(f"{name}: {sample_nodes[0]}, {sample_nodes[1]}, ..., {sample_nodes[-1]}", file=sys.stderr)
        glike_samples[name] = sample_nodes

    def report_bad_trees(trees, demo, samples):
        fail = False
        for i, t in enumerate(trees):
            try:
                glike.glike_trees([t], demo, samples=samples)
            except Exception as e:
                print(f"Bad tree: i={i}, t={t} ({e})")
                fail = True
        assert not fail

    STATUS_EVERY = 100
    calls = 0
    def fun(c_afr_e5, c_afr_e04, c_eur_e3, c_eur_e2, c_asia_e2, c_eur_e01, c_asia_e01, c_nat_e01, c_admix_e0,
            g_eur_e1, g_asia_e1,
            t_e0, t_e1, t_e2, t_e3, t_e4,
            a_afr_admix, a_eur_admix):
        demo = aa5_demography(c_afr_e5, c_afr_e04, c_eur_e3, c_eur_e2, c_asia_e2, c_eur_e01, c_asia_e01, c_nat_e01, c_admix_e0,
                              g_eur_e1, g_asia_e1,
                              t_e0, t_e1, t_e2, t_e3, t_e4,
                              a_afr_admix, a_eur_admix)
        # Uncomment this out to report information about local trees that cause gLike to throw
        # errors, then add them to SKIP_TSINFER_TREES.
        #report_bad_trees(trees, demo, glike_samples)
        result = glike.glike_trees(trees, demo, samples=glike_samples)
        global calls
        calls += 1
        if calls % STATUS_EVERY == 1:
            print(f"After {calls} evaluations, best result is {result}")
        return result

    param_order = ["c_afr_e5", "c_afr_e04", "c_eur_e3", "c_eur_e2", "c_asia_e2", "c_eur_e01", "c_asia_e01",
                   "c_nat_e01", "c_admix_e0", "g_eur_e1", "g_asia_e1", "t_e0", "t_e1", "t_e2",
                   "t_e3", "t_e4", "a_afr_admix", "a_eur_admix", ]
    bounds = [PARAMETERS[p]["bounds"] for p in param_order]
    x0 = {p: random.uniform(*PARAMETERS[p]["bounds"]) for p in param_order}
    # Give these even admixture parameters to start, just so we don't accidentally start
    # somewhere near the "edge"
    x0["a_eur_admix"] = 0.4
    x0["a_afr_admix"] = 0.4

    if "INIT_TRUTH" in os.environ:
        print("Starting from the known correct values for x0")
        x0 = {p: PARAMETERS[p]["truth"] for p in param_order}

    if USE_SYMBOLIC_BOUNDS:
        idx = param_order.index("a_afr_admix")
        bounds[idx] = (ADMIX_LOWER_BOUND, "1-a_eur_admix")
        idx = param_order.index("a_eur_admix")
        bounds[idx] = (ADMIX_LOWER_BOUND, "1-a_afr_admix")

    print(f"========== STARTING CONDITIONS ==========")
    print(f"PARAMETER ORDER: {param_order}")
    print(f"BOUNDS: {bounds}")
    print(f"X0: {x0}")
    print(f"=========================================")

    x, logp = glike.maximize(fun, x0, bounds = bounds)
    print(f"FINAL x = {x}")
    print(f"FINAL logp = {logp}")
