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

## Parameters: ground truth and bounds values.
PARAMETERS = {
 "c_afr_e2": {"truth": 6.839945280437756e-05, "bounds": (6.839945280437756e-06, 0.0006839945280437756)},
 "c_afr_e1": {"truth": 3.45447008428907e-05, "bounds": (3.45447008428907e-06, 0.000345447008428907)},
 "c_afr_e0": {"truth": 1.1570737191765288e-06, "bounds": (1.1570737191765288e-07, 1.1570737191765288e-05)},
 "g_afr_e0": {"truth": 0.0166, "bounds": (0.01, 0.025)},
 "t_e0": {"truth": 204.6, "bounds": (100.0, 400.0)},
 "t_e1": {"truth": 5920.0, "bounds": (3062.0, 7451.0)},
}

# Crank this up from the default of 100, because this is a pretty simple model (so its fast)
# and it might help improve results?
DISCRETIZE = 500

def CR(coal_rate):
    return coal_rate

def iicr_demography(c_afr_e2, c_afr_e1, c_afr_e0,
                    g_afr_e0,
                    t_e0, t_e1):
    demo = glike.Demo()
    # Epoch 0
    demo.add_phase(glike.Phase(
       0, t_e0, # Simulated Truth: 0 - 204.6
       [CR(c_afr_e0)],
       grs = [g_afr_e0],
       populations = ["AFR"]),
       discretize=DISCRETIZE)

    # Epoch 1
    c_afr_e1 = CR(c_afr_e0) * math.exp(t_e0 * g_afr_e0)
    demo.add_phase(glike.Phase(
       t_e0, t_e1, # Simulated Truth: 204.6 - 5920.0
       [CR(c_afr_e1)],
       grs = [0],
       populations = ["AFR"]),
       discretize=DISCRETIZE)
    # Epoch 2
    demo.add_phase(glike.Phase(
      t_e1, math.inf,
      [CR(c_afr_e2)],
      grs = [0],
      populations = ["AFR"]),
       discretize=DISCRETIZE)
    return demo


def iicr_demography_exp(c_afr_e2, c_afr_e0,
                        g_afr_e0,
                        t_e0, t_e1):
    demo = glike.Demo()
    # Epoch 0
    demo.add_phase(glike.Phase(
       0, t_e0, # Simulated Truth: 0 - 204.6
       [CR(c_afr_e0)],
       grs = [g_afr_e0],
       populations = ["AFR"]),
       discretize=DISCRETIZE)

    # Epoch 1
    c_afr_e1 = CR(c_afr_e0) * math.exp(t_e0 * g_afr_e0)
    demo.add_phase(glike.Phase(
       t_e0, t_e1, # Simulated Truth: 204.6 - 5920.0
       [CR(c_afr_e1)],
       grs = [0],
       populations = ["AFR"]),
       discretize=DISCRETIZE)
    # Epoch 2
    demo.add_phase(glike.Phase(
      t_e1, math.inf,
      [CR(c_afr_e2)],
      grs = [0],
      populations = ["AFR"]),
      discretize=DISCRETIZE)
    return demo

if __name__ == "__main__":
    if sys.argv[1] == "DEMO_DUMP":
        truth = {k: v["truth"] for k, v in PARAMETERS.items()}
        demo = iicr_demography(**truth)
        #demography = glike.demo_to_demography(demo)
        #print("##### msprime Demography (not accurate)")
        #print(demography.debug())
        #print()
        print("##### Raw Demography")
        demo.print()
        exit(0)
    tree_file = sys.argv[1]
    ts = tskit.load(tree_file)
    if "TREES10" in os.environ:
        SPACING = 2975000  # 10 trees (on simulated 100MB data)
    else:
        assert "TREES100" in os.environ
        SPACING = 297500  # 100 trees (on simulated 100MB data)
    trees = [ts.at(pos).copy() for pos in range(SPACING, 30000000, SPACING)]
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

    STATUS_EVERY = 100
    calls = 0
    if "USE_EXP" in os.environ:
        print("Using exponential to determine Ne after growth")
        def fun(c_afr_e2, c_afr_e0, g_afr_e0, t_e0, t_e1):
            demo = iicr_demography_exp(c_afr_e2, c_afr_e0, g_afr_e0, t_e0, t_e1)
            result = glike.glike_trees(trees, demo, samples=glike_samples)
            global calls
            calls += 1
            if calls % STATUS_EVERY == 1:
                print(f"After {calls} evaluations, best result is {result}")
            return result
        param_order = ["c_afr_e2", "c_afr_e0", "g_afr_e0", "t_e0", "t_e1"]
    else:
        def fun(c_afr_e2, c_afr_e1, c_afr_e0, g_afr_e0, t_e0, t_e1):
            demo = iicr_demography(c_afr_e2, c_afr_e1, c_afr_e0, g_afr_e0, t_e0, t_e1)
            result = glike.glike_trees(trees, demo, samples=glike_samples)
            global calls
            calls += 1
            if calls % STATUS_EVERY == 1:
                print(f"After {calls} evaluations, best result is {result}")
            return result
        param_order = ["c_afr_e2", "c_afr_e1", "c_afr_e0", "g_afr_e0", "t_e0", "t_e1"]

    bounds = [PARAMETERS[p]["bounds"] for p in param_order]
    random.seed(42)
    x0 = {p: random.uniform(*PARAMETERS[p]["bounds"]) for p in param_order}

    print(f"========== STARTING CONDITIONS ==========")
    print(f"PARAMETER ORDER: {param_order}")
    print(f"BOUNDS: {bounds}")
    print(f"X0: {x0}")
    print(f"=========================================")

    x, logp = glike.maximize(fun, x0, bounds = bounds)
    print(f"FINAL x = {x}")
    print(f"FINAL logp = {logp}")
