# Randomize the initialization values for a concrete model. If you already have a
# JSON file for solver input (produced by "mrpast process") then you can randomize
# it via "python randomize.py solver_in.json > new_solver_in.json"
import sys
import json
import random

if __name__ == "__main__":
    model_file = sys.argv[1]

    with open(model_file) as f:
        model = json.load(f)

    model_epochs = model["epoch_times_gen"] or []
    model_params = model["smatrix_values_ne__gen"]

    for me in model_epochs:
        me["init"] = random.uniform(me["lb"], me["ub"])
    for mp in model_params:
        mp["init"] = random.uniform(mp["lb"], mp["ub"])

    print(json.dumps(model, indent=2))
