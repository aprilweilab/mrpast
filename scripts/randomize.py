# Randomize the initialization values for a concrete model. If you already have a
# JSON file for solver input (produced by "mrpast process") then you can randomize
# it via "python randomize.py solver_in.json > new_solver_in.json"
import sys
import json
import random
import itertools

if __name__ == "__main__":
    model_file = sys.argv[1]

    with open(model_file) as f:
        model = json.load(f)

    parameter_keys = (
        "epoch_times_gen",
        "smatrix_values_ne__gen",
        "amatrix_parameters",
    )
    all_params = itertools.chain.from_iterable(
        map(lambda k: model.get(k, []) or [], parameter_keys)
    )

    for p in all_params:
        if p["lb"] == p["ub"]:
            p["init"] = p["lb"]
        else:
            p["init"] = random.uniform(p["lb"], p["ub"])

    print(json.dumps(model, indent=2))
