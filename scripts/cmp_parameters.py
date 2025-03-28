# Rank parameters from best to worst, based on their percentage error. Display
# both percent and absolute error.
import json
import sys
from tabulate import tabulate

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: cmp_params.py <JSON file1> <JSON file2>")
        exit(1)
    with open(sys.argv[1]) as f:
        output1 = json.load(f)
    with open(sys.argv[2]) as f:
        output2 = json.load(f)

    all_params1 = (output1.get("epoch_times_gen", []) or []) + (
        output1.get("smatrix_values_ne__gen", []) or {}
    )
    all_params2 = (output2.get("epoch_times_gen", []) or []) + (
        output2.get("smatrix_values_ne__gen", []) or {}
    )
    assert len(all_params1) == len(all_params2)
    total_err1 = 0.0
    total_err2 = 0.0
    results = []
    for param_idx, (param1, param2) in enumerate(zip(all_params1, all_params2)):
        gt = param1["ground_truth"]
        assert param2["ground_truth"] == gt
        final1 = param1["final"]
        final2 = param2["final"]
        relerr1 = abs(gt - final1) / gt
        total_err1 += relerr1
        relerr2 = abs(gt - final2) / gt
        total_err2 += relerr2
        epochs = set()
        for app in param1["apply_to"]:
            epochs.add(app["epoch"])
        results.append(
            (
                param_idx,
                param1["description"],
                relerr1,
                relerr2,
                gt,
                final1,
                final2,
                list(sorted(epochs)),
            )
        )
    results = sorted(results, key=lambda x: x[3])

    print(
        tabulate(
            results,
            headers=[
                "Index",
                "Description",
                "Relative Error (1)",
                "Relative Error (2)",
                "Truth",
                "Final (1)",
                "Final (2)",
                "Epochs",
            ],
        )
    )
    print()
    print(f"Error1 = {total_err1} vs. Error2 = {total_err2}")
