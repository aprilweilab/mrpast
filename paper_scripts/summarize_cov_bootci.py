import json
import sys
import pandas as pd
from mrpast.result import summarize_bootstrap_data
import os


def get_name(filename):
    if "6deme2epoch.split" in filename:
        return "6D2E"
    elif "6deme2epoch" in filename:
        return "6D2E (old)"
    elif "6deme1epoch" in filename:
        return "6D1E"
    elif "5deme2epoch.diffcoal" in filename:
        return "5D2E (B)"
    elif "5deme2epoch" in filename:
        return "5D2E (A)"
    elif "5deme1epoch" in filename:
        return "5D1E"
    elif "20deme1epoch" in filename:
        return "20D1E"
    elif "ooa_2" in filename:
        return "OOA2"
    elif "ooa_3" in filename:
        return "OOA3"
    elif "ooa_4" in filename:
        return "OOA4"
    return "UNKNOWN"


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(
            "Usage: summarize_cov_bootci.py <input type> <... multiple JSON files with GIM output>"
        )
        exit(1)

    generate_data = [
        ("BOOTCI", 0.99),
        ("BOOTCI", 0.95),
        ("BOOTCI", 0.9),
        ("BOOTCI", 0.75),
    ]
    result = []
    for method, interval_conf in generate_data:
        within = 0
        total_params = 0
        by_model = {}
        # Relative to ground truth.
        excess = {}
        deficiency = {}
        size = {}
        for fn in sys.argv[1:]:
            df = pd.read_csv(fn)
            df = summarize_bootstrap_data(df, interval_conf=interval_conf)

            within_model = 0
            within_excess = []
            within_deficiency = []
            within_size_rel = []
            num_params = len(df)
            total_params += num_params
            for i, row in df.iterrows():
                gt = row["Ground Truth"]
                value = row["Optimized Value"]
                ci = (value - row["err_low"], value + row["err_low"])
                if row["covered"]:
                    within_model += 1
                    within += 1
                    d = max(ci[1] - gt, gt - ci[0])
                    assert d >= 0
                    within_excess.append(d / gt)
                else:
                    d = max(ci[0] - gt, gt - ci[1])
                    assert d >= 0
                    within_deficiency.append(d / gt)
                within_size_rel.append((ci[1] - ci[0]) / gt)
            name = get_name(fn)
            result.append(
                {
                    "Method": method,
                    "Interval": interval_conf,
                    "Model": name,
                    "CI Coverage": within_model / num_params,
                }
            )
            for e in within_excess:
                result.append(
                    {
                        "Method": method,
                        "Interval": interval_conf,
                        "Model": name,
                        "CI Relative Excess": e,
                    }
                )
            for d in within_deficiency:
                result.append(
                    {
                        "Method": method,
                        "Interval": interval_conf,
                        "Model": name,
                        "CI Relative Deficiency": d,
                    }
                )
            for s in within_size_rel:
                result.append(
                    {
                        "Method": method,
                        "Interval": interval_conf,
                        "Model": name,
                        "CI Relative Size": s,
                    }
                )
        print(f"Method: {method} @ {interval_conf*100}%", file=sys.stderr)
        print(f"Total parameters: {total_params}", file=sys.stderr)
        print(f"Within CI: {within}", file=sys.stderr)
        print(f"Coverage: {(within/total_params)*100.0}%", file=sys.stderr)
        print("", file=sys.stderr)
    print(json.dumps(result, indent=2))
