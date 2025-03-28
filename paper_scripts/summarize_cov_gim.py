import json
import sys


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
            "Usage: summarize_cov_gim.py <input type> <... multiple JSON files with GIM output>"
        )
        exit(1)

    def get_ci(param, method, interval):
        mult = {0.99: 2.576, 0.95: 1.96, 0.9: 1.645, 0.75: 1.150}[interval]
        se = param["gim_se"]
        return (param["final"] - (se * mult), param["final"] + (se * mult))

    generate_data = [("GIM", 0.99), ("GIM", 0.95), ("GIM", 0.9), ("GIM", 0.75)]

    result = []
    for method, interval in generate_data:
        within = 0
        total_params = 0
        by_model = {}
        # Relative to ground truth.
        excess = {}
        deficiency = {}
        size = {}
        for fn in sys.argv[1:]:
            within_model = 0
            within_excess = []
            within_deficiency = []
            within_size_rel = []
            with open(fn) as f:
                data = json.load(f)
            all_params = (data.get("epoch_times_gen", []) or []) + (
                data.get("smatrix_values_ne__gen", []) or {}
            )
            total_params += len(all_params)
            for param in all_params:
                ci = get_ci(param, method, interval)
                if param["ground_truth"] >= ci[0] and param["ground_truth"] <= ci[1]:
                    within_model += 1
                    within += 1
                    d = max(
                        ci[1] - param["ground_truth"], param["ground_truth"] - ci[0]
                    )
                    assert d >= 0
                    within_excess.append(d / param["ground_truth"])
                else:
                    d = max(
                        ci[0] - param["ground_truth"], param["ground_truth"] - ci[1]
                    )
                    assert d >= 0
                    within_deficiency.append(d / param["ground_truth"])
                within_size_rel.append((ci[1] - ci[0]) / param["ground_truth"])
            name = get_name(fn)
            result.append(
                {
                    "Method": method,
                    "Interval": interval,
                    "Model": name,
                    "CI Coverage": within_model / len(all_params),
                }
            )
            for e in within_excess:
                result.append(
                    {
                        "Method": method,
                        "Interval": interval,
                        "Model": name,
                        "CI Relative Excess": e,
                    }
                )
            for d in within_deficiency:
                result.append(
                    {
                        "Method": method,
                        "Interval": interval,
                        "Model": name,
                        "CI Relative Deficiency": d,
                    }
                )
            for s in within_size_rel:
                result.append(
                    {
                        "Method": method,
                        "Interval": interval,
                        "Model": name,
                        "CI Relative Size": s,
                    }
                )
        print(f"Method: {method} @ {interval*100}%", file=sys.stderr)
        print(f"Total parameters: {total_params}", file=sys.stderr)
        print(f"Within CI: {within}", file=sys.stderr)
        print(f"Coverage: {(within/total_params)*100.0}%", file=sys.stderr)
        print("", file=sys.stderr)
    print(json.dumps(result, indent=2))
