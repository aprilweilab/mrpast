import json
import sys
import os


def get_type(param):
    if "migration rate" in param["description"].lower():
        return "migration"
    elif "coalescence rate" in param["description"].lower():
        return "Ne"
    elif "growth rate" in param["description"].lower():
        return "growth"
    else:
        return "epoch"


def get_name(filename):
    if "6deme2epoch.split" in filename:
        return "6D2E"
    elif "6deme1epoch" in filename:
        return "6D1E"
    elif "6deme2epoch" in filename:
        return "6D2E (old)"
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


PLOIDY = 2


def cr2ne(rate):
    return 1 / (PLOIDY * rate)


DO_ABS = True


def diff(a, b):
    if DO_ABS:
        return abs(a - b)
    return a - b


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: relerr.py <output files>")
        exit(1)
    print("WARNING: Assumes ploidy=2", file=sys.stderr)

    output = []
    for fn in sys.argv[1:]:
        error_list = []
        ar_err = 0.0
        ar_err_mig = 0.0
        ar_err_coal = 0.0
        ar_err_ne = 0.0
        with open(fn) as f:
            data = json.load(f)
        all_params = (data.get("epoch_times_gen", []) or []) + (
            data.get("smatrix_values_ne__gen", []) or {}
        )
        num_migs = 0
        num_coals = 0
        for param in all_params:
            rel_err = (
                diff(param["final"], param["ground_truth"]) / param["ground_truth"]
            )
            ar_err += rel_err
            if "migration rate" in param["description"].lower():
                error_list.append(rel_err)
                ar_err_mig += rel_err
                num_migs += 1
            elif "coalescence rate" in param["description"].lower():
                ar_err_coal += rel_err
                ne_err = diff(
                    cr2ne(param["final"]), cr2ne(param["ground_truth"])
                ) / cr2ne(param["ground_truth"])
                ar_err_ne += ne_err
                error_list.append(ne_err)
                num_coals += 1
            else:
                error_list.append(rel_err)
        for i in range(len(all_params)):
            output.append(
                {
                    "file": os.path.basename(fn),
                    "parameter": i,
                    "relative_error": error_list[i],
                    "type": get_type(all_params[i]),
                    "name": get_name(fn),
                    "avg_rel_err": ar_err / len(all_params),
                    "mig_rel_err": ar_err_mig / num_migs,
                    "coal_rel_err": ar_err_coal / num_coals,
                    "ne_rel_err": ar_err_ne / num_coals,
                }
            )
    print(json.dumps(output, indent=2))
