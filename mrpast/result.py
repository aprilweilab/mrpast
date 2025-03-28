# Migration Rate and Population Size Across Space and Time (mrpast)
# Copyright (C) 2025 April Wei
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# with this program.  If not, see <https://www.gnu.org/licenses/>.
from tabulate import tabulate
from typing import Optional, Dict, Any, List
import json
import math
import mrpast.model
import numpy as np
import numpy as np
import pandas as pd

try:
    import networkx as nx
    import matplotlib.pyplot as plt
except ImportError:
    nx = None  # type: ignore
    plt = None  # type: ignore


def load_json_pandas(
    filename: str, interval_field: Optional[str] = None
) -> pd.DataFrame:
    """
    Load a solver output JSON file as a Pandas DataFrame.

    :param filename: The JSON filename.
    :param interval_field: Optionally, the name of a field in the JSON file (on each
        parameter) to use for computing the parameter confidence intervals. Typically
        this is "gim_ci", which is present on an output file if the "mrpast confidence"
        command was used to generate the JSON.
    :return: Pandas DataDrame, where coalescent rates have been converted into effective
        population sizes (Ne).
    """
    result = []
    with open(filename) as f:
        data = json.load(f)
    ploidy = data["ploidy"]

    def clamp(param, value):
        if value < param["lb"]:
            return param["lb"]
        if value > param["ub"]:
            return param["ub"]
        return value

    def get_interval(param, idx):
        if interval_field is not None:
            if isinstance(interval_field, str):
                return param[interval_field][idx]
            return param[interval_field[idx]]
        return float("NaN")

    epochs = data.get("epoch_times_gen")
    for i, p in enumerate(epochs if epochs else []):
        v = clamp(p, p["final"])
        del p["apply_to"]
        p.update(
            {
                "label": f"E{i}",
                "Ground Truth": clamp(p, p["ground_truth"]),
                "err_low": v - clamp(p, get_interval(p, 0)),
                "err_hi": clamp(p, get_interval(p, 1)) - v,
                "Optimized Value": v,
                "Parameter Type": "Epoch time",
            }
        )
        result.append(p)
    mcounter = 0
    ncounter = 0
    gcounter = 0
    for p in data["smatrix_values_ne__gen"]:
        if p["description"].startswith("Migration rate"):
            v = clamp(p, p["final"])
            epochs = list(sorted(set([a.get("epoch") for a in p["apply_to"]])))
            del p["apply_to"]
            p.update(
                {
                    "label": f"M{mcounter}",
                    "Ground Truth": clamp(p, p["ground_truth"]),
                    "err_low": v - clamp(p, get_interval(p, 0)),
                    "err_hi": clamp(p, get_interval(p, 1)) - v,
                    "Optimized Value": v,
                    "Parameter Type": "Migration rate",
                    "Epochs": epochs,
                }
            )
            result.append(p)
            mcounter += 1
        elif p["description"].startswith("Growth rate"):
            v = clamp(p, p["final"])
            epochs = list(sorted(set([a.get("epoch") for a in p["apply_to"]])))
            del p["apply_to"]
            p.update(
                {
                    "label": f"G{gcounter}",
                    "Ground Truth": clamp(p, p["ground_truth"]),
                    "err_low": v - clamp(p, get_interval(p, 0)),
                    "err_hi": clamp(p, get_interval(p, 1)) - v,
                    "Optimized Value": v,
                    "Parameter Type": "Growth rate",
                    "Epochs": epochs,
                }
            )
            result.append(p)
            gcounter += 1
        elif p["description"].startswith("Coalescence rate"):

            def coal2ne(rate):
                return 1 / (ploidy * clamp(p, rate))

            epochs = list(sorted(set([a.get("epoch") for a in p["apply_to"]])))
            del p["apply_to"]
            v = coal2ne(p["final"])
            # These are flipped because of coal2ne...
            lower_ci = coal2ne(get_interval(p, 1))
            upper_ci = coal2ne(get_interval(p, 0))
            p.update(
                {
                    "label": f"P{ncounter}",
                    "Ground Truth": coal2ne(p["ground_truth"]),
                    "err_low": v - lower_ci,
                    "err_hi": upper_ci - v,
                    "Optimized Value": v,
                    "Parameter Type": "Effective popsize",
                    "Epochs": epochs,
                }
            )
            result.append(p)
            ncounter += 1
    return pd.DataFrame.from_dict(result)


def summarize_bootstrap_data(
    bootstrap_df: pd.DataFrame,
    use_median: bool = True,
    interval_conf: float = 0.95,
) -> pd.DataFrame:
    """
    Given a Pandas DataFrame loaded from a bootstrap CSV file, produce a new DataFrame
    that summarizes the data. Confidence intervals are calculated, and the resulting
    parameter estimates are the mean (or median) of the value over all bootstrap samples.

    :param bootstrap_df: A DataFrame as loaded via pandas.read_csv() with the CSV that is
        generated by "mrpast confidence --bootstrap".
    :param use_median: Set to False if you want to use the mean instead of the median for
        summarizing parameter values over all bootstrap samples. The median is more robust
        to parameter estimates that hit the lower or upper bounds during maximum likelihood
        estimation.
    :param interval_conf: Defaults to 0.95. Set to one of 0.99, 0.95, 0.9, 0.75 or 1.0.
        1.0 means use the entire range of the bootstrap values instead of the standard
        deviation plus a confidence interval. The other values are the confidence for the
        normal distribution confidence intervals based on sample standard deviation.
    :return: DataFrame with one row per parameter, summarizing the value and confidence
        interval.
    """
    ci_mult = {0.99: 2.576, 0.95: 1.96, 0.9: 1.645, 0.75: 1.150}.get(interval_conf)
    assert (
        ci_mult is not None or interval_conf == 1.0
    ), f"Unsupported confidence {interval_conf}; try 1.0, 0.99, 0.95, 0.9, or 0.75"

    def get_singular(df, label, field):
        value = set(df[df["label"] == label][field])
        assert len(value) == 1
        return list(value)[0]

    new_data = []
    for label in set(bootstrap_df["label"]):
        truth = get_singular(bootstrap_df, label, "Ground Truth")
        values = bootstrap_df[bootstrap_df["label"] == label]["Optimized Value"]
        median = np.median(values)
        mean = np.average(values)
        value = median if use_median else mean
        if ci_mult is None:
            c95_low = value - min(values)
            c95_hi = max(values) - value
        else:
            # numpy defaults to population stddev, so set degrees of freedom to 1
            std_err = np.std(values, ddof=1)
            c95_low = max(0, (ci_mult * std_err))
            c95_hi = max(0, (ci_mult * std_err))
        new_data.append(
            {
                "label": label,
                "kind": get_singular(bootstrap_df, label, "kind"),
                "description": get_singular(bootstrap_df, label, "description"),
                "Parameter Type": get_singular(bootstrap_df, label, "Parameter Type"),
                "Ground Truth": truth,
                "Optimized Value": value,
                "err_low": c95_low,
                "err_hi": c95_hi,
                "min": min(values),
                "max": max(values),
                "covered": (truth >= (value - c95_low) and truth <= (value + c95_hi)),
                "param_index": get_singular(bootstrap_df, label, "Unnamed: 0"),
            }
        )
    return pd.DataFrame.from_dict(new_data).sort_values("param_index")


def draw_graphs(
    model_file: str,
    ax,
    grid_cols: Optional[int] = None,
    epoch_spacing: float = 0.5,
    epoch_label_spacing: float = 0.15,
    max_node_size: int = 800,
    migrate_color: Optional[str] = None,
    popsize_color: Optional[str] = None,
    x_offset: float = 0.25,
):
    """
    Draw the topology of the given input mrpast model file on the given matplotlib axis.

    :param model_file: The mrpast model filename.
    :param ax: The matplotlib axis object.
    :param grid_cols: Default is None. When set to an integer, layout the graphs in a grid
        with the given number of columns. For example, if you have a 6-deme model then you
        might want to set grid_cols=2 or grid_cols=3 to lay the graph out as 3x2 or 2x3.
    :param epoch_spacing: Spacing between each epoch in the figure.
    :param epoch_label_spacing: Spacing between the epoch label and the epoch graph.
    :param max_node_size: The maximum size that a particular node can be.
    :param migrate_color: The color to use for migration edges. By default, a spectrum of
        colors is used which indicates a higher (darker color) or lower (lighter color) rate.
    :param popsize_color: The color to use for deme nodes. By default, the
        matplotlib.pyplot.cm.Dark2 colormap is used.
    :param x_offset: The offset from the X-axis to start drawing.
    """
    assert (
        nx is not None and plt is not None
    ), "Plotting requires networkx and matplotlib; run 'pip install networkx matplotlib'"
    G = nx.DiGraph()
    model_config = mrpast.model.load_model_config(model_file)
    ploidy = model_config["ploidy"]
    epochs = len(model_config["coalescence"]["vectors"])
    base_node_id = 0
    node_sizes = []
    node_colors = []
    y_offset = 0.0
    pos: Optional[Dict[Any, Any]] = None

    max_popsize = 0
    for epoch in reversed(range(epochs)):
        for param_idx in model_config["coalescence"]["vectors"][epoch]:
            if param_idx > 0:
                param = model_config["coalescence"]["parameters"][param_idx - 1]
                pop_size = 1 / (ploidy * param["ground_truth"])
                if pop_size > max_popsize:
                    max_popsize = pop_size

    for epoch in reversed(range(epochs)):
        coal_vector = model_config["coalescence"]["vectors"][epoch]
        mig_matrix = model_config["migration"]["matrices"][epoch]

        num_demes = len(coal_vector)
        pop_sizes = [0 for _ in range(num_demes)]

        for i in range(num_demes):
            param_idx = coal_vector[i]
            if param_idx > 0:
                param = model_config["coalescence"]["parameters"][param_idx - 1]
                pop_sizes[i] = 1 / (ploidy * param["ground_truth"])
        nodes = [i for i in range(num_demes) if pop_sizes[i] > 0]
        node_sizes.extend(
            [max(0, (pop_sizes[i] / max_popsize) * max_node_size) for i in nodes]
        )
        node_colors.extend(nodes)

        for i in nodes:
            G.add_node(base_node_id + i)
            for j in range(num_demes):
                param_idx = mig_matrix[i][j]
                if param_idx > 0:
                    param = model_config["migration"]["parameters"][param_idx - 1]
                    G.add_edge(
                        base_node_id + i,
                        base_node_id + j,
                        weight=param["ground_truth"] * 5,
                    )

        if grid_cols is not None:
            if pos is None:
                pos = {}
            pos.update(
                {
                    base_node_id
                    + node: (x_offset + (i % grid_cols), y_offset - (i // grid_cols))
                    for i, node in enumerate(nodes)
                }
            )
            ax.text(x_offset, y_offset + epoch_label_spacing, f"Epoch {epoch}")
            y_offset -= (len(nodes) // grid_cols) + epoch_spacing
        else:
            pos = None
        base_node_id += num_demes

    weights = [G[u][v]["weight"] for u, v in G.edges()]
    weights = [
        math.log10(w) for w in weights
    ]  # Log-scale the weights, to differentiate better
    if popsize_color is None:
        node_options = {
            "node_size": node_sizes,
            "node_color": node_colors,
            "cmap": plt.cm.Dark2,  # type: ignore
            "linewidths": 2,
        }
    else:
        node_options = {
            "node_size": node_sizes,
            "node_color": popsize_color,
            "linewidths": 2,
        }
    if migrate_color is None:
        edge_options = {
            "edge_color": weights,
            "width": 2,
            "edge_vmin": min(weights)
            * 1.25,  # These are log scale, so this makes a smaller (more negative) value
            "edge_vmax": max(weights),
            "edge_cmap": plt.cm.Oranges,  # type: ignore
            "connectionstyle": "arc3,rad=0.1",
        }
    else:
        edge_options = {
            "edge_color": migrate_color,
            "width": 2,
            "connectionstyle": "arc3,rad=0.1",
        }
    print(f"Node sizes: {node_sizes}")
    print(f"Edge weights: {weights}")
    nx.draw_networkx_nodes(G, pos=pos, **node_options, ax=ax)
    nx.draw_networkx_edges(G, pos=pos, **edge_options, ax=ax)
    _ = plt.axis("off")


def get_matching_colors(num_demes, demes=[]):
    if demes:
        return {
            f"{demes[i]}": plt.cm.Dark2(c)
            for i, c in enumerate(np.linspace(0, 1, num_demes))
        }
    return {
        f"pop_{i}": plt.cm.Dark2(c) for i, c in enumerate(np.linspace(0, 1, num_demes))
    }


def tab_show(filename: str, sort_by: str = "Index"):
    """
    Print an ASCII table showing the parameter values and their error from ground truth, for the
    given JSON output from the solver.
    """
    with open(filename) as f:
        output = json.load(f)

    all_params: List[Dict[str, Any]] = (output.get("epoch_times_gen", []) or []) + (
        output.get("smatrix_values_ne__gen", []) or []
    )
    results = []
    total_rel = 0.0
    total_abs = 0.0
    for param_idx, param in enumerate(all_params):
        gt = param["ground_truth"]
        final = param["final"]
        abserr = abs(gt - final)
        total_abs += abserr
        relerr = abserr / gt
        total_rel += relerr
        epochs = set()
        for app in param["apply_to"]:
            epochs.add(app["epoch"])
        results.append(
            (
                param_idx,
                param["description"],
                relerr,
                abserr,
                gt,
                final,
                list(sorted(epochs)),
            )
        )

    headers = [
        "Index",
        "Description",
        "Relative Error",
        "Absolute Error",
        "Truth",
        "Final",
        "Epochs",
    ]

    try:
        sort_key = headers.index(sort_by)
    except ValueError:
        raise RuntimeError(f"Unexpected sort_by key: {sort_by}. Try one of {headers}.")
    results = sorted(results, key=lambda x: x[sort_key])
    print(tabulate(results, headers=headers))
    print()
    print(f"Total absolute error: {total_abs}")
    print(f"Total relative error: {total_rel}")

    return total_rel
