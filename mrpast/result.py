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
from mrpast.model import UserModel
from tabulate import tabulate
from typing import Optional, Dict, Any, List, Iterable, Tuple
import copy
import itertools
import json
import math
import mrpast.model
import numpy
import os
import pandas
import sys

try:
    import networkx as nx
except ImportError:
    nx = None  # type: ignore

try:
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import seaborn as sns
except ImportError:
    plt = None  # type: ignore
    mpl = None  # type: ignore
    sns = None  # type: ignore


# b is within 0.1% of a
def _is_nearly(a, b, relerr=0.001):
    return abs(a - b) / a <= relerr


# A fixed parameter is a constant value specified by the user, and not used by the solver.
def _param_is_fixed(parameter: Dict[str, Any]) -> bool:
    return parameter["lb"] == parameter["ub"]


# A synthetic parameter is either fixed or constructed, and is invisible to the user.
def _param_is_synthetic(parameter: Dict[str, Any]) -> bool:
    return parameter["kind_index"] < 0


# A constructed parameter is a user-visible parameter that is completely determined by the
# value(s) of other parameters.
def _param_is_constructed(parameter: Dict[str, Any]) -> bool:
    return len(parameter.get("one_minus", [])) > 0


def _clamp(param, value):
    if value < param["lb"]:
        return param["lb"]
    if value > param["ub"]:
        return param["ub"]
    return value


def load_json_pandas(
    filename: str, interval_field: Optional[str] = None, skip_fixed: bool = True
) -> pandas.DataFrame:
    """
    Load a solver output JSON file as a Pandas DataFrame.

    :param filename: The JSON filename.
    :param interval_field: Optionally, the name of a field in the JSON file (on each
        parameter) to use for computing the parameter confidence intervals. Typically
        this is "gim_ci", which is present on an output file if the "mrpast confidence"
        command was used to generate the JSON.
    :param skip_fixed: Set to False to keep the fixed values that were part of the solution, otherwise
        only the parameters will be returned.
    :return: Pandas DataDrame, where coalescent rates have been converted into effective
        population sizes (Ne).
    """
    result = []
    with open(filename) as f:
        data = json.load(f)
    ploidy = data["ploidy"]

    def get_interval(param, idx):
        if interval_field is not None and not (
            _param_is_fixed(param) or _param_is_constructed(param)
        ):
            if isinstance(interval_field, str):
                return param[interval_field][idx]
            return param[interval_field[idx]]
        return float("NaN")

    epochs = data.get("epoch_times_gen")
    for i, p in enumerate(epochs if epochs else []):
        v = _clamp(p, p["final"])
        del p["apply_to"]
        is_fixed = _param_is_fixed(p)
        if is_fixed and skip_fixed:
            continue
        truth = _clamp(p, p["ground_truth"])
        err_low = v - _clamp(p, get_interval(p, 0))
        err_hi = _clamp(p, get_interval(p, 1)) - v
        p.update(
            {
                "label": f"E{i}",
                "Ground Truth": truth,
                "err_low": err_low,
                "err_hi": err_hi,
                "Optimized Value": v,
                "Parameter Type": "Epoch time",
                "Fixed": is_fixed,
                "Lower Bound": p["lb"],
                "Upper Bound": p["ub"],
                "covered": (truth >= (v - err_low) and truth <= (v + err_hi)),
                "Lower Bound": p["lb"],
                "Upper Bound": p["ub"],
            }
        )
        result.append(p)
    mcounter = 0
    ncounter = 0
    gcounter = 0
    for p in data["smatrix_values_ne__gen"]:
        is_fixed = _param_is_fixed(p)
        if is_fixed and skip_fixed:
            continue
        if p["kind"] == "migration":
            v = _clamp(p, p["final"])
            epochs = list(sorted(set([a.get("epoch") for a in p["apply_to"]])))
            del p["apply_to"]
            truth = _clamp(p, p["ground_truth"])
            err_low = v - _clamp(p, get_interval(p, 0))
            err_hi = _clamp(p, get_interval(p, 1)) - v
            p.update(
                {
                    "label": f"M{mcounter}",
                    "Ground Truth": truth,
                    "err_low": err_low,
                    "err_hi": err_hi,
                    "Optimized Value": v,
                    "Parameter Type": "Migration rate",
                    "Epochs": epochs,
                    "Fixed": is_fixed,
                    "Lower Bound": p["lb"],
                    "Upper Bound": p["ub"],
                    "covered": (truth >= (v - err_low) and truth <= (v + err_hi)),
                    "Lower Bound": p["lb"],
                    "Upper Bound": p["ub"],
                }
            )
            result.append(p)
            mcounter += 1
        elif p["kind"] == "growth":
            v = _clamp(p, p["final"])
            epochs = list(sorted(set([a.get("epoch") for a in p["apply_to"]])))
            del p["apply_to"]
            truth = _clamp(p, p["ground_truth"])
            err_low = v - _clamp(p, get_interval(p, 0))
            err_hi = _clamp(p, get_interval(p, 1)) - v
            p.update(
                {
                    "label": f"G{gcounter}",
                    "Ground Truth": truth,
                    "err_low": err_low,
                    "err_hi": err_hi,
                    "Optimized Value": v,
                    "Parameter Type": "Growth rate",
                    "Epochs": epochs,
                    "Fixed": is_fixed,
                    "Lower Bound": p["lb"],
                    "Upper Bound": p["ub"],
                    "covered": (truth >= (v - err_low) and truth <= (v + err_hi)),
                    "Lower Bound": p["lb"],
                    "Upper Bound": p["ub"],
                }
            )
            result.append(p)
            gcounter += 1
        elif p["kind"] == "coalescence":

            def coal2ne(rate):
                return 1 / (ploidy * _clamp(p, rate))

            epochs = list(sorted(set([a.get("epoch") for a in p["apply_to"]])))
            del p["apply_to"]
            v = coal2ne(p["final"])
            # These are flipped because of coal2ne...
            lower_ci = coal2ne(get_interval(p, 1))
            upper_ci = coal2ne(get_interval(p, 0))
            truth = coal2ne(p["ground_truth"])
            err_low = v - lower_ci
            err_hi = upper_ci - v
            p.update(
                {
                    "label": f"P{ncounter}",
                    "Ground Truth": truth,
                    "err_low": err_low,
                    "err_hi": err_hi,
                    "Optimized Value": v,
                    "Parameter Type": "Effective popsize",
                    "Epochs": epochs,
                    "Fixed": is_fixed,
                    "Lower Bound": coal2ne(p["ub"]),
                    "Upper Bound": coal2ne(p["lb"]),
                    "covered": (truth >= (v - err_low) and truth <= (v + err_hi)),
                    "Lower Bound": coal2ne(p["ub"]),
                    "Upper Bound": coal2ne(p["lb"]),
                }
            )
            result.append(p)
            ncounter += 1
    for acounter, p in enumerate(data.get("amatrix_parameters", []) or []):
        v = _clamp(p, p["final"])
        del p["apply_to"]
        is_fixed = _param_is_fixed(p)
        if is_fixed and skip_fixed:
            continue
        truth = _clamp(p, p["ground_truth"]) if not is_fixed else p["init"]
        err_low = v - _clamp(p, get_interval(p, 0))
        err_hi = _clamp(p, get_interval(p, 1)) - v
        p.update(
            {
                "label": f"A{acounter}",
                "Ground Truth": truth,
                "err_low": err_low,
                "err_hi": err_hi,
                "Optimized Value": v,
                "Parameter Type": "Admixture proportion",
                "Epochs": [],  # TODO
                "Fixed": is_fixed,
                "Lower Bound": p["lb"],
                "Upper Bound": p["ub"],
                "covered": (truth >= (v - err_low) and truth <= (v + err_hi)),
                "Lower Bound": p["lb"],
                "Upper Bound": p["ub"],
            }
        )
        if "one_minus" in p:
            del p["one_minus"]
        result.append(p)

    return pandas.DataFrame.from_dict(result)


def summarize_bootstrap_data(
    bootstrap_df: pandas.DataFrame,
    use_median: bool = True,
    interval_conf: float = 0.95,
    use_percentile: bool = False,
) -> pandas.DataFrame:
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
    :param use_percentile: Instead of computing the sample stddev, just use the percentile
        from the bootstraps for each parameter.
    :return: DataFrame with one row per parameter, summarizing the value and confidence
        interval.
    """
    if not use_percentile:
        ci_mult = {0.99: 2.576, 0.95: 1.96, 0.9: 1.645, 0.75: 1.150}.get(interval_conf)
        assert (
            ci_mult is not None or interval_conf == 1.0
        ), f"Unsupported confidence {interval_conf}; try 1.0, 0.99, 0.95, 0.9, or 0.75"

    def get_singular(df, label, field):
        items = df[df["label"] == label][field]
        # We string-ify lists so that we can uniquify them
        if items.dtype == object:
            items = items.astype(str)
        value = set(items)
        assert len(value) == 1
        return list(value)[0]

    new_data = []
    for label in set(bootstrap_df["label"]):
        truth = get_singular(bootstrap_df, label, "Ground Truth")
        values = bootstrap_df[bootstrap_df["label"] == label]["Optimized Value"]
        median = numpy.median(values)
        mean = numpy.average(values)
        value = median if use_median else mean
        if use_percentile:
            err_low, err_hi = numpy.percentile(
                values, [1 - interval_conf, interval_conf]
            )
            err_low = value - err_low
            err_hi = err_hi - value
            std_err = None
        elif ci_mult is None:
            err_low = value - min(values)
            err_hi = max(values) - value
            std_err = None
        else:
            # numpy defaults to population stddev, so set degrees of freedom to 1
            std_err = numpy.std(values, ddof=1)
            err_low = max(0, (ci_mult * std_err))
            err_hi = max(0, (ci_mult * std_err))
        new_data.append(
            {
                "label": label,
                "kind": get_singular(bootstrap_df, label, "kind"),
                "description": get_singular(bootstrap_df, label, "description"),
                "Parameter Type": get_singular(bootstrap_df, label, "Parameter Type"),
                "Ground Truth": truth,
                "Optimized Value": value,
                "err_low": err_low,
                "err_hi": err_hi,
                "min": min(values),
                "max": max(values),
                "std": std_err,
                "Epochs": get_singular(bootstrap_df, label, "Epochs"),
                "covered": (truth >= (value - err_low) and truth <= (value + err_hi)),
            }
        )
        if "Upper Bound" in bootstrap_df.columns:
            new_data[-1]["Upper Bound"] = get_singular(
                bootstrap_df, label, "Upper Bound"
            )
        if "Lower Bound" in bootstrap_df.columns:
            new_data[-1]["Lower Bound"] = get_singular(
                bootstrap_df, label, "Lower Bound"
            )
    return pandas.DataFrame.from_dict(new_data).sort_values("label")


def draw_graphs(
    model_file: str,
    ax,
    grid_cols: Optional[int] = None,
    epoch_spacing: float = 0.5,
    epoch_label_spacing: Optional[float] = 0.15,
    max_node_size: int = 800,
    migrate_color: Optional[str] = None,
    popsize_color: Optional[str] = None,
    x_offset: float = 0.25,
    coal_values: Optional[List[float]] = None,
    mig_values: Optional[List[float]] = None,
    cax=None,
    cmap=None,
    min_max_migrate: Optional[Tuple[float, float]] = None,
):
    """
    Draw the topology of the given input mrpast model file on the given matplotlib axis.

    :param model_file: The mrpast model filename.
    :param ax: The matplotlib axis object.
    :param grid_cols: Default is None. When set to an integer, layout the graphs in a grid
        with the given number of columns. For example, if you have a 6-deme model then you
        might want to set grid_cols=2 or grid_cols=3 to lay the graph out as 3x2 or 2x3.
    :param epoch_spacing: Spacing between each epoch in the figure.
    :param epoch_label_spacing: Spacing between the epoch label and the epoch graph. Set to
        None to disable epoch labels.
    :param max_node_size: The maximum size that a particular node can be.
    :param migrate_color: The color to use for migration edges. By default, a spectrum of
        colors is used which indicates a higher (darker color) or lower (lighter color) rate.
    :param popsize_color: The color to use for deme nodes. By default, the
        matplotlib.pyplot.cm.Dark2 colormap is used.
    :param x_offset: The offset from the X-axis to start drawing.
    :param coal_values: If non-None, use this list of coalescence rate values instead of the
        ground-truth values from the model. This only works if the model has densely packed
        parameters, with no parameter index gaps.
    :param mig_values: If non-None, use this list of migration rate values instead of the
        ground-truth values from the model. This only works if the model has densely packed
        parameters, with no parameter index gaps.
    :param min_max_migrate: Optional tuple of (min, max) float values, which are the minimum
        and maximum possible migration rate values for the purposes of coloring the edges.
    """
    assert (
        nx is not None and plt is not None
    ), "Plotting requires networkx and matplotlib; run 'pip install networkx matplotlib'"
    if cmap is None:
        cmap = plt.cm.RdYlBu  # type: ignore
    G = nx.DiGraph()
    model = mrpast.model.UserModel.from_file(model_file)
    base_node_id = 0
    node_sizes = []
    node_colors = []
    y_offset = 0.0
    pos: Optional[Dict[Any, Any]] = None

    # FIXME: neither of these functions work if there are gaps in the parameter indexing.
    def get_coal_value(coal_param_idx):
        if coal_values is not None:
            return coal_values[coal_param_idx - 1]
        return model.coalescence.get_parameter(coal_param_idx).ground_truth

    def get_mig_value(mig_param_idx):
        if mig_values is not None:
            return mig_values[mig_param_idx - 1]
        return model.migration.get_parameter(mig_param_idx).ground_truth

    # Use the average of the ground truth values as our normalizing factor.
    avg_mig = 0.0
    for p in model.migration.parameters:
        avg_mig += p.ground_truth
    avg_mig /= len(model.migration.parameters)

    def norm_mig(mr):
        return math.log10(mr / avg_mig)

    max_popsize = 0
    for entry in model.coalescence.entries:
        if isinstance(entry.rate, mrpast.model.ParamRef):
            cr = get_coal_value(entry.rate.param)
        else:
            cr = entry.rate
        pop_size = 1 / (model.ploidy * cr)
        if pop_size > max_popsize:
            max_popsize = pop_size

    for epoch in reversed(range(model.num_epochs)):
        pop_sizes = [0 for _ in range(model.num_demes)]
        for i in range(model.num_demes):
            entry = model.coalescence.get_entry(epoch, i)
            if entry is not None:
                if isinstance(entry.rate, mrpast.model.ParamRef):
                    cr = get_coal_value(entry.rate.param)
                else:
                    cr = entry.rate
                pop_sizes[i] = 1 / (model.ploidy * cr)

        nodes = [i for i in range(model.num_demes) if pop_sizes[i] > 0]
        node_sizes.extend(
            [max(0, (pop_sizes[i] / max_popsize) * max_node_size) for i in nodes]
        )
        node_colors.extend(nodes)

        for i in nodes:
            G.add_node(base_node_id + i)
            for j in range(model.num_demes):
                entry = model.migration.get_entry(epoch, i, j)
                if entry is not None:
                    if isinstance(entry.rate, mrpast.model.ParamRef):
                        w = norm_mig(get_mig_value(entry.rate.param))
                    else:
                        w = norm_mig(entry.rate)
                    G.add_edge(
                        base_node_id + i,
                        base_node_id + j,
                        weight=w,
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
            if epoch_label_spacing is not None:
                ax.text(x_offset, y_offset + epoch_label_spacing, f"Epoch {epoch}")
            y_offset -= (len(nodes) // grid_cols) + epoch_spacing
        else:
            pos = None
        base_node_id += model.num_demes

    weights = [G[u][v]["weight"] for u, v in G.edges()]
    min_edge_color = (
        norm_mig(min_max_migrate[0]) if min_max_migrate else min(weights)
    ) * 1.25
    max_edge_color = norm_mig(min_max_migrate[1]) if min_max_migrate else max(weights)
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
            "edge_vmin": min_edge_color,
            "edge_vmax": max_edge_color,
            "edge_cmap": cmap,
            "connectionstyle": "arc3,rad=0.1",
        }
    else:
        edge_options = {
            "edge_color": migrate_color,
            "width": 2,
            "connectionstyle": "arc3,rad=0.1",
        }
    # print(f"Node sizes: {node_sizes}")
    # print(f"Edge weights: {weights}")
    nx.draw_networkx_nodes(G, pos=pos, **node_options, ax=ax)
    nx.draw_networkx_edges(G, pos=pos, **edge_options, ax=ax)
    ax.axis("off")

    # Colorbar "legend"
    norm_weights = mpl.colors.Normalize(vmin=min_edge_color, vmax=max_edge_color)
    if cax is not None:
        plt.colorbar(
            plt.cm.ScalarMappable(cmap=cmap, norm=norm_weights),
            orientation="vertical",
            cax=cax,
        )


def get_matching_colors(num_demes, demes=[]):
    if demes:
        return {
            f"{demes[i]}": plt.cm.Dark2(c)
            for i, c in enumerate(numpy.linspace(0, 1, num_demes))
        }
    return {
        f"pop_{i}": plt.cm.Dark2(c)
        for i, c in enumerate(numpy.linspace(0, 1, num_demes))
    }


def tab_show(
    filename: str,
    sort_by: str = "Index",
    show_popsize: bool = False,
    bound_column: bool = False,
):
    """
    Print an ASCII table showing the parameter values and their error from ground truth, for the
    given JSON output from the solver.
    """
    with open(filename) as f:
        output = json.load(f)

    parameter_keys = (
        "epoch_times_gen",
        "smatrix_values_ne__gen",
        "amatrix_parameters",
    )

    ploidy = output["ploidy"]

    def coal2ne(param, rate):
        return 1 / (ploidy * _clamp(param, rate))

    all_params: Iterable[Dict[str, Any]] = itertools.chain.from_iterable(
        map(lambda k: output.get(k, []) or [], parameter_keys)
    )
    results = []
    total_rel = 0.0
    total_abs = 0.0
    for param_idx, param in enumerate(all_params):
        if _param_is_synthetic(param) or _param_is_fixed(param):
            continue
        gt = param["ground_truth"]
        final = param["final"]
        description = param["description"]
        lower = param["lb"]
        upper = param["ub"]

        if show_popsize and param["kind"] == "coalescence":
            gt = coal2ne(param, gt)
            final = coal2ne(param, final)
            lower = coal2ne(param, lower)
            upper = coal2ne(param, upper)
            description = description.replace("Coalescence rate", "Ne")
        abserr = abs(gt - final)
        total_abs += abserr
        relerr = abserr / gt
        total_rel += relerr
        epochs = set()
        for app in param["apply_to"]:
            epochs.add(app["epoch"])
        results.append(
            [
                param_idx,
                description,
                relerr,
                abserr,
                gt,
                final,
                list(sorted(epochs)),
            ]
        )
        if bound_column:
            results[-1].append(_is_nearly(lower, final) or _is_nearly(upper, final))

    headers = [
        "Index",
        "Description",
        "Relative Error",
        "Absolute Error",
        "Truth",
        "Final",
        "Epochs",
    ]
    if bound_column:
        headers.append("On Bounds")

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


# FIXME This is gross. Would be better to store the number of demes (and other model info)
# in the JSON file directory.
# nstates = D*(D+1) / 2, below is the positive root
def _demes_from_states(nstates: int) -> int:
    return int(math.sqrt(8 * nstates + 1) / 2 - 1 / 2)


def _verify_timeslices(time_slice_lists: List[numpy.typing.NDArray], labels: List[str]):
    """
    Print warnings if there are properties of the time slices that may cause artifacts when comparing
    the corresponding models/data.
    """
    mints = min(map(len, time_slice_lists))
    maxts = max(map(len, time_slice_lists))
    if mints != maxts:
        print(
            f"WARNING: Models are using different numbers of time slices. Min time slices={mints}, max={maxts}",
            file=sys.stderr,
        )

    print("Median generations between time slices:", file=sys.stderr)
    for l, ts_list in enumerate(time_slice_lists):
        print(
            f"  {labels[l]}: {numpy.median(ts_list[1:] - ts_list[:-1])}",
            file=sys.stderr,
        )
    print("  A large discrepancy can affect the comparison plots.", file=sys.stderr)


def coal_dist_compare(
    result_files: List[str],
    model: str,
    labels: List[str] = [],
    max_generation: int = 200_000,
    do_cdf: bool = True,
    plot: Optional[str] = "DISPLAY",
) -> pandas.DataFrame:
    """
    Compare the pairwise coalescence distributions of the given mrpast JSON files. All inputs MUST have
    the same number of demes! This method does not check the names of demes, just the count.

    Note: the timeslices used by each result (coalescence matrix) may differ, which can make comparison
    difficult. The ideal scenario is that all timeslice boundaries are the same. If this is not practical,
    then there should at least be a large number of time slices (100+) to make all the curves smoother and
    more comparable.

    :param result_files: A list of filenames of mrpast solver inputs or solver outputs. Solver inputs and
        outputs have the same file format, and only differ in whether there is output information like
        the parameter inferred values. For this method, we only need the coalescence matrix which is present
        in both input and output.
    :type result_files: List[str]
    :param model: The filename of the model to use for extracting deme information.
    :type model: str
    :param labels: Optional list of labels corresponding to each input filename. By default, the comparison
        will just use the filename (basename), unless this list is provided. Empty list means just use the
        default (filenames).
    :type labels: List[str]
    :param max_generation: The maximum number of generations to display. Default: 200,000.
    :type max_generation: int
    :param do_cdf: Plot the CDF of the coalescences. Default: True.
    :type do_cdf: bool
    :param plot: If special value "DISPLAY", then display the plot via IPython. If None, then do not attempt
        to render a plot at all. Otherwise, the string value is a filename to write the figure to (using
        matplotlib.pyplot.savefig()). Default: "DISPLAY".
    :param plot: Optional[str]
    """
    out_matrices = []
    out_times = []
    all_epochs = []
    for fn in result_files:
        with open(fn) as f:
            data = json.load(f)
        out_matrices.append(data["coal_count_matrices"])
        maxgen = max(data["time_slices_gen"][-1] + 1, max_generation)
        out_times.append(numpy.array(data["time_slices_gen"] + [maxgen]))
        all_epochs.extend([e["ground_truth"] for e in data["epoch_times_gen"]])
    all_epochs = list(sorted(set(all_epochs)))
    nstates = len(out_matrices[0][0])

    user_model = UserModel.from_file(model)
    ndemes = user_model.num_demes
    assert _demes_from_states(nstates) == ndemes
    deme_names = user_model.pop_names

    if not labels:
        labels = list(map(os.path.basename, result_files))
    assert len(labels) == len(
        result_files
    ), "Number of labels does not match number of input files"
    _verify_timeslices(out_times, labels)

    deme_pair_map = {}
    ct = 0
    for i in range(ndemes):
        for j in range(i, ndemes):
            deme_pair_map[(i, j)] = ct
            ct += 1
    assert ct == len(out_matrices[0][0]), (ndemes, ct, len(out_matrices[0][0]))

    def norm_row(row):
        return numpy.array(row) / sum(row)

    def norm_matrix(matrix):
        sumv = 0
        for row in matrix:
            sumv += sum(row)
        return numpy.array(matrix) / sumv

    def double_norm(matrix):
        new_matrix = [norm_row(r) for r in matrix]
        return norm_matrix(new_matrix)

    def to_cdf(row):
        new_row = copy.copy(row)
        for i in range(0, len(new_row)):
            if i > 0:
                new_row[i] += new_row[i - 1]
        return new_row

    df_rows = []
    for d1 in range(ndemes):
        for d2 in range(d1, ndemes):
            state = deme_pair_map[d1, d2]

            for i, matrix_list in enumerate(out_matrices):
                for m in matrix_list:
                    m = double_norm(m)
                    row_value = m[state]
                    if do_cdf:
                        row_value = to_cdf(row_value)
                    for j in range(len(row_value)):
                        df_rows.append(
                            {
                                "ARGs": labels[i],
                                "state": state,
                                "time": out_times[i][j],
                                "coals": row_value[j],
                            }
                        )
    data = pandas.DataFrame.from_dict(df_rows)

    if plot is not None:
        assert (
            plt is not None
        ), "Plotting requires matplotlib and seaborn (pip install them)"
        plt.rc("font", **{"size": 12})
        num_cols = 3
        num_rows = 2
        fig, axs = plt.subplots(
            num_rows, num_cols, figsize=(num_cols * 6, num_rows * 5)
        )

        row = 0
        col = 0
        for d1 in range(ndemes):
            for d2 in range(d1, ndemes):
                state = deme_pair_map[d1, d2]

                data_subset = data[data["state"] == state]

                if do_cdf:
                    sns.scatterplot(
                        data=data_subset,
                        x="time",
                        y="coals",
                        hue="ARGs",
                        ax=axs[row][col],
                        lw=0,
                        alpha=0.75,
                    )
                else:
                    # The faded area is 95% confidence interval
                    sns.lineplot(
                        data=data_subset,
                        x="time",
                        y="coals",
                        hue="ARGs",
                        ax=axs[row][col],
                        errorbar=("pi", 100),
                    )
                axs[row][col].set_xlabel("Time (Generations)")
                axs[row][col].set_ylabel("Normalized coalescences")
                if d1 == d2:
                    axs[row][col].set_title(f"Within {deme_names[d1]}")
                else:
                    axs[row][col].set_title(f"Across {deme_names[d1]},{deme_names[d2]}")
                axs[row][col].set_xscale("log")
                if col > 0:
                    axs[row][col].set_ylabel(None)
                if row == 0:
                    axs[row][col].set_xlabel(None)
                if (row, col) != (0, 0):
                    axs[row][col].get_legend().remove()
                axs[row][col].set_yticks([])

                col += 1
                if col >= num_cols:
                    col = 0
                    row += 1

        fig.tight_layout()
        fig.subplots_adjust(hspace=0.25)
        if plot == "DISPLAY":
            pass
        else:
            fig.savefig(plot)
    return data
