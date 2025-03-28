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
import math
import copy
from typing import Dict, Any, Optional

try:
    import demes
except ImportError:
    demes = None  # type: ignore

TIME_UNIT_GENS = "generations"


def convert_from_demes(demes_file: str) -> Dict[str, Any]:
    assert demes is not None, f"'demes' package not found; try 'pip install demes'"

    output_model: Dict[str, Any] = {
        "ploidy": 2,
        "pop_names": [],
        "coalescence": {
            "vectors": [],
            "parameters": [],
        },
        "migration": {
            "matrices": [],
            "parameters": [],
        },
        "growth": {
            "vectors": [],
            "parameters": [],
        },
        "epochTimeSplit": [],
        "populationConversion": [],
    }
    name2deme = {}

    def make_matrix(N: int):
        return [copy.copy([0 for _ in range(N)]) for _ in range(N)]

    def make_vector(N: int):
        return copy.copy([0 for _ in range(N)])

    with open(demes_file) as f:
        demes_model = demes.load(f)
    assert (
        not demes_model.pulses
    ), "Admixture pulses are not supported in MrPast models yet."
    assert (
        demes_model.time_units == TIME_UNIT_GENS
    ), "Only time_units supported is generations"
    num_demes = len(demes_model.demes)

    # List of epochs, ordered by most recent time first.
    epoch_set = set()
    for i, d in enumerate(demes_model.demes):
        output_model["pop_names"].append(d.name)
        name2deme[d.name] = i
        for e in d.epochs:
            epoch_set.add(e.start_time)
            epoch_set.add(e.end_time)
    epoch_delimiters = list(sorted(epoch_set))
    assert epoch_delimiters[0] == 0
    assert math.isinf(epoch_delimiters[-1])
    del epoch_delimiters[-1]
    num_epochs = len(epoch_delimiters)
    epochs_by_start = {e: i for i, e in enumerate(epoch_delimiters)}
    for _ in range(num_epochs - 1):
        output_model["populationConversion"].append(list(range(len(demes_model.demes))))
    coal_vects = output_model["coalescence"]["vectors"]
    mig_mats = output_model["migration"]["matrices"]
    growth_vects = output_model["growth"]["vectors"]

    def add_param(area, gt_value, sub_field=True, lb=1e-5, ub=0.01):
        if sub_field:
            dest = output_model[area]["parameters"]
        else:
            dest = output_model[area]
        index = len(dest) + 1
        dest.append({"ground_truth": gt_value, "lb": lb, "ub": ub, "index": index})
        return index

    # When the same value is used across multiple epochs for the "same" parameter, we remember that and don't
    # create another model parameter. The one kind of annoying case that we still don't cover is that of symmetric
    # migration: we will create a parameter for both directions even if the rate value is identicial.
    remembered_params: Dict[Any, Any] = {}

    seen_growth = False
    for _ in range(num_epochs):
        coal_vects.append(make_vector(len(demes_model.demes)))
        mig_mats.append(make_matrix(len(demes_model.demes)))
        growth_vects.append(make_vector(len(demes_model.demes)))
    for d_idx, d in enumerate(demes_model.demes):
        last_epoch = 0
        for e in d.epochs:
            e_end = epochs_by_start.get(e.start_time, None)
            e_idx: Optional[int] = epochs_by_start[e.end_time]
            # The next epoch and the end of the demes' epoch are not the same, so this demes' epoch spans
            # multiple epochs in the output model.
            while e_idx != e_end:
                assert e_idx is not None
                growth_rate = -math.log(e.start_size / e.end_size) / e.time_span
                if growth_rate > 0.0:
                    seen_growth = True
                    if ("growth", growth_rate, d_idx) in remembered_params:
                        growth_vects[e_idx][d_idx] = remembered_params[
                            ("growth", growth_rate, d_idx)
                        ]
                    else:
                        param_idx = add_param("growth", growth_rate)
                        remembered_params[("growth", growth_rate, d_idx)] = param_idx
                        growth_vects[e_idx][d_idx] = param_idx
                coal_rate = 1 / (2 * e.end_size)
                if ("coalescence", coal_rate, d_idx) in remembered_params:
                    coal_vects[e_idx][d_idx] = remembered_params[
                        ("coalescence", coal_rate, d_idx)
                    ]
                else:
                    param_idx = add_param("coalescence", coal_rate)
                    remembered_params[("coalescence", coal_rate, d_idx)] = param_idx
                    coal_vects[e_idx][d_idx] = param_idx
                last_epoch = max(last_epoch, e_idx)
                e_idx = (e_idx + 1) if (e_idx + 1) < len(epochs_by_start) else None
        # Set the population splits.
        for a in d.ancestors:
            output_model["populationConversion"][last_epoch][d_idx] = name2deme[a]
        # Mark the population as dead in the relevant epochs.
        for e_idx in range(last_epoch + 1, num_epochs - 1):
            output_model["populationConversion"][e_idx][d_idx] = float("NaN")

    prev_row = list(range(num_demes))
    for row in output_model["populationConversion"]:
        for i in range(len(row)):
            if math.isnan(row[i]):
                last_dest = prev_row[i]
                current_dest = row[last_dest]
                assert not math.isnan(current_dest)
                row[i] = current_dest
        prev_row = row

    for m in demes_model.migrations:
        e_end = epochs_by_start.get(m.start_time, None)
        e_idx = epochs_by_start[m.end_time]
        while e_idx != e_end:
            assert e_idx is not None
            # Our model reverses destination and source (backwards vs. forward in time)
            dest_idx = name2deme[m.dest]
            src_idx = name2deme[m.source]
            if ("migration", m.rate, dest_idx, src_idx) in remembered_params:
                mig_mats[e_idx][dest_idx][src_idx] = remembered_params[
                    ("migration", m.rate, dest_idx, src_idx)
                ]
            else:
                param_idx = add_param("migration", m.rate)
                remembered_params[("migration", m.rate, dest_idx, src_idx)] = param_idx
                mig_mats[e_idx][dest_idx][src_idx] = param_idx
            e_idx = (e_idx + 1) if (e_idx + 1) < len(epochs_by_start) else None

    next_lb = 100.0
    epoch_times = list(epochs_by_start.keys())
    for i in range(len(epoch_times)):
        if i == 0:
            continue
        timeval = epoch_times[i]
        next_val = (
            epoch_times[i + 1] if (i + 1) < len(epoch_times) else timeval + next_lb
        )
        next_ub = (timeval + next_val) // 2
        add_param("epochTimeSplit", timeval, sub_field=False, lb=next_lb, ub=next_ub)
        next_lb = next_ub

    if not seen_growth:
        del output_model["growth"]

    return output_model
