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
from mrpast.model import (
    UserModel,
    DEFAULT_PLOIDY,
    DemeDemeRates,
    DemeRates,
    DemeDemeEntry,
    DemeRateEntry,
    SymbolicEpochs,
    FloatParameter,
    ParamRef,
)

try:
    import demes
except ImportError:
    demes = None  # type: ignore

TIME_UNIT_GENS = "generations"


def convert_from_demes(demes_file: str) -> Dict[str, Any]:
    assert demes is not None, f"'demes' package not found; try 'pip install demes'"

    output_model = UserModel(
        ploidy=DEFAULT_PLOIDY,
        pop_count=0,
        pop_names=[],
        migration=DemeDemeRates(entries=[], parameters=[]),
        coalescence=DemeRates(entries=[], parameters=[]),
        growth=DemeRates(entries=[], parameters=[]),
        epochs=SymbolicEpochs([]),
        pop_convert=[],
    )
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
        output_model.pop_names.append(d.name)
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
        output_model.pop_convert.append(list(range(len(demes_model.demes))))

    def add_param(param_list, gt_value, lb=1e-5, ub=0.01):
        index = len(param_list) + 1
        param_list.append(
            FloatParameter(ground_truth=gt_value, lb=lb, ub=ub, index=index)
        )
        return index

    # When the same value is used across multiple epochs for the "same" parameter, we remember that and don't
    # create another model parameter. The one kind of annoying case that we still don't cover is that of symmetric
    # migration: we will create a parameter for both directions even if the rate value is identicial.
    remembered_params: Dict[Any, Any] = {}

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
                    gkey = ("growth", growth_rate, d_idx)
                    if gkey not in remembered_params:
                        param_idx = add_param(
                            output_model.growth.parameters, growth_rate
                        )
                        remembered_params[gkey] = param_idx
                    param_idx = remembered_params[gkey]
                    output_model.growth.entries.append(
                        DemeRateEntry(
                            epoch=e_idx, deme=d_idx, rate=ParamRef(param=param_idx)
                        )
                    )

                coal_rate = 1 / (2 * e.end_size)
                ckey = ("coalescence", coal_rate, d_idx)
                if ckey not in remembered_params:
                    param_idx = add_param(
                        output_model.coalescence.parameters, coal_rate
                    )
                    remembered_params[ckey] = param_idx
                param_idx = remembered_params[ckey]
                output_model.coalescence.entries.append(
                    DemeRateEntry(
                        epoch=e_idx, deme=d_idx, rate=ParamRef(param=param_idx)
                    )
                )

                last_epoch = max(last_epoch, e_idx)
                e_idx = (e_idx + 1) if (e_idx + 1) < len(epochs_by_start) else None
        if d.ancestors:
            assert (
                len(d.ancestors) == 1
            ), "mrpast does not support admixture; only 1-to-1 population split"
            # Set the population splits.
            for a in d.ancestors:
                output_model.pop_convert[last_epoch][d_idx] = name2deme[a]
        # Mark the population as dead in the relevant epochs.
        for e_idx in range(last_epoch + 1, num_epochs - 1):
            output_model.pop_convert[e_idx][d_idx] = float("NaN")

    prev_row = list(range(num_demes))
    for row in output_model.pop_convert:
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
            dest_idx = name2deme[m.dest]
            src_idx = name2deme[m.source]
            mkey = ("migration", m.rate, dest_idx, src_idx)
            if mkey not in remembered_params:
                param_idx = add_param(output_model.migration.parameters, m.rate)
                remembered_params[mkey] = param_idx
            param_idx = remembered_params[mkey]
            # Our model reverses destination and source (backwards vs. forward in time)
            output_model.migration.entries.append(
                DemeDemeEntry(
                    epoch=e_idx,
                    source=dest_idx,
                    dest=src_idx,
                    rate=ParamRef(param=param_idx),
                )
            )
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
        add_param(output_model.epochs.epoch_times, timeval, lb=next_lb, ub=next_ub)
        next_lb = next_ub

    output_model.pop_count = len(output_model.pop_names)
    return output_model
