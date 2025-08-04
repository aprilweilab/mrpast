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
from collections import defaultdict
from tabulate import tabulate
from tqdm import tqdm
from multiprocessing import Pool
from typing import Tuple, List, Union
import msprime
import os

from mrpast.helpers import (
    load_ratemap,
)
from mrpast.model import UserModel, ParamRef


def build_demography(user_model: UserModel) -> Tuple[msprime.Demography, List[int]]:
    """
    Create the msprime Demography object for the given mrpast model filename.
    """

    # ASSUMPTION: Epoch_0 (the most recent) lists all the populations. Any population can be "unused"
    # in any epoch. If Epoch_i "introduces" a new population (backwards in time) then that population
    # would just have no coalescence+migration in Epoch_i-1.

    # For a given epoch, if there is no coalescence parameter then that population must be INACTIVE
    # in that epoch, which means:
    # 1. It is either in the msprime INACTIVE or PREVIOUSLY_ACTIVE state, i.e. it is involved in an
    #    ancestral population split.
    # 2. There can be no migration to/from it during that epoch.
    #
    # For example, if you have an ancestral population you can start out with it being active or inactive:
    # * Giving the population a coalescence parameter will ensure it starts out active.
    # * Otherwise it will start inactive, and become active during the first epoch where it has a
    #   coalescence parameter.
    # * You cannot give the population a coalescence parameter in an Epoch other than 0, unless the
    #   population has become active by that Epoch due to a migration event.

    def Ne_from_coal_rate(coal_rate: float) -> float:
        return 1 / (user_model.ploidy * coal_rate)

    # Migration is thought of backwards here, and in the MrPast model. So if we
    # have non-zero migration from A->B, it means that in forward-time people migrated
    # from B->A
    num_epochs = user_model.num_epochs
    num_pops = user_model.pop_count
    demography = msprime.Demography()
    active_pops = []
    # Pass 1: Add all of the populations and determine their initial size.
    for i in range(num_pops):
        initially_active = False
        size = None
        for epoch in range(num_epochs):
            # ne = 1 / (lambda * ploidy)
            rate_entry = user_model.coalescence.get_entry(epoch, i)
            if rate_entry is not None:
                assert isinstance(rate_entry.rate, ParamRef)
                rate_param = user_model.coalescence.get_parameter(rate_entry.rate.param)
                if epoch == 0:
                    initially_active = True
                    active_pops.append(i)
                size = Ne_from_coal_rate(rate_param.ground_truth)
                break
        assert (
            rate_param is not None
        ), "Every population must have at least one epoch with a coalescent rate"
        growth_rate = None
        grow_entry = user_model.growth.get_entry(epoch, i)
        if grow_entry is not None:
            assert isinstance(grow_entry.rate, ParamRef)
            grate_param = user_model.growth.get_parameter(grow_entry.rate.param)
            growth_rate = grate_param.ground_truth
        demography.add_population(
            initial_size=size,
            initially_active=initially_active,
            growth_rate=growth_rate,
            name=user_model.pop_names[i],
        )
    del rate_param

    # Pass 2: Find all population splits by identifying changes in populationConversion
    dead_pops = {}
    for dest_epoch in range(1, num_epochs):
        source_epoch = dest_epoch - 1
        derived = defaultdict(list)
        for src_pop in range(num_pops):
            dest_pop = user_model.pop_convert[source_epoch][src_pop]
            if dest_pop != src_pop and src_pop not in dead_pops:
                derived[dest_pop].append(src_pop)
                dead_pops[src_pop] = dest_epoch
        epoch_time = user_model.epochs.epoch_times[source_epoch].ground_truth
        for dest_pop, derived_list in derived.items():
            demography.add_population_split(
                time=epoch_time, derived=list(derived_list), ancestral=dest_pop
            )

    # Pass 3: setup the initial migration rates and rate change events.
    for epoch in range(num_epochs):
        epoch_time = None
        if epoch > 0:
            epoch_time = user_model.epochs.epoch_times[epoch - 1].ground_truth
        for i in range(num_pops):
            # We skip any dead populations. The simulator will yell at us for trying to make changes to them.
            if dead_pops.get(i, 2**32) <= epoch:
                continue

            # Handle coalescence rate (effective population size) changes
            if epoch > 0:
                crate_entry = user_model.coalescence.get_entry(epoch, i)
                prev_crate_entry = user_model.coalescence.get_entry(epoch - 1, i)
                grate_entry = user_model.growth.get_entry(epoch, i)
                prev_grate_entry = user_model.growth.get_entry(epoch - 1, i)
                if (
                    prev_crate_entry != crate_entry
                    and prev_crate_entry is not None
                    and crate_entry is not None
                ) or (grate_entry != prev_grate_entry):
                    assert isinstance(crate_entry.rate, ParamRef)
                    crate_param = user_model.coalescence.get_parameter(
                        crate_entry.rate.param
                    )
                    if crate_param is not None:
                        size = Ne_from_coal_rate(crate_param.ground_truth)
                    else:
                        size = 0

                    growth_rate = 0
                    if grate_entry is not None:
                        assert isinstance(grate_entry.rate, ParamRef)
                        grate_param = user_model.growth.get_parameter(
                            grate_entry.rate.param
                        )
                        growth_rate = grate_param.ground_truth
                        print(f"Changing growth rate to {growth_rate}")

                    demography.add_population_parameters_change(
                        epoch_time,
                        initial_size=size,
                        population=i,
                        growth_rate=growth_rate,
                    )

            # Handle migration with all other populations.
            for j in range(user_model.pop_count):
                mrate_entry = user_model.migration.get_entry(epoch, i, j)
                prev_mrate_entry = None
                if epoch > 0:
                    prev_mrate_entry = user_model.migration.get_entry(epoch - 1, i, j)

                # "The entry of [migration rate matrix] is the expected number of migrants moving from population i
                #    to population j per generation, divided by the size of population j."
                if epoch == 0:
                    if mrate_entry is not None:
                        assert isinstance(mrate_entry.rate, ParamRef)
                        mrate_param = user_model.migration.get_parameter(
                            mrate_entry.rate.param
                        )
                        demography.set_migration_rate(i, j, mrate_param.ground_truth)
                else:
                    # We have a change in rate.
                    if mrate_entry != prev_mrate_entry:
                        if mrate_entry is None:
                            demography.add_migration_rate_change(
                                time=epoch_time, rate=0, source=i, dest=j
                            )
                        else:
                            assert isinstance(mrate_entry.rate, ParamRef)
                            mrate_param = user_model.migration.get_parameter(
                                mrate_entry.rate.param
                            )
                            demography.add_migration_rate_change(
                                time=epoch_time,
                                rate=mrate_param.ground_truth,
                                source=i,
                                dest=j,
                            )
    demography.sort_events()
    return demography, active_pops


def _run_simulation(
    model_file: str,
    arg_prefix: str,
    seq_len: int,
    num_replicates: int,
    ident: int,
    recomb_rate: Union[float, msprime.RateMap] = 1e-8,
    samples_per_pop: int = 10,
    debug_demo: bool = True,
    seed: int = 42,
) -> int:

    model = UserModel.from_file(model_file)
    demography, active_pops = build_demography(model)

    table = [
        ["Sequence Length", seq_len],
        ["Recombination rate", recomb_rate],
        ["Samples/population", samples_per_pop],
        ["Ploidy", model.ploidy],
        ["Epochs", model.num_epochs],
        ["Populations", model.pop_count],
    ]
    print("Preparing simulation with parameters:")
    print(tabulate(table, headers=["Parameter", "Value"]))
    print()

    if debug_demo:
        print(demography.debug())

    replicates = msprime.sim_ancestry(
        samples={i: samples_per_pop for i in active_pops},
        demography=demography,
        recombination_rate=recomb_rate if recomb_rate != 0 else None,
        sequence_length=seq_len,
        random_seed=seed,
        num_replicates=num_replicates,
    )

    total_trees = 0
    for i, tree_sequence in tqdm(enumerate(replicates)):
        total_trees += tree_sequence.num_trees
        tree_sequence.dump(f"{arg_prefix}_{ident}-{i}.trees")
    return total_trees


def run_simulation(
    model: str,
    arg_prefix: str,
    seq_len: int,
    num_replicates: int,
    recomb_rate: Union[float, str] = 1e-8,
    samples_per_pop: int = 10,
    debug_demo: bool = True,
    jobs: int = 1,
    seed: int = 42,
) -> int:
    # If user provided a filename, load it as recombination map
    if isinstance(recomb_rate, str):
        assert os.path.isfile(
            recomb_rate
        ), "Invalid recombination rate: provide a filename or constant rate"
        recomb_rate_or_map = load_ratemap(recomb_rate)
        # msprime is really picky about the ratemap length exactly matching the
        # length of the simulated sequence.
        recomb_rate_or_map = recomb_rate_or_map.slice(left=0, right=seq_len, trim=True)
    else:
        recomb_rate_or_map = float(recomb_rate)

    if num_replicates % jobs == 0:
        reps_per_job = num_replicates // jobs
        work = [
            (
                model,
                arg_prefix,
                seq_len,
                reps_per_job,
                ident,
                recomb_rate_or_map,
                samples_per_pop,
                debug_demo,
                seed + ident,
            )
            for ident in range(jobs)
        ]
    else:
        work = [
            (
                model,
                arg_prefix,
                seq_len,
                1,
                ident,
                recomb_rate_or_map,
                samples_per_pop,
                debug_demo,
                seed + ident,
            )
            for ident in range(num_replicates)
        ]
    if jobs == 1:
        tree_counts = [_run_simulation(*work[0])]
    else:
        with Pool(jobs) as p:
            tree_counts = p.starmap(_run_simulation, work)
    return sum(tree_counts)
