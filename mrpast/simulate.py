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


MAX_EPOCH = 2**32


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
        rate_value = None
        for epoch in range(num_epochs):
            # ne = 1 / (lambda * ploidy)
            rate_value = user_model.coalescence.get_sim_value(epoch, i)
            if rate_value is not None:
                if epoch == 0:
                    initially_active = True
                    active_pops.append(i)
                size = Ne_from_coal_rate(rate_value)
                break
        assert (
            rate_value is not None
        ), "Every population must have at least one epoch with a coalescent rate"
        growth_rate = user_model.growth.get_sim_value(epoch, i)
        demography.add_population(
            initial_size=size,
            initially_active=initially_active,
            growth_rate=growth_rate,  # May be None
            name=user_model.pop_names[i],
        )

    # Pass 2: Find all population splits and admixture events.
    def proportionAsFloat(proportion: Union[float, ParamRef]) -> float:
        if isinstance(proportion, ParamRef):
            return user_model.admixture.get_parameter(proportion.param).ground_truth
        return float(proportion)

    dead_pops = {}
    for epoch in range(1, user_model.num_epochs):
        # Collect all entries by their derived population.
        by_derived = defaultdict(list)
        for i, entry in enumerate(user_model.admixture.entries):
            if entry.epoch == epoch:
                by_derived[entry.derived].append(
                    (entry.ancestral, proportionAsFloat(entry.proportion))
                )
        epoch_start = user_model.epochs.epoch_times[epoch - 1].ground_truth
        for derived_deme in range(user_model.num_demes):
            # Population split if we have a 1-to-1 mapping.
            if len(by_derived[derived_deme]) == 1:
                ancestral, proportion = by_derived[derived_deme][0]
                assert abs(1.0 - proportion) < 1e6
                demography.add_population_split(
                    time=epoch_start, derived=[derived_deme], ancestral=ancestral
                )
                dead_pops[derived_deme] = epoch
            # Otherwise it is admixture.
            elif len(by_derived[derived_deme]) > 1:
                demography.add_admixture(
                    epoch_start,
                    derived=derived_deme,
                    ancestral=list(map(lambda t: t[0], by_derived[derived_deme])),
                    proportions=list(map(lambda t: t[1], by_derived[derived_deme])),
                )
                dead_pops[derived_deme] = epoch

    # Pass 3: setup the initial migration rates and rate change events.
    for epoch in range(num_epochs):
        epoch_time = None
        if epoch > 0:
            epoch_time = user_model.epochs.epoch_times[epoch - 1].ground_truth
        for i in range(num_pops):
            # We skip any dead populations. The simulator will yell at us for trying to make changes to them.
            if dead_pops.get(i, MAX_EPOCH) <= epoch:
                continue

            # Handle coalescence rate (effective population size) changes
            if epoch > 0:
                coal_rate = user_model.coalescence.get_sim_value(epoch, i)
                prev_coal_rate = user_model.coalescence.get_sim_value(epoch - 1, i)
                grow_rate = user_model.growth.get_sim_value(epoch, i)
                prev_grow_rate = user_model.growth.get_sim_value(epoch - 1, i)
                if (
                    prev_coal_rate != coal_rate
                    and prev_coal_rate is not None
                    and coal_rate is not None
                ) or (grow_rate != prev_grow_rate):
                    if coal_rate is not None:
                        size = Ne_from_coal_rate(coal_rate)
                    else:
                        size = 0

                    demography.add_population_parameters_change(
                        epoch_time,
                        initial_size=size,
                        population=i,
                        growth_rate=grow_rate if grow_rate is not None else 0,
                    )

            # Handle migration with all other populations.
            for j in range(user_model.pop_count):
                mig_rate = user_model.migration.get_sim_value(epoch, i, j)
                if epoch > 0:
                    prev_mig_rate = user_model.migration.get_sim_value(epoch - 1, i, j)

                # "The entry of [migration rate matrix] is the expected number of migrants moving from population i
                #    to population j per generation, divided by the size of population j."
                if epoch == 0:
                    if mig_rate is not None:
                        demography.set_migration_rate(i, j, mig_rate)
                else:
                    # We have a change in rate.
                    if mig_rate != prev_mig_rate:
                        if mig_rate is None:
                            demography.add_migration_rate_change(
                                time=epoch_time, rate=0, source=i, dest=j
                            )
                        else:
                            demography.add_migration_rate_change(
                                time=epoch_time,
                                rate=mig_rate,
                                source=i,
                                dest=j,
                            )

    # Add any mass migration (pulse) events
    for entry in user_model.pulse.entries:
        demography.add_mass_migration(
            time=entry.get_sim_time(user_model.pulse),
            source=entry.source,
            dest=entry.dest,
            proportion=entry.get_sim_proportion(user_model.pulse),
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
