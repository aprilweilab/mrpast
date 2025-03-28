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
from mrpast.model import (
    load_model_config,
    SymbolicMatrices,
    SymbolicEpochs,
)


def build_demography(model: str) -> Tuple[msprime.Demography, List[int]]:
    """
    Create the msprime Demography object for the given mrpast model filename.
    """
    config = load_model_config(model)
    ploidy = config["ploidy"]
    M_symbolic = SymbolicMatrices.from_config(config["migration"])
    q_symbolic = SymbolicMatrices.from_config(config["coalescence"], is_vector=True)
    if "growth" in config:
        g_symbolic = SymbolicMatrices.from_config(config["growth"], is_vector=True)
    else:
        g_symbolic = None
    epochs_symbolic = SymbolicEpochs.from_config(config.get("epochTimeSplit", []))
    assert M_symbolic.num_epochs >= 1
    popConvert = config.get("populationConversion", [])
    if popConvert is None:
        popConvert = []
    assert len(popConvert) == M_symbolic.num_epochs - 1

    npops = [M_symbolic.num_demes(i) for i in range(M_symbolic.num_epochs)]
    assert all(
        map(lambda n: n == npops[0], npops)
    ), "All matrices must have the same population size"
    for row in popConvert:
        assert len(row) == npops[0]

    pop_names = config.get("pop_names", [f"pop_{i}" for i in range(npops[0])])
    assert len(pop_names) == npops[0], f"pop_names has wrong length, must be {npops[0]}"

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
        return 1 / (ploidy * coal_rate)

    # Migration is thought of backwards here, and in the MrPast model. So if we
    # have non-zero migration from A->B, it means that in forward-time people migrated
    # from B->A
    num_epochs = len(npops)
    num_pops = npops[0]
    demography = msprime.Demography()
    active_pops = []
    # Pass 1: Add all of the populations and determine their initial size.
    for i in range(num_pops):
        initially_active = False
        size = None
        for epoch in range(num_epochs):
            # ne = 1 / (lambda * ploidy)
            rate_param = q_symbolic.get_parameter(q_symbolic.matrices[epoch][i])
            if rate_param is not None:
                if epoch == 0:
                    initially_active = True
                    active_pops.append(i)
                size = Ne_from_coal_rate(rate_param.ground_truth)
                break
        assert (
            rate_param is not None
        ), "Every population must have at least one epoch with a coalescent rate"
        growth_rate = None
        if g_symbolic is not None:
            grate_param = g_symbolic.get_parameter(g_symbolic.matrices[epoch][i])
            if grate_param is not None:
                growth_rate = grate_param.ground_truth
        demography.add_population(
            initial_size=size,
            initially_active=initially_active,
            growth_rate=growth_rate,
            name=pop_names[i],
        )
    # Pass 2: Find all population splits by identifying changes in populationConversion
    dead_pops = {}
    for dest_epoch in range(1, num_epochs):
        source_epoch = dest_epoch - 1
        derived = defaultdict(list)
        for src_pop in range(num_pops):
            dest_pop = popConvert[source_epoch][src_pop]
            if dest_pop != src_pop and src_pop not in dead_pops:
                derived[dest_pop].append(src_pop)
                dead_pops[src_pop] = dest_epoch
        epoch_time = epochs_symbolic.epoch_times[source_epoch].ground_truth
        for dest_pop, derived_list in derived.items():
            demography.add_population_split(
                time=epoch_time, derived=list(derived_list), ancestral=dest_pop
            )
    del rate_param

    # Pass 3: setup the initial migration rates and rate change events.
    for epoch in range(num_epochs):
        epoch_time = None
        if epoch > 0:
            epoch_time = epochs_symbolic.epoch_times[epoch - 1].ground_truth
        for i in range(num_pops):
            # We skip any dead populations. The simulator will yell at us for trying to make changes to them.
            if dead_pops.get(i, 2**32) <= epoch:
                continue

            # Handle coalescence rate (effective population size) changes
            if epoch > 0:
                crate_idx = q_symbolic.matrices[epoch][i]
                prev_crate_idx = q_symbolic.matrices[epoch - 1][i]
                prev_grate_param_idx = None
                grate_param_idx = None
                if g_symbolic is not None:
                    grate_param_idx = g_symbolic.matrices[epoch][i]
                    prev_grate_param_idx = g_symbolic.matrices[epoch - 1][i]
                if (
                    prev_crate_idx != crate_idx
                    and prev_crate_idx != 0
                    and crate_idx != 0
                ) or (grate_param_idx != prev_grate_param_idx):
                    crate_param = q_symbolic.get_parameter(crate_idx)
                    if crate_param is not None:
                        size = Ne_from_coal_rate(crate_param.ground_truth)
                    else:
                        size = 0

                    growth_rate = 0
                    if g_symbolic is not None:
                        grate_param = g_symbolic.get_parameter(grate_param_idx)
                        if grate_param is not None:
                            growth_rate = grate_param.ground_truth
                            print(f"Changing growth rate to {growth_rate}")

                    demography.add_population_parameters_change(
                        epoch_time,
                        initial_size=size,
                        population=i,
                        growth_rate=growth_rate,
                    )

            # Handle migration with all other populations.
            for j in range(num_pops):
                mrate_idx = M_symbolic.matrices[epoch][i, j]
                mrate_param = M_symbolic.get_parameter(mrate_idx)
                prev_mrate_idx = None
                if epoch > 0:
                    prev_mrate_idx = M_symbolic.matrices[epoch - 1][i, j]

                # "The entry of [migration rate matrix] is the expected number of migrants moving from population i
                #    to population j per generation, divided by the size of population j."
                if epoch == 0:
                    if mrate_param is not None:
                        demography.set_migration_rate(i, j, mrate_param.ground_truth)
                else:
                    assert prev_mrate_idx is not None
                    # We have a change in rate.
                    if mrate_idx != prev_mrate_idx:
                        if mrate_param is None:
                            demography.add_migration_rate_change(
                                time=epoch_time, rate=0, source=i, dest=j
                            )
                        else:
                            demography.add_migration_rate_change(
                                time=epoch_time,
                                rate=mrate_param.ground_truth,
                                source=i,
                                dest=j,
                            )
    demography.sort_events()
    return demography, active_pops


def _run_simulation(
    model: str,
    arg_prefix: str,
    seq_len: int,
    num_replicates: int,
    ident: int,
    recomb_rate: Union[float, msprime.RateMap] = 1e-8,
    samples_per_pop: int = 10,
    debug_demo: bool = True,
    seed: int = 42,
) -> int:

    demography, active_pops = build_demography(model)

    config = load_model_config(model)
    ploidy = config["ploidy"]
    M_symbolic = SymbolicMatrices.from_config(config["migration"])
    npops = [M_symbolic.num_demes(i) for i in range(M_symbolic.num_epochs)]
    assert all(
        map(lambda n: n == npops[0], npops)
    ), "All matrices must have the same population size"

    table = [
        ["Sequence Length", seq_len],
        ["Recombination rate", recomb_rate],
        ["Samples/population", samples_per_pop],
        ["Ploidy", ploidy],
        ["Epochs", M_symbolic.num_epochs],
    ]
    for i in range(len(npops)):
        table.append([f"Population Count (E{i})", npops[i]])
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
