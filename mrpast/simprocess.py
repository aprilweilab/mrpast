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
from copy import deepcopy
from functools import lru_cache
from multiprocessing import Pool
from numpy.typing import NDArray
from tqdm import tqdm
from typing import Dict, List, Tuple, Optional, BinaryIO, Set, Iterable
import hashlib
import itertools
import msprime
import numpy as np
import orjson
import tskit

from mrpast.helpers import count_lines
from mrpast.model import PopMap

DemePairIndex = Dict[Tuple[int, int], int]
CoalescenceList = List[Tuple[int, int, float, float]]


# Load a tskit/SINGER-style rate map. The values in the tuple are (start position,
# end position, rate for that region)
def load_rate_map(rate_map_file: str) -> List[Tuple[float, float, float]]:
    result = []
    with open(rate_map_file) as f:
        for line in f:
            parts = list(map(str.strip, line.split()))
            result.append((float(parts[0]), float(parts[1]), float(parts[2])))
    return result


def update_coalescence_map(
    tree_seq_file: str,
    output_file: BinaryIO,
    min_time_precision: float = 1.0,
    sample_every_bp: Optional[int] = 125_000,
    leave_out_pops: Set[int] = set(),
    recomb_map_file: Optional[str] = None,
    recomb_rate_threshold: Optional[float] = 1e-09,
):
    """
    Given a tree-seq, get the list of all (popI, popJ, coalescence time, coalescence weight)
    """
    rate_map = None
    if recomb_map_file is not None:
        rate_map = load_rate_map(recomb_map_file)

    def sample_pop(sample_id: int, tree: tskit.Tree) -> int:
        return tree.population(sample_id)

    # Recursively collect the vector [p_0, p_1, ..., p_k] where there are k populations in the
    # ARG, and p_i is the count of samples below the current node that are in population i.
    # Invoke the given callback for every pair (i, j) of population IDs, with the coalescence
    # time and number of samples that coalesced between those two populations.
    # For any given node, all the samples to the left already coalesced with each other, same
    # with the right. So we're adding the coalescence BETWEEN all the ones on the left and all
    # the ones on the right. When computing coalescence this way we will only get one order
    # of sample nodes (a, b). We have two nodes "l" and "r" that are the children of "n". For
    # a given population "i", PC_i(l) and PC_i(r) are the number of samples for population "i"
    # beneath "l" and "r" respectively. Note that the sets of samples in population "i" beneath
    # "l" and "r" are DISJOINT, so we have to count them both. This is why we need to count
    # (i, j) as well as (j, i) for all populations i*j.
    def pops_below(tree: tskit.Tree, node: int, num_pop: int, callback) -> List[int]:
        if tree.is_sample(node):
            result = [0] * num_pop
            pop = sample_pop(node, tree)
            if pop not in leave_out_pops:
                result[pop] = 1
            return result
        ntime = int(tree.time(node) / min_time_precision) * min_time_precision
        child_data = []
        for c in tree.children(node):
            child_data.append(pops_below(tree, c, num_pop, callback))
        while len(child_data) > 1:
            left, right = child_data[-2:]
            for i, j in itertools.product(range(num_pop), repeat=2):
                num_coalescences = left[i] * right[j]
                if num_coalescences == 0:
                    continue
                if i <= j:
                    callback(i, j, ntime, num_coalescences)
                else:
                    callback(j, i, ntime, num_coalescences)
            new_result = [a + b for a, b in zip(left, right)]
            child_data = child_data[:-2] + [new_result]
        return child_data[0]

    tree_sequence = tskit.load(tree_seq_file)
    num_populations = tree_sequence.num_populations

    orj_flags = (
        orjson.OPT_NAIVE_UTC | orjson.OPT_SERIALIZE_NUMPY | orjson.OPT_APPEND_NEWLINE
    )
    last_map_bp = -1.0
    rate_idx = 0
    # Any region outside of provided rate map is assumed to have this rate.
    current_rate = 1e-8
    next_sample_bp = float(sample_every_bp) if sample_every_bp is not None else 0.0
    for tree in tree_sequence.trees():
        tree_start = tree.interval[0]
        maybe_sample = rate_map is None
        if rate_map is not None:
            assert recomb_rate_threshold is not None
            while tree_start > last_map_bp and rate_idx < len(rate_map):
                current_bp_start = rate_map[rate_idx][0]
                last_map_bp = rate_map[rate_idx][1]
                current_rate = rate_map[rate_idx][2]
                rate_idx += 1
            maybe_sample = current_rate <= recomb_rate_threshold
            # We want to keep moving the goal post, because we have these criteria:
            # 1. The locations that we sample have to be deterministic given a RateMap and a
            #    specific chromosome.
            # 2. We must have at _least_ sample_every_bp base-pairs between each tree.
            # 3. Related to #1, if we have different ARG samples then the trees themselves can
            #    span different regions, so we cannot select a region based on how the trees are
            #    laid out.
            next_sample_bp = max(next_sample_bp, current_bp_start)
        if not maybe_sample:
            continue
        if sample_every_bp is None:
            tree_length = tree.interval[1] - tree.interval[0]
            weight = tree_length / tree_sequence.sequence_length
        else:
            if tree_start >= next_sample_bp:
                next_sample_bp += sample_every_bp
            else:
                continue
            weight = 1.0
        result_map: Dict[Tuple[int, int, float], float] = defaultdict(float)

        def add_result(i, j, node_time, num_coals):
            result_map[(i, j, float(node_time))] += weight * num_coals

        for r in tree.roots:
            pops_below(tree, r, num_populations, add_result)
        # Each line in the output represents a single tree's worth of coalescence information
        data_str = orjson.dumps(
            [(p1, p2, t, w) for (p1, p2, t), w in result_map.items()],
            option=orj_flags,
        )
        output_file.write(data_str)


def get_coalescence_list(
    tree_seq_file: str,
    output_file_name: str,
    tree_sample_rate: int,
    min_time_precision=1.0,
    leave_out_pops=[],
    rate_map: Optional[str] = None,
    rate_map_threshold: float = 1.0,
) -> str:
    with open(output_file_name, "wb") as fout:
        update_coalescence_map(
            tree_seq_file,
            fout,
            min_time_precision=min_time_precision,
            sample_every_bp=tree_sample_rate,
            leave_out_pops=set(leave_out_pops),
            recomb_map_file=rate_map,
            recomb_rate_threshold=rate_map_threshold,
        )
    return output_file_name


def get_time_slices(
    coal_filenames: List[str],
    ntimes: int,
    max_generation: float,
    verbose: bool = True,
    left_skewed: bool = True,
):
    binned_counts: Dict[float, float] = defaultdict(float)
    total_coal = 0
    for fn in coal_filenames:
        with open(fn) as f:
            for line in f:
                coal_list = orjson.loads(line)
                for _, _, t, w in coal_list:
                    if t < max_generation:
                        binned_counts[t] += w
                        total_coal += w

    if left_skewed:
        coal_per_time = [total_coal / 2]
        for _ in range(ntimes - 1):
            coal_per_time.append(coal_per_time[-1] / 2)
    else:
        uniform_part = total_coal / ntimes
        coal_per_time = [uniform_part for _ in range(ntimes)]

    # Sorted by time ascending already.
    input_list = sorted([(t, c) for t, c in binned_counts.items() if c > 0])

    # We want an equal number of coalescence events in each of our time slices.
    time_slices = []
    last_coal_time = -1.0
    current_coal = 0.0
    next_mark = coal_per_time.pop()
    for coal_time, coal_count in input_list:
        current_coal += coal_count
        if current_coal >= next_mark and coal_time != last_coal_time:
            current_coal = 0
            time_slices.append(coal_time)
            last_coal_time = coal_time
            next_mark = coal_per_time.pop()

    @lru_cache(maxsize=None)
    def discretize(timeval: float) -> int:
        for k, t in enumerate(time_slices):
            if timeval <= t:
                return k
        return len(time_slices)

    # Histogram of coalescence events, index is time slice.
    hist: Dict[int, float] = defaultdict(float)
    for coal_time, coal_count in input_list:
        hist[discretize(coal_time)] += coal_count
    if verbose:
        print("================ HIST ================")
        print(
            orjson.dumps(
                {str(k): v for k, v in hist.items()},
                option=orjson.OPT_INDENT_2 | orjson.OPT_SORT_KEYS,
            ).decode("utf-8")
        )
        print("======================================")
    return time_slices


class CoalSampler:
    """
    Class that abstracts away the sampling of coalescences between pairs of populations
    at a particular (discretized) time.
    """

    ALL_COAL_SUFFIX = "all_coal"

    def __init__(self, num_samples: int = 1):
        self.coal_matrices: List[NDArray] = []
        self.num_samples = num_samples

    def init(self, nstates: int, ntimes: int, num_trees: List[int], seed: int):
        raise NotImplementedError("Derived class must implement")

    def add_sample(self, tree_number: int, state: int, tau: float, weight: float):
        raise NotImplementedError("Derived class must implement")

    def merge(self, other_samplers: List["CoalSampler"]):
        raise NotImplementedError("Derived class must implement")

    def finalize(self, jobs: int):
        raise NotImplementedError("Derived class must implement")

    def copy(self):
        return deepcopy(self)

    @staticmethod
    def hash_tree_list(tree_list: Iterable[int]) -> str:
        m = hashlib.md5()
        m.update(str(sorted(tree_list)).encode("utf-8"))
        return m.hexdigest()

    def sample_hashes(self) -> List[str]:
        """
        A one-way hash of the trees used by each of the bootstrap samples, for validation against
        another run (to ensure equivalence).
        """
        raise NotImplementedError("Child class must implement.")


class StandardCoalSampler(CoalSampler):
    """
    Sample all coalescences and convert into a single coalescence matrix.
    """

    def init(self, nstates: int, ntimes: int, num_trees: List[int], seed: int):
        self.num_trees = num_trees
        self.coal_matrices.append(np.zeros((nstates, ntimes)))

    def add_sample(self, tree_number: int, state: int, tau: float, weight: float):
        self.coal_matrices[0][state, int(tau)] += weight

    def merge(self, other_samplers: List["CoalSampler"]):
        for s in other_samplers:
            self.coal_matrices[0] += s.coal_matrices[0]

    def sample_hashes(self) -> List[str]:
        return [CoalSampler.hash_tree_list(range(sum(self.num_trees)))]

    def finalize(self, jobs: int):
        pass


def random_cmatrix(
    data_points: List[Tuple[int, int, float]], nstates: int, ntimes: int
):
    N = len(data_points)
    cmatrix = np.zeros((nstates, ntimes))
    for idx in np.random.randint(low=0, high=N - 1, size=N):
        state, tau, weight = data_points[idx]
        cmatrix[state, tau] += weight
    return cmatrix


class BootstrapCoalSampler(CoalSampler):
    """
    Same as StandardCoalSampler, but also does K bootstrap samples of size N.
    """

    def __init__(self, num_sample_sets: int):
        super().__init__(num_sample_sets)
        self.num_sample_sets = num_sample_sets
        self.coal_matrices = []

    def init(self, nstates: int, ntimes: int, num_trees: List[int], seed: int):
        self.nstates = nstates
        self.ntimes = ntimes
        for i in range(self.num_sample_sets):
            self.coal_matrices.append(np.zeros((nstates, ntimes)))

        # Compute the list of indexes that will be used for each.
        total_trees = sum(num_trees)
        self.num_copies = [
            [0 for _ in range(total_trees)] for _ in range(self.num_sample_sets)
        ]
        np.random.seed(seed)
        for i in range(self.num_sample_sets):
            for j in np.random.randint(low=0, high=total_trees - 1, size=total_trees):
                self.num_copies[i][j] += 1

    def add_sample(self, tree_number: int, state: int, tau: float, weight: float):
        for i in range(self.num_sample_sets):
            copies = self.num_copies[i][tree_number]
            if copies > 0:
                self.coal_matrices[i][state, int(tau)] += copies * weight

    def merge(self, other_samplers: List["CoalSampler"]):
        for s in other_samplers:
            for i in range(0, self.num_sample_sets):
                self.coal_matrices[i] += s.coal_matrices[i]

    def finalize(self, jobs: int):
        pass

    def sample_hashes(self) -> List[str]:
        result = []
        for i in range(self.num_sample_sets):
            tree_list = []
            for j in range(len(self.num_copies[i])):
                for _ in range(self.num_copies[i][j]):
                    tree_list.append(j)
            result.append(CoalSampler.hash_tree_list(tree_list))
        return result

    def copy(self):
        rv = BootstrapCoalSampler(self.num_sample_sets)
        rv.coal_matrices = deepcopy(self.coal_matrices)
        rv.nstates = self.nstates
        rv.ntimes = self.ntimes
        # INTENTIONAL REFERENCE COPY! This is huge, and immutable.
        rv.num_copies = self.num_copies


class JackknifeCoalSampler(CoalSampler):
    """
    Same as StandardCoalSampler, but also does K jackknife samples of N-1 trees
    each, where N is the block count (each block contains >= 1 tree).
    """

    def __init__(self, num_blocks: int):
        super().__init__(num_blocks)
        self.num_blocks = num_blocks

    def init(self, nstates: int, ntimes: int, num_trees: List[int], seed: int):
        self.num_trees = sum(num_trees)
        self.trees_per_block = (
            self.num_trees + (self.num_blocks - 1)
        ) // self.num_blocks
        for i in range(self.num_blocks):
            self.coal_matrices.append(np.zeros((nstates, ntimes)))
        self.next_tree_at = self.trees_per_block
        self.current_leave_out = 1
        print(
            f"Tree jackknife (trees={self.num_trees}, per_block={self.trees_per_block})"
        )

    def add_sample(self, tree_number: int, state: int, tau: float, weight: float):
        while tree_number >= self.next_tree_at:
            self.current_leave_out += 1
            self.next_tree_at += self.trees_per_block
        for i in range(0, self.num_blocks):
            if i != self.current_leave_out:
                self.coal_matrices[i][state, int(tau)] += weight

    def merge(self, other_samplers: List["CoalSampler"]):
        for s in other_samplers:
            for i in range(0, self.num_blocks):
                self.coal_matrices[i] += s.coal_matrices[i]

    def finalize(self, jobs: int):
        pass

    def sample_hashes(self) -> List[str]:
        result: List[str] = []
        # TODO: implement the hashing of trees for this class so we can validate the result.
        # The Jackknife is not officially supported by mrpast, so may alternatively remove this functionality.
        print(f"WARNING: JackknifeCoalSampler does not properly support sample_hashes")
        return result


# Merge a bunch of coal files (assumed to be in the same order, for the same number of trees)
# by summing the coalescence counts
def merge_coals(coal_filenames: List[str], new_coal_filename: str):
    file_objs = [open(fn) for fn in coal_filenames]
    with open(new_coal_filename, "wb") as fout:
        for lines in zip(*file_objs):
            combined: Dict[Tuple[int, int, float], float] = defaultdict(float)
            for line in lines:
                for pop1, pop2, ctime, cweight in orjson.loads(line):
                    combined[(pop1, pop2, ctime)] += cweight
            fout.write(
                orjson.dumps(
                    [(p1, p2, ct, cw) for (p1, p2, ct), cw in combined.items()],
                    option=orjson.OPT_NAIVE_UTC | orjson.OPT_APPEND_NEWLINE,
                )
            )


def process_coal_file(
    sampler: CoalSampler,
    deme_pair_to_index: DemePairIndex,
    filename: str,
    time_slices: List[float],
    trees_before: int,
):
    @lru_cache(maxsize=None)
    def discretize(timeval: float) -> int:
        for k, t in enumerate(time_slices):
            if timeval <= t:
                return k
        return len(time_slices)

    tree_number = trees_before
    with open(filename) as f:
        for line in tqdm(f):
            coal_list = orjson.loads(line)
            for pop_i, pop_j, t, w in coal_list:
                pop_i = int(pop_i)
                pop_j = int(pop_j)
                a_b = deme_pair_to_index[pop_i, pop_j]
                # Add a sample for a particular tree, Markov state a_b, discretized
                # time tau, and coalescence "weight" w.
                sampler.add_sample(tree_number, a_b, discretize(t), w)
            tree_number += 1
    return sampler


def sample_coal_matrices(
    coal_filenames: List[str],
    time_slices: List[float],
    deme_pair_to_index: DemePairIndex,
    sampler: CoalSampler,
    jobs: int = 1,
    seed: int = 42,
):
    assert len(coal_filenames) > 0
    # Produce a matrix with nstates rows and ntime_slices columns
    nstates = len(set(deme_pair_to_index.values()))
    ntimes = len(time_slices) + 1
    num_trees = [count_lines(fn) for fn in coal_filenames]

    # Initialize the sampler. This can be slightly expensive if the sampler needs to count
    # all the trees (each line in the coal_filenames).
    sampler.init(nstates, ntimes, num_trees, seed=seed)

    # This is a large amount of data, and iterating the files this way is actually much more
    # effecient than collecting by (non-discretized) time first.
    print(f"Processing {len(coal_filenames)} coalescence files")
    trees_before_list = [0] * len(num_trees)
    for i, t in enumerate(num_trees):
        if i > 0:
            trees_before_list[i] = num_trees[i - 1] + trees_before_list[i - 1]

    args = [
        (deepcopy(sampler), deme_pair_to_index, fn, time_slices, trees_before)
        for i, (fn, trees_before) in enumerate(zip(coal_filenames, trees_before_list))
    ]
    with Pool(jobs) as p:
        forked_samplers = p.starmap(process_coal_file, args)

    print(f"Merging and finalizing {len(forked_samplers)} samplers...")
    sampler.merge(forked_samplers)
    sampler.finalize(jobs)


def simulate_muts_and_export(
    ts_file: str,
    vcf_file_out: str,
    mut_rate: float = 1e-8,
    leave_out_pops: Set[int] = set(),
    seed: int = 42,
    ploidy: int = 2,
) -> PopMap:
    ts = tskit.load(ts_file)
    if mut_rate > 0.0:
        mts = msprime.sim_mutations(ts, rate=mut_rate, random_seed=seed)
    else:
        mts = ts
    ts = None

    # This is just for sanity checking; we expect and use the sample-based ordering.
    node2indiv = {}
    for indiv in mts.individuals():
        for n in indiv.nodes:
            node2indiv[n] = indiv.id

    def pop_name(pop: tskit.Population):
        if pop.metadata:
            meta = pop.metadata
        else:
            meta = {}
        return meta.get("name", f"pop_{pop.id}")

    pop_map = {}
    ordered_pops = []
    individuals = set()
    for pop in mts.populations():
        p = pop.id
        skip = p in leave_out_pops
        ordered_pops.append((int(p), pop_name(pop), skip))
        if not skip:
            for sample_id in mts.samples(population=p):
                individual_idx = sample_id // ploidy
                pop_map[individual_idx] = int(p)
                individuals.add(individual_idx)
                assert not node2indiv or (
                    node2indiv[sample_id] == individual_idx
                ), "Tree sequence does not follow expected individual/sample ordering"
    print(
        f"Converting {len(list(filter(lambda p: not p[2], ordered_pops)))} haploid samples to VCF"
    )
    ordered_pops = list(sorted(set(ordered_pops)))

    emit_indivs = list(sorted(individuals))
    print(f"Emitting individuals={list(map(lambda i: int(i), emit_indivs))}")
    with open(vcf_file_out, "w") as f:
        mts.write_vcf(
            f,
            individuals=emit_indivs,
            individual_names=list(map(lambda i: f"tsk_{i}", emit_indivs)),
        )

    # This assumes that tskit is emitting the samples in the order that we requested
    new_pop_list: List[List[int]] = [[] for _ in range(len(ordered_pops))]
    for new_indiv_idx, old_indiv_idx in enumerate(emit_indivs):
        pop_idx = pop_map[old_indiv_idx]
        new_pop_list[pop_idx].append(new_indiv_idx)

    print("Ordered populations:")
    pop_names = []
    for pop_id, pop_name, skip in sorted(set(ordered_pops)):
        print(f"  {pop_id}: {pop_name}{' SKIPPING' if skip else ''}")
        pop_names.append(pop_name)
    assert len(new_pop_list) == len(pop_names)
    return PopMap(mapping=new_pop_list, names=pop_names)
