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

try:
    from yaml import CLoader as Loader, CDumper as Dumper  # type: ignore
except ImportError:
    from yaml import Loader, Dumper  # type: ignore
from yaml import load, dump
from typing import Dict, Any, List, Tuple, Optional, Union
from dataclasses import dataclass, field
from dataclasses_json import dataclass_json
from enum import Enum
import random
import json
import os


DEFAULT_PLOIDY = 2


def load_model_config(filename: str) -> Dict[str, Any]:
    with open(filename) as f:
        config = load(f, Loader=Loader)
    assert (
        "migration" in config
    ), "Model configuration must contain at least one migration matrix"
    assert (
        "coalescence" in config
    ), "Model configuration must contain at least one coalescence matrix"
    return config


class ParameterKind(str, Enum):
    PARAM_KIND_MIGRATION = "migration"
    PARAM_KIND_COAL = "coalescence"
    PARAM_KIND_GROWTH = "growth"
    PARAM_KIND_ADMIXTURE = "admixture"


class TimeSliceAdjustment(str, Enum):
    GROWTH_RATE = "growth_rate"


@dataclass_json
@dataclass
class FloatParameter:
    ground_truth: float
    lb: float
    ub: float
    index: int


DemePairIndex = Dict[Tuple[int, int], int]


@dataclass_json
@dataclass
class SymbolicEpochs:
    # Transition time from epoch n to epoch n+1
    epoch_times: List[FloatParameter]

    @staticmethod
    def from_config(config_entry: Optional[List[Dict[str, Any]]]) -> "SymbolicEpochs":
        entry = [] if config_entry is None else config_entry
        epoch_times = [FloatParameter.from_dict(p) for p in entry]  # type: ignore
        return SymbolicEpochs(epoch_times=epoch_times)

    @property
    def num_epochs(self):
        return len(self.epoch_times) + 1


def resolve_name(name: str, name2index: Dict[str, int]) -> int:
    if name not in name2index:
        raise RuntimeError(f"Unknown deme name: {name}")
    return name2index[name]


@dataclass_json
@dataclass
class ParamRef:
    param: int


class ResolveableEntry:
    def resolve_names(self, name2index: Dict[str, int]):
        raise NotImplementedError("Derived class must implement")


@dataclass_json
@dataclass
class DemeDemeEntry(ResolveableEntry):
    epoch: int
    source: Union[int, str]
    dest: Union[int, str]
    rate: Union[float, ParamRef]

    def resolve_names(self, name2index: Dict[str, int]):
        if isinstance(self.source, str):
            self.source = resolve_name(self.source, name2index)
        if isinstance(self.dest, str):
            self.dest = resolve_name(self.dest, name2index)

    def unresolve_names(self, names: List[str]):
        if isinstance(self.source, int):
            self.source = names[self.source]
        if isinstance(self.dest, int):
            self.dest = names[self.dest]


@dataclass_json
@dataclass
class DemeRateEntry(ResolveableEntry):
    epoch: int
    deme: Union[int, str]
    rate: Union[float, ParamRef]

    def resolve_names(self, name2index: Dict[str, int]):
        if isinstance(self.deme, str):
            self.deme = resolve_name(self.deme, name2index)

    def unresolve_names(self, names: List[str]):
        if isinstance(self.deme, int):
            self.deme = names[self.deme]


@dataclass
class ParamContainer:
    def get_parameter(self, param_idx: int) -> Optional[FloatParameter]:
        assert hasattr(self, "parameters")
        assert param_idx > 0
        for p in self.parameters:
            if p.index == param_idx:
                return p
        return None

    def resolve_names(self, name2index: Dict[str, int]):
        assert hasattr(self, "entries")
        [e.resolve_names(name2index) for e in self.entries]

    def unresolve_names(self, names: List[str]):
        assert hasattr(self, "entries")
        [e.unresolve_names(names) for e in self.entries]

    @staticmethod
    def _from_config(config_entry: Dict[str, Any], entry_type, base_type):
        parameters = []
        for p in config_entry.get("parameters", []):
            parameters.append(FloatParameter.from_dict(p))  # type: ignore
        assert len(parameters) == len(
            set(map(lambda p: p.index, parameters))
        ), "Duplicate parameter index"
        entries = []
        for entry in config_entry.get("entries", []):
            entries.append(entry_type.from_dict(entry))
        result = base_type(entries=entries, parameters=parameters)
        for e in result.entries:
            if isinstance(e, ParamRef):
                param = result.get_parameter(e.param)
                assert param is not None, f"Invalid parameter reference: {e.param}"
                assert (
                    param.lb <= param.ub
                ), "Invalid parameter bound: lower must be <= upper"
        return result


@dataclass
class DemeDemeRates(ParamContainer):
    entries: List[DemeDemeEntry]
    parameters: List[FloatParameter]

    # Oldest epoch number referenced.
    @property
    def max_epoch(self) -> int:
        return max(map(lambda e: e.epoch, self.entries))

    # Largest deme number referenced.
    @property
    def max_deme(self) -> int:
        return max(
            max(map(lambda e: int(e.source), self.entries)),
            max(map(lambda e: int(e.dest), self.entries)),
        )

    def get_entry(self, epoch: int, source: int, dest: int) -> Optional[DemeDemeEntry]:
        for e in self.entries:
            if e.epoch == epoch and e.source == source and e.dest == dest:
                return e
        return None

    @staticmethod
    def from_config(config_entry: Dict[str, Any]) -> "DemeDemeRates":
        return ParamContainer._from_config(config_entry, DemeDemeEntry, DemeDemeRates)


@dataclass
class DemeRates(ParamContainer):
    entries: List[DemeRateEntry]
    parameters: List[FloatParameter]

    # Oldest epoch number referenced.
    @property
    def max_epoch(self) -> int:
        return max(list(map(lambda e: e.epoch, self.entries)) + [0])

    # Largest deme number referenced.
    @property
    def max_deme(self) -> int:
        return max(list(map(lambda e: int(e.deme), self.entries)) + [0])

    def get_entry(self, epoch: int, deme: int) -> Optional[DemeRateEntry]:
        for e in self.entries:
            if e.epoch == epoch and e.deme == deme:
                return e
        return None

    @staticmethod
    def from_config(config_entry: Dict[str, Any]) -> "DemeRates":
        return ParamContainer._from_config(config_entry, DemeRateEntry, DemeRates)


@dataclass
class UserModel:
    ploidy: int
    pop_count: int
    pop_names: List[str]
    migration: DemeDemeRates
    coalescence: DemeRates
    epochs: SymbolicEpochs
    growth: DemeRates
    pop_convert: List[List[int]]

    # TODO: validation
    # 1. admixture proportions should sum to 1
    # 2. admixture epochs >= 1 (cannot do "after" epoch 0)
    # 3. admixture source != dest
    # 4. epoch/demes counts < num pops / epochs

    @staticmethod
    def from_file(filename: str) -> "UserModel":
        # Load the configuration and create the symbolic model from it.
        config = load_model_config(filename)
        mig = DemeDemeRates.from_config(config.get("migration", {}))
        coal = DemeRates.from_config(config.get("coalescence", {}))
        assert (
            len(coal.entries) > 0
        ), "Must have 'coalescence' key and at least one entry"
        epochs = SymbolicEpochs.from_config(config.get("epochTimeSplit", []))
        grow = DemeRates.from_config(config.get("growth", {}))
        pop_names = config.get("pop_names", [])
        pop_count = config.get("pop_count", len(pop_names))
        assert (
            pop_count > 0
        ), f"Model with no populations: specify 'pop_names' or 'pop_count' in your model"
        if not pop_names:
            pop_names = [f"pop_{i}" for i in range(pop_count)]
        ploidy = int(config.get("ploidy", DEFAULT_PLOIDY))
        pop_convert = config.get("populationConversion", []) or []

        result = UserModel(
            ploidy=ploidy,
            pop_count=pop_count,
            pop_names=pop_names,
            migration=mig,
            coalescence=coal,
            epochs=epochs,
            growth=grow,
            pop_convert=pop_convert,
        )
        result.resolve_names()

        # Validation
        nepoch = result.epochs.num_epochs
        assert (
            result.ploidy >= 1 and result.ploidy <= 8
        ), f"Unexpected ploidy value of {result.ploidy}"
        assert (
            nepoch == len(result.pop_convert) + 1
        ), "For k epochs, there must be k-1 populationConversion lists"
        assert (
            result.migration.max_deme < result.pop_count
        ), f"Migration entries reference a deme {result.migration.max_deme} that exceeds number of populations {pop_count}"
        assert (
            result.migration.max_epoch < nepoch
        ), f"Migration entries reference an epoch {result.migration.max_epoch} that exceeds number of populations {nepoch}"
        assert (
            result.coalescence.max_deme < pop_count
        ), f"Coalescence entries reference a deme {result.coalescence.max_deme} that exceeds number of populations {pop_count}"
        assert (
            result.coalescence.max_epoch < nepoch
        ), f"Coalescence entries reference an epoch {result.coalescence.max_epoch} that exceeds number of populations {nepoch}"
        assert (
            result.growth.max_deme < pop_count
        ), f"Growth entries reference a deme {result.growth.max_deme} that exceeds number of populations {pop_count}"
        assert (
            result.growth.max_epoch < nepoch
        ), f"Growth entries reference an epoch {result.growth.max_epoch} that exceeds number of populations {nepoch}"

        return result

    def resolve_names(self):
        name2index = {k: v for v, k in enumerate(self.pop_names)}
        self.migration.resolve_names(name2index)
        self.coalescence.resolve_names(name2index)
        self.growth.resolve_names(name2index)

    def unresolve_names(self):
        self.migration.unresolve_names(self.pop_names)
        self.coalescence.unresolve_names(self.pop_names)
        self.growth.unresolve_names(self.pop_names)

    @property
    def num_epochs(self) -> int:
        return self.epochs.num_epochs

    @property
    def num_demes(self) -> int:
        return self.pop_count

    def get_pair_ordering(self) -> DemePairIndex:
        """
        Returns a map from (i,j) indexes to a fixed index that orders all such pairs.
        (i,j) and (j,i) map to the same index.
        """
        ordered = {}
        counter = 0
        for i in range(self.num_demes):
            for j in range(i, self.num_demes):
                ordered[i, j] = counter
                ordered[j, i] = counter
                counter += 1
        return ordered

    def to_solver_model(self, generate_ground_truth: bool = False):
        return ModelSolverInput.construct(
            self,
            init_from_ground_truth=generate_ground_truth,
            pop_names=self.pop_names,
            ploidy=self.ploidy,
        )


@dataclass_json
@dataclass
class VariableApplication:
    """
    An application of a variable to a matrix. This is interpreted as:
      M[i, j] += (coeff * variable)
    In this way variables can be reused across a matrix when there is dependence between cells.
    """

    epoch: int
    i: int
    j: int
    coeff: float
    adjustment: Optional[TimeSliceAdjustment] = None


@dataclass_json
@dataclass
class BoundedVariable:
    """
    Initial value with boundary constraints
    """

    init: float
    lb: float
    ub: float
    kind: str
    kind_index: int
    apply_to: List[VariableApplication] = field(default_factory=list)
    description: str = ""
    ground_truth: Optional[float] = None
    final: Optional[float] = None

    @staticmethod
    def from_float_parameter(
        param: FloatParameter,
        kind: str,
        kind_index: int,
        init_random: bool = True,
        desc: str = "",
    ) -> "BoundedVariable":
        if init_random:
            init = random.uniform(param.lb, param.ub)
        else:
            init = param.ground_truth
        return BoundedVariable(
            init=init,
            lb=param.lb,
            ub=param.ub,
            kind=kind,
            kind_index=kind_index,
            description=desc,
            ground_truth=param.ground_truth,
        )

    def randomize(self) -> "BoundedVariable":
        return BoundedVariable(
            init=random.uniform(self.lb, self.ub),
            lb=self.lb,
            ub=self.ub,
            kind=self.kind,
            kind_index=self.kind_index,
            apply_to=self.apply_to,
            description=self.description,
            ground_truth=self.ground_truth,
        )


# Instead of enforcing that here with parameter applications, the solver enforces that later.
def construct_stoch_matrix(
    model: UserModel,
    epoch: int,
    M_parameters: Dict[int, BoundedVariable],
    Q_parameters: Dict[int, BoundedVariable],
    G_parameters: Dict[int, BoundedVariable],
    init_from_ground_truth: bool = False,
):
    """
    The output is populated in M_parameters and Q_parameters which contain all of the model parameters
    that were bound to at least one use in the model.
    """
    deme_pair_index = model.get_pair_ordering()
    nstates = len(set(deme_pair_index.values()))

    # Add symbolic state transitions for moving between demes.
    # We are moving an individual from i->j
    for i in range(model.num_demes):
        for j in range(model.num_demes):
            m_param_entry = model.migration.get_entry(epoch, i, j)
            if m_param_entry is not None:
                assert isinstance(m_param_entry.rate, ParamRef)
                migration_param_idx = m_param_entry.rate.param
                assert i != j, f"Migration to self is not allowed (deme {i})"
                parameter = model.migration.get_parameter(migration_param_idx)
                assert parameter is not None
                # Migration rate is a variable to be optimized. The actual stochastic state transition rates
                # are derived from these by a linear combination.
                if migration_param_idx not in M_parameters:
                    M_parameters[migration_param_idx] = (
                        BoundedVariable.from_float_parameter(
                            parameter,
                            ParameterKind.PARAM_KIND_MIGRATION,
                            int(migration_param_idx),
                            not init_from_ground_truth,
                            desc=f"Migration rate from {i}->{j}",
                        )
                    )
                migrate_boundedvar = M_parameters[migration_param_idx]

                # The second individual stays at k, for all possible k values.
                for k in range(model.num_demes):
                    from_idx = deme_pair_index[i, k]  # Moving from state
                    to_idx = deme_pair_index[j, k]  # Moving to state
                    assert from_idx != to_idx
                    # Both nodes started in i (i==k), and one of them is moving to j, but we don't know which one, so
                    # we multiply the coefficient *2 to cover both cases.
                    if (k == i) and (j != i):
                        multiplier = 2
                    else:
                        multiplier = 1

                    migrate_boundedvar.apply_to.append(
                        VariableApplication(
                            i=from_idx, j=to_idx, coeff=multiplier, epoch=epoch
                        )
                    )

    # Add symbolic state transitions for the coalescent states
    for i in range(model.num_demes):
        c_param_entry = model.coalescence.get_entry(epoch, i)
        if c_param_entry is not None:
            assert isinstance(c_param_entry.rate, ParamRef)
            coal_param_idx = c_param_entry.rate.param
            parameter = model.coalescence.get_parameter(coal_param_idx)
            assert parameter is not None

            # Coalescence rate is a variable to be optimized. The actual stochastic state transition rates
            # are derived from these by a linear combination.
            if coal_param_idx not in Q_parameters:
                Q_parameters[coal_param_idx] = BoundedVariable.from_float_parameter(
                    parameter,
                    ParameterKind.PARAM_KIND_COAL,
                    int(coal_param_idx),
                    not init_from_ground_truth,
                    desc=f"Coalescence rate for deme {i}",
                )
            coalesce_boundedvar = Q_parameters[coal_param_idx]

            # Coalescence state can only be reached from the state where both individuals are in the same deme.
            from_idx = deme_pair_index[i, i]
            # Case: they do coalesce, which moves to the "final" state (last column)
            coalesce_boundedvar.apply_to.append(
                VariableApplication(i=from_idx, j=nstates, coeff=1, epoch=epoch)
            )

            # Optional growth-rate adjustment the coalescent rate. At discretized time t, the adjustment
            # to the coalescent rate is 1/e^{-alpha*t}, where alpha is defined per-deme.
            g_param_entry = model.growth.get_entry(epoch, i)
            if g_param_entry is not None:
                assert isinstance(g_param_entry.rate, ParamRef)
                growth_param_idx = g_param_entry.rate.param
                if growth_param_idx not in G_parameters:
                    growth_param = model.growth.get_parameter(growth_param_idx)
                    assert growth_param is not None
                    G_parameters[growth_param_idx] = (
                        BoundedVariable.from_float_parameter(
                            growth_param,
                            ParameterKind.PARAM_KIND_GROWTH,
                            int(growth_param_idx),
                            not init_from_ground_truth,
                            desc=f"Growth rate for deme {i}",
                        )
                    )
                growth_boundedvar = G_parameters[growth_param_idx]
                # We apply the growth to both the "coalescence" and "stay put" cases, because they should be
                # symmetric for the Q-matrix to work properly.
                growth_boundedvar.apply_to.append(
                    VariableApplication(
                        i=from_idx,
                        j=nstates,
                        coeff=1,
                        epoch=epoch,
                        adjustment=TimeSliceAdjustment.GROWTH_RATE,
                    )
                )


# This assumes that all the population conversion matrices are written in terms of epoch0, instead
# of being sequential compositions epoch0 -> epoch1 -> epoch2 -> ...
def pop_conv_to_state_space(model: UserModel, from_epoch: int) -> List[int]:
    pop_convert_single = model.pop_convert[from_epoch]
    deme_pair_index = model.get_pair_ordering()
    nstates = len(set(deme_pair_index.values()))
    result = [-1] * nstates
    for i in range(len(pop_convert_single)):
        for j in range(len(pop_convert_single)):
            from_state = deme_pair_index[i, j]
            dest_i = pop_convert_single[i]
            dest_j = pop_convert_single[j]
            to_state = deme_pair_index[dest_i, dest_j]
            assert result[from_state] == -1 or result[from_state] == to_state
            result[from_state] = to_state
    return result


CMatrix = List[List[float]]


@dataclass_json
@dataclass
class ModelSolverInput:
    # The suffix of each field represents its units. A suffix of "foo_bar" represents Foo*Bar units. A suffix
    # of "foo__bar"  is Foo/Bar ("Foos per Bar"). The abbreviations for units are:
    #  ne  == Number of effective_popsize individuals. So 1 Ne == effective_popsize individuals
    #  gen == Generations
    #  idv == Individuals
    # Unitless fields have no suffix.

    # (Symbolic) time when transitions between epochs occurs (nepochs - 1)
    epoch_times_gen: List[BoundedVariable]
    # (Symbolic) paramters that affect transition probabilities between states (the Q-matrix) (depends on migration rate topologies)
    smatrix_values_ne__gen: List[BoundedVariable]
    # (Concrete) coalescence counts matrix, ((nstates-1) x ntimeslices)
    coal_count_matrices: List[CMatrix]
    # (Concrete) times when transitions between time slices/buckets occur (ntimeslices-1)
    time_slices_gen: List[float]
    # (Concrete) population conversion matrix, by epoch. E.g., the 0th vector (row) describes how populations
    #    move from epoch0->epoch1.
    pop_convert: List[List[int]]
    # How to make use of coal_count_matrices
    observation_mode: str = "average"
    # How to interpret the samples in coal_count_matrices
    sampling_description: Optional[str] = None
    # Sampling hashes: a hash of the trees involved in each sample
    sampling_hashes: Optional[List[str]] = None

    # Extra data: unused by the solver, but useful to have as pass-through.
    pop_names: Optional[List[str]] = None
    ploidy: int = 2

    @staticmethod
    def construct(
        model: UserModel,
        init_from_ground_truth: bool = False,
        pop_names: Optional[List[str]] = None,
        ploidy: int = 2,
    ):
        # Construct the symbolic stochastic matrix values based on the inputs
        Q_parameters: Dict[int, BoundedVariable] = {}
        M_parameters: Dict[int, BoundedVariable] = {}
        G_parameters: Dict[int, BoundedVariable] = {}
        [
            construct_stoch_matrix(
                model,
                e,
                M_parameters,
                Q_parameters,
                G_parameters,
                init_from_ground_truth=init_from_ground_truth,
            )
            for e in range(model.epochs.num_epochs)
        ]
        # The order of smatrix parameters MUST stay as:
        #   all migration rate parameters
        #   all coalescent rate parameters
        smatrix_values = []
        for mig_idx, param in sorted(M_parameters.items(), key=lambda t: t[0]):
            smatrix_values.append(param)
        for coal_idx, param in sorted(Q_parameters.items(), key=lambda t: t[0]):
            smatrix_values.append(param)
        for growth_idx, param in sorted(G_parameters.items(), key=lambda t: t[0]):
            smatrix_values.append(param)

        # Generate random initial values for the epoch times
        epoch_times = [
            BoundedVariable.from_float_parameter(
                t, "epoch", int(i), not init_from_ground_truth, desc=f"Epoch {i}->{i+1}"
            )
            for i, t in enumerate(model.epochs.epoch_times)
        ]
        # Convert the population conversion matrices for single individuals to paired, to match our state-space
        pop_convert = [
            pop_conv_to_state_space(model, e) for e in range(model.num_epochs - 1)
        ]
        # Construct the (JSON-ifiable) input object that will be sent to the solver component
        return ModelSolverInput(
            epoch_times_gen=epoch_times,
            smatrix_values_ne__gen=smatrix_values,
            coal_count_matrices=[],
            time_slices_gen=[],
            pop_convert=pop_convert,
            pop_names=pop_names,
            ploidy=ploidy,
        )

    def randomize(self) -> "ModelSolverInput":
        return ModelSolverInput(
            epoch_times_gen=[e.randomize() for e in self.epoch_times_gen],
            smatrix_values_ne__gen=[v.randomize() for v in self.smatrix_values_ne__gen],
            coal_count_matrices=self.coal_count_matrices,
            time_slices_gen=self.time_slices_gen,
            pop_convert=self.pop_convert,
            sampling_description=self.sampling_description,
            sampling_hashes=self.sampling_hashes,
        )

    @property
    def num_epochs(self):
        return len(self.epoch_times_gen) + 1


# TODO: create class that represents a solver output, and put the below methods on it.


def count_output_params_from_obj(solver_out_obj: Dict[str, Any]) -> int:
    return len(solver_out_obj["smatrix_values"]) + len(
        solver_out_obj.get("epoch_times", []) or []
    )


def count_output_params(solver_output: str):
    with open(solver_output) as f:
        output_json = json.load(f)
    return count_output_params_from_obj(output_json)


def get_negLL_from_obj(solver_out_obj: Dict[str, Any]) -> float:
    return solver_out_obj["negLL"]


def get_negLL(solver_output: str) -> float:
    with open(solver_output) as f:
        output_json = json.load(f)
    return get_negLL_from_obj(output_json)


def copy_model_with_outputs(
    input_model: str,
    solver_output: str,
    output_model: str,
    check_descriptions: bool = True,
    errors: Dict[int, float] = {},
):
    config = load_model_config(input_model)
    with open(solver_output) as f:
        solver_out_obj = json.load(f)

    # The ordering of parameters is as follows:
    #   all migration rates
    #   all coalescence rates
    non_epoch_outparams = solver_out_obj["smatrix_values"]
    out_param_idx = 0
    applied_errors = 0
    for model_param in config["migration"]["parameters"]:
        output_param = non_epoch_outparams[out_param_idx]
        lb, ub = model_param["lb"], model_param["ub"]
        model_param["ground_truth"] = max(
            lb, min(ub, output_param["final"] + errors.get(out_param_idx, 0.0))
        )
        applied_errors += 1 if out_param_idx in errors else 0
        out_param_idx += 1
        assert not check_descriptions or "Migration rate" in output_param["description"]
    for model_param in config["coalescence"]["parameters"]:
        output_param = non_epoch_outparams[out_param_idx]
        lb, ub = model_param["lb"], model_param["ub"]
        model_param["ground_truth"] = max(
            lb, min(ub, output_param["final"] + errors.get(out_param_idx, 0.0))
        )
        applied_errors += 1 if out_param_idx in errors else 0
        out_param_idx += 1
        assert (
            not check_descriptions or "Coalescence rate" in output_param["description"]
        )
    if "epochTimeSplit" in config:
        out_param_idx = 0
        for model_param in config["epochTimeSplit"] or []:
            output_param = solver_out_obj["epoch_times"][out_param_idx]
            lb, ub = model_param["lb"], model_param["ub"]
            model_param["ground_truth"] = max(
                lb, min(ub, output_param["final"] + errors.get(out_param_idx, 0.0))
            )
            applied_errors += 1 if out_param_idx in errors else 0
            out_param_idx += 1

    assert applied_errors == len(errors), "Invalid parameter index in errors?"

    if os.path.exists(output_model):
        raise RuntimeError(
            f"Will not replace existing {output_model}; remove file and try again"
        )
    with open(output_model, "w") as f:
        f.write(dump(config, Dumper=Dumper))


def print_model_warnings(model_filename: str):
    model = UserModel.from_file(model_filename)
    last_ub = 0.0
    for i, et in enumerate(model.epochs.epoch_times):
        if et.lb < last_ub:
            print(
                f"WARNING: Time between epochs {i} and {i+1} has overlapping bounds with the previous epoch. "
                "This can cause problems where the solver picks out of order epoch times. Consider setting "
                "bounds so that the time axis is partitioned among the epoch parameters."
            )
        last_ub = et.ub
    min_suggested_ne = 250
    for i, param in enumerate(model.coalescence.parameters):
        if param.ground_truth >= (1 / (model.ploidy * min_suggested_ne)):
            print(
                f"WARNING: Coalescence rate parameter {i} results in an Ne less than the suggested minimum ({min_suggested_ne})"
            )


@dataclass_json
@dataclass
class PopMap:
    """
    Maps populations to the individuals contained within, and the names of the populations.
    """

    mapping: List[List[int]]
    names: List[str]

    @property
    def num_pops(self):
        assert len(self.mapping) == len(self.names)
        return len(self.mapping)
