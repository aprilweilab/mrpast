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
import numpy as np

try:
    from yaml import CLoader as Loader, CDumper as Dumper  # type: ignore
except ImportError:
    from yaml import Loader, Dumper  # type: ignore
from yaml import load, dump
from typing import Dict, Any, List, Tuple, Optional
from numpy.typing import NDArray
from dataclasses import dataclass, field
from dataclasses_json import dataclass_json
from enum import Enum
import random
import itertools
import json
import os


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


class TimeSliceAdjustment(str, Enum):
    GROWTH_RATE = "growth_rate"


@dataclass
class FloatParameter:
    ground_truth: float
    lower_bound: float
    upper_bound: float
    index: int

    @staticmethod
    def from_config(config_entry: Dict[str, Any]) -> "FloatParameter":
        assert int(config_entry["index"]) > 0, "Parameter index must be positive"
        return FloatParameter(
            ground_truth=float(config_entry["ground_truth"]),
            lower_bound=float(config_entry["lb"]),
            upper_bound=float(config_entry["ub"]),
            index=int(config_entry["index"]),
        )


@dataclass
class SymbolicEpochs:
    # Transition time from epoch n to epoch n+1
    epoch_times: List[FloatParameter]

    @staticmethod
    def from_config(config_entry: Optional[List[Dict[str, Any]]]) -> "SymbolicEpochs":
        entry = [] if config_entry is None else config_entry
        epoch_times = [FloatParameter.from_config(p) for p in entry]
        return SymbolicEpochs(epoch_times=epoch_times)

    @property
    def num_epochs(self):
        return len(self.epoch_times) + 1


DemePairIndex = Dict[Tuple[int, int], int]


@dataclass
class SymbolicMatrices:
    parameters: List[FloatParameter]
    matrices: List[NDArray[Any]]

    @property
    def num_epochs(self):
        return len(self.matrices)

    def num_demes(self, epoch: int):
        return len(self.matrices[epoch])

    def assert_square(self, name: str):
        for m in self.matrices:
            assert m.shape[0] == m.shape[1], f"Matrix {name} must be square"

    def assert_vector(self, name: str):
        for m in self.matrices:
            assert len(m.shape) == 1, f"{name} must be vectors (not matrices)"

    def get_parameter(self, param_idx: Optional[int]) -> Optional[FloatParameter]:
        if param_idx == 0 or param_idx is None:
            return None
        assert param_idx > 0
        for p in self.parameters:
            if p.index == param_idx:
                return p
        return None

    def get_pair_ordering(self, epoch: int) -> DemePairIndex:
        """
        Returns a map from (i,j) indexes to a fixed index that orders all such pairs.
        (i,j) and (j,i) map to the same index.
        """
        ordered = {}
        counter = 0
        for i in range(self.num_demes(epoch)):
            for j in range(i, self.num_demes(epoch)):
                ordered[i, j] = counter
                ordered[j, i] = counter
                counter += 1
        return ordered

    @staticmethod
    def from_config(
        config_entry: Dict[str, Any], is_vector: bool = False
    ) -> "SymbolicMatrices":
        assert (
            "parameters" in config_entry
        ), 'Cannot convert configuration to SymbolicMatrices: need "parameters" field'
        if is_vector:
            assert (
                "vectors" in config_entry
            ), 'Cannot convert configuration to SymbolicMatrices: need "vectors" field'
        else:
            assert (
                "matrices" in config_entry
            ), 'Cannot convert configuration to SymbolicMatrices: need "matrices" field'

        parameters = []
        for p in config_entry["parameters"]:
            parameters.append(FloatParameter.from_config(p))
        assert len(parameters) == len(
            set(map(lambda p: p.index, parameters))
        ), "Duplicate parameter index"
        matrices = []
        if is_vector:
            data = config_entry["vectors"]
        else:
            data = config_entry["matrices"]
        for m in data:
            m = np.array(m)
            if is_vector:
                assert len(m.shape) == 1, "Expected D-length vector, not DxD matrix"
            else:
                assert len(m.shape) == 2, "Expected DxD matrix"
            matrices.append(m)

        return SymbolicMatrices(parameters=parameters, matrices=matrices)


# The objective function takes as input values for parameters P_0...P_L and returns the negative log likelihood. The parameters have to be mapped
# into the square stochastic matrix (concretely) based on a mapping from symbolic->concrete. This mapping is _fixed_ at configuration time, we can
# compute it once and then just reuse it.


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
        lb = param.lower_bound
        ub = param.upper_bound
        if init_random:
            init = random.uniform(lb, ub)
        else:
            init = param.ground_truth
        return BoundedVariable(
            init=init,
            lb=lb,
            ub=ub,
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


# The Q-matrix requires that the diagonal be the negative sum of the off-diagonal for each row.
# Instead of enforcing that here with parameter applications, the solver enforces that later.
def construct_stoch_matrix(
    M_symbolic: SymbolicMatrices,
    q_symbolic: SymbolicMatrices,
    g_symbolic: Optional[SymbolicMatrices],
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
    deme_pair_index = M_symbolic.get_pair_ordering(epoch)
    nstates = len(set(deme_pair_index.values()))
    assert M_symbolic.num_demes(epoch) == q_symbolic.num_demes(epoch)
    ndemes = M_symbolic.num_demes(epoch)
    M_symbolic.assert_square("migration")
    q_symbolic.assert_vector("coalescence")
    if g_symbolic is not None:
        g_symbolic.assert_vector("growth")

    # Add symbolic state transitions for moving between demes.
    # We are moving an individual from i->j
    for i in range(ndemes):
        for j in range(ndemes):
            migration_param_idx = M_symbolic.matrices[epoch][i, j]
            can_migration_from_i_to_j = migration_param_idx != 0
            if can_migration_from_i_to_j:
                assert i != j, f"Migration to self is not allowed (deme {i})"
                parameter = M_symbolic.get_parameter(migration_param_idx)
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

                growth_boundedvar = None
                if g_symbolic is not None:
                    growth_param_idx = g_symbolic.matrices[epoch][i]
                    if growth_param_idx > 0:
                        if growth_param_idx not in G_parameters:
                            growth_param = g_symbolic.get_parameter(growth_param_idx)
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

                # The second individual stays at k, for all possible k values.
                for k in range(ndemes):
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
    for i in range(ndemes):
        coal_param_idx = q_symbolic.matrices[epoch][i]
        if coal_param_idx > 0:
            parameter = q_symbolic.get_parameter(coal_param_idx)
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
            if g_symbolic is not None:
                growth_param_idx = g_symbolic.matrices[epoch][i]
                if growth_param_idx > 0:
                    if growth_param_idx not in G_parameters:
                        growth_param = g_symbolic.get_parameter(growth_param_idx)
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
def pop_conv_to_state_space(
    M_symbolic: SymbolicMatrices, pop_convert_single: List[int], from_epoch: int
) -> List[int]:
    deme_pair_index0 = M_symbolic.get_pair_ordering(0)
    nstates0 = len(set(deme_pair_index0.values()))
    result = [-1] * nstates0
    deme_pair_index_current = M_symbolic.get_pair_ordering(from_epoch + 1)
    for i in range(len(pop_convert_single)):
        for j in range(len(pop_convert_single)):
            from_state = deme_pair_index0[i, j]
            dest_i = pop_convert_single[i]
            dest_j = pop_convert_single[j]
            to_state = deme_pair_index_current[dest_i, dest_j]
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
        M_symbolic: SymbolicMatrices,
        q_symbolic: SymbolicMatrices,
        g_symbolic: Optional[SymbolicMatrices],
        epochs_symbolic: SymbolicEpochs,
        pop_convert_single: List[List[int]],
        init_from_ground_truth: bool = False,
        pop_names: Optional[List[str]] = None,
        ploidy: int = 2,
    ):
        assert (
            M_symbolic.num_epochs == q_symbolic.num_epochs == epochs_symbolic.num_epochs
        ), "There must be an equal number of coalescence and migration matrices (one for each epoch)"
        # Construct the symbolic stochastic matrix values based on the inputs
        Q_parameters: Dict[int, BoundedVariable] = {}
        M_parameters: Dict[int, BoundedVariable] = {}
        G_parameters: Dict[int, BoundedVariable] = {}
        [
            construct_stoch_matrix(
                M_symbolic,
                q_symbolic,
                g_symbolic,
                e,
                M_parameters,
                Q_parameters,
                G_parameters,
                init_from_ground_truth=init_from_ground_truth,
            )
            for e in range(M_symbolic.num_epochs)
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
            for i, t in enumerate(epochs_symbolic.epoch_times)
        ]
        # Convert the population conversion matrices for single individuals to paired, to match our state-space
        pop_convert = [
            pop_conv_to_state_space(M_symbolic, pop_convert_single[e], e)
            for e in range(M_symbolic.num_epochs - 1)
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
    config = load_model_config(model_filename)
    last_ub = 0
    for i, et in enumerate(config.get("epochTimeSplit", []) or []):
        if et["lb"] < last_ub:
            print(
                f"WARNING: Time between epochs {i} and {i+1} has overlapping bounds with the previous epoch. "
                "This can cause problems where the solver picks out of order epoch times. Consider setting "
                "bounds so that the time axis is partitioned among the epoch parameters."
            )
        last_ub = et["ub"]
    q_symbolic = SymbolicMatrices.from_config(config["coalescence"], is_vector=True)
    min_suggested_ne = 250
    for i, param in enumerate(q_symbolic.parameters):
        if param.ground_truth >= (1 / (config["ploidy"] * min_suggested_ne)):
            print(
                f"WARNING: Coalescence rate parameter {i} results in an Ne less than the suggested minimum ({min_suggested_ne})"
            )


def validate_model(model_filename: str):
    def _must_have_keys(object: Dict[str, Any], *keys):
        for k in keys:
            assert k in object, f"Required key missing: {k}"

    config = load_model_config(model_filename)
    _must_have_keys(
        config,
        "migration",
        "coalescence",
        "epochTimeSplit",
        "populationConversion",
        "ploidy",
    )
    ploidy = config["ploidy"]
    assert ploidy >= 1 and ploidy < 8, f"Unexpected ploidy value of {ploidy}"
    M_symbolic = SymbolicMatrices.from_config(config["migration"])
    q_symbolic = SymbolicMatrices.from_config(config["coalescence"], is_vector=True)
    if "growth" in config:
        G_symbolic = SymbolicMatrices.from_config(config["growth"], is_vector=True)
    else:
        G_symbolic = None
    assert (
        M_symbolic.num_epochs == q_symbolic.num_epochs
    ), "Coalescence and migration matrices must have same number of epochs"
    assert M_symbolic.num_epochs >= 1, "There must be at least 1 epoch"
    epochs_symbolic = SymbolicEpochs.from_config(config.get("epochTimeSplit", []))
    assert (
        M_symbolic.num_epochs == len(epochs_symbolic.epoch_times) + 1
    ), "For k epochs, there must be k-1 epochTimeSplit parameters"
    popConvert = config.get("populationConversion", [])
    if popConvert is None:
        popConvert = []
    assert (
        M_symbolic.num_epochs == len(popConvert) + 1
    ), "For k epochs, there must be k-1 populationConversion lists"

    M_symbolic.assert_square("migration")
    q_symbolic.assert_vector("coalescence")
    if G_symbolic is not None:
        G_symbolic.assert_vector("growth")

    npops = M_symbolic.num_demes(0)
    for e in range(M_symbolic.num_epochs):
        assert (
            M_symbolic.num_demes(e) == npops
        ), f"Migration rate matrix {e} has wrong dimensions (number of demes)"
        assert (
            q_symbolic.num_demes(e) == npops
        ), f"Coalescence rate matrix {e} has wrong dimensions (number of demes)"
        for i, j in itertools.product(range(npops), repeat=2):
            paramIdx = M_symbolic.matrices[e][i][j]
            if paramIdx != 0:
                assert (
                    M_symbolic.get_parameter(paramIdx) is not None
                ), f"Invalid parameter index {paramIdx} at epoch {e} in migration rates"
            paramIdx = q_symbolic.matrices[e][i]
            if paramIdx != 0:
                assert (
                    q_symbolic.get_parameter(paramIdx) is not None
                ), f"Invalid parameter index {paramIdx} at epoch {e} in coalescence rates"
        if e > 0:
            assert (
                len(popConvert[e - 1]) == npops
            ), f"populationConversion between epoch {e-1} and {e} has wrong number of demes (should be {npops})"
            for i, j in enumerate(popConvert[e - 1]):
                assert (
                    j < npops
                ), f"populationConversion values must be between 0...{npops-1}, saw value {j} for epoch {e-1}->{e}"

    for param in M_symbolic.parameters:
        assert (param.ground_truth >= param.lower_bound) and (
            param.ground_truth <= param.upper_bound
        ), f"Migration parameter {param.index} has a ground truth outside of its bounds"
    for param in q_symbolic.parameters:
        assert (param.ground_truth >= param.lower_bound) and (
            param.ground_truth <= param.upper_bound
        ), f"Coalescence parameter {param.index} has a ground truth outside of its bounds"
    if G_symbolic is not None:
        for param in G_symbolic.parameters:
            assert (param.ground_truth >= param.lower_bound) and (
                param.ground_truth <= param.upper_bound
            ), f"Growth parameter {param.index} has a ground truth outside of its bounds"


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
