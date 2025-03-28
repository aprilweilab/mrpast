/*
 * Migration Rate and Population Size Across Space and Time (mrpast)
 * Copyright (C) 2025 April Wei
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef MRP_SOLVER_COMMON_H
#define MRP_SOLVER_COMMON_H

#include <cmath>
#include <iostream>
#include <limits>

extern "C" {
#include <nlopt.h>
}

#define CHECK_CASTS 0
// Randomly init every time instead of using the "init" value from the input
// JSON
#define RANDOM_INIT 0

#if CHECK_CASTS
inline Eigen::Index SIZE_T_TO_INDEX(const size_t szt) {
    assert((szt) < std::numeric_limits<Eigen::Index>::max());
    return static_cast<Eigen::Index>(szt);
}
#else
#define SIZE_T_TO_INDEX(szt) static_cast<Eigen::Index>(szt)
#endif

class ModelAssertFailure : public std::runtime_error {
public:
    explicit ModelAssertFailure(char const* const message)
        : std::runtime_error(message) {}
};

#define RELEASE_ASSERT(condition)                                                                                      \
    do {                                                                                                               \
        if (!(condition)) {                                                                                            \
            std::cerr << "Release assert failed: " #condition << std::endl;                                            \
            abort();                                                                                                   \
        }                                                                                                              \
    } while (0)

#define MODEL_ASSERT(condition)                                                                                        \
    do {                                                                                                               \
        if (!(condition)) {                                                                                            \
            throw ModelAssertFailure("Model assert failed: " #condition);                                              \
        }                                                                                                              \
    } while (0)

#define MODEL_ASSERT_MSG(condition, msg)                                                                               \
    do {                                                                                                               \
        if (!(condition)) {                                                                                            \
            std::cerr << msg << std::endl;                                                                             \
            throw ModelAssertFailure("Model assert failed: " #condition);                                              \
        }                                                                                                              \
    } while (0)

// Outer-most keys for the JSON input file.
#define UNIT_MUL              "_"
#define UNIT_PER              "__"
#define UNIT_EFFECTIVEPOPSIZE "ne"
#define UNIT_GENERATIONS      "gen"
#define UNIT_INDIVIDUALS      "idv"
constexpr const char* EPOCH_TIMES_KEY = "epoch_times_" UNIT_GENERATIONS;
constexpr const char* SMATRIX_VALS_KEY = "smatrix_values_" UNIT_EFFECTIVEPOPSIZE UNIT_PER UNIT_GENERATIONS;
constexpr const char* COAL_COUNTS_KEY = "coal_count_matrices";
constexpr const char* TIME_SLICES_KEY = "time_slices_" UNIT_GENERATIONS;
constexpr const char* POP_CONVERT_KEY = "pop_convert";
constexpr const char* EFFECTIVE_POPSIZE_KEY = "effective_popsize_" UNIT_INDIVIDUALS;
constexpr const char* OBSERVATION_MODE_KEY = "observation_mode";
constexpr const char* SAMPLING_HASHES_KEY = "sampling_hashes";

// The minimum value that we use in any computations, the square-root of the
// minimum representable floating point value.
static const double FLOAT_EPS = std::numeric_limits<double>::epsilon();
static const double MIN_FLOAT_VAL = std::sqrt(FLOAT_EPS);

// DERIVATIVE-FREE: Generally significantly faster than the gradient-based
// approaches (which I am using numeric differentiation for). The best ones also tend to find lower
// minima than the gradient algorithms. I suspect there are a lot of local minima and we need to explore
// that more:
// - Is it a result of having under-constrained model setup?
// - Is it just a result of having only coalescence times as input (too many deg of freedom)?
// - Is it a result of the numerical properties of the expm() function we're using, or similar
// floating point precision problems?

// constexpr auto MRP_OPT_ALG = NLOPT_LN_COBYLA;
// constexpr auto MRP_OPT_ALG = NLOPT_LN_BOBYQA;
// constexpr auto MRP_OPT_ALG = NLOPT_LN_PRAXIS;
constexpr auto MRP_OPT_ALG = NLOPT_LN_SBPLX; // Best in initial testing
// GRADIENT-BASED:
// constexpr auto MRP_OPT_ALG = NLOPT_LD_LBFGS; // Best of gradient methods in initial testing

#endif /* MRP_SOLVER_COMMON_H */