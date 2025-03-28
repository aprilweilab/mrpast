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

/******************************** WARNING ************************************
 * This file is not thread-safe. The way we use NLOPT needs some global state
 * so you will get wackiness if you call any of the functions in this file
 * from multiple threads.
 */
#include <cassert>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <signal.h>
#include <stdexcept>
#include <vector>

#include "objective.h"
#include "solve.h"

constexpr double numericDiffRelStepSize = 1e-12;
// The solver terminates when the optimum improves by less than this amount.
constexpr double solverMinChange = 1e-12;
// The solver terminates when there are more objective function evaluations than
// this.
constexpr int maxObjectiveExecs = 5000000;

constexpr size_t SHOW_STATUS_EVERY = 10000;

double nlopt_objective(unsigned n, const double* params, double* grad, void* my_func_data) {
    NegLogLikelihoodCostFunctor* functor = static_cast<NegLogLikelihoodCostFunctor*>(my_func_data);
    assert(n == functor->m_schema.totalParams());
    const double result = (*functor)(params);
    functor->m_minValue = std::min<double>(result, functor->m_minValue);
    if (functor->m_numObjCalls == 1 || functor->m_numObjCalls % SHOW_STATUS_EVERY == 0) {
        std::cerr << "calls: " << functor->m_numObjCalls << ", f() = " << result
                  << ", min(f()) = " << functor->m_minValue << std::endl;
    }
    functor->m_numObjCalls++;
    // If the nlopt method requires a derivative then we use numeric
    // differentiation via the central difference method. Currently
    // derivative-free approaches perform better. If we need to use a
    // derivative-based method we should use something better than central
    // difference.
    if (nullptr != grad) {
        std::vector<double> xModified(n);
        for (size_t i = 0; i < n; i++) {
            xModified[i] = params[i];
        }
        for (size_t i = 0; i < n; i++) {
            const double origX = params[i];
            double stepSize = std::max(numericDiffRelStepSize * origX, MIN_FLOAT_VAL);
            xModified[i] = origX + stepSize;
            const double y1 = (*functor)(xModified.data());
            xModified[i] = origX - stepSize;
            const double y2 = (*functor)(xModified.data());
            xModified[i] = origX;
            grad[i] = (y1 - y2) / stepSize;
        }
    }
    return result;
}

nlopt_opt* global_opt = nullptr;

// Handles ctrl+c interrupts gracefully so we can get the currently-best output
// parameters when we terminate a solver run early.
void interrupt_handler(int signal) {
    (void)signal;
    if (nullptr != global_opt) {
        nlopt_force_stop(*global_opt);
    } else {
        exit(1);
    }
}

double solve(NegLogLikelihoodCostFunctor& cost,
             std::vector<double>& paramVector,
             nlopt_algorithm optAlg,
             const bool verbose,
             double timeoutSeconds) {
    nlopt_opt opt = nlopt_create(optAlg, cost.m_schema.totalParams());

    RELEASE_ASSERT(NLOPT_SUCCESS == nlopt_set_maxtime(opt, timeoutSeconds));

    std::vector<double> lb(cost.m_schema.totalParams());
    std::vector<double> ub(cost.m_schema.totalParams());
    for (size_t i = 0; i < cost.m_schema.totalParams(); i++) {
        cost.m_schema.getBounds(i, lb[i], ub[i]);
    }
    RELEASE_ASSERT(NLOPT_SUCCESS == nlopt_set_lower_bounds(opt, lb.data()));
    RELEASE_ASSERT(NLOPT_SUCCESS == nlopt_set_upper_bounds(opt, ub.data()));
    RELEASE_ASSERT(NLOPT_SUCCESS == nlopt_set_min_objective(opt, nlopt_objective, &cost));
    RELEASE_ASSERT(NLOPT_SUCCESS == nlopt_set_ftol_abs(opt, solverMinChange));
    RELEASE_ASSERT(NLOPT_SUCCESS == nlopt_set_maxeval(opt, maxObjectiveExecs));
    double minf = 0.0;

    global_opt = &opt;
    struct sigaction handler;
    handler.sa_handler = interrupt_handler;
    sigemptyset(&handler.sa_mask);
    handler.sa_flags = 0;
    sigaction(SIGINT, &handler, nullptr);

    cost.resetSolveStats();
    try {
        nlopt_result result = nlopt_optimize(opt, paramVector.data(), &minf);
        if (verbose) {
            std::cerr << "total calls: " << cost.m_numObjCalls << std::endl;
            std::cerr << "found minimum at f(";
            cost.m_schema.dumpParameters(paramVector.data(), paramVector.size());
            std::cerr << ") = " << minf << std::endl;
        }
    } catch (std::exception& e) {
        std::cerr << "nlopt failed: " << e.what() << std::endl;
    }

    global_opt = nullptr;

    return minf;
}
