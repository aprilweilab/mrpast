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
#ifndef MRPAST_SOLVE_H
#define MRPAST_SOLVE_H

#include <limits>
extern "C" {
#include <nlopt.h>
}

#include "objective.h"

/**
 * Solve (minimize) the cost function. NOT THREAD SAFE.
 */
double solve(NegLogLikelihoodCostFunctor& cost,
             std::vector<double>& paramVector,
             nlopt_algorithm optAlg,
             bool verbose,
             double timeoutSeconds = std::numeric_limits<double>::max());

#endif /* MRPAST_SOLVE_H */