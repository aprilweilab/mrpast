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
#include <cassert>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <signal.h>
#include <stdexcept>
#include <vector>

#include "args.hxx"

extern "C" {
#include <nlopt.h>
}

#include "common.h"
#include "objective.h"
#include "solve.h"

// Don't allow the solver to run for more than 24 hours by default. User can override with --timeout
constexpr double DEFAULT_TIMEOUT_SECONDS = 24 * 60 * 60;

int main(int argc, char** argv) {
    std::cerr << std::setprecision(100);

    args::ArgumentParser parser("mrpast solver.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::Positional<std::string> infile(parser, "infile", "The input file: an internal JSON format.");
    args::Positional<std::string> outfile(
        parser, "outfile", "The output file: same as input file, but with inferred parameter values.");
    args::ValueFlag<double> timeout(
        parser,
        "timeout",
        "Timeout in seconds. Solver will stop running and emit the best output after this amount of time.",
        {'t', "timeout"});
    args::Flag randomInit(
        parser, "randomInit", "Ignore the initial values in the input file and randomize them", {'r', "random-init"});
    args::ValueFlag<size_t> randomRestarts(
        parser,
        "randomRestart",
        "If the init value results in negLL=NaN, restart with a different init value this many times (deterministic).",
        {"random-restarts"});
    args::ValueFlag<size_t> selectMatrix(parser,
                                         "selectMatrix",
                                         "Select a specific coal matrix from the list of matrices in the input JSON.",
                                         {'m', "select-matrix"});
    args::Flag normCoals(
        parser, "normCoals", "Normalize coal matrix prior to computing likelihoods", {'n', "norm-coals"});
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help&) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    } catch (args::ValidationError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (!infile || !outfile) {
        std::cout << parser;
        return 1;
    }

    const bool logParams = true;
    NegLogLikelihoodCostFunctor cost(*infile, logParams, false, OBS_MODE_UNSPECIFIED, normCoals);
    if (selectMatrix) {
        const size_t selected = *selectMatrix;
        std::cerr << "Using only matrix " << selected << std::endl;
        cost.restrictToCMatrix(selected);
    }

    std::vector<double> paramVector(cost.m_schema.totalParams());
    if (randomInit) {
        setRandSeedToTime();
        cost.m_schema.randomParamVector(paramVector.data());
    } else {
        cost.m_schema.initParamVector(paramVector.data());
    }

    std::cerr << "Starting with X = ";
    cost.m_schema.dumpParameters(paramVector.data(), paramVector.size());
    std::cerr << std::endl;

    const double timeoutSec = timeout ? *timeout : DEFAULT_TIMEOUT_SECONDS;
    const size_t maxRestarts = randomRestarts ? *randomRestarts : DEFAULT_RESTARTS;
    double minf = std::numeric_limits<double>::quiet_NaN();
    size_t restarts = 0;
    do {
        minf = solve(cost, paramVector, MRP_OPT_ALG, true, timeoutSec);
        if (std::isnan(minf)) {
            // For our first restart, initialize the seed to the integer sum of all parameters.
            if (restarts == 0) {
                uint64_t seed = 0;
                for (const double v : paramVector) {
                    seed += (uint64_t)(v * 1001001001); // Multiply by a big number so we don't just get 0 for int
                }
                setRandSeed(seed);
            }
            cost.m_schema.randomParamVector(paramVector.data());
            restarts++;
            std::cerr << "Restarting solver due to NaN (" << restarts << "/" << maxRestarts << ")" << std::endl;
        }
    } while (std::isnan(minf) && restarts <= maxRestarts);

    // Don't let the solver silently return NaN to the Python code.
    user_exc_check(!std::isnan(minf),
                   "Likelihood is NaN after "
                       << restarts << " restarts; your model may be problematic (try restricting the bounds more)");

    if (*outfile != std::string("-")) {
        std::ofstream outputData(*outfile);
        std::cerr << "Writing results to " << *outfile << std::endl;
        outputData << cost.m_schema.toJsonOutput(paramVector.data(), minf) << std::endl;
    } else {
        std::cout << cost.m_schema.toJsonOutput(paramVector.data(), minf) << std::endl;
    }
    return 0;
}
