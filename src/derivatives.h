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
#ifndef MRPAST_DERIVATIVES_H
#define MRPAST_DERIVATIVES_H

#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

#include <Eigen/Dense>

#include "common.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// When values get small, do one-sided derivatives.
#define ALLOW_ONE_SIDED 1

// Warn about too few coal matrices if less than this.
constexpr size_t NUM_COAL_MATRICES_WARN_THRESHOLD = 10;

// Default value for when we stop looking for a better solution.
constexpr double RIDDERS_DEFAULT_ERROR = 1e-16;
// The Ridders-based derivatives start at a step factor of 10% of the size of the values.
constexpr double RIDDERS_INIT_STEP_FACTOR = 0.1;

// If the input value x is less than this, only perform the one-sided (upward) derivative on f(x).
constexpr double MIN_TWO_SIDED_VAL = 1e-8;

enum DerivativeDirection {
    DERIV_DIR_UP = 0x1U,
    DERIV_DIR_DOWN = 0x2U,
    DERIV_DIR_BOTH = (DERIV_DIR_UP | DERIV_DIR_DOWN),
};

// Clamp a value to a given lower-bound.
inline double lclamp(double value, double lb) { return std::max(value, lb); }

template <typename CostFunc>
static inline double
costVary1Param(CostFunc& cost, std::vector<double>& params, const size_t index, const double delta) {
    const double originalValue = params.at(index);
    params.at(index) += delta;
    const double result = cost(params.data());
    params.at(index) = originalValue;
    return result;
}

template <typename CostFunc>
static inline double costVary2Param(CostFunc& cost,
                                    std::vector<double>& params,
                                    const size_t index1,
                                    const size_t index2,
                                    const double delta1,
                                    const double delta2) {
    const double originalValue1 = params.at(index1);
    const double originalValue2 = params.at(index2);
    params.at(index1) += delta1;
    params.at(index2) += delta2;
    const double result = cost(params.data());
    params.at(index1) = originalValue1;
    params.at(index2) = originalValue2;
    return result;
}

template <typename CostFunc>
inline double calc2ndOrder3pt(CostFunc& cost,
                              std::vector<double>& workingParams,
                              const std::vector<DerivativeDirection>& directions,
                              const size_t xIndex,
                              const size_t yIndex,
                              const double h1,
                              const double h2) {
    double negH1 = 0.0;
    double posH1 = 0.0;
    double negH2 = 0.0;
    double posH2 = 0.0;
    switch (directions[xIndex]) {
    case DERIV_DIR_DOWN:
        negH1 = -2 * h1;
        posH1 = 0;
        break;
    case DERIV_DIR_UP:
        negH1 = 0;
        posH1 = +2 * h1;
        break;
    case DERIV_DIR_BOTH:
        negH1 = -h1;
        posH1 = +h1;
        break;
    }
    switch (directions[yIndex]) {
    case DERIV_DIR_DOWN:
        negH2 = -2 * h2;
        posH2 = 0;
        break;
    case DERIV_DIR_UP:
        negH2 = 0;
        posH2 = +2 * h2;
        break;
    case DERIV_DIR_BOTH:
        negH2 = -h2;
        posH2 = +h2;
        break;
    }
    const double A = costVary2Param(cost, workingParams, xIndex, yIndex, posH1, posH2);
    const double B = costVary2Param(cost, workingParams, xIndex, yIndex, negH1, negH2);
    const double C = costVary2Param(cost, workingParams, xIndex, yIndex, posH1, negH2);
    const double D = costVary2Param(cost, workingParams, xIndex, yIndex, negH1, posH2);
    return (A + B - C - D) / (4 * h1 * h2);
}

#define CHEAP_DIFFERENTIATION 0

#if CHEAP_DIFFERENTIATION
constexpr double numericDiffRelStepSize = 0.01; // 1%

template <typename CostFunc> double cheapF_x(CostFunc& cost, std::vector<double>& workingParams, size_t index) {
    const double step = workingParams[index] * numericDiffRelStepSize;
    const double A = costVary1Param(cost, workingParams, index, +step);
    const double B = costVary1Param(cost, workingParams, index, -step);
    return (A - B) / (2 * step);
}

template <typename CostFunc> double cheapF_xx(CostFunc& cost, std::vector<double>& workingParams, size_t index) {
    const double step = workingParams[index] * numericDiffRelStepSize;
    const double orig = cost(params.data());
    const double A = costVary1Param(cost, workingParams, index, +step);
    const double B = costVary1Param(cost, workingParams, index, -step);
    return (A - 2 * orig + B) / (step * step);
}

template <typename CostFunc>
double cheapF_xy(CostFunc& cost, std::vector<double>& workingParams, size_t index1, size_t index2) {
    const double step1 = workingParams[index1] * numericDiffRelStepSize;
    const double step2 = workingParams[index2] * numericDiffRelStepSize;
    return calc2ndOrder3pt(cost, workingParams, index1, index2, step1, step2);
}
#endif

// Ridder's method for more accurate central differences. A good explanation is
// at http://ceres-solver.org/numerical_derivatives.html#ridders-method
template <typename CostFunc>
inline double riddersD1(CostFunc& cost,
                        std::vector<double>& workingParams,
                        const std::vector<DerivativeDirection>& directions,
                        const size_t xIndex,
                        const double h,
                        const double errThreshold = RIDDERS_DEFAULT_ERROR,
                        const bool secondOrder = false,
                        // const size_t yIndex = std::numeric_limits<size_t>::max(),
                        // const double h2 = 0,
                        double* finalStepSize = nullptr,
                        double* estimatedError = nullptr) {
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    constexpr size_t MAX_COLUMNS = 25; // Somewhat arbitrary.
    double bestValue = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> column1;
    std::vector<double> column2;
    double current = 0.0;
    size_t col = 1;
    size_t row = 1;
    size_t iterations = 0;
    const double riddersFactor = 2;
    const double valueAtX = cost(workingParams.data());
    double errEst = std::numeric_limits<double>::max();
    for (size_t col = 0; col < MAX_COLUMNS; col++) {
        const double denominator = std::pow(riddersFactor, col);
        const double hFactor = h / denominator;
        double left = nan;
        double right = nan;
        double center = nan;
        switch (directions[xIndex]) {
        case DERIV_DIR_DOWN:
            left = costVary1Param(cost, workingParams, xIndex, -2 * hFactor);
            center = secondOrder ? costVary1Param(cost, workingParams, xIndex, -hFactor) : nan;
            right = valueAtX;
            break;
        case DERIV_DIR_UP:
            left = valueAtX;
            center = secondOrder ? costVary1Param(cost, workingParams, xIndex, +hFactor) : nan;
            right = costVary1Param(cost, workingParams, xIndex, +2 * hFactor);
            break;
        case DERIV_DIR_BOTH:
            left = costVary1Param(cost, workingParams, xIndex, +hFactor);
            center = secondOrder ? valueAtX : nan;
            right = costVary1Param(cost, workingParams, xIndex, -hFactor);
            break;
        }
        if (secondOrder) {
            current = ((left + right) - (2 * center)) / (hFactor * hFactor);
        } else {
            current = (left - right) / (2 * hFactor);
        }

        column1 = std::move(column2);
        column2 = {current};

        for (size_t row = 1; row <= col; row++) {
            const double coeff4 = std::pow((riddersFactor * riddersFactor), row);
            current = ((coeff4 * column2.back()) - column1.at(row - 1)) / (coeff4 - 1);

            const double errorA = std::fabs(current - column2.back());
            const double errorB = std::fabs(current - column1.at(row - 1));
            const double currentErr = std::max(errorA, errorB);
            if (currentErr < errEst) {
                errEst = currentErr;
                bestValue = current;
                if (finalStepSize != nullptr) {
                    *finalStepSize = hFactor;
                }
                if (errEst <= errThreshold) {
                    break;
                }
            }
            column2.push_back(current);
        }
        // Numerical instability test (swiped this approach from Ceres-solver)
        if (col > 0 && (std::fabs(current - column1.back()) >= 2 * errEst)) {
            break;
        }

        if (errEst <= errThreshold) {
            break;
        }
    }

    if (estimatedError != nullptr) {
        *estimatedError = errEst;
    }
    return bestValue;
}

template <typename CostFunc>
inline double calc2ndOrder5pt(CostFunc& cost,
                              std::vector<double>& workingParams,
                              const std::vector<DerivativeDirection>& directions,
                              const size_t xIndex,
                              const size_t yIndex,
                              const double h1,
                              const double h2) {
    double negH1 = 0.0;
    double posH1 = 0.0;
    double negH2 = 0.0;
    double posH2 = 0.0;
    switch (directions[xIndex]) {
    case DERIV_DIR_DOWN:
        negH1 = -2 * h1;
        posH1 = 0;
        break;
    case DERIV_DIR_UP:
        negH1 = 0;
        posH1 = +2 * h1;
        break;
    case DERIV_DIR_BOTH:
        negH1 = -h1;
        posH1 = +h1;
        break;
    }
    switch (directions[yIndex]) {
    case DERIV_DIR_DOWN:
        negH2 = -2 * h2;
        posH2 = 0;
        break;
    case DERIV_DIR_UP:
        negH2 = 0;
        posH2 = +2 * h2;
        break;
    case DERIV_DIR_BOTH:
        negH2 = -h2;
        posH2 = +h2;
        break;
    }
    const double A = costVary2Param(cost, workingParams, xIndex, yIndex, posH1, posH2);
    const double B = costVary2Param(cost, workingParams, xIndex, yIndex, negH1, negH2);
    const double C = costVary2Param(cost, workingParams, xIndex, yIndex, posH1, negH2);
    const double D = costVary2Param(cost, workingParams, xIndex, yIndex, negH1, posH2);
    double deriv1 = (A + B - C - D) / (4 * h1 * h2);

    const double E = costVary2Param(cost, workingParams, xIndex, yIndex, posH1 / 2, posH2 / 2);
    const double F = costVary2Param(cost, workingParams, xIndex, yIndex, negH1 / 2, negH2 / 2);
    const double G = costVary2Param(cost, workingParams, xIndex, yIndex, posH1 / 2, negH2 / 2);
    const double H = costVary2Param(cost, workingParams, xIndex, yIndex, negH1 / 2, posH2 / 2);
    double deriv2 = (E + F - G - H) / (h1 * h2);

    return (4 * deriv2 - deriv1) / 3;
}

template <typename CostFunc>
inline double calc2ndOrder10pt(CostFunc& cost,
                               std::vector<double>& workingParams,
                               const std::vector<DerivativeDirection>& directions,
                               const size_t xIndex,
                               const size_t yIndex,
                               const double h1,
                               const double h2) {
    const double d1 = calc2ndOrder5pt(cost, workingParams, directions, xIndex, yIndex, h1, h2);
    const double d2 = calc2ndOrder5pt(cost, workingParams, directions, xIndex, yIndex, h1 / 2, h2 / 2);
    return (4 * d2 - d1) / 3;
}

template <typename CostFunc>
inline double calc2ndOrder20pt(CostFunc& cost,
                               std::vector<double>& workingParams,
                               const std::vector<DerivativeDirection>& directions,
                               const size_t xIndex,
                               const size_t yIndex,
                               const double h1,
                               const double h2) {
    const double d1 = calc2ndOrder10pt(cost, workingParams, directions, xIndex, yIndex, h1, h2);
    const double d2 = calc2ndOrder10pt(cost, workingParams, directions, xIndex, yIndex, h1 / 2, h2 / 2);
    return (4 * d2 - d1) / 3;
}

// Our MLE solution is based on a bounded constraint solution, and we require
// that the lower bound be not too close to 0. Hence (unlike δaδi) we always use
// the central difference approximation.
template <typename CostFunc>
inline MatrixXd calculateSensitivity(CostFunc& cost, const std::vector<double>& optParams) {
    // Make a copy of the parameters that we can change for each component. Should
    // be much faster than allocating/copying every time.
    std::vector<double> paramsWorking = optParams;

    // We are computing an approximation to the expectation of the second order
    // gradeints, for which we just use the optimal parameters (theta-hat) and the
    // "true" observed data (coalescence matrix).
    cost.resetCMatrixCache();

    std::vector<DerivativeDirection> directions;

    const size_t N = optParams.size();
    MatrixXd H = MatrixXd::Zero(N, N);
    std::vector<double> stepSizes(N);
    for (Eigen::Index i = 0; i < N; i++) {
        const double origI = paramsWorking.at(i);
#if ALLOW_ONE_SIDED
        if (origI <= MIN_TWO_SIDED_VAL) {
            directions.push_back(DERIV_DIR_UP);
        } else {
#endif
            directions.push_back(DERIV_DIR_BOTH);
#if ALLOW_ONE_SIDED
        }
#endif
#if CHEAP_DIFFERENTIATION
        const double f_xx = cheapF_xx(cost, paramsWorking, i);
#else
        double stepSize = RIDDERS_INIT_STEP_FACTOR * origI;
        const double f_xx =
            riddersD1(cost, paramsWorking, directions, i, stepSize, RIDDERS_DEFAULT_ERROR, true, &stepSize);
        stepSizes[i] = stepSize;
#endif
        H(i, i) = -f_xx;
    }

    for (Eigen::Index i = 0; i < N; i++) {
        for (Eigen::Index j = i + 1; j < N; j++) {
#if CHEAP_DIFFERENTIATION
            H(j, i) = H(i, j) = -cheapF_xy(cost, paramsWorking, i, j);
#else
            H(j, i) = H(i, j) = -calc2ndOrder20pt(cost, paramsWorking, directions, i, j, stepSizes[i], stepSizes[j]);
#endif
        }
    }
    return std::move(H);
}

// The "score function" of the composite likelihood.
template <typename CostFunc>
inline MatrixXd
calculateVariability(CostFunc& cost, const std::vector<double>& optParams, const bool bootstrap = true) {
    // To calculate the expectation of the jacobian we need multiple inputs, hence
    // we iterate over all provided coalescence matrices. These matrices could
    // come from sampling via the ARG inference method (taking many MCMC samples,
    // sufficiently distanced from each other) or they could come from
    // bootstrapping the data in some way. The main thing that matters is that we
    // have enough of them and that they are all coming from f(y) which we can
    // think of as being "possible ARGs that explain the genetic data".

    // Make a copy of the parameters that we can change for each component. Should
    // be much faster than allocating/copying every time.
    std::vector<double> paramsWorking = optParams;

    size_t numSamples = 0;
    // If we are bootstrapping, we will use each individual CMatrix that was
    // provided to us, regardless of what the MLE configuration was (average, only
    // use 1 matrix, etc).
    if (bootstrap) {
        numSamples = cost.numCoalMatrices();
        // Otherwise we reset to whatever the configuration was, and treat it as a
        // single sample.
    } else {
        numSamples = 1;
        cost.resetCMatrixCache();
    }
    const size_t N = optParams.size();

    std::vector<DerivativeDirection> directions;
    for (Eigen::Index j = 0; j < N; j++) {
        const double value = optParams[j];
#if ALLOW_ONE_SIDED
        if (value <= MIN_TWO_SIDED_VAL) {
            directions.push_back(DERIV_DIR_UP);
        } else {
#endif
            directions.push_back(DERIV_DIR_BOTH);
#if ALLOW_ONE_SIDED
        }
#endif
    }

    MatrixXd expectedJ = MatrixXd::Zero(N, N);
    VectorXd expectedGrad = VectorXd::Zero(N);
    for (size_t sample = 0; sample < numSamples; sample++) {
        if (bootstrap) {
            // Only use the ith coal matrix
            cost.restrictToCMatrix(sample);
        }

        // f(x + h) − f(x − h) / 2h
        VectorXd gradient(N);
        for (Eigen::Index j = 0; j < N; j++) {
            const double value = optParams[j];
#if CHEAP_DIFFERENTIATION
            gradient[j] = cheapF_x(cost, paramsWorking, j);
#else
            const double step = lclamp(MIN_FLOAT_VAL, RIDDERS_INIT_STEP_FACTOR * value);
            gradient[j] = riddersD1(cost, paramsWorking, directions, j, step);
#endif
        }
        expectedJ += (gradient * gradient.transpose());
    }
    expectedJ /= (double)numSamples;

    return std::move(expectedJ);
}

#endif /* MRPAST_DERIVATIVES_H */