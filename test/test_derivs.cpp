#include <gtest/gtest.h>

// We use autodiff to check the quality of our numerical differentiation implementation.
// autodiff produces numerical evaluations of _analytic_ derivatives, via Eigen and the
// std:: math library. In contrast, our implementation produces numerical approximations
// of a derivative at a point -- autodiff should always be more accurate for the set of
// primitives that it supports.
#include <autodiff/reverse/var.hpp>

#include <string>
#include <vector>

#include "derivatives.h"

// Set to 1 to emit CSV data for plotting derivative results.
#define EMIT_CSV_DATA 1

inline void ENSURE_UNCHANGED(
        const std::vector<double>& origParams,
        const std::vector<double>& finalParams) {
    RELEASE_ASSERT(origParams.size() == finalParams.size());
    for (size_t i = 0; i < origParams.size(); i++) {
        RELEASE_ASSERT(origParams[i] == finalParams[i]);
    }
}

template <typename CostFunc>
double oldF_x(CostFunc& cost,
              const std::vector<double>& params,
              size_t index) {
    std::vector<double> workingParams = params;
    const double step = params[index] * 0.01; // 1%
    const double A = costVary1Param(cost, workingParams, index, +step);
    const double B = costVary1Param(cost, workingParams, index, -step);
    ENSURE_UNCHANGED(params, workingParams);
    return (A - B) / (2*step);
}

template <typename CostFunc>
double oldF_xx(CostFunc& cost,
               const std::vector<double>& params,
               size_t index) {
    std::vector<double> workingParams = params;
    const double step = params[index] * 0.01; // 1%
    const double orig = cost(params.data());
    const double A = costVary1Param(cost, workingParams, index, +step);
    const double B = costVary1Param(cost, workingParams, index, -step);
    ENSURE_UNCHANGED(params, workingParams);
    return (A - 2*orig + B) / (step*step);
}

template <typename CostFunc>
double oldF_xy(CostFunc& cost,
               const std::vector<double>& params,
               size_t index1,
               size_t index2) {
    std::vector<DerivativeDirection> directions(params.size(), DERIV_DIR_BOTH);
    std::vector<double> workingParams = params;
    const double step1 = params[index1] * 0.01; // 1%
    const double step2 = params[index2] * 0.01; // 1%
    auto result = calc2ndOrder3pt(cost, workingParams, directions, index1, index2, step1, step2);
    ENSURE_UNCHANGED(params, workingParams);
    return result;
}

template <typename CostFunc>
double newF_x(CostFunc& cost,
              const std::vector<double>& params,
              size_t index) {
    std::vector<DerivativeDirection> directions(params.size(), DERIV_DIR_BOTH);
    std::vector<double> workingParams = params;
    auto result = riddersD1(cost, workingParams, directions, index, 0.5*params[index]);
    ENSURE_UNCHANGED(params, workingParams);
    return result;
}

template <typename CostFunc>
double newF_xx(CostFunc& cost,
               const std::vector<double>& params,
               size_t index) {
    std::vector<DerivativeDirection> directions(params.size(), DERIV_DIR_BOTH);
    std::vector<double> workingParams = params;
    auto result = riddersD1(cost, workingParams, directions, index, 0.5*params[index], RIDDERS_DEFAULT_ERROR, true);
    ENSURE_UNCHANGED(params, workingParams);
    return result;
}

template <typename CostFunc>
double newF_xy(CostFunc& cost,
               const std::vector<double>& params,
               size_t index1,
               size_t index2) {
    std::vector<DerivativeDirection> directions(params.size(), DERIV_DIR_BOTH);
    std::vector<double> workingParams = params;
    double step1 = 0.5*params[index1];
    const double bestValue1 = riddersD1(cost, workingParams, directions, index1, step1,
                                        RIDDERS_DEFAULT_ERROR, true, &step1);
    double step2 = 0.5*params[index2];
    const double bestValue2 = riddersD1(cost, workingParams, directions, index2, step2,
                                        RIDDERS_DEFAULT_ERROR, true, &step2);

    auto result = calc2ndOrder10pt(cost, workingParams, directions, index1, index2, step1, step2);
    ENSURE_UNCHANGED(params, workingParams);
    return result;
}

template <typename T>
void emitCsvLine(std::vector<T> items) {
#if EMIT_CSV_DATA
    for (size_t i = 0; i < items.size(); i++) {
        if (i != 0) {
            std::cout << ",";
        }
        std::cout << items[i];
    }
    std::cout << std::endl;
#endif
}

inline std::string TS(double value) {
    std::stringstream ss;
    ss << std::setprecision(100) << value;
    return std::move(ss.str());
}

TEST(FirstOrder, ThreeD_Test1) {
    std::string experiment = "x*log(x*y)*exp(z)";
    const double yVal = 0.5;
    const double zVal = 2.0;

    emitCsvLine<std::string>({"Experiment", "Derivative", "x", "y", "z", "NewD", "OldD"});

    auto costFunc = [&](const double* parameters) {
        const double x = parameters[0];
        const double y = parameters[1];
        const double z = parameters[2];
        return x * std::log(x*y) * std::exp(z);
    };

    for (double xVal = 0.05; xVal <= 1.0; xVal += 0.05) {
        autodiff::var xAD = xVal;
        autodiff::var yAD = yVal;
        autodiff::var zAD = zVal;
        autodiff::var uAD =
            xAD * autodiff::reverse::detail::log(xAD*yAD) * autodiff::reverse::detail::exp(zAD);
        // First and second order analytic partial derivatives.
        auto [ux, uy, uz] = derivativesx(uAD, wrt(xAD, yAD, zAD));
        auto [uxx, uxy, uxz] = derivativesx(ux, wrt(xAD, yAD, zAD));

        // First order gradient w.r.t X
        std::vector<double> params = {xVal, yVal, zVal};
        const double oldDx = oldF_x(costFunc, params, 0);
        double oldErr = std::fabs((double)ux - oldDx);
        const double newDx = newF_x(costFunc, params, 0);
        double newErr = std::fabs((double)ux - newDx);
        emitCsvLine<std::string>({experiment, "f_x", TS(xVal), TS(yVal), TS(zVal), TS(newDx), TS(oldDx)});

        const double oldDxx = oldF_xx(costFunc, params, 0);
        oldErr = std::fabs((double)uxx - oldDxx);
        const double newDxx = newF_xx(costFunc, params, 0);
        newErr = std::fabs((double)uxx - newDxx);
        emitCsvLine<std::string>({experiment, "f_xx", TS(xVal), TS(yVal), TS(zVal), TS(newDxx), TS(oldDxx)});
   }

    for (double yVal = 0.05; yVal <= 1.0; yVal += 0.05) {
        const double xVal = 0.5;
        autodiff::var xAD = xVal;
        autodiff::var yAD = yVal;
        autodiff::var zAD = zVal;
        autodiff::var uAD =
            xAD * autodiff::reverse::detail::log(xAD*yAD) * autodiff::reverse::detail::exp(zAD);
        // First and second order analytic partial derivatives.
        auto [ux, uy, uz] = derivativesx(uAD, wrt(xAD, yAD, zAD));
        auto [uxx, uxy, uxz] = derivativesx(ux, wrt(xAD, yAD, zAD));

        std::vector<double> params = {xVal, yVal, zVal};
        const double oldDxy = oldF_xy(costFunc, params, 0, 1);
        double oldErr = std::fabs((double)uxy - oldDxy);
        const double newDxy = newF_xy(costFunc, params, 0, 1);
        double newErr = std::fabs((double)uxy - newDxy);
        emitCsvLine<std::string>({experiment, "f_xy", TS(xVal), TS(yVal), TS(zVal), TS(newDxy), TS(oldDxy)});
   }
}

TEST(FirstOrder, TwoD_Test1) {
    std::string experiment = "log(x^2)*sin(y)";
    const double yVal = 0.25;

    std::cout << std::setprecision(100);
    emitCsvLine<std::string>({"Experiment", "Derivative", "x", "y", "z", "NewD", "OldD"});

    auto costFunc = [&](const double* parameters) {
        const double x = parameters[0];
        const double y = parameters[1];
        return std::log(x*x) * std::sin(y);
    };

    for (double xVal = 0.05; xVal <= 1.0; xVal += 0.05) {
        autodiff::var xAD = xVal;
        autodiff::var yAD = yVal;
        autodiff::var uAD =
            autodiff::reverse::detail::log(xAD*xAD) * autodiff::reverse::detail::sin(yAD);
        // First and second order analytic partial derivatives.
        auto [ux, uy] = derivativesx(uAD, wrt(xAD, yAD));
        auto [uxx, uxy] = derivativesx(ux, wrt(xAD, yAD));

        // First order gradient w.r.t X
        std::vector<double> params = {xVal, yVal};
        const double oldDx = oldF_x(costFunc, params, 0);
        double oldErr = std::fabs((double)ux - oldDx);
        const double newDx = newF_x(costFunc, params, 0);
        double newErr = std::fabs((double)ux - newDx);
        emitCsvLine<std::string>({experiment, "f_x", TS(xVal), TS(yVal), TS(0.0), TS(newDx), TS(oldDx)});

        const double oldDxx = oldF_xx(costFunc, params, 0);
        oldErr = std::fabs((double)uxx - oldDxx);
        const double newDxx = newF_xx(costFunc, params, 0);
        newErr = std::fabs((double)uxx - newDxx);
        emitCsvLine<std::string>({experiment, "f_xx", TS(xVal), TS(yVal), TS(0.0), TS(newDxx), TS(oldDxx)});
   }

    for (double yVal = 0.05; yVal <= 1.0; yVal += 0.05) {
        const double xVal = 0.5;
        autodiff::var xAD = xVal;
        autodiff::var yAD = yVal;
        autodiff::var uAD =
            autodiff::reverse::detail::log(xAD*xAD) * autodiff::reverse::detail::sin(yAD);
        // First and second order analytic partial derivatives.
        auto [ux, uy] = derivativesx(uAD, wrt(xAD, yAD));
        auto [uxx, uxy] = derivativesx(ux, wrt(xAD, yAD));

        std::vector<double> params = {xVal, yVal};
        const double oldDxy = oldF_xy(costFunc, params, 0, 1);
        double oldErr = std::fabs((double)uxy - oldDxy);
        const double newDxy = newF_xy(costFunc, params, 0, 1);
        double newErr = std::fabs((double)uxy - newDxy);
        emitCsvLine<std::string>({experiment, "f_xy", TS(xVal), TS(yVal), TS(0.0), TS(newDxy), TS(oldDxy)});
   }
}

TEST(Helpers, Vary1Param) {
    auto costFunc = [&](const double* parameters) {
        const double x = parameters[0];
        const double y = parameters[1];
        return std::log(x*x) * std::sin(y);
    };

	// Just verify that this function properly restores the original values for params.
    std::vector<double> parameters = {1, 2};
    costVary1Param(costFunc, parameters, 0, 5);
    ASSERT_EQ(parameters[0], 1);
    ASSERT_EQ(parameters[1], 2);
    costVary1Param(costFunc, parameters, 1, 10);
    ASSERT_EQ(parameters[0], 1);
    ASSERT_EQ(parameters[1], 2);
}

TEST(Helpers, Vary2Param) {
    auto costFunc = [&](const double* parameters) {
        const double x = parameters[0];
        const double y = parameters[1];
        const double z = parameters[2];
        return std::log(x*x) * std::sin(y) + z;
    };

	// Just verify that this function properly restores the original values for params.
    std::vector<double> parameters = {1, 2, 3};
    costVary2Param(costFunc, parameters, 0, 1, 10, 10);
    ASSERT_EQ(parameters[0], 1);
    ASSERT_EQ(parameters[1], 2);
    ASSERT_EQ(parameters[2], 3);
    costVary2Param(costFunc, parameters, 1, 2, 10, 10);
    ASSERT_EQ(parameters[0], 1);
    ASSERT_EQ(parameters[1], 2);
    ASSERT_EQ(parameters[2], 3);
}