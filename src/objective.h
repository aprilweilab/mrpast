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
#ifndef MRP_SOLVER_OBJECTIVE_H
#define MRP_SOLVER_OBJECTIVE_H

#include <cmath>
#include <iostream>
#include <list>
#include <string>
#include <vector>

#include "json.hpp"

#include "common.h"

using json = nlohmann::json;

// We rescale epoch times to "millions of generations". This is just for improved matrix
// operations on the Hessian, since all our rates are between 0.0 and 1.0, and having a
// hugely different scale produces problems with the inverse of the Hessian.
static constexpr double EPOCH_TIME_RESCALE = 1e-6;

// How we should treat the observed data (coal matrices).
enum ObservationMode {
    OBS_MODE_UNSPECIFIED = 0,
    OBS_MODE_AVERAGE = 1, // Average each C[i, j] across the coal matrices
    OBS_MODE_FIRST = 2,   // Only evaluate using the first matrix in the list
};

static inline ObservationMode parse_obs_mode(const std::string& obsModeStr) {
    if (obsModeStr == "average") {
        return OBS_MODE_AVERAGE;
    }
    if (obsModeStr == "first") {
        return OBS_MODE_FIRST;
    }
    std::cerr << "Invalid " << OBSERVATION_MODE_KEY << " value: " << obsModeStr << std::endl;
    abort();
}

enum Adjustment {
    ADJUST_NONE = 0,
    ADJUST_GROWTH_RATE = 1,
    ADJUST_INV_GROWTH_RATE = 2,
};

static inline Adjustment parse_adjust(const std::string& adjustStr) {
    if (adjustStr == "growth_rate") {
        return ADJUST_GROWTH_RATE;
    }
    if (adjustStr == "inverse_growth_rate") {
        return ADJUST_INV_GROWTH_RATE;
    }
    abort();
}

// An application of a variable to a matrix. This is interpreted as:
//  M[i, j] += (coeff * variable)
// In this way variables can be reused across a matrix when there is dependence
// between cells.
struct VariableApplication {
    double coeff;
    size_t i;
    size_t j;
    size_t epoch;
    Adjustment adjust;
};

// An optimizable parameter to the model.
struct BoundedVariable {
    double init;
    double lb;
    double ub;
    std::list<VariableApplication> applications;
    std::string description;
    double ground_truth;
};

/**
 * The schema describes the parameters, their bounds, and how they apply to the
 * transition matrix of the Markov model.
 */
class ParameterSchema {
public:
    using VarList = std::vector<BoundedVariable>;

    explicit ParameterSchema(const bool logSpace = false)
        : m_logSpace(logSpace) {}

    /**
     * Are the parameters being manipulated in log10-space?
     */
    inline bool isLogSpace() const { return m_logSpace; }

    /**
     * Convert parameter values to the parameter space which is seen by the
     * solver. The non-solver view is always normal space, but solver space is
     * sometimes log10.
     */
    inline double toParam(double value, const size_t index) const {
        assert(index < m_paramRescale.size());
        value *= m_paramRescale[index];
        if (m_logSpace) {
            return std::log10(value);
        }
        return value;
    }

    /**
     * Convert parameter values from the parameter space which is seen by the
     * solver. The non-solver view is always normal space, but solver space is
     * sometimes log10.
     */
    inline double fromParam(const double value, const size_t index) const {
        double rv = value;
        if (m_logSpace) {
            rv = std::pow(10.0, rv);
        }
        rv /= m_paramRescale.at(index);
        RELEASE_ASSERT(!std::isnan(rv));
        return rv;
    }

    /**
     * See toParam(). This just does multiple parameters at once.
     */
    void toParams(std::vector<double>& parameters) const {
        for (size_t i = 0; i < parameters.size(); i++) {
            parameters[i] = toParam(parameters[i], i);
        }
    }

    /**
     * See fromParam(). This just does multiple parameters at once.
     */
    void fromParams(std::vector<double>& parameters) const {
        for (size_t i = 0; i < parameters.size(); i++) {
            parameters[i] = fromParam(parameters[i], i);
        }
    }

    inline void dumpParameters(const double* parameters, const size_t numParams) const {
        std::cerr << "[";
        for (size_t i = 0; i < numParams; i++) {
            if (i != 0) {
                std::cerr << ", ";
            }
            std::cerr << fromParam(parameters[i], i);
        }
        std::cerr << "]";
    }

    /**
     * Construct from the JSON data sent to the solver.
     */
    void load(const json& inputData);

    size_t numEpochs() const { return m_numEpochs; }

    size_t totalParams() const { return m_eParams.size() + m_sParams.size(); }

    size_t numStates(size_t epoch) const { return m_sStates.at(epoch); }

    /**
     * Populate the parameters[] vector with random values based on the
     * lower/upper bounds from the parameter schema.
     */
    void randomParamVector(double* parameters) const;

    /**
     * Populate the parameters[] vector with the initial value provided in the
     * parameter schema.
     */
    void initParamVector(double* parameters) const;

    void getBounds(size_t paramIdx, double& lowerBound, double& upperBound) const;

    double getEpochStartTime(double const* parameters, size_t epoch) const;

    template <typename T> void addToParamOutput(json& outputData, const std::string& field, const std::vector<T>& data);

    json toJsonOutput(const double* parameters, double negLL) const;

    void fromJsonOutput(const json& jsonOutput, double* parameters, std::string key = "final") const;

    bool samplesAreJackknifed() const { return m_inputJson["sampling_description"] == "jackknife"; }

    // Parameters for epoch times
    VarList m_eParams;
    // Fixed values for epoch times.
    VarList m_eFixed;
    // Parameters for stochastic matrix (by epoch)
    VarList m_sParams;
    // Fixed values for stochastic matrix. These are just parameters that are
    // marked "fixed".
    VarList m_sFixed;

private:
    // We save the actual JSON object, because the output is identical to the
    // input except for
    // some additional fields. This way we can just copy this JSON object for the
    // output.
    json m_inputJson;
    // Total number of epochs
    size_t m_numEpochs{};
    // Number of states per stochastic matrix (by epoch)
    std::vector<size_t> m_sStates;

    // Rescaling of the parameter that should occur from the solver/derivative perspective.
    // This scaling is linear, as non-linear scaling will screw up the confidence intervals.
    std::vector<double> m_paramRescale;

    // Should the parameters be transformed into log space?
    bool m_logSpace;
};

template <typename T>
void ParameterSchema::addToParamOutput(json& outputData, const std::string& field, const std::vector<T>& data) {
    size_t p = 0;
    json& epochTimes = outputData[EPOCH_TIMES_KEY];
    for (size_t i = 0; i < m_eParams.size(); i++) {
        epochTimes[i][field] = data[p];
        p++;
    }
    json& sMatrixValues = outputData[SMATRIX_VALS_KEY];
    for (size_t i = 0; i < m_sParams.size(); i++) {
        sMatrixValues[i][field] = data[p];
        p++;
    }
}

class ObservedData;
class CachedCMatrix;

/**
 * The negative log-likelihood function to be minimized.
 */
class NegLogLikelihoodCostFunctor {
public:
    /**
     * The negative log-likelihood function to be minimized.
     * @param[in] jsonFile The file containing the concrete model JSON.
     * @param[in] logParam If true, the solver will see a log10 view of the
     * parameters while the actual computations will be done in the normal scale.
     * This helps with numerical stability of the solver as the log10 scale is
     * what the derivatives will be seen through (helps with derivative-free
     * solvers also).
     * @param[in] maximization By default this is negative log-likelihood, but
     * setting this to true will make it just log-likelihood.
     * @param[in] mode Deprecated. Do not use.
     * @param[in] normalize Set to true to normalize the coal matrix so that each
     *  state sums to 1.
     */
    explicit NegLogLikelihoodCostFunctor(const std::string& jsonFile,
                                         bool logParams = false,
                                         bool maximization = false,
                                         ObservationMode mode = OBS_MODE_UNSPECIFIED,
                                         bool normalize = false);
    virtual ~NegLogLikelihoodCostFunctor();

    // This class is not safe for copying because of the PIMPL idiom for the
    // observed data and the matrix cache. Deleting all of these prevents misuse.
    NegLogLikelihoodCostFunctor& operator=(const NegLogLikelihoodCostFunctor& rhs) = delete;
    NegLogLikelihoodCostFunctor& operator=(NegLogLikelihoodCostFunctor&& rhs) = delete;
    NegLogLikelihoodCostFunctor(const NegLogLikelihoodCostFunctor& rhs) = delete;
    NegLogLikelihoodCostFunctor(const NegLogLikelihoodCostFunctor&& rhs) = delete;

    double operator()(double const* parameters) const;

    /**
     * How many coalescent count matrices were loaded from the JSON file.
     */
    size_t numCoalMatrices() const;

    /**
     * Only use the coalescence count matrix with this index (0 <= index <
     * numCoalMatrices()).
     */
    void restrictToCMatrix(size_t index);

    /**
     * Use whatever count matrix(es) that the schema requests.
     */
    void resetCMatrixCache();

    const std::vector<double>& getTimeSlices() const;

    /**
     * Emit the model JSON, replacing the input CoalMatrix (presumably derived from ARGs)
     * with the theoretical CoalMatrix (derived from the given parameters and the Markov model)
     */
    void outputTheoreticalCoalMatrix(double const* parameters, std::ostream& out) const;

    ParameterSchema m_schema;
    ObservedData* m_observed; // PIMPL
    CachedCMatrix* m_cache;   // PIMPL

    void resetSolveStats() {
        m_numObjCalls = 1;
        m_minValue = std::numeric_limits<double>::max();
        m_callsAtMinValue = 0;
    }

    size_t m_numObjCalls;
    double m_minValue;
    size_t m_callsAtMinValue;

private:
    const bool m_maximization;
    ObservationMode m_mode;
};

#endif /* MRP_SOLVER_OBJECTIVE_H */