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
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include "Eigen/Dense"
#include "json.hpp"
#include "unsupported/Eigen/MatrixFunctions"

#include "common.h"
#include "objective.h"

using json = nlohmann::json;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Turn this on and recompile, to track down numerical problems.
#define TRACE_MATRICES 0

// When computing probabilities in multi-epoch models, this uses the additive probability rule
// P(A or B) = P(A) + P(B) - P(A and B), instead of negating all of the probabilities to only
// use P(not(A) and not(B)). The former results in fewer subtractions, which anecdotally appears
// to produce slightly better point estimates with definitely better confidence intervals.
#define FEWER_SUBTRACTIONS_PER_EPOCH 1

#define DUMP_MATRIX(m, desc)                                                                                           \
    do {                                                                                                               \
        std::cerr << "% " << desc << ": " << std::endl;                                                                \
        std::cerr << "[";                                                                                              \
        for (Eigen::Index i = 0; i < (m).rows(); i++) {                                                                \
            if (i > 0) {                                                                                               \
                std::cerr << "; ";                                                                                     \
            }                                                                                                          \
            for (Eigen::Index j = 0; j < (m).cols(); j++) {                                                            \
                if (j > 0) {                                                                                           \
                    std::cerr << ", ";                                                                                 \
                }                                                                                                      \
                std::cerr << (m).coeff(i, j);                                                                          \
            }                                                                                                          \
        }                                                                                                              \
        std::cerr << "]" << std::endl;                                                                                 \
    } while (0)

#if TRACE_MATRICES
#define TRACE_MATRIX(m, desc) DUMP_MATRIX
#define TRACELN(msg)          std::cerr << msg << std::endl;
#else
#define TRACE_MATRIX(m, desc)
#define TRACELN(msg)
#endif

double randDouble(const double lower, const double upper) {
    static unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    static std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> sampler(lower, upper);
    return sampler(generator);
}

inline bool nearlyEqual(const double val1, const double val2, const double epsilon) {
    return ((val1 + epsilon >= val2) && (val1 - epsilon <= val2));
}

MatrixXd loadJsonMatrix(json& inputJson, const std::string& key) {
    const auto& inputMatrix = inputJson.at(key);
    const Eigen::Index nRows = SIZE_T_TO_INDEX(inputMatrix.size());
    if (nRows == 0) {
        return {};
    }
    const Eigen::Index nCols = SIZE_T_TO_INDEX(inputMatrix[0].size());
    // We have counts for all the states that involve demeA x demeB.
    MatrixXd result = MatrixXd::Zero(nRows, nCols);
    for (size_t i = 0; i < nRows; i++) {
        for (size_t j = 0; j < nCols; j++) {
            result(i, j) = inputMatrix.at(i).at(j);
        }
    }
    return std::move(result);
}

std::vector<MatrixXd> loadCMatrices(json& inputListOfMatrices) {
    std::vector<MatrixXd> result;
    for (const auto& inputMatrix : inputListOfMatrices) {
        const Eigen::Index nRows = SIZE_T_TO_INDEX(inputMatrix.size());
        if (nRows == 0) {
            result.emplace_back();
            continue;
        }
        const Eigen::Index nCols = SIZE_T_TO_INDEX(inputMatrix[0].size());
        // We have counts for all the states that involve demeA x demeB.
        MatrixXd curMatrix = MatrixXd::Zero(nRows, nCols);
        for (Eigen::Index i = 0; i < nRows; i++) {
            for (Eigen::Index j = 0; j < nCols; j++) {
                curMatrix(i, j) = inputMatrix.at(i).at(j);
            }
        }
        result.push_back(std::move(curMatrix));
    }
    return std::move(result);
}

/**
 * Both input and output values are in generations.
 */
std::vector<double> getTimeSlices(json& inputJson) {
    std::vector<double> result;
    for (const auto& slice : inputJson.at(TIME_SLICES_KEY)) {
        result.push_back((double)slice);
    }
    return std::move(result);
}

#ifndef NDEBUG
void verifyQMatrix(const MatrixXd& sMatrix) {
    const double nearlyZero = 0.000001;
    for (Eigen::Index i = 0; i < sMatrix.rows(); i++) {
        MODEL_ASSERT_MSG(sMatrix.row(i).sum() <= nearlyZero, sMatrix);
        for (Eigen::Index j = 0; j < sMatrix.cols(); j++) {
            MODEL_ASSERT(i == j || sMatrix(i, j) >= 0);
        }
    }
}
#endif

double getGroundTruth(const json& parameter) {
    double ground_truth = std::numeric_limits<float>::quiet_NaN();
    if (parameter.contains("ground_truth") && !parameter["ground_truth"].is_null()) {
        ground_truth = parameter["ground_truth"];
    }
    return ground_truth;
}

void loadParamList(const json& parameterList,
                   ParameterSchema::VarList& allParams,
                   std::vector<size_t>& freeParamIdx,
                   std::vector<size_t>& fixedParamIdx,
                   std::vector<size_t>& oneMinusIdx,
                   std::vector<double>& rescaling) {
    for (const auto& parameter : parameterList) {
        std::vector<size_t> oneMinus;
        if (parameter.contains("one_minus")) {
            for (const auto& om : parameter["one_minus"]) {
                oneMinus.push_back((size_t)om);
            }
        }
        const double ground_truth = getGroundTruth(parameter);
        if (isJsonParamFixed(parameter)) {
            fixedParamIdx.push_back(allParams.size());
        } else if (!oneMinus.empty()) {
            oneMinusIdx.push_back(allParams.size());
        } else {
            freeParamIdx.push_back(allParams.size());
            rescaling.emplace_back(1.0);
        }
        BoundedVariable bv = {parameter["init"],
                              parameter["lb"],
                              parameter["ub"],
                              {},
                              parameter["description"],
                              ground_truth,
                              parameter["kind_index"],
                              std::move(oneMinus)};
        allParams.push_back(std::move(bv));
    }
}

void ParameterSchema::load(const json& inputData) {
    m_paramRescale.resize(0);
    m_inputJson = inputData;
    std::vector<size_t> ignored;
    // EPOCH PARAMETERS
    loadParamList(inputData.at(EPOCH_TIMES_KEY), m_eParams, m_eParamIdx, m_eFixedIdx, ignored, m_paramRescale);
    m_numEpochs = m_eParams.size() + 1;
    m_sStates = std::vector<size_t>(m_numEpochs);
    // STOCHASTIC MATRIX (Q-Matrix) PARAMETERS
    for (const auto& parameter : inputData.at(SMATRIX_VALS_KEY)) {
        std::vector<VariableApplication> applications;
        const auto& applyTo = parameter["apply_to"];
        assert(!applyTo.empty());
        for (const auto& app : applyTo) {
            const size_t i = app["i"];
            const size_t j = app["j"];
            const size_t epoch = app["epoch"];
            if (j >= m_sStates.at(epoch)) {
                m_sStates[epoch] = j + 1;
            }
            auto adjustmentJson = app.value("adjustment", json());
            Adjustment adjust = ADJUST_NONE;
            if (!adjustmentJson.is_null()) {
                adjust = parse_adjust(adjustmentJson);
            }
            applications.push_back({app["coeff"], i, j, epoch, adjust});
        }
        const double ground_truth = getGroundTruth(parameter);
        BoundedVariable bv = {parameter["init"],
                              parameter["lb"],
                              parameter["ub"],
                              std::move(applications),
                              parameter["description"],
                              ground_truth};
        if (isJsonParamFixed(parameter)) {
            m_sFixedIdx.push_back(m_sParams.size());
        } else {
            m_sParamIdx.push_back(m_sParams.size());
            m_paramRescale.emplace_back(1.0);
        }
        m_sParams.push_back(std::move(bv));
    }
    // ADMIXTURE STATE MATRIX PARAMETERS
    loadParamList(
        inputData.at(AMATRIX_PARAMS_KEY), m_aParams, m_aParamIdx, m_aFixedIdx, m_aOneMinusIdx, m_paramRescale);
    for (const auto& application : inputData.at(AMATRIX_APPS_KEY)) {
        m_admixtureApps.push_back({
            application["coeff"],
            application["i"],
            application["j"],
            application["epoch"],
            application["vars"].at(0),
            application["vars"].at(1),
        });
    }
    const size_t paramCount = (m_sParamIdx.size() + m_eParamIdx.size() + m_aParamIdx.size());
    RELEASE_ASSERT(paramCount == m_paramRescale.size());
}

/**
 * Populate the parameters[] vector with random values based on the lower/upper
 * bounds from the parameter schema.
 */
void ParameterSchema::randomParamVector(double* parameters) const {
    size_t p = 0;
    for (const size_t index : m_eParamIdx) {
        const auto& parameter = m_eParams.at(index);
        RELEASE_ASSERT(p < totalParams());
        parameters[p] = toParam(randDouble(parameter.lb, parameter.ub), p);
        p++;
    }
    for (const size_t index : m_sParamIdx) {
        const auto& parameter = m_sParams.at(index);
        RELEASE_ASSERT(p < totalParams());
        parameters[p] = toParam(randDouble(parameter.lb, parameter.ub), p);
        p++;
    }
    for (const size_t index : m_aParamIdx) {
        const auto& parameter = m_aParams.at(index);
        RELEASE_ASSERT(p < totalParams());
        parameters[p] = toParam(randDouble(parameter.lb, parameter.ub), p);
        p++;
    }
}

void ParameterSchema::initParamsViaList(double* parameters,
                                        size_t& i,
                                        const ParameterSchema::VarList& paramVars,
                                        const std::vector<size_t>& paramIdx) const {
    for (const size_t index : paramIdx) {
        const auto& parameter = paramVars.at(index);
        RELEASE_ASSERT(index < totalParams());
        parameters[i] = toParam(parameter.init, i);
        i++;
    }
}

/**
 * Populate the parameters[] vector with the initial value provided in the
 * parameter schema.
 */
void ParameterSchema::initParamVector(double* parameters) const {
    size_t p = 0;
    initParamsViaList(parameters, p, m_eParams, m_eParamIdx);
    initParamsViaList(parameters, p, m_sParams, m_sParamIdx);
    initParamsViaList(parameters, p, m_aParams, m_aParamIdx);
}

void ParameterSchema::getBounds(size_t paramIdx, double& lowerBound, double& upperBound) const {
    if (paramIdx >= paramStartAdmix()) {
        const size_t index = m_aParamIdx.at(paramIdx - paramStartAdmix());
        lowerBound = toParam(m_aParams.at(index).lb, paramIdx);
        upperBound = toParam(m_aParams.at(index).ub, paramIdx);
    } else if (paramIdx >= paramStartSMatrix()) {
        const size_t index = m_sParamIdx.at(paramIdx - paramStartSMatrix());
        lowerBound = toParam(m_sParams.at(index).lb, paramIdx);
        upperBound = toParam(m_sParams.at(index).ub, paramIdx);
    } else {
        const size_t index = m_eParamIdx.at(paramIdx);
        lowerBound = toParam(m_eParams.at(index).lb, paramIdx);
        upperBound = toParam(m_eParams.at(index).ub, paramIdx);
    }
}

std::vector<double> ParameterSchema::getEpochStartTimes(double const* parameters) const {
    std::vector<double> result(m_eParams.size());
    for (size_t i = 0; i < m_eParamIdx.size(); i++) {
        const size_t epochIndex = m_eParamIdx[i];
        result.at(epochIndex) = fromParam(parameters[i], i);
    }
    for (size_t epochIndex : m_eFixedIdx) {
        result.at(epochIndex) = m_eParams[epochIndex].init;
    }
    return std::move(result);
}

std::vector<double> getAdmixtureValues(const ParameterSchema& schema, double const* parameters) {
    // Copy all of the parameter and fixed values into the admixture value vector.
    const size_t firstParamIdx = schema.m_eParamIdx.size() + schema.m_sParamIdx.size();
    const size_t numAdmixtureValues = schema.m_aParams.size();
    std::vector<double> admixtureValues(numAdmixtureValues);
    // Copy parameters from the solver input
    size_t paramIdx = firstParamIdx;
    for (const size_t index : schema.m_aParamIdx) {
        RELEASE_ASSERT(paramIdx < schema.totalParams());
        admixtureValues.at(index) = schema.fromParam(parameters[paramIdx], paramIdx);
        paramIdx++;
    }
    // Populate the fixed values  (TODO: could prepopulate this during parsing)
    for (const size_t index : schema.m_aFixedIdx) {
        admixtureValues.at(index) = schema.m_aParams.at(index).init;
    }
    // Populate the "one-minus" values (parameters determined by other parameters)
    for (const size_t index : schema.m_aOneMinusIdx) {
        double sum = 0.0;
        for (const size_t otherIndex : schema.m_aParams.at(index).oneMinus) {
            sum += admixtureValues.at(otherIndex);
        }
        admixtureValues.at(index) = 1.0 - sum;
    }
    return std::move(admixtureValues);
}

json ParameterSchema::toJsonOutput(const double* parameters, const double negLL) const {
    json output = m_inputJson; // Make a copy of the input.
    size_t i = 0;
    size_t p = 0;
    RELEASE_ASSERT(output.at(EPOCH_TIMES_KEY).size() == m_eParams.size());
    for (auto& parameter : output.at(EPOCH_TIMES_KEY)) {
        if (isJsonParamFixed(parameter)) {
            parameter["final"] = m_eParams.at(i).init;
        } else {
            RELEASE_ASSERT(p < totalParams());
            parameter["final"] = fromParam(parameters[p], p);
            p++;
        }
        i++;
    }
    RELEASE_ASSERT(i == m_eParams.size());
    RELEASE_ASSERT(p == m_eParamIdx.size());

    RELEASE_ASSERT(output.at(SMATRIX_VALS_KEY).size() == m_sParams.size());
    i = 0;
    for (auto& parameter : output.at(SMATRIX_VALS_KEY)) {
        if (isJsonParamFixed(parameter)) {
            parameter["final"] = m_sParams.at(i).init;
        } else {
            RELEASE_ASSERT(p < totalParams());
            parameter["final"] = fromParam(parameters[p], p);
            p++;
        }
        i++;
    }
    RELEASE_ASSERT(i == m_sParams.size());
    RELEASE_ASSERT((p - m_eParamIdx.size()) == m_sParamIdx.size());

    RELEASE_ASSERT(output.at(AMATRIX_PARAMS_KEY).size() == m_aParams.size());
    i = 0;
    std::vector<double> admixtureValues = getAdmixtureValues(*this, parameters);
    for (auto& parameter : output.at(AMATRIX_PARAMS_KEY)) {
        parameter["final"] = admixtureValues.at(i);
        i++;
    }
    RELEASE_ASSERT(i == m_aParams.size());

    output["negLL"] = negLL;
    return std::move(output);
}

void ParameterSchema::fromJsonOutput(const json& jsonOutput, double* parameters, std::string key) const {
    RELEASE_ASSERT(false); // FIXME
    size_t p = 0;
    for (const auto& paramVal : jsonOutput[EPOCH_TIMES_KEY]) {
        if (!isJsonParamFixed(paramVal)) {
            parameters[p] = toParam((double)paramVal[key], p);
            p++;
        }
    }
    RELEASE_ASSERT(m_eParams.size() == p);
    for (const auto& paramVal : jsonOutput[SMATRIX_VALS_KEY]) {
        if (!isJsonParamFixed(paramVal)) {
            parameters[p] = toParam((double)paramVal[key], p);
            p++;
        }
    }
    RELEASE_ASSERT(m_sParams.size() + m_eParams.size() == p);
}

/**
 * Given the vector of concrete parameters values (parameters[]) and the current
 * epoch, produce the concrete infinitesimal rate matrix based on the mappings in this
 * schema.
 */
MatrixXd
createQMatrix(const ParameterSchema& schema, double const* parameters, size_t epoch, const double timeSinceEpochStart) {
    assert(epoch < schema.numEpochs());

    const Eigen::Index nStates = SIZE_T_TO_INDEX(schema.numStates(epoch));
    // Skip this many epoch transition times.
    const size_t firstParamIdx = schema.m_eParamIdx.size();
    MatrixXd qMatrix = MatrixXd::Zero(nStates, nStates);

    auto applyParameter =
        [&](const std::vector<VariableApplication>& applications, const double parameterValue, bool& anyAdjusted) {
            for (const auto& application : applications) {
                if (application.epoch == epoch) {
                    const Eigen::Index i = SIZE_T_TO_INDEX(application.i);
                    const Eigen::Index j = SIZE_T_TO_INDEX(application.j);
                    // Adjustments always come last. See model.py
                    if (application.adjust == ADJUST_GROWTH_RATE) {
                        const double origRate = qMatrix(i, j);
                        // The integral on interval {0, u} is (1/a - exp(-a*u)/a)
                        // we can divide by alpha to get the average rate over time period lower:upper.
                        const double integratedRate = (std::exp(parameterValue * timeSinceEpochStart) - 1) /
                                                      (parameterValue * timeSinceEpochStart);
                        qMatrix(i, j) = origRate * integratedRate * application.coeff;
                        anyAdjusted = true;
                    } else if (application.adjust == ADJUST_INV_GROWTH_RATE) {
                        // For backwards compatibility with JSON input files. Just ignore.
                    } else {
                        qMatrix(i, j) += parameterValue * application.coeff;
                        assert(!anyAdjusted);
                    }
                }
            }
        };

    // Apply the non-parameters first.
    bool anyAdjusted = false;
    for (const size_t index : schema.m_sFixedIdx) {
        const auto& param = schema.m_sParams.at(index);
        applyParameter(param.applications, param.init, anyAdjusted);
    }
    for (size_t i = 0; i < schema.m_sParamIdx.size(); i++) {
        const size_t paramIdx = firstParamIdx + i;
        RELEASE_ASSERT(paramIdx < schema.totalParams());
        double pVal = schema.fromParam(parameters[paramIdx], paramIdx);
        const size_t index = schema.m_sParamIdx[i];
        applyParameter(schema.m_sParams.at(index).applications, pVal, anyAdjusted);
    }
    // Sum the off-diagonal and set the diagonal to the negative sum. This makes a valid
    // Q-matrix.
    for (Eigen::Index i = 0; i < qMatrix.rows(); i++) {
        double rowSum = 0.0;
        for (Eigen::Index j = 0; j < qMatrix.cols(); j++) {
            if (i != j) {
                rowSum += qMatrix(i, j);
            }
        }
        qMatrix(i, i) = -rowSum;
    }
#ifndef NDEBUG
    verifyQMatrix(qMatrix);
#endif
    return std::move(qMatrix);
}

#ifndef NDEBUG
void verifyASMatrix(const MatrixXd& ASMatrix) {
    const double nearlyZero = 0.000001;
    for (Eigen::Index i = 0; i < ASMatrix.rows(); i++) {
        MODEL_ASSERT_MSG(std::abs(1.0 - ASMatrix.row(i).sum()) <= nearlyZero, ASMatrix);
        for (Eigen::Index j = 0; j < ASMatrix.cols(); j++) {
            MODEL_ASSERT(ASMatrix(i, j) >= 0);
        }
    }
}
#endif

/**
 * Given the vector of concrete parameters values (parameters[]) and the current
 * epoch, produce the concrete admixture proportion matrix.
 */
MatrixXd createASMatrix(const ParameterSchema& schema, double const* parameters, size_t epoch) {
    assert(epoch < schema.numEpochs());

    std::vector<double> admixtureValues = getAdmixtureValues(schema, parameters);

    // Now apply these values to construct the matrix, we leave off the coalescence state.
    const Eigen::Index nStates = SIZE_T_TO_INDEX(schema.numStates(epoch) - 1);
    MatrixXd ASMatrix = MatrixXd::Zero(nStates, nStates);
    for (const auto& app : schema.m_admixtureApps) {
        if (app.epoch == epoch) {
            ASMatrix(SIZE_T_TO_INDEX(app.i), SIZE_T_TO_INDEX(app.j)) +=
                app.coeff * admixtureValues.at(app.v1) * admixtureValues.at(app.v2);
        }
    }
#ifndef NDEBUG
    verifyASMatrix(ASMatrix);
#endif
    return std::move(ASMatrix);
}

inline size_t numTimeSlices(const std::vector<double>& timeSlices) { return timeSlices.size() + 1; }

struct ModelProbabilities {
    MatrixXd locations;
    MatrixXd coalescence;
};

/**
 * Computes the model probability from t=0 to t=k (cumulative).
 */
ModelProbabilities probabilitiesUpToTime(const MatrixXd& sMatrix, const double timeK) {
    // XXX this solution should be interpretable in the following way:
    // - Row is the starting state at t=0
    // - Col is the ending state at t=K
    // The solution is a single matrix that can be broken into two:
    //  1. P_m[i, j] == "the probability that we started at state i and ended up
    //  at state j".
    //  2. P_c[i] == "the probability that we started at state i and ended up
    //  coalescing"
    // For both of these, state i = (a, b) implying one individual in deme a and
    // one individual in deme b Each row in the full matrix (P_m concatenated with
    // P_c) should:
    // - Sum to 1 (approximately, there will be floating pt error)
    // - Have no cell < 0
    // - Have no cell > 1
    const MatrixXd intermediate = sMatrix * timeK;
    TRACE_MATRIX(intermediate, "exponential parameter @ t=" << timeK);
    const MatrixXd stateProbabilities = intermediate.exp();
    TRACE_MATRIX(stateProbabilities, "State probability @ t=" << timeK);
    const size_t N = stateProbabilities.cols(); // NxN matrix
    return {stateProbabilities.block(0, 0, N - 1, N - 1), stateProbabilities.topRightCorner(N - 1, 1)};
}

struct TimeMarker {
    double time;
    size_t epoch;
    bool isSlice;
};

// Takes two vectors, one for time slice and one for epochs, and outputs a
// vector of objects indicating what kind of time marker it is (a time slice, or
// epoch boundary?) and which Epoch it belongs to.
inline std::vector<TimeMarker> combineTimeVectors(const std::vector<double>& timeSlices,
                                                  const std::vector<double>& epochTimes) {
    std::vector<TimeMarker> tsAndEpoch;
    size_t epochCounter = 0;
    {
        size_t posT = 0;
        size_t posE = 0;
        while (posT < timeSlices.size() || posE < epochTimes.size()) {
            if (posT >= timeSlices.size()) {
                tsAndEpoch.push_back({epochTimes[posE], ++epochCounter, false});
                posE++;
            } else if (posE >= epochTimes.size()) {
                tsAndEpoch.push_back({timeSlices[posT], epochCounter, true});
                posT++;
            } else if (timeSlices[posT] < epochTimes[posE]) {
                tsAndEpoch.push_back({timeSlices[posT], epochCounter, true});
                posT++;
            } else if (timeSlices[posT] == epochTimes[posE]) {
                tsAndEpoch.push_back({timeSlices[posT], ++epochCounter, true});
                posT++;
                posE++;
            } else {
                tsAndEpoch.push_back({epochTimes[posE], ++epochCounter, false});
                posE++;
            }
        }
    }
    return std::move(tsAndEpoch);
}

static inline void rowNorm(MatrixXd& matrix) {
    for (Eigen::Index j = 0; j < matrix.rows(); j++) {
        matrix.row(j).array() /= matrix.row(j).sum();
    }
}

/**
 * Given concrete matrices and vectors (so the symbolic ones have been populated
 * by the current parameter values), compute the probabilities.
 */
MatrixXd modelPMFByTimeWithEpochs(const ParameterSchema& schema,
                                  double const* parameters,
                                  const std::vector<double>& timeSlices,
                                  const std::vector<double>& epochTimes,
                                  std::vector<std::pair<double, MatrixXd>>* qMatricesByTime = nullptr) {
    const size_t numEpochs = schema.numEpochs();
    assert(numEpochs > 0);
    assert(numEpochs == epochTimes.size() + 1);

    // Combine all the times into a single timeline just to simplify things.
    std::vector<TimeMarker> allTimes = combineTimeVectors(timeSlices, epochTimes);

    // Epoch0 is special. The number of states in this epoch corresponds with the
    // concrete data that we have, which we need _ALL_ model probabilities to fit
    // (dimension-wise) for the likelihood calculation. I.e., epochK may have
    // fewer states, but after we calculate the state probabilities against those
    // fewer states we always have to map them back to Epoch0's states.
    const Eigen::Index nStatesEpoch0 = SIZE_T_TO_INDEX(schema.numStates(0)) - 1;

    // Maps the current epoch's states back to Epoch0, in proportions.
    MatrixXd currentStateMap = MatrixXd::Identity(nStatesEpoch0, nStatesEpoch0);
    // Resulting coalescence probabilities per time slice
    const Eigen::Index nTimeBins = SIZE_T_TO_INDEX(timeSlices.size()) + 1;
    MatrixXd probsByTime = MatrixXd::Ones(nStatesEpoch0, nTimeBins);
    // The probability of non-coalescence states at the end of the previous epoch.
    // This is the "starting state" of the current epoch, and is updated at the end of each epoch.
    MatrixXd currentStateProbs = MatrixXd::Identity(nStatesEpoch0, nStatesEpoch0);
#if FEWER_SUBTRACTIONS_PER_EPOCH
    // The probability that a state has coalesced by the end of the last epoch.
    VectorXd probCoalByLastEpoch = VectorXd::Zero(nStatesEpoch0);
#else
    // The probability that a state has NOT coalesced by the end of the last
    // epoch.
    VectorXd probNotCoalByLastEpoch = VectorXd::Ones(nStatesEpoch0);
#endif

    size_t currentEpoch = 0;
    Eigen::Index currentSlice = 0;
    double epochStart = 0;
    for (size_t i = 0; i < allTimes.size(); i++) {
        const double time = allTimes[i].time;
        const size_t newEpoch = allTimes[i].epoch;
        const double deltaT = time - epochStart;

        const MatrixXd curMatrix = std::move(createQMatrix(schema, parameters, currentEpoch, deltaT));
        if (qMatricesByTime != nullptr) {
            qMatricesByTime->emplace_back(time, curMatrix);
        }

        ModelProbabilities probabilities = probabilitiesUpToTime(curMatrix, deltaT);
        probabilities.coalescence = currentStateProbs * probabilities.coalescence;
        TRACE_MATRIX(probabilities.coalescence,
                     "coalescence probabilities in epoch " << currentEpoch << " up to time " << time);

        // We only record probabilities at time slice boundaries. All other
        // calculations update our running values for epoch calculations, which
        // _implicitly_ affect the probabilities.
        if (allTimes[i].isSlice) {
#if FEWER_SUBTRACTIONS_PER_EPOCH
            probsByTime.col(currentSlice) = ((probCoalByLastEpoch.array() + probabilities.coalescence.array()) -
                                             (probCoalByLastEpoch.array() * probabilities.coalescence.array()))
                                                .matrix();
#else
            // ProbabilityWeDidCoalesce = !(DidntCoalesceLastEpoch &&
            // !DidCoalesceByNow)
            probsByTime.col(currentSlice) =
                (1 - (probNotCoalByLastEpoch.array() * (1 - probabilities.coalescence.array()))).matrix();
#endif
            TRACE_MATRIX(probsByTime.col(currentSlice), "CDF coalescence for time " << time);
            currentSlice++;
        }

        // The epoch has changed, update all the intermediate mappings/values
        if (newEpoch != currentEpoch) {
            assert(newEpoch == currentEpoch + 1);

            // Map our current location probabilities back to the epoch0 locations.
            MatrixXd locProbMappedToEpoch0 = currentStateMap * probabilities.locations;
            TRACE_MATRIX(locProbMappedToEpoch0, "End of epoch " << currentEpoch << " locs");
            // locProbMappedToEpoch0 contains the location probabilities in terms of
            // the initial states, but only for the time period between the end of the
            // previous epoch (E_k) and the start of the current time slice (which is
            // the start of a new epoch, E_k+2). This reweights all the migration
            // states based on the location probabilities at the end of E_k.
            currentStateProbs = currentStateProbs * locProbMappedToEpoch0;
            TRACE_MATRIX(currentStateProbs, " ... multiplied by end of previous epoch");
            rowNorm(currentStateProbs); // Normalize each row to be a probability.
            TRACE_MATRIX(currentStateProbs, " ... normalized and converted to EOPE");

            const auto ASMatrix = createASMatrix(schema, parameters, newEpoch);
            TRACE_MATRIX(ASMatrix, "ASMatrix(" << newEpoch << ")");
            currentStateMap = currentStateMap * ASMatrix;
            TRACE_MATRIX(currentStateMap, "currentStateMap(" << newEpoch << ")");
            currentStateProbs = currentStateProbs * currentStateMap;

#if FEWER_SUBTRACTIONS_PER_EPOCH
            probCoalByLastEpoch = ((probCoalByLastEpoch.array() + probabilities.coalescence.array()) -
                                   (probCoalByLastEpoch.array() * probabilities.coalescence.array()))
                                      .matrix();
            TRACE_MATRIX(probCoalByLastEpoch, "probCoalByLastEpoch");
#else
            probNotCoalByLastEpoch =
                (probNotCoalByLastEpoch.array() * (1 - probabilities.coalescence.array())).matrix();
            TRACE_MATRIX(probNotCoalByLastEpoch, "probNotHavingCoalLastEpoch");
#endif

            epochStart = time;
            currentEpoch = newEpoch;
        }
    }
    // If the solver pushes an epoch time into the last time slice, it will just be 1-p, and
    // not directly affect the result (but from a degrees-of-freedom perspective, it should).

    TRACE_MATRIX(probsByTime, "CDF of state probabilities");
    // Convert the CDF --> PMF. Also adds one more column for the "infinity"
    // (backwards in time) bucket, which for the CDF will be probability 1, and
    // for the PMF will get the residual probability.
    static double minPmfVal = FLOAT_EPS;
    for (Eigen::Index i = nTimeBins; i > 1; i--) {
        // We linearize the last bit if it becomes non-monotonic. This can happen
        // for certain random initialization values that are far away from the
        // optimum.
        probsByTime.col(i - 1) = (probsByTime.col(i - 1) - probsByTime.col(i - 2)).array().max(minPmfVal).matrix();
    }
    TRACELN("-----------------------");

    return std::move(probsByTime);
}

class ObservedData {
public:
    std::vector<MatrixXd> countMatrices;
    std::vector<double> timeSlices;
};

class CachedCMatrix {
public:
    void set(const MatrixXd& mat) {
        this->matrix = mat;
        this->isSet = true;
    }

    MatrixXd matrix;
    bool isSet = false;
};

NegLogLikelihoodCostFunctor::NegLogLikelihoodCostFunctor(const std::string& jsonFile,
                                                         bool logParams,
                                                         const bool maximization,
                                                         const ObservationMode mode,
                                                         const bool normalize)
    : m_schema(logParams),
      m_observed(new ObservedData),
      m_cache(new CachedCMatrix),
      m_maximization(maximization) {
    resetSolveStats();
    std::ifstream inputData(jsonFile);
    json inputJson = json::parse(inputData);

    // Parse all the inputs.
    auto& coalCountsJson = inputJson.at(COAL_COUNTS_KEY);
    m_observed->countMatrices = loadCMatrices(coalCountsJson);
    RELEASE_ASSERT(!m_observed->countMatrices.empty());

    // Normalize the matrices so that each row sums to 1
    if (normalize) {
        for (size_t i = 0; i < m_observed->countMatrices.size(); i++) {
            auto& matrix = m_observed->countMatrices[i];
            for (Eigen::Index j = 0; j < matrix.rows(); j++) {
                matrix.row(j).array() /= matrix.row(j).sum();
            }
        }
    }
    m_observed->timeSlices = ::getTimeSlices(inputJson);
    m_schema.load(inputJson);

#if DEBUG_OUTPUT
    std::cerr << "Loaded " << m_observed->countMatrices.size() << " count matrices" << std::endl;
    std::cerr << "... with " << m_observed->countMatrices.at(0).rows() << "x" << m_observed->countMatrices.at(0).cols()
              << " matrices" << std::endl;
    std::cerr << "Loaded " << m_observed->timeSlices.size() << " timeslices" << std::endl;
    std::cerr << "Loaded " << m_schema.totalParams() << " solver parameters: " << std::endl;
#endif

    if (mode == OBS_MODE_UNSPECIFIED) {
        m_mode = parse_obs_mode(inputJson.at(OBSERVATION_MODE_KEY));
    } else {
        m_mode = mode;
    }

    resetCMatrixCache();
}

NegLogLikelihoodCostFunctor::~NegLogLikelihoodCostFunctor() {
    delete m_observed;
    delete m_cache;
}

size_t NegLogLikelihoodCostFunctor::numCoalMatrices() const { return m_observed->countMatrices.size(); }

const std::vector<double>& NegLogLikelihoodCostFunctor::getTimeSlices() const { return m_observed->timeSlices; }

/**
 * Force the cached CMatrix to be a specific CMatrix.
 */
void NegLogLikelihoodCostFunctor::restrictToCMatrix(const size_t index) {
    assert(!m_observed->countMatrices.empty());
    m_cache->set(m_observed->countMatrices.at(index));
}

/**
 * Set the cached CMatrix back to whatever the schema requires (average, first,
 * all).
 */
void NegLogLikelihoodCostFunctor::resetCMatrixCache() {
    assert(!m_observed->countMatrices.empty());
    switch (m_mode) {
    case OBS_MODE_AVERAGE: {
        const size_t numMats = m_observed->countMatrices.size();
        assert(numMats > 0);
        MatrixXd avg = m_observed->countMatrices.at(0);
        for (size_t i = 1; i < numMats; i++) {
            avg += m_observed->countMatrices[i];
        }
        avg /= (double)numMats;
        m_cache->set(avg);
    } break;
    case OBS_MODE_FIRST: {
        m_cache->set(m_observed->countMatrices.at(0));
    } break;
    default: break;
    }
}

inline double computeLLForTime(Eigen::Index timeK, const MatrixXd& probabilityAtTime, const MatrixXd& countMatrix) {
    MatrixXd timeKProbVector = probabilityAtTime.col(timeK);
    TRACELN("=====================================");
    TRACE_MATRIX(timeKProbVector.transpose(), "timeKProbVector");
    MatrixXd countsForK = countMatrix.col(timeK).transpose();
    TRACE_MATRIX(countsForK, "countsForK");
    // max(FLOAT_EPS) to avoid zeroes in the logarithm
    MatrixXd logModelProbForK = (timeKProbVector.array().max(FLOAT_EPS).log()).matrix();
    TRACE_MATRIX(logModelProbForK.transpose(), "logModelProbForK");
    const auto product = (countsForK * logModelProbForK);
    return product(0, 0);
}

double NegLogLikelihoodCostFunctor::operator()(double const* parameters) const {
    // One for each transition between epochs
    std::vector<double> epochTimes = m_schema.getEpochStartTimes(parameters);
    double cost = 0.0;
    try {
        auto probabilityAtTime = modelPMFByTimeWithEpochs(m_schema, parameters, m_observed->timeSlices, epochTimes);
        RELEASE_ASSERT(m_cache->isSet);
        for (Eigen::Index timeK = 0; timeK < probabilityAtTime.cols(); timeK++) {
            const double ll = computeLLForTime(timeK, probabilityAtTime, m_cache->matrix);
            if (m_maximization) {
                cost += ll;
            } else {
                cost -= ll;
            }
            TRACELN("Cost after timeK=" << timeK << ": " << cost);
        }
    } catch (const ModelAssertFailure& e) {
        std::cerr << "FAILURE PARAMETERS: " << std::endl;
        for (size_t i = 0; i < m_schema.totalParams(); i++) {
            std::cerr << m_schema.fromParam(parameters[i], i) << ", ";
        }
        std::cerr << std::endl;
        throw e;
    }
#if 0
    // Debugging code for NaN issues. Generally NaNs will happen just due to edge cases with things like
    // growth rates, where a really large growth rate can cause the Q-matrix to have huge rates that cannot
    // be solved for.
    if (std::isnan(cost)) {
        std::cerr << "WARNING: NaN cost" << std::endl;
        std::cerr << "FAILURE PARAMETERS: " << std::endl;
        for (size_t i = 0; i < m_schema.totalParams(); i++) {
            std::cerr << m_schema.fromParam(parameters[i], i) << ", ";
        }
        std::cerr << std::endl;
        exit(1);
    }
#endif
    return cost;
}

void NegLogLikelihoodCostFunctor::outputTheoreticalCoalMatrix(double const* parameters, std::ostream& out) const {
    std::vector<double> epochTimes = m_schema.getEpochStartTimes(parameters);

    auto probabilityAtTime = modelPMFByTimeWithEpochs(m_schema, parameters, m_observed->timeSlices, epochTimes);

    json outputJson = m_schema.toJsonOutput(parameters, std::numeric_limits<double>::quiet_NaN());
    std::vector<std::vector<double>> matrix;
    for (Eigen::Index i = 0; i < probabilityAtTime.rows(); i++) {
        std::vector<double> row;
        for (Eigen::Index j = 0; j < probabilityAtTime.cols(); j++) {
            row.push_back(probabilityAtTime(i, j));
        }
        matrix.push_back(row);
    }
    std::vector<std::vector<std::vector<double>>> matrices = {matrix};
    outputJson[COAL_COUNTS_KEY] = matrices;
    out << outputJson;
}
