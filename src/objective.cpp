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

// Flags that affect how we approximate the application of the growth rate.
// Use the average coalescence rate from start of epoch to time T_k. The smaller T_k is, the
// more accurate this will be.
#define USE_AVG_COAL_RATE 1
// Use the the actual coalescence rate at time T_k, which will be the largest rate. This over-
// approximates the rate and under-approximates the Ne over T_k. When T_k is large, the coal
// rate is large (Ne is small), and this will _really_ over-approximate.
#define USE_OVER_APPROX_COAL_RATE 0

// When computing probabilities in multi-epoch models, this uses the additive probability rule
// P(A or B) = P(A) + P(B) - P(A and B), instead of negating all of the probabilities to only
// use P(not(A) and not(B)). The former results in fewer subtractions, which anecdotally appears
// to produce slightly better point estimates with definitely better confidence intervals.
#define FEWER_SUBTRACTIONS_PER_EPOCH 1

#define DUMP_MATRIX(m, desc)                                                                                           \
    do {                                                                                                               \
        std::cerr << "% " << desc << ":" << std::endl;                                                                 \
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

MatrixXd loadPopConversion(json& inputJson) { return loadJsonMatrix(inputJson, POP_CONVERT_KEY); }

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

void ParameterSchema::load(const json& inputData) {
    m_paramRescale.resize(0);
    m_inputJson = inputData;
    for (const auto& parameter : inputData.at(EPOCH_TIMES_KEY)) {
        double ground_truth = std::numeric_limits<float>::quiet_NaN();
        if (parameter.contains("ground_truth")) {
            ground_truth = parameter["ground_truth"];
        }
        BoundedVariable bv = {
            parameter["init"], parameter["lb"], parameter["ub"], {}, parameter["description"], ground_truth};
        const bool isFixed = parameter.contains("fixed") && (bool)parameter["fixed"];
        if (isFixed) {
            m_eFixed.push_back(std::move(bv));
        } else {
            m_eParams.push_back(std::move(bv));
            m_paramRescale.emplace_back(EPOCH_TIME_RESCALE);
        }
    }
    RELEASE_ASSERT(m_eParams.empty() || m_eFixed.empty()); // Only all fixed or all params
    m_numEpochs = m_eParams.size() + m_eFixed.size() + 1;
    m_sStates = std::vector<size_t>(m_numEpochs);
    for (const auto& parameter : inputData.at(SMATRIX_VALS_KEY)) {
        std::list<VariableApplication> applications;
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
        double ground_truth = std::numeric_limits<float>::quiet_NaN();
        if (parameter.contains("ground_truth")) {
            ground_truth = parameter["ground_truth"];
        }
        BoundedVariable bv = {parameter["init"],
                              parameter["lb"],
                              parameter["ub"],
                              std::move(applications),
                              parameter["description"],
                              ground_truth};
        const bool isFixed = parameter.contains("fixed") && (bool)parameter["fixed"];
        if (isFixed) {
            m_sFixed.push_back(std::move(bv));
        } else {
            m_sParams.push_back(std::move(bv));
            m_paramRescale.emplace_back(1.0);
        }
    }
    RELEASE_ASSERT(m_sParams.size() + m_eParams.size() == m_paramRescale.size());
}

/**
 * Populate the parameters[] vector with random values based on the lower/upper
 * bounds from the parameter schema.
 */
void ParameterSchema::randomParamVector(double* parameters) const {
    size_t p = 0;
    for (const auto& parameter : m_eParams) {
        parameters[p] = toParam(randDouble(parameter.lb, parameter.ub), p);
        p++;
    }
    for (const auto& parameter : m_sParams) {
        parameters[p] = toParam(randDouble(parameter.lb, parameter.ub), p);
        p++;
    }
}

/**
 * Populate the parameters[] vector with the initial value provided in the
 * parameter schema.
 */
void ParameterSchema::initParamVector(double* parameters) const {
    size_t p = 0;
    for (const auto& parameter : m_eParams) {
        parameters[p] = toParam(parameter.init, p);
        p++;
    }
    for (const auto& parameter : m_sParams) {
        parameters[p] = toParam(parameter.init, p);
        p++;
    }
}

void ParameterSchema::getBounds(size_t paramIdx, double& lowerBound, double& upperBound) const {
    if (paramIdx < m_eParams.size()) {
        lowerBound = toParam(m_eParams.at(paramIdx).lb, paramIdx);
        upperBound = toParam(m_eParams.at(paramIdx).ub, paramIdx);
    } else {
        const size_t sIdx = paramIdx - m_eParams.size();
        lowerBound = toParam(m_sParams.at(sIdx).lb, paramIdx);
        upperBound = toParam(m_sParams.at(sIdx).ub, paramIdx);
    }
}

double ParameterSchema::getEpochStartTime(double const* parameters, size_t epoch) const {
    assert(epoch < m_numEpochs);
    double epochTime = 0.0;
    if (epoch > 0) {
        const size_t paramIdx = epoch - 1;
        if (!m_eParams.empty()) {
            epochTime = fromParam(parameters[paramIdx], paramIdx);
        } else {
            epochTime = m_eFixed.at(paramIdx).init;
        }
    }
    return epochTime;
}

json ParameterSchema::toJsonOutput(const double* parameters, const double negLL) const {
    json output = m_inputJson; // Make a copy of the input.
    size_t e = 0;
    size_t p = 0;
    RELEASE_ASSERT(output.at(EPOCH_TIMES_KEY).size() == m_eParams.size() + m_eFixed.size());
    for (auto& parameter : output.at(EPOCH_TIMES_KEY)) {
        if (!m_eParams.empty()) {
            parameter["final"] = fromParam(parameters[p], p);
            p++;
        } else {
            parameter["final"] = m_eFixed.at(e).init;
            e++;
        }
    }
    RELEASE_ASSERT(output.at(SMATRIX_VALS_KEY).size() == m_sParams.size());
    for (auto& parameter : output.at(SMATRIX_VALS_KEY)) {
        parameter["final"] = fromParam(parameters[p], p);
        p++;
    }
    output["negLL"] = negLL;
    return std::move(output);
}

void ParameterSchema::fromJsonOutput(const json& jsonOutput, double* parameters, std::string key) const {
    size_t p = 0;
    for (const auto& paramVal : jsonOutput[EPOCH_TIMES_KEY]) {
        const bool isFixed = paramVal.contains("fixed") && (bool)paramVal["fixed"];
        if (!isFixed) {
            parameters[p] = toParam((double)paramVal[key], p);
            p++;
        }
    }
    RELEASE_ASSERT(m_eParams.size() == p);
    for (const auto& paramVal : jsonOutput[SMATRIX_VALS_KEY]) {
        const bool isFixed = paramVal.contains("fixed") && (bool)paramVal["fixed"];
        if (!isFixed) {
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

    bool anyAdjusted = false;
    const Eigen::Index nStates = SIZE_T_TO_INDEX(schema.numStates(epoch));
    // Skip this many epoch transition times.
    size_t firstParamIdx = schema.m_eParams.size();
    MatrixXd qMatrix = MatrixXd::Zero(nStates, nStates);
    for (size_t p = 0; p < schema.m_sParams.size(); p++) {
        const size_t paramIdx = firstParamIdx + p;
        double pVal = schema.fromParam(parameters[paramIdx], paramIdx);
        for (const auto& application : schema.m_sParams[p].applications) {
            if (application.epoch == epoch) {
                const Eigen::Index i = SIZE_T_TO_INDEX(application.i);
                const Eigen::Index j = SIZE_T_TO_INDEX(application.j);
                // Adjustments always come last. See model.py
                if (application.adjust == ADJUST_GROWTH_RATE) {
                    const double origRate = qMatrix(i, j);
#if USE_AVG_COAL_RATE
                    // The integral on interval {0, u} is (1/a - exp(-a*u)/a)
                    // we can divide by alpha to get the average rate over time period lower:upper.
                    const double integratedRate =
                        (std::exp(pVal * timeSinceEpochStart) - 1) / (pVal * timeSinceEpochStart);
                    qMatrix(i, j) = origRate * integratedRate * application.coeff;
#elif USE_OVER_APPROX_COAL_RATE
                    // rateB = base_rate * math.exp(r_EU * i)
                    qMatrix(i, j) = origRate * std::exp(pVal * timeSinceEpochStart) * application.coeff;
                    // std::cout << "OrigRate = " << origRate << ", NewRate = " << qMatrix(i, j) << "\n";
#else
                    static_assert(false, "Invalid preprocessor definition combination");
#endif
                    anyAdjusted = true;
                } else if (application.adjust == ADJUST_INV_GROWTH_RATE) {
                    // For backwards compatibility with JSON input files. Just ignore.
                } else {
                    qMatrix(i, j) += pVal * application.coeff;
                    assert(!anyAdjusted);
                }
            }
        }
    }
    for (const auto& param : schema.m_sFixed) {
        for (const auto& application : param.applications) {
            if (application.epoch == epoch) {
                const Eigen::Index i = SIZE_T_TO_INDEX(application.i);
                const Eigen::Index j = SIZE_T_TO_INDEX(application.j);
                qMatrix(i, j) += param.init * application.coeff;
            }
        }
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
    auto intermediate = sMatrix * timeK;
    TRACE_MATRIX(intermediate, "exponential parameter @ t=" << timeK);
    MatrixXd stateProbabilities = intermediate.exp();
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
    std::vector<TimeMarker> result;
    size_t posT = 0;
    size_t posE = 0;
    size_t epochCounter = 0;
    while (posT < timeSlices.size() || posE < epochTimes.size()) {
        if (posT >= timeSlices.size()) {
            result.push_back({epochTimes[posE], ++epochCounter, false});
            posE++;
        } else if (posE >= epochTimes.size()) {
            result.push_back({timeSlices[posT], epochCounter, true});
            posT++;
        } else if (timeSlices[posT] < epochTimes[posE]) {
            result.push_back({timeSlices[posT], epochCounter, true});
            posT++;
        } else if (timeSlices[posT] == epochTimes[posE]) {
            result.push_back({timeSlices[posT], ++epochCounter, true});
            posT++;
            posE++;
        } else {
            result.push_back({epochTimes[posE], ++epochCounter, false});
            posE++;
        }
    }
    return std::move(result);
}

/**
 * Given concrete matrices and vectors (so the symbolic ones have been populated
 * by the current parameter values), compute the probabilities.
 */
MatrixXd modelPMFByTimeWithEpochs(const ParameterSchema& schema,
                                  double const* parameters,
                                  const std::vector<double>& timeSlices,
                                  const std::vector<double>& epochTimes,
                                  const MatrixXd& popConvert,
                                  std::vector<std::pair<double, MatrixXd>>* qMatricesByTime = nullptr) {
    const size_t numEpochs = schema.numEpochs();
    assert(numEpochs > 0);
    assert(numEpochs == epochTimes.size() + 1);
    // Each row of popConvert is for the transition between epoch e and epoch e+1
    assert(numEpochs == popConvert.rows() + 1);

    // Combine all the times into a single timeline just to simplify things.
    std::vector<TimeMarker> allTimes = combineTimeVectors(timeSlices, epochTimes);

    // Epoch0 is special. The number of states in this epoch corresponds with the
    // concrete data that we have, which we need _ALL_ model probabilities to fit
    // (dimension-wise) for the likelihood calculation. I.e., epochK may have
    // fewer states, but after we calculate the state probabilities against those
    // fewer states we always have to map them back to Epoch0's states.
    const Eigen::Index nStatesEpoch0 = SIZE_T_TO_INDEX(schema.numStates(0)) - 1;

    // Maps the current epoch's states back to Epoch0. E.g., if position "i" has
    // value "k", it means that state "k" in the current epoch maps back to state
    // "i" in Epoch0.
    const double minPop = 0;
    const double maxPop = (double)nStatesEpoch0 - 1;
    VectorXd currentStateMap = VectorXd::LinSpaced(nStatesEpoch0, minPop, maxPop);
    // Resulting coalescence probabilities per time slice
    const Eigen::Index nTimeBins = SIZE_T_TO_INDEX(timeSlices.size()) + 1;
    MatrixXd probsByTime = MatrixXd::Ones(nStatesEpoch0, nTimeBins);
    // The probability of non-coalescence states at the end of the previous epoch.
    // This is the "starting state" of the current epoch. EOPE = "end of previous
    // epoch"
    MatrixXd migrationStateProbsEOPE = MatrixXd::Identity(nStatesEpoch0, nStatesEpoch0);
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
        probabilities.coalescence = migrationStateProbsEOPE * probabilities.coalescence(currentStateMap, 0);
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
            MatrixXd locProbMappedToEpoch0 = probabilities.locations(currentStateMap, currentStateMap);
            TRACE_MATRIX(locProbMappedToEpoch0, "End of epoch " << currentEpoch << " locs");
            // locProbMappedToEpoch0 contains the location probabilities in terms of
            // the initial states, but only for the time period between the end of the
            // previous epoch (E_k) and the start of the current time slice (which is
            // the start of a new epoch, E_k+2). This reweights all the migration
            // states based on the location probabilities at the end of E_k.
            migrationStateProbsEOPE = migrationStateProbsEOPE * locProbMappedToEpoch0;
            TRACE_MATRIX(migrationStateProbsEOPE, " ... multiplied by end of previous epoch");
            // Normalize each row to be a probability.
            for (Eigen::Index j = 0; j < migrationStateProbsEOPE.rows(); j++) {
                migrationStateProbsEOPE.row(j).array() /= migrationStateProbsEOPE.row(j).sum();
            }
            TRACE_MATRIX(migrationStateProbsEOPE, " ... normalized and converted to EOPE");

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
            currentStateMap = popConvert.row(SIZE_T_TO_INDEX(currentEpoch));
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
    MatrixXd popConvert;
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
    m_observed->popConvert = loadPopConversion(inputJson);

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
    std::vector<double> epochTimes; // One for each transition between epochs
    epochTimes.reserve(m_schema.numEpochs() - 1);
    for (size_t epoch = 0; epoch < m_schema.numEpochs(); epoch++) {
        if (epoch > 0) {
            epochTimes.emplace_back(m_schema.getEpochStartTime(parameters, epoch));
        }
    }
    double cost = 0.0;
    try {
        auto probabilityAtTime =
            modelPMFByTimeWithEpochs(m_schema, parameters, m_observed->timeSlices, epochTimes, m_observed->popConvert);
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
    std::vector<double> epochTimes;
    epochTimes.reserve(m_schema.numEpochs() - 1);
    for (size_t epoch = 0; epoch < m_schema.numEpochs(); epoch++) {
        if (epoch > 0) {
            epochTimes.emplace_back(m_schema.getEpochStartTime(parameters, epoch));
        }
    }

    auto probabilityAtTime =
        modelPMFByTimeWithEpochs(m_schema, parameters, m_observed->timeSlices, epochTimes, m_observed->popConvert);

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