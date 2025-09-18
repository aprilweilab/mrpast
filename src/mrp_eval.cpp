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
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include <Eigen/Dense>
#include <args.hxx>
#include <json.hpp>

extern "C" {
#include <nlopt.h>
}

#include "common.h"
#include "derivatives.h"
#include "objective.h"
#include "solve.h"

using json = nlohmann::json;

// Compute the parameter CIs with a single inverse calculation. This is much
// more numerically stable than taking the inverse of J and GIM.
#define SINGLE_INVERSE_CALCULATION 1

#if TRACE_MATRICES
#define TRACE_MATRIX(m, desc)                                                                                          \
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

#define TRACELN(msg) std::cerr << msg << std::endl;
#else
#define TRACE_MATRIX(m, desc)
#define TRACELN(msg)
#endif

template <typename Out> inline void split(const std::string& s, char delim, Out result) {
    std::istringstream iss(s);
    std::string item;
    while (std::getline(iss, item, delim)) {
        *result++ = item;
    }
}

inline std::vector<std::string> split(const std::string& s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return std::move(elems);
}

inline void ltrim(std::string& s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) { return !std::isspace(ch); }));
}

inline void rtrim(std::string& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
}

inline void trim(std::string& s) {
    ltrim(s);
    rtrim(s);
}

inline bool ends_with(std::string const& string1, std::string const& string2) {
    if (string1.length() < string2.length()) {
        return false;
    }
    return (string2 == string1.substr(string1.length() - string2.length()));
}

template <typename T> T convert(const std::string& strValue);

template <> double convert(const std::string& strValue) {
    char* endPtr = nullptr;
    double result = std::strtod(strValue.c_str(), &endPtr);
    if (endPtr != (strValue.c_str() + strValue.size())) {
        std::stringstream ssErr;
        ssErr << "Cannot parse as double: " << strValue;
        throw std::runtime_error(ssErr.str());
    }
    return result;
}

// This only works for "simple" CSV files, where there are no commas (escaped or otherwise)
// in the actual data, so we can just split by ",".
std::vector<double> loadNegLL(std::ifstream& inStream) {
    const std::string colName = "negLL";
    const std::string onlyRowIndex = "0"; // Only take the 0th parameter row
    std::string line;
    size_t column = std::numeric_limits<size_t>::max();
    bool first = true;
    std::vector<double> result;
    while (std::getline(inStream, line)) {
        auto tokens = split(line, ',');
        if (first) {
            for (column = 0; column < tokens.size(); column++) {
                trim(tokens[column]);
                if (tokens[column] == colName) {
                    break;
                }
            }
            first = false;
        } else {
            if (column >= tokens.size()) {
                std::stringstream ssErr;
                ssErr << "Cannot find the \"" << colName << "\" column in CSV file";
                throw std::runtime_error(ssErr.str());
            }
            trim(tokens[0]); // This is the index column
            // The CSV file is setup like a dataframe, so we have a row for each parameter that was
            // inferred. We only take the 0th parameter, because the likelihood is the same for all
            // parameters from the same bootstrap sample.
            if (tokens[0] == onlyRowIndex) {
                trim(tokens[column]);
                result.emplace_back(convert<double>(tokens[column]));
            }
        }
    }
    return std::move(result);
}

// Return the inverse of the given matrix.
inline MatrixXd invert(const MatrixXd& m) {
    Eigen::FullPivLU<MatrixXd> lu(m);
    RELEASE_ASSERT(lu.isInvertible());
    return lu.inverse();
}

inline bool isInvertible(const MatrixXd& m) {
    Eigen::FullPivLU<MatrixXd> lu(m);
    return lu.isInvertible();
}

inline MatrixXd copyWithoutRowCol(const MatrixXd& m, const size_t holdOut) {
    if (holdOut >= m.rows()) {
        return m;
    }
    MatrixXd result = MatrixXd::Zero(m.rows() - 1, m.cols() - 1);
    for (Eigen::Index j = 0; j < m.rows(); j++) {
        if (j == holdOut) {
            continue;
        }
        const size_t jDiff = j > holdOut ? 1 : 0;
        for (Eigen::Index k = 0; k < m.cols(); k++) {
            if (k == holdOut) {
                continue;
            }
            const size_t kDiff = k > holdOut ? 1 : 0;
            result(j - jDiff, k - kDiff) = m(j, k);
        }
    }
    return std::move(result);
}

// NOTE: This is not really required any longer, as all the parameters are scaled such that the
// H matrix is always invertible. However, if someone produces an H that is not invertible this
// will hold one parameter constant (usually an epoch time) in order to make the rest of the
// matrix invertible.
inline MatrixXd findFirstInvertible(const MatrixXd& input, size_t& heldOut) {
    if (isInvertible(input)) {
        heldOut = std::numeric_limits<size_t>::max();
        return input;
    }
    std::cerr << "WARNING: Original matrix is NOT invertible, attempting to solve this by leaving out "
                 " a single parameter (fixing it to the theta-hat value)"
              << std::endl;
    MatrixXd result;
    for (heldOut = 0; heldOut < input.rows(); heldOut++) {
        result = copyWithoutRowCol(input, heldOut);
        if (isInvertible(result)) {
            std::cerr << "Leaving out parameter " << heldOut << " produces an invertible matrix" << std::endl;
            return result;
        }
    }
    heldOut = std::numeric_limits<size_t>::max();
    return input;
}

// Prevent parameters with lower bounds that are negative.
bool checkCostFunc(const NegLogLikelihoodCostFunctor& cost) {
    const double lbThreshold = 0.0;
    for (size_t i = 0; i < cost.m_schema.totalParams(); i++) {
        double lowerBound = 0.0;
        double upperBound = 0.0;
        cost.m_schema.getBounds(i, lowerBound, upperBound);
        if (lowerBound < lbThreshold) {
            std::cerr << "ERROR: Lowerbound of parameter " << i << " is " << lowerBound << " < " << lbThreshold
                      << std::endl;
            return false;
        }
    }
    return true;
}

/**
 * Get the J and H matrices for the likelihood function.
 *
 * @param[in] fixNonInvertible When set to true, attempt to fix non-invertible H
 *  matrices by holding one of the parameters at a constant value (and otherwise
 *  ignoring it in all all derivative calculations). Rarely needed now that we
 *  scale the epoch times. USE WITH CAUTION.
 */
std::pair<MatrixXd, MatrixXd> getJandH(NegLogLikelihoodCostFunctor& cost,
                                       const std::vector<double>& paramVector,
                                       size_t& leftOut,
                                       bool bootstrap = true,
                                       bool fixNonInvertible = false) {
    Eigen::MatrixXd sensitivityH = calculateSensitivity(cost, paramVector);
    TRACE_MATRIX(sensitivityH, "H(theta)");
    if (fixNonInvertible) {
        sensitivityH = findFirstInvertible(sensitivityH, leftOut);
    } else {
        leftOut = std::numeric_limits<size_t>::max();
    }

    auto variabilityJ = calculateVariability(cost, paramVector, bootstrap);
    TRACE_MATRIX(variabilityJ, "J(theta)");
    if (fixNonInvertible) {
        variabilityJ = copyWithoutRowCol(variabilityJ, leftOut);
    }
    return {variabilityJ, sensitivityH};
}

args::Group arguments("arguments");
args::HelpFlag h(arguments, "help", "help", {'h', "help"});
args::Flag
    normCoals(arguments, "normCoals", "Normalize coal matrix prior to computing likelihoods", {'n', "norm-coals"});

void intervalsCommand(args::Subparser& parser) {
    args::Positional<std::string> inputFile(
        parser, "inputFile", "Input file: a JSON file as emitted by mrpast's solver");
    args::Flag fixNonInvertible(parser,
                                "fixNonInvertible",
                                "When the H matrix is non-invertible, try holding a single parameter constant (usually "
                                "an epoch time) to see if it becomes invertible",
                                {'f', "fix-non-invertible"});
    parser.Parse();
    if (!inputFile) {
        std::cerr << parser << std::endl;
        exit(1);
    }

    // We do not transform parameters into the log space when evaluating results.
    const bool logParams = false;
    // The normal mode is minimization (negative log likelihood), but here we want
    // just the log likelihood.
    const bool maximization = true;

    NegLogLikelihoodCostFunctor cost(*inputFile, logParams, maximization, OBS_MODE_UNSPECIFIED, normCoals);
    if (!checkCostFunc(cost)) {
        exit(2);
    }

    std::vector<double> paramVector(cost.m_schema.totalParams());
    std::ifstream resultData(*inputFile);
    json prevResultJson = json::parse(resultData);
    cost.m_schema.fromJsonOutput(prevResultJson, paramVector.data());

    std::cerr << "MLE THETA = ";
    cost.m_schema.dumpParameters(paramVector.data(), paramVector.size());
    std::cerr << std::endl;

    if (cost.numCoalMatrices() < NUM_COAL_MATRICES_WARN_THRESHOLD) {
        std::cerr << "#########################" << std::endl;
        std::cerr << "WARNING: only found " << cost.numCoalMatrices() << " coalescence matrices, and these are what "
                  << "will be used to compute the jacobian. Your gradients may "
                     "be low quality."
                  << std::endl;
        std::cerr << "#########################" << std::endl;
    }

    // Construct a new copy of the JSON output, because we are going to add
    // statistics to it.
    json outputResult = cost.m_schema.toJsonOutput(paramVector.data(), prevResultJson.at("negLL"));
    VectorXd uncertaintySE = VectorXd::Zero(SIZE_T_TO_INDEX(paramVector.size()));

    size_t leftOut = std::numeric_limits<size_t>::max();
    auto [variabilityJ, sensitivityH] =
        getJandH(cost, paramVector, leftOut, /*bootstrap=*/true, /*fixNonInvertible=*/fixNonInvertible);

#if SINGLE_INVERSE_CALCULATION
    MatrixXd H_inverse = invert(sensitivityH);
    MatrixXd GIM_inverse = H_inverse * variabilityJ * H_inverse;
    TRACE_MATRIX(GIM_inverse, "GIM_inverse");
#else
    MatrixXd GIM = (sensitivityH * invert(variabilityJ)) * sensitivityH;
    TRACE_MATRIX(GIM, "GIM");
    MatrixXd GIM_inverse = invert(GIM);
#endif
    const VectorXd gimDiag = GIM_inverse.diagonal().array().sqrt();
    RELEASE_ASSERT(gimDiag.size() == paramVector.size() ||
                   gimDiag.size() == (paramVector.size() - 1) && leftOut < paramVector.size());

    for (size_t i = 0; i < gimDiag.size(); i++) {
        const double se = cost.m_schema.fromParam(gimDiag(i), i);
        if (i < leftOut) {
            uncertaintySE(i) = se;
        } else {
            uncertaintySE(i + 1) = cost.m_schema.fromParam(gimDiag(i), i);
        }
    }
    std::cerr << "Computed GIM-based uncertainty." << std::endl;

    std::vector<double> standardErrsGIM;
    std::vector<std::pair<double, double>> confIntervalsGIM;
    for (size_t i = 0; i < cost.m_schema.totalParams(); i++) {
        const double value = cost.m_schema.fromParam(paramVector.at(i), i);
        const double uncert = uncertaintySE(SIZE_T_TO_INDEX(i));
        const double vPlusUncert = value + uncert;
        standardErrsGIM.push_back(vPlusUncert - value);
        const double ciPlus = (value + 1.96 * uncert);
        const double ciMinus = (value - 1.96 * uncert);
        confIntervalsGIM.emplace_back(lclamp(ciMinus, 0), ciPlus);
    }

    cost.m_schema.addToParamOutput(outputResult, "gim_se", standardErrsGIM);
    cost.m_schema.addToParamOutput(outputResult, "gim_ci", confIntervalsGIM);

    std::cout << outputResult;
}

void selectCommand(args::Subparser& parser) {
    args::PositionalList<std::string> inputFiles(
        parser, "inputFiles", "Input file: two or more mrpast result files (JSON from the solver)");
    args::Flag bootstrap(
        parser,
        "bootstrap",
        "Compare all bootstrapped likelihood values from the corresponding CSV (from 'mrpast confidence --bootstrap'), "
        "to get a distribution of AIC values instead of a single value",
        {'b', "bootstrap"});
    args::Flag fixNonInvertible(parser,
                                "fixNonInvertible",
                                "When the H matrix is non-invertible, try holding a single parameter constant (usually "
                                "an epoch time) to see if it becomes invertible",
                                {'f', "fix-non-invertible"});
    args::Flag failOnNegTrace(parser,
                              "failOnNegTrace",
                              "When the trace (for adjusted AIC) is negative, fail the comparison",
                              {'n', "fail-bad-trace"});
    parser.Parse();
    if (!inputFiles || inputFiles->size() < 2) {
        std::cerr << parser << std::endl;
        exit(1);
    }
    std::vector<std::string> resultFiles = *inputFiles;

    // We do not transform parameters into the log space when evaluating results (derivatives would be log-scaled).
    const bool logParams = false;
    // The normal mode is minimization (negative log likelihood), but here we want
    // just the log likelihood.
    const bool maximization = true;

    // Validation of the inputs: they should all come from the same underlying data (we can't validate that), which
    // means that the sampled trees, number of coalescence matrices, etc., should all match.
    size_t numMatrices = 0;
    std::vector<std::string> prevHashes;
    for (size_t i = 0; i < resultFiles.size(); i++) {
        if (!ends_with(resultFiles[i], ".json")) {
            throw std::runtime_error("ERROR: Input files should have .json extension");
        }
        NegLogLikelihoodCostFunctor cost(resultFiles[i], logParams, maximization, OBS_MODE_UNSPECIFIED, normCoals);
        if (numMatrices == 0) {
            numMatrices = cost.numCoalMatrices();
        } else if (numMatrices != cost.numCoalMatrices()) {
            throw std::runtime_error("ERROR: Input files have different number of samples (coal. matrices)");
        }

        std::ifstream resultStream(resultFiles[i]);
        const json resultJson = json::parse(resultStream);
        // Validate that the hashes (unique combination of trees sampled) match for all bootstrap/MCMC samples.
        if (prevHashes.empty()) {
            prevHashes = resultJson[SAMPLING_HASHES_KEY].get<std::vector<std::string>>();
        } else if (prevHashes != resultJson[SAMPLING_HASHES_KEY].get<std::vector<std::string>>()) {
            throw std::runtime_error("ERROR: sampling_hashes mismatch between compared outputs");
        }
    }

    const size_t numSamples = bootstrap ? numMatrices : 1;

    json aicOutput;
    aicOutput["aic_values"] = json();
    std::vector<std::vector<double>> bySampleAIC(numSamples);
    std::vector<std::vector<double>> bySampleAIC_CL(numSamples);
    std::vector<std::vector<double>> bySampleCL(numSamples);
    for (size_t i = 0; i < resultFiles.size(); i++) {
        const auto& resultFile = resultFiles[i];
        std::cerr << "Calculating AIC for " << resultFile << std::endl;

        NegLogLikelihoodCostFunctor cost(resultFile, logParams, maximization, OBS_MODE_UNSPECIFIED, normCoals);
        std::vector<double> paramVector(cost.m_schema.totalParams());

        std::ifstream resultStream(resultFile);
        const json resultJson = json::parse(resultStream);
        cost.m_schema.fromJsonOutput(resultJson, paramVector.data());

        size_t leftOut = std::numeric_limits<size_t>::max();
        auto [variabilityJ, sensitivityH] =
            getJandH(cost, paramVector, leftOut, /*bootstrap=*/true, /*fixNonInvertible=*/fixNonInvertible);
        const double currentTrace = (variabilityJ * invert(sensitivityH)).trace();
        if (currentTrace <= 0) {
            std::cerr << "WARNING: Non-positive trace value of " << currentTrace << std::endl;
            if (failOnNegTrace) {
                throw std::runtime_error("Bad trace value: failing");
            }
        }

        std::vector<double> likelihoods;
        if (!bootstrap) {
            const double logLikelihood = cost(paramVector.data());
            likelihoods.push_back(logLikelihood);
        } else {
            // Load the likelihood values from the CSV (take the 0'th parameter)
            // The number of likelihood values should match the # of coal matrices in the JSON
            // Compute all metrics and add to a list
            std::stringstream ssCsvFile;
            ssCsvFile << resultFile.substr(0, resultFile.size() - 4) << "bootstrap.csv";
            std::ifstream csvFile(ssCsvFile.str());
            likelihoods = loadNegLL(csvFile);
            // We have negative log-likelihoods, convert to log-likelihoods
            for (size_t i = 0; i < likelihoods.size(); i++) {
                likelihoods[i] = -likelihoods[i];
            }
        }
        RELEASE_ASSERT(likelihoods.size() == numSamples);

        for (size_t i = 0; i < likelihoods.size(); i++) {
            const double numParams = (double)paramVector.size();
            const double logLikelihood = likelihoods[i];
            const double AIC = (2 * numParams) - (2 * logLikelihood);
            const double AIC_cl = (2 * currentTrace) - (2 * logLikelihood);

            json aicOutVal;
            aicOutVal["sample"] = i;
            aicOutVal["cL"] = logLikelihood;
            aicOutVal["AIC"] = AIC;
            aicOutVal["AIC_cl"] = AIC_cl;
            aicOutVal["file"] = resultFile;
            aicOutVal["param_count"] = paramVector.size();
            aicOutVal["param_trace"] = currentTrace;
            aicOutput["aic_values"].push_back(aicOutVal);

            bySampleAIC.at(i).emplace_back(AIC);
            bySampleAIC_CL.at(i).emplace_back(AIC_cl);
            bySampleCL.at(i).emplace_back(logLikelihood);
        }
    }

    std::vector<size_t> winsAIC(resultFiles.size());
    std::vector<size_t> winsAIC_CL(resultFiles.size());
    std::vector<size_t> winsCL(resultFiles.size());
    for (size_t i = 0; i < numSamples; i++) {
        size_t bestIndexAIC = 0;
        double bestAIC = std::numeric_limits<double>::max();
        size_t bestIndexAIC_CL = 0;
        double bestAIC_CL = std::numeric_limits<double>::max();
        size_t bestIndexCL = 0;
        double bestCL = std::numeric_limits<double>::lowest();
        for (size_t j = 0; j < resultFiles.size(); j++) {
            const double AIC = bySampleAIC.at(i).at(j);
            const double AIC_CL = bySampleAIC_CL.at(i).at(j);
            const double CL = bySampleCL.at(i).at(j);
            if (AIC < bestAIC) {
                bestAIC = AIC;
                bestIndexAIC = j;
            }
            if (AIC_CL < bestAIC_CL) {
                bestAIC_CL = AIC_CL;
                bestIndexAIC_CL = j;
            }
            if (CL > bestCL) {
                bestCL = CL;
                bestIndexCL = j;
            }
        }
        winsAIC.at(bestIndexAIC)++;
        winsAIC_CL.at(bestIndexAIC_CL)++;
        winsCL.at(bestIndexCL)++;
    }

    aicOutput["wins"] = json();
    for (size_t i = 0; i < resultFiles.size(); i++) {
        std::cerr << "Filename: " << resultFiles[i] << ", Wins(AIC): " << winsAIC.at(i)
                  << ", Wins(AIC_CL): " << winsAIC_CL.at(i) << ", Wins(CL): " << winsCL.at(i) << std::endl;
        json winRow;
        winRow["AIC"] = winsAIC.at(i);
        winRow["AIC_CL"] = winsAIC_CL.at(i);
        winRow["CL"] = winsCL.at(i);
        winRow["file"] = resultFiles.at(i);
        aicOutput["wins"].push_back(winRow);
    }

    std::cout << aicOutput << std::endl;
}

void contourCommand(args::Subparser& parser) {
    args::Positional<std::string> inputFile(
        parser, "inputFile", "Input file: a JSON file as emitted by mrpast's solver");
    args::Positional<size_t> param1(parser, "param1", "Index of the first parameter to use.");
    args::Positional<size_t> param2(parser, "param2", "Index of the second parameter to use.");
    parser.Parse();
    if (!inputFile || !param1 || !param2) {
        std::cerr << parser << std::endl;
        exit(1);
    }

    // Do 100 values between lower-upper in each dimension (parameter)
    const size_t numDimValues = 100;

    // We do not transform parameters into the log space when evaluating results.
    const bool logParams = false;
    // The normal mode is minimization (negative log likelihood), but here we want
    // just the log likelihood.
    const bool maximization = true;

    NegLogLikelihoodCostFunctor cost(*inputFile, logParams, maximization);
    std::vector<double> paramVector(cost.m_schema.totalParams());
    std::ifstream resultData(*inputFile);
    json prevResultJson = json::parse(resultData);
    cost.m_schema.fromJsonOutput(prevResultJson, paramVector.data());

    const size_t p1Idx = *param1;
    size_t p2Idx = *param2;
    double p1Lower = 0.0;
    double p1Upper = 0.0;
    double p2Lower = 0.0;
    double p2Upper = 0.0;
    cost.m_schema.getBounds(p1Idx, p1Lower, p1Upper);
    cost.m_schema.getBounds(p2Idx, p2Lower, p2Upper);
    std::cerr << "Varying parameter " << p1Idx << " between " << p1Lower << " and " << p1Upper << std::endl;
    std::cerr << "Varying parameter " << p2Idx << " between " << p2Lower << " and " << p2Upper << std::endl;

    std::vector<double> X;
    std::vector<double> Y;
    std::vector<std::vector<double>> Z;
    const double p1Step = (p1Upper - p1Lower) / (double)numDimValues;
    const double p2Step = (p2Upper - p2Lower) / (double)numDimValues;
    for (double p1Val = p1Lower; p1Val < p1Upper; p1Val += p1Step) {
        X.push_back(p1Val);
        Z.emplace_back();
        paramVector.at(p1Idx) = p1Val;
        for (double p2Val = p2Lower; p2Val < p2Upper; p2Val += p2Step) {
            if (p1Val == p1Lower) {
                Y.push_back(p2Val);
            }
            paramVector.at(p2Idx) = p2Val;
            const double z = cost(paramVector.data());
            Z.back().push_back(z);
        }
    }
    json result;
    result["X"] = X;
    result["Y"] = Y;
    result["Z"] = Z;
    std::cout << result;
}

void emitCoalCommand(args::Subparser& parser) {
    args::Positional<std::string> inputFile(
        parser, "inputFile", "Input file: a JSON file as emitted by mrpast's solver");
    args::Flag useGroundTruth(parser,
                              "useGroundTruth",
                              "Use the ground_truth from the input file, instead of the final value",
                              {'g', "ground-truth"});
    parser.Parse();
    if (!inputFile) {
        std::cerr << parser << std::endl;
        exit(1);
    }

    NegLogLikelihoodCostFunctor cost(*inputFile, true);
    std::vector<double> paramVector(cost.m_schema.totalParams());
    std::ifstream resultData(*inputFile);
    const json prevResultJson = json::parse(resultData);

    std::string field = "final";
    if (useGroundTruth) {
        field = "ground_truth";
    }
    cost.m_schema.fromJsonOutput(prevResultJson, paramVector.data(), field);

    cost.outputTheoreticalCoalMatrix(paramVector.data(), std::cout);
}

int main(int argc, char** argv) {
    std::cerr << std::setprecision(100);

    const std::string CMD_INTERVALS = "intervals";
    const std::string CMD_SELECT = "select";
    const std::string CMD_CONTOUR = "contour";
    const std::string CMD_EMIT_COAL = "emit-coal";

    args::ArgumentParser parser("Evaluate the results of the mrpast solver");
    args::Group commands(parser, "commands");
    args::Command intervals(commands,
                            CMD_INTERVALS,
                            "Construct confidence intervals for parameters using the GIM-based method.",
                            &intervalsCommand);
    args::Command select(
        commands, CMD_SELECT, "Run AIC-based model selection on a list of solver results.", &selectCommand);
    args::Command contour(
        commands, CMD_CONTOUR, "Construct 3D contour data along a grid for the given solver result.", &contourCommand);
    args::Command emitCoal(
        commands,
        CMD_EMIT_COAL,
        "Emit a new output containing the theoretical coalescence matrix for the given solver result.",
        &emitCoalCommand);
    args::GlobalOptions globals(parser, arguments);

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help&) {
        std::cout << parser;
    } catch (args::Error& e) {
        std::cerr << e.what() << std::endl << parser;
        return 1;
    }
    return 0;
}