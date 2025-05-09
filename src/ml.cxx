#ifndef GUARD_ML_H
#define GUARD_ML_H

#include "../include/ml.hxx"
#include "../include/defaults.hxx"
#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/OnnxSessionManager.hxx"
#include "../include/utility/utility.hxx"
#include "../include/vectoroperations.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>

#include "TInterpreter.h"
#include "TMVA/RModel.hxx"
#include "TMVA/RModelParser_ONNX.hxx"
#include "TSystem.h"
#include <assert.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <onnxruntime_cxx_api.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>

using json = nlohmann::json;

namespace ml {
/**
 * @brief Function to perform a standard transformation of input variables for
 * NN evaluation.
 *
 * @param df The input dataframe
 * @param inputname name of the variable which should be transformed
 * @param outputname name of the output column
 * @param param_file path to a json file with a dictionary of mean and std
 * values
 * @param var_type variable data type for correct processing e.g. "i" for
 * integer or "f" for float
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode StandardTransformer(ROOT::RDF::RNode df,
                                     const std::string &inputname,
                                     const std::string &outputname,
                                     const std::string &param_file,
                                     const std::string &var_type) {
    // read params from file
    Logger::get("StandardTransformer")->debug("reading file {}", param_file);
    Logger::get("StandardTransformer")
        ->warn("Not using CorrectionManager is deprecated, this function will "
               "be removed in a future release, switch to the new function "
               "using CorrectionManager");
    std::string replace_str = std::string("EVTID");
    std::string odd_file_path =
        std::string(param_file)
            .replace(param_file.find(replace_str), replace_str.length(),
                     std::string("odd"));
    std::string even_file_path =
        std::string(param_file)
            .replace(param_file.find(replace_str), replace_str.length(),
                     std::string("even"));

    std::ifstream odd_file(odd_file_path);
    json odd_info = json::parse(odd_file);
    std::ifstream even_file(even_file_path);
    json even_info = json::parse(even_file);

    // odd or even files mean that they are trained on odd or even events, so it
    // has to be applied on the opposite
    auto transform_int = [odd_info, even_info,
                          inputname](const unsigned long long event_id,
                                     const int input_var) {
        float shifted = -10;
        if (int(event_id) % 2 == 0) {
            shifted = (float(input_var) - float(odd_info[inputname]["mean"])) /
                      float(odd_info[inputname]["std"]);
        } else if (int(event_id) % 2 == 1) {
            shifted = (float(input_var) - float(even_info[inputname]["mean"])) /
                      float(even_info[inputname]["std"]);
        }
        Logger::get("StandardTransformer")
            ->debug("transforming var {} from {} to {}", inputname, input_var,
                    shifted);
        return shifted;
    };
    auto transform_float = [odd_info, even_info,
                            inputname](const unsigned long long event_id,
                                       const float input_var) {
        float shifted = -10;
        if (int(event_id) % 2 == 0) {
            shifted = (float(input_var) - float(odd_info[inputname]["mean"])) /
                      float(odd_info[inputname]["std"]);
        } else if (int(event_id) % 2 == 1) {
            shifted = (float(input_var) - float(even_info[inputname]["mean"])) /
                      float(even_info[inputname]["std"]);
        }
        Logger::get("StandardTransformer")
            ->debug("transforming var {} from {} to {}", inputname, input_var,
                    shifted);
        return shifted;
    };

    const std::string event_id = std::string("event");
    if (var_type.rfind("i", 0) == 0) {
        auto df1 = df.Define(outputname, transform_int, {event_id, inputname});
        return df1;
    } else if (var_type.rfind("f", 0) == 0) {
        auto df1 =
            df.Define(outputname, transform_float, {event_id, inputname});
        return df1;
    } else {
        Logger::get("StandardTransformer")
            ->debug("transformation failed due to wrong variable type: {}",
                    var_type);
        return df;
    }
}
/**
 * @brief Function to perform a standard transformation of input variables for
 * NN evaluation.
 *
 * @param df The input dataframe
 * @param correctionManager correction manager
 * @param inputname name of the variable which should be transformed
 * @param outputname name of the output column
 * @param param_file path to a json file with a dictionary of mean and std
 * values
 * @param var_type variable data type for correct processing e.g. "i" for
 * integer or "f" for float
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode
StandardTransformer(ROOT::RDF::RNode df,
                    correctionManager::CorrectionManager &correctionManager,
                    const std::string &inputname, const std::string &outputname,
                    const std::string &param_file,
                    const std::string &var_type) {
    // read params from file
    Logger::get("StandardTransformer")->debug("reading file {}", param_file);
    std::string replace_str = std::string("EVTID");
    std::string odd_file_path =
        std::string(param_file)
            .replace(param_file.find(replace_str), replace_str.length(),
                     std::string("odd"));
    std::string even_file_path =
        std::string(param_file)
            .replace(param_file.find(replace_str), replace_str.length(),
                     std::string("even"));

    json odd_info = *correctionManager.loadjson(odd_file_path);
    json even_info = *correctionManager.loadjson(even_file_path);

    // odd or even files mean that they are trained on odd or even events, so it
    // has to be applied on the opposite
    auto transform_int = [odd_info, even_info,
                          inputname](const unsigned long long event_id,
                                     const int input_var) {
        float shifted = -10;
        if (int(event_id) % 2 == 0) {
            shifted = (float(input_var) - float(odd_info[inputname]["mean"])) /
                      float(odd_info[inputname]["std"]);
        } else if (int(event_id) % 2 == 1) {
            shifted = (float(input_var) - float(even_info[inputname]["mean"])) /
                      float(even_info[inputname]["std"]);
        }
        Logger::get("StandardTransformer")
            ->debug("transforming var {} from {} to {}", inputname, input_var,
                    shifted);
        return shifted;
    };
    auto transform_float = [odd_info, even_info,
                            inputname](const unsigned long long event_id,
                                       const float input_var) {
        float shifted = -10;
        if (int(event_id) % 2 == 0) {
            shifted = (float(input_var) - float(odd_info[inputname]["mean"])) /
                      float(odd_info[inputname]["std"]);
        } else if (int(event_id) % 2 == 1) {
            shifted = (float(input_var) - float(even_info[inputname]["mean"])) /
                      float(even_info[inputname]["std"]);
        }
        Logger::get("StandardTransformer")
            ->debug("transforming var {} from {} to {}", inputname, input_var,
                    shifted);
        return shifted;
    };

    const std::string event_id = std::string("event");
    if (var_type.rfind("i", 0) == 0) {
        auto df1 = df.Define(outputname, transform_int, {event_id, inputname});
        return df1;
    } else if (var_type.rfind("f", 0) == 0) {
        auto df1 =
            df.Define(outputname, transform_float, {event_id, inputname});
        return df1;
    } else {
        Logger::get("StandardTransformer")
            ->debug("transformation failed due to wrong variable type: {}",
                    var_type);
        return df;
    }
}

} // end namespace ml
#endif /* GUARD_ML_H */