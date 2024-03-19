#ifndef GUARD_ML_H
#define GUARD_ML_H

#include "../include/basefunctions.hxx"
#include "../include/defaults.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/vectoroperations.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>

#include "TInterpreter.h"
#include "TMVA/RModel.hxx"
#include "TMVA/RModelParser_ONNX.hxx"
#include "TSystem.h"
#include <onnxruntime_cxx_api.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

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
 * @brief Function to produce columns with mass values for the evaluation of
 * parametric neural networks for the NMSSM analysis.
 *
 * @param df The input dataframe
 * @param outputname name of the output column
 * @param mass_paramfile path to a json file with a dictionary for mass value
 * transformations
 * @param mass_parameter actual mass value to be transformed
 * @param boson one of the two boson with unknown masses "X" or "Y"
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode DefineMassColumns(ROOT::RDF::RNode df,
                                   const std::string &outputname,
                                   const std::string &mass_paramfile,
                                   const std::string &mass_parameter,
                                   const std::string &boson) {
    std::ifstream mass_file(mass_paramfile);
    json mass_transform = json::parse(mass_file);

    std::string mass_str;
    if (boson == "X") {
        mass_str = "massX";
    } else if (boson == "Y") {
        mass_str = "massY";
    } else {
        Logger::get("DefineMassColumns")->debug("boson {} not defined", boson);
    }

    auto get_mass = [mass_transform, mass_str, mass_parameter]() {
        float mass = float(mass_transform[mass_str][mass_parameter]);
        return mass;
    };
    Logger::get("DefineMassColumns")
        ->debug("output column produced for {}", outputname);
    auto df1 = df.Define(outputname, get_mass, {});
    return df1;
}

namespace sofie {
/**
 * @brief Function which includes the model header file and defines a functor
 * for the evaluation in RDFs.
 *
 * @param headerModelFile path to a header file with the model
 * @param ninputs number of input variable to the model
 * @param nslots  number of data processing slots
 */
void CompileModelForRDF(const std::string &headerModelFile,
                        unsigned int ninputs, unsigned int nslots = 0) {

    std::string modelName =
        headerModelFile.substr(0, headerModelFile.find(".hxx"));
    std::string cmd =
        std::string("#include \"") + headerModelFile + std::string("\"");
    auto ret = gInterpreter->Declare(cmd.c_str());
    if (!ret) {
        throw std::runtime_error("Error compiling: " + cmd);
    }
    Logger::get("SOFIEGenerator")->debug("compiled: {}", cmd);

    cmd = std::string("#include \"data/ml/SofieHelpers.hxx\"");
    ret = gInterpreter->Declare(cmd.c_str());
    if (!ret) {
        throw std::runtime_error("Error compiling: " + cmd);
    }
    Logger::get("SOFIEGenerator")
        ->debug("Compiled custom SOFIE helpers: {}", cmd);

    cmd = "auto sofie_functor_" + modelName + " = SofieFunctor<" +
          std::to_string(ninputs) + ",TMVA_SOFIE_" + modelName + "::Session>(" +
          std::to_string(nslots) + ");";
    ret = gInterpreter->Declare(cmd.c_str());
    if (!ret) {
        throw std::runtime_error("Error compiling : " + cmd);
    }
    Logger::get("SOFIEGenerator")->debug("compiled: {}", cmd);
    Logger::get("SOFIEGenerator")->debug("Model is ready to be evaluated");

    return;
}
/**
 * @brief Function which generates a header file of a model from an onnx file.
 *
 * @param input_vec vector of input variables to the model
 * @param model_file path to the onnx file with the saved model
 */
std::string SOFIEGenerator(std::vector<std::string> input_vec,
                           const std::string &model_file) {

    std::string modelFileName = std::filesystem::path(model_file).filename();
    std::string::size_type const p(modelFileName.find_last_of('.'));
    std::string modelName = modelFileName.substr(0, p);

    std::string modelHeaderFile = modelName + std::string(".hxx");

    if (!std::filesystem::exists(modelHeaderFile)) {
        Logger::get("NMSSMEvaluate")
            ->debug("loading model file {}", model_file);
        bool verboseParser = true;
        TMVA::Experimental::SOFIE::RModelParser_ONNX parser;
        TMVA::Experimental::SOFIE::RModel model =
            parser.Parse(model_file, verboseParser);
        Logger::get("NMSSMEvaluate")->debug("finished loading model");

        Logger::get("SOFIEGenerator")->debug("generating model code");
        model.Generate();
        Logger::get("SOFIEGenerator")
            ->debug("dumping model code to {}", modelHeaderFile);
        model.OutputGenerated(modelHeaderFile);

        Logger::get("SOFIEGenerator")->debug("compiling model code");
        CompileModelForRDF(modelHeaderFile, input_vec.size());
    } else {
        Logger::get("SOFIEGenerator")
            ->debug("model already compiled, skipping");
    }

    std::string sofie_func_str = "sofie_functor_" + modelName + "(rdfslot_, ";
    sofie_func_str += input_vec[0];
    for (unsigned i = 1; i < input_vec.size(); i++) {
        sofie_func_str += ", " + input_vec[i];
    }
    sofie_func_str += ")";

    return sofie_func_str;
}
/**
 * @brief Function for the evaluation of parametric neural networks for the
 * NMSSM analysis.
 *
 * @param df The input dataframe
 * @param input_vec vector of input variables to the model
 * @param outputname_vector name of the output column with the full prediction
 * vector
 * @param outputname_class name of the output column with the class with the
 * highest value
 * @param outputname_max_value name of the output column with the highest value
 * @param model_file path to the onnx file with the saved model
 * @param massX_parameter mass value of X the model should be evaluated with
 * @param massY_parameter mass value of Y the model should be evaluated with
 * @return a new dataframe containing the new columns
 */
ROOT::RDF::RNode NMSSMEvaluate_ONNX(
    ROOT::RDF::RNode df, const std::vector<std::string> &input_vec,
    const std::string &outputname_vector, const std::string &outputname_class,
    const std::string &outputname_max_value, const std::string &model_file,
    const std::string &massX_parameter, const std::string &massY_parameter) {

    using namespace TMVA::Experimental::SOFIE;
    Logger::get("NMSSMEvaluate")
        ->debug("mass hypothesis: mX={}, mY={}", massX_parameter,
                massY_parameter);

    std::vector<std::string> inputs = input_vec;
    inputs.push_back(std::string("massX_") + massX_parameter);
    inputs.push_back(std::string("massY_") + massY_parameter);

    std::string replace_str = std::string("EVTID");
    std::string model_file_path_odd =
        std::string(model_file)
            .replace(model_file.find(replace_str), replace_str.length(),
                     std::string("odd"));
    std::string model_file_path_even =
        std::string(model_file)
            .replace(model_file.find(replace_str), replace_str.length(),
                     std::string("even"));

    auto eval_func_odd = SOFIEGenerator(inputs, model_file_path_odd);
    auto eval_func_even = SOFIEGenerator(inputs, model_file_path_even);
    Logger::get("NMSSMEvaluate")
        ->debug("evaluating model code with: {} and {}", eval_func_even,
                eval_func_odd);

    std::string eval_func =
        "(event%2==0) ? " + eval_func_odd + " : " + eval_func_even + ";";

    auto df1 = df.Define(outputname_vector, eval_func);
    auto df2 = df1.Define(
        outputname_class,
        [](const std::vector<float> &prediction) {
            int max_idx = std::distance(
                prediction.begin(),
                std::max_element(prediction.begin(), prediction.end()));
            return max_idx;
        },
        {outputname_vector});
    auto df3 = df2.Define(
        outputname_max_value,
        [](const std::vector<float> &prediction, const int &max_idx) {
            float max_value = prediction[max_idx];
            return max_value;
        },
        {outputname_vector, outputname_class});

    return df3;
}

} // end namespace sofie
} // end namespace ml
#endif /* GUARD_ML_H */