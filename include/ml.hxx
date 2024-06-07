#ifndef GUARD_ML_H
#define GUARD_ML_H

#include "../include/utility/OnnxSessionManager.hxx"
#include "TMVA/RModel.hxx"
#include "TMVA/RModelParser_ONNX.hxx"
#include "utility/utility.hxx"
#include <cstddef>

namespace ml {

ROOT::RDF::RNode StandardTransformer(ROOT::RDF::RNode df,
                                     const std::string &inputname,
                                     const std::string &outputname,
                                     const std::string &paramfile,
                                     const std::string &var_type);

/// Generic Function to evaluate an ONNX model using the ONNX Runtime
/// Due to unknowns reasons, this function must be implemented inline in the
/// header file, otherwise the linker will complain about undefined references.
/// Moving the implementation to the source file will result in a linker error.
/// Why, I don't know...
/// This generic implementation currenty supports only NNs with one input tensor
/// and one output tensor
///
/// \param df the dataframe to add the quantity to
/// \param OnnxSessionManager The OnnxSessionManager object to handle the
/// runtime session. By default this is called onnxSessionManager and created in
/// the main function
/// \param outputname Name of the output column
/// \param model_file_path Path to the ONNX model file
/// \param input_vec Vector of input variable names,
/// the order of the variables must match the order of the input
/// nodes in the ONNX model
///
/// \returns a dataframe with the filter applied

template <std::size_t nParameters>
inline ROOT::RDF::RNode GenericOnnxEvaluator(
    ROOT::RDF::RNode df, OnnxSessionManager &onnxSessionManager,
    const std::string &outputname_vector, const std::string &outputname_class,
    const std::string &outputname_max_value, const std::string &model_file_path,
    const std::vector<std::string> &input_vec) {

    std::vector<std::string> InputList;
    for (auto i = 0; i < input_vec.size(); i++) {
        InputList.push_back(std::string(input_vec[i]));
    }

    // print content of InputList
    for (auto i = 0; i < InputList.size(); ++i) {
        Logger::get("OnnxEvaluate")
            ->debug("input: {} ( {} / {} )", InputList[i], i + 1, nParameters);
    }

    if (nParameters != InputList.size()) {
        Logger::get("OnnxEvaluate")
            ->error("Number of input parameters does not match the number of "
                    "input variables: {} vs {}",
                    nParameters, InputList.size());
        throw std::runtime_error("Number of input parameters does not match");
    }

    std::string replace_str = std::string("EVTID");
    std::string model_file_path_odd =
        std::string(model_file_path)
            .replace(model_file_path.find(replace_str), replace_str.length(),
                     std::string("odd"));
    std::string model_file_path_even =
        std::string(model_file_path)
            .replace(model_file_path.find(replace_str), replace_str.length(),
                     std::string("even"));

    // Load the model and create InferenceSession
    std::vector<int64_t> input_node_dims;
    std::vector<int64_t> output_node_dims;
    int num_input_nodes;
    int num_output_nodes;
    Ort::AllocatorWithDefaultOptions allocator;

    auto session_odd = onnxSessionManager.getSession(model_file_path_odd);
    auto session_even = onnxSessionManager.getSession(model_file_path_even);

    // for even/odd splitting
    auto col_names = df.GetColumnNames();
    bool found_col = (std::find(col_names.begin(), col_names.end(),
                                "event_id") != col_names.end());
    auto df1 = df;
    if (!found_col) {
        df1 = df.Define("event_id",
                        [](unsigned long long event) { return float(event); },
                        {"event"});
    }

    InputList.push_back(std::string("event_id"));
    // should be the same for even and odd
    onnxhelper::prepare_model(session_odd, allocator, input_node_dims, output_node_dims,
                              num_input_nodes, num_output_nodes);

    
    auto NNEvaluator = [session_odd, session_even, allocator, input_node_dims,
                        output_node_dims, num_input_nodes,
                        num_output_nodes, InputList](std::vector<float> inputs) {
        int event_id = (int)inputs.back();
        inputs.pop_back();

        TStopwatch timer;
        timer.Start();
        std::vector<float> output;
        if (event_id % 2 == 0) {
            output = onnxhelper::run_interference(
                session_odd, allocator, inputs, input_node_dims,
                output_node_dims, num_input_nodes, num_output_nodes);
        } else {
            output = onnxhelper::run_interference(
                session_even, allocator, inputs, input_node_dims,
                output_node_dims, num_input_nodes, num_output_nodes);
        }

        timer.Stop();
        Logger::get("OnnxEvaluate")
            ->debug("Inference time: {} mus", timer.RealTime() * 1000 * 1000);
        return output;
    };
    auto df2 = df1.Define(
        outputname_vector,
        utility::PassAsVec<(nParameters + 1), float>(NNEvaluator), InputList);
    auto df3 = df2.Define(
        outputname_class,
        [](const std::vector<float> &prediction) {
            int max_idx = std::distance(
                prediction.begin(),
                std::max_element(prediction.begin(), prediction.end()));
            return max_idx;
        },
        {outputname_vector});
    auto df4 = df3.Define(
        outputname_max_value,
        [](const std::vector<float> &prediction, const int &max_idx) {
            float max_value = prediction[max_idx];
            return max_value;
        },
        {outputname_vector, outputname_class});

    return df4;
}
} // end namespace ml
#endif /* GUARD_ML_H */