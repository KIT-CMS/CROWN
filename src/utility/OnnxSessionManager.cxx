#include "../../include/utility/Logger.hxx"
#include "../../include/utility/utility.hxx"
#include <memory>
#include <numeric>
#include <onnxruntime_cxx_api.h>
#include <string>
#include <unordered_map>

namespace onnxhelper {

template <typename T> T vectorProduct(const std::vector<T> &v) {
    return accumulate(v.begin(), v.end(), 1, std::multiplies<T>());
}

std::vector<float> run_interference(Ort::Session *session,
                                    Ort::AllocatorWithDefaultOptions allocator,
                                    std::vector<float> &evt_input,
                                    std::vector<int64_t> input_node_dims,
                                    std::vector<int64_t> output_node_dims,
                                    const int num_input_nodes,
                                    const int num_output_nodes) {

    Ort::MemoryInfo memory_info = Ort::MemoryInfo::CreateCpu(
        OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault);

    std::vector<Ort::Value> inputTensors;
    std::vector<Ort::Value> outputTensors;

    size_t inputTensorSize = vectorProduct(input_node_dims);
    size_t outputTensorSize = vectorProduct(output_node_dims);

    for (auto i = 0; i < input_node_dims.size(); ++i) {
        Logger::get("OnnxInterference")
            ->debug("input_node_dims[{}]: {}", i, input_node_dims[i]);
    }
    for (auto i = 0; i < output_node_dims.size(); ++i) {
        Logger::get("OnnxInterference")
            ->debug("output_node_dims[{}]: {}", i, output_node_dims[i]);
    }
    Logger::get("OnnxInterference")
        ->debug("num_input_nodes: {}", num_input_nodes);
    Logger::get("OnnxInterference")
        ->debug("num_output_nodes: {}", num_output_nodes);
    Logger::get("OnnxInterference")
        ->debug("inputTensorSize: {}", inputTensorSize);
    Logger::get("OnnxInterference")
        ->debug("outputTensorSize: {}", outputTensorSize);

    std::vector<float> outputTensorValues(outputTensorSize);

    inputTensors.push_back(Ort::Value::CreateTensor<float>(
        memory_info, evt_input.data(), inputTensorSize, input_node_dims.data(),
        input_node_dims.size()));
    outputTensors.push_back(Ort::Value::CreateTensor<float>(
        memory_info, outputTensorValues.data(), outputTensorSize,
        output_node_dims.data(), output_node_dims.size()));
    assert(inputTensors[0].IsTensor());
    // print content of inputTensors
    Logger::get("OnnxInterference")->debug("Checking input tensors");
    const float *input_ptr = inputTensors[0].GetTensorMutableData<float>();
    for (auto i = 0; i < evt_input.size(); ++i) {
        Logger::get("OnnxInterference")->debug("input: {}", input_ptr[i]);
    }
    // convert input_node_names and output_node_names to const char*
    std::vector<const char *> input_node_names_cstr;
    std::vector<const char *> output_node_names_cstr;

    auto input_name_tmp = session->GetInputNameAllocated(0, allocator);
    auto output_name_tmp = session->GetOutputNameAllocated(0, allocator);
    std::string input_name_str(input_name_tmp.get());
    std::string output_name_str(output_name_tmp.get());
    input_node_names_cstr.push_back(input_name_str.c_str());
    output_node_names_cstr.push_back(output_name_str.c_str());

    Logger::get("OnnxInterference")->debug("Input name: {} ", input_name_str);
    Logger::get("OnnxInterference")->debug("Output name: {} ", output_name_str);

    session->Run(Ort::RunOptions{nullptr}, input_node_names_cstr.data(),
                 inputTensors.data(), input_node_names_cstr.size(),
                 output_node_names_cstr.data(), outputTensors.data(),
                 output_node_names_cstr.size());

    const float *output_ptr = outputTensors[0].GetTensorMutableData<float>();
    std::vector<float> output;
    for (size_t i = 0; i < outputTensorSize; ++i) {
        output.push_back(output_ptr[i]);
        Logger::get("OnnxInterference")->debug("output: {}", output[i]);
    }
    return output;
}

void prepare_model(Ort::Session *session,
                   Ort::AllocatorWithDefaultOptions allocator,
                   std::vector<int64_t> &input_node_dims,
                   std::vector<int64_t> &output_node_dims, int &num_input_nodes,
                   int &num_output_nodes) {

    num_input_nodes = session->GetInputCount();
    num_output_nodes = session->GetOutputCount();

    if (num_output_nodes != 1 or num_input_nodes != 1) {
        Logger::get("OnnxPrepareModel")
            ->error("Number of input and output nodes is not 1 Input: {} "
                    "Output: {}",
                    num_input_nodes, num_output_nodes);
        throw std::runtime_error("Number of nodes is not 1");
    }

    Logger::get("OnnxPrepareModel")
        ->debug("Number of input nodes: {}", num_input_nodes);
    Logger::get("OnnxPrepareModel")
        ->debug("Number of output nodes: {}", num_output_nodes);

    auto input_type_info = session->GetInputTypeInfo(0);
    auto input_tensor_info = input_type_info.GetTensorTypeAndShapeInfo();
    input_node_dims = input_tensor_info.GetShape();

    auto output_type_info = session->GetOutputTypeInfo(0);
    auto output_tensor_info = output_type_info.GetTensorTypeAndShapeInfo();
    output_node_dims = output_tensor_info.GetShape();
}
} // namespace onnxhelper