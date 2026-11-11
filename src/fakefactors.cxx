#ifndef GUARDCROWNFAKEFACTORS_CXX
#define GUARDCROWNFAKEFACTORS_CXX

#include "../include/fakefactors.hxx"
#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "correction.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace fakefactors {
namespace generic {

/**
 * @brief Constructor of the helper that maps the packed input column onto the
 * inputs of a correction. Each input of the correction is looked up by name in
 * the layout; the inputs "process" and "syst" are filled at evaluation time.
 *
 * @param correction_name name of the correction, used in the error messages
 * @param inputs inputs declared by the correction
 * @param layout names of the quantities in the packed input column, in order
 * @param proc process name passed to the "process" input of the correction
 * @throws std::runtime_error if an input is not in the layout or the
 * correction has no "syst" input
 */
ArgBuilder::ArgBuilder(const std::string &correction_name,
                       const std::vector<correction::Variable> &inputs,
                       const std::vector<std::string> &layout,
                       const std::string &proc)
    : process(proc) {
    bool has_syst = false;
    for (const auto &var : inputs) {
        const std::string name = var.name();
        if (name == "process") {
            slots.push_back({Kind::Process, 0, var.type()});
        } else if (name == "syst") {
            syst_pos = slots.size();
            has_syst = true;
            slots.push_back({Kind::Syst, 0, var.type()});
        } else {
            auto it = std::find(layout.begin(), layout.end(), name);
            if (it == layout.end()) {
                throw std::runtime_error("fakefactors::generic: input '" +
                                         name + "' of correction '" +
                                         correction_name +
                                         "' is not in the input layout");
            }
            slots.push_back({Kind::Value,
                             static_cast<size_t>(it - layout.begin()),
                             var.type()});
        }
    }
    if (!has_syst) {
        throw std::runtime_error("fakefactors::generic: correction '" +
                                 correction_name + "' has no 'syst' input");
    }
}

/**
 * @brief Builds the arguments for the evaluation of the correction for one
 * event. Integer inputs of the correction are cast to int, all other ones to
 * double.
 *
 * @param values the packed input column of the event
 * @param syst name of the variation, or nominal, for the "syst" input
 * @return the arguments in the order of the inputs of the correction
 */
std::vector<correction::Variable::Type>
ArgBuilder::build(const std::vector<float> &values,
                  const std::string &syst) const {
    std::vector<correction::Variable::Type> args;
    args.reserve(slots.size());
    for (const auto &slot : slots) {
        switch (slot.kind) {
        case Kind::Process:
            args.emplace_back(process);
            break;
        case Kind::Syst:
            args.emplace_back(syst);
            break;
        case Kind::Value:
            if (slot.type == correction::Variable::VarType::integer) {
                args.emplace_back(static_cast<int>(values.at(slot.index)));
            } else {
                args.emplace_back(static_cast<double>(values.at(slot.index)));
            }
        }
    }
    return args;
}

/**
 * @brief Constructor of the handler for the compound non-closure correction of
 * a process.
 *
 * @param p part of the variation names before the variable, e.g.
 * "CMS_fake_t_QCD_non_closure_"
 * @param vars variables stacked in the compound correction
 * @param c the compound correction
 * @param pos position of the "syst" input in the arguments of the correction
 */
NonClosureHandler::NonClosureHandler(const std::string &p,
                                     const std::vector<std::string> &vars,
                                     const correction::CompoundCorrection *c,
                                     size_t pos)
    : prefix(p), compound(c), syst_pos(pos), variables(vars) {}

/**
 * @brief Evaluates the compound non-closure correction. A coarse variation
 * (prefix followed by non_closure_Corr...) is the nominal value shifted by the
 * quadratic sum over the variables of the difference of each per variable
 * variation to the nominal; Up adds and Down subtracts it. Any other variation
 * (nominal, per variable) is evaluated directly.
 *
 * @param systematic name of the variation or nominal
 * @param args the arguments of the event built by the ArgBuilder of the
 * correction
 * @return the value of the correction
 */
float NonClosureHandler::evaluate(
    const std::string &systematic,
    std::vector<correction::Variable::Type> args) const {
    args[syst_pos] = "nominal";
    float nominal_value = compound->evaluate(args);

    bool is_target_variation =
        (systematic != "nominal") && (systematic.find(prefix) == 0) &&
        (systematic.find("non_closure_Corr") != std::string::npos);

    if (!is_target_variation) {
        if (systematic == "nominal") {
            return nominal_value;
        }
        args[syst_pos] = systematic;
        return compound->evaluate(args);
    }

    std::string coarse_correction = systematic.substr(prefix.length());

    double sum_sq_diff = 0.0;
    for (const auto &variable : variables) {
        args[syst_pos] = prefix + variable + "_" + coarse_correction;
        float diff =
            nominal_value - static_cast<float>(compound->evaluate(args));
        sum_sq_diff += (diff * diff);
    }

    float total_uncertainty = std::sqrt(sum_sq_diff);
    bool is_up = (systematic.compare(systematic.length() - 2, 2, "Up") == 0);

    return is_up ? (nominal_value + total_uncertainty)
                 : (nominal_value - total_uncertainty);
}

/**
 * @brief Packs columns into one std::vector<float> column, so that the terms
 * read their inputs by position. Each column is cast to float, whatever its
 * type.
 *
 * @param df the input dataframe
 * @param outputname name of the packed column
 * @param input_columns names of the columns to pack, in the order of the layout
 * @return a new dataframe containing the packed column
 */
ROOT::RDF::RNode build_inputs(ROOT::RDF::RNode df,
                              const std::string &outputname,
                              const std::vector<std::string> &input_columns) {
    if (input_columns.empty()) {
        throw std::runtime_error("fakefactors::generic::build_inputs requires "
                                 "at least one input column");
    }
    std::string expression = "std::vector<float>{";
    for (size_t i = 0; i < input_columns.size(); ++i) {
        expression += "static_cast<float>(" + input_columns[i] + ")";
        if (i + 1 < input_columns.size())
            expression += ", ";
    }
    expression += "}";
    return df.Define(outputname, expression);
}

/**
 * @brief Evaluates a correction of a process (fake factor or DR->SR) with
 * correctionlib. The term is 0 for events with a negative guard.
 *
 * @param df the input dataframe
 * @param correctionManager the correction manager to load corrections
 * @param outputname name of the output column
 * @param guard pt of the hadronic tau, events with a negative value are not
 * evaluated
 * @param inputs_column the packed input column
 * @param layout names of the quantities in the packed input column, in order
 * @param process process name passed to the correction
 * @param correction name of the correction in the file
 * @param variation name of the uncertainty variation or nominal
 * @param file correctionlib json file with the correction
 * @return a new dataframe containing the output column
 */
ROOT::RDF::RNode
correction_term(ROOT::RDF::RNode df,
                correctionManager::CorrectionManager &correctionManager,
                const std::string &outputname, const std::string &guard,
                const std::string &inputs_column,
                const std::vector<std::string> &layout,
                const std::string &process, const std::string &correction,
                const std::string &variation, const std::string &file) {

    Logger::get("GenericFakeFactor")
        ->debug("Setting up {} for {}", correction, process);

    auto corr = correctionManager.loadCorrection(file, correction);
    auto args = std::make_shared<ArgBuilder>(correction, corr->inputs(), layout,
                                             process);

    auto calc = [corr, args, variation](float guard,
                                        const std::vector<float> &in) {
        return guard >= 0.0f ? static_cast<float>(
                                   corr->evaluate(args->build(in, variation)))
                             : 0.0f;
    };
    return df.Define(outputname, calc, {guard, inputs_column});
}

/**
 * @brief Evaluates the compound non-closure correction of a process with
 * correctionlib, see NonClosureHandler. The term is 0 for events with a
 * negative guard.
 *
 * @param df the input dataframe
 * @param correctionManager the correction manager to load corrections
 * @param outputname name of the output column
 * @param guard pt of the hadronic tau, events with a negative value are not
 * evaluated
 * @param inputs_column the packed input column
 * @param layout names of the quantities in the packed input column, in order
 * @param process process name passed to the correction
 * @param compound_correction name of the compound correction in the file
 * @param non_closure_prefix part of the variation names before the variable
 * @param non_closure_variables variables stacked in the compound correction
 * @param variation name of the uncertainty variation or nominal
 * @param file correctionlib json file with the corrections
 * @return a new dataframe containing the output column
 */
ROOT::RDF::RNode non_closure_term(
    ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correctionManager,
    const std::string &outputname, const std::string &guard,
    const std::string &inputs_column, const std::vector<std::string> &layout,
    const std::string &process, const std::string &compound_correction,
    const std::string &non_closure_prefix,
    const std::vector<std::string> &non_closure_variables,
    const std::string &variation, const std::string &file) {

    Logger::get("GenericFakeFactor")
        ->debug("Setting up {} for {}", compound_correction, process);

    auto compound =
        correctionManager.loadCompoundCorrection(file, compound_correction);
    auto args = std::make_shared<ArgBuilder>(
        compound_correction, compound->inputs(), layout, process);
    auto handler = std::make_shared<NonClosureHandler>(
        non_closure_prefix, non_closure_variables, compound, args->syst_pos);

    auto calc = [args, handler, variation](float guard,
                                           const std::vector<float> &in) {
        return guard >= 0.0f
                   ? handler->evaluate(variation, args->build(in, "nominal"))
                   : 0.0f;
    };
    return df.Define(outputname, calc, {guard, inputs_column});
}

/**
 * @brief Evaluates the process fractions of a leg with correctionlib. The
 * output has one entry per process, all 0 for events with a negative guard.
 *
 * @param df the input dataframe
 * @param correctionManager the correction manager to load corrections
 * @param outputname name of the output column
 * @param guard pt of the hadronic tau, events with a negative value are not
 * evaluated
 * @param inputs_column the packed input column
 * @param layout names of the quantities in the packed input column, in order
 * @param correction name of the fractions correction in the file
 * @param variation name of the uncertainty variation or nominal
 * @param processes processes of the leg, in the order of the output
 * @param file correctionlib json file with the correction
 * @return a new dataframe containing the output column
 */
ROOT::RDF::RNode
fractions(ROOT::RDF::RNode df,
          correctionManager::CorrectionManager &correctionManager,
          const std::string &outputname, const std::string &guard,
          const std::string &inputs_column,
          const std::vector<std::string> &layout, const std::string &correction,
          const std::string &variation,
          const std::vector<std::string> &processes, const std::string &file) {

    Logger::get("GenericFakeFactor")
        ->debug("Setting up {} for {} processes", correction, processes.size());

    auto corr = correctionManager.loadCorrection(file, correction);
    std::vector<ArgBuilder> args;
    for (const auto &process : processes)
        args.emplace_back(correction, corr->inputs(), layout, process);

    auto calc = [corr, args, variation](float guard,
                                        const std::vector<float> &in) {
        std::vector<float> result(args.size(), 0.0f);
        if (guard >= 0.0f) {
            for (size_t i = 0; i < args.size(); ++i)
                result[i] = corr->evaluate(args[i].build(in, variation));
        }
        return result;
    };
    return df.Define(outputname, calc, {guard, inputs_column});
}

namespace {

/**
 * @brief Wrapper of an ONNX model of the session manager. It keeps the sizes of
 * the input and output tensors read at setup and checks the input size at every
 * call.
 */
struct OnnxModel {
    Ort::Session *session;
    Ort::AllocatorWithDefaultOptions allocator;
    std::vector<int64_t> input_dims, output_dims;
    int n_inputs, n_outputs;
    size_t input_size = 1, output_size = 1;

    OnnxModel(OnnxSessionManager &manager, const std::string &path)
        : session(manager.getSession(path)) {
        onnxhelper::prepare_model(session, allocator, input_dims, output_dims,
                                  n_inputs, n_outputs);
        for (auto dim : input_dims)
            input_size *= dim;
        for (auto dim : output_dims)
            output_size *= dim;
    }

    std::vector<float> operator()(std::vector<float> in) const {
        if (in.size() != input_size) {
            throw std::runtime_error(
                "fakefactors::generic: " + std::to_string(in.size()) +
                " inputs for a model with " + std::to_string(input_size));
        }
        return onnxhelper::run_interference(session, allocator, in, input_dims,
                                            output_dims, n_inputs, n_outputs);
    }
};

} // namespace

/**
 * @brief Evaluates an ONNX model with one output (fake factor or DR->SR). The
 * term is 0 for events with a negative guard.
 *
 * @param df the input dataframe
 * @param onnxSessionManager the manager that holds the session of the model
 * @param outputname name of the output column
 * @param guard pt of the hadronic tau, events with a negative value are not
 * evaluated
 * @param inputs_column the packed input column, in the feature order of the
 * model
 * @param model_file_path path of the ONNX model
 * @return a new dataframe containing the output column
 * @throws std::runtime_error if the model has not exactly one output
 */
ROOT::RDF::RNode onnx_term(ROOT::RDF::RNode df,
                           OnnxSessionManager &onnxSessionManager,
                           const std::string &outputname,
                           const std::string &guard,
                           const std::string &inputs_column,
                           const std::string &model_file_path) {

    auto model =
        std::make_shared<OnnxModel>(onnxSessionManager, model_file_path);
    if (model->output_size != 1) {
        throw std::runtime_error(
            "fakefactors::generic: model " + model_file_path + " has " +
            std::to_string(model->output_size) + " outputs, expected 1");
    }
    auto calc = [model](float guard, const std::vector<float> &in) {
        return guard >= 0.0f ? (*model)(in)[0] : 0.0f;
    };
    return df.Define(outputname, calc, {guard, inputs_column});
}

/**
 * @brief Evaluates an ONNX model with one output per process (process
 * fractions). The output has one entry per process, all 0 for events with a
 * negative guard.
 *
 * @param df the input dataframe
 * @param onnxSessionManager the manager that holds the session of the model
 * @param outputname name of the output column
 * @param guard pt of the hadronic tau, events with a negative value are not
 * evaluated
 * @param inputs_column the packed input column, in the feature order of the
 * model
 * @param model_file_path path of the ONNX model
 * @param n_processes number of processes, the model must have one output each
 * @return a new dataframe containing the output column
 * @throws std::runtime_error if the number of outputs differs from n_processes
 */
ROOT::RDF::RNode
onnx_fractions(ROOT::RDF::RNode df, OnnxSessionManager &onnxSessionManager,
               const std::string &outputname, const std::string &guard,
               const std::string &inputs_column,
               const std::string &model_file_path, const size_t n_processes) {

    auto model =
        std::make_shared<OnnxModel>(onnxSessionManager, model_file_path);
    if (model->output_size != n_processes) {
        throw std::runtime_error(
            "fakefactors::generic: model " + model_file_path + " has " +
            std::to_string(model->output_size) + " outputs for " +
            std::to_string(n_processes) + " processes");
    }
    auto calc = [model, n_processes](float guard,
                                     const std::vector<float> &in) {
        return guard >= 0.0f ? (*model)(in)
                             : std::vector<float>(n_processes, 0.0f);
    };
    return df.Define(outputname, calc, {guard, inputs_column});
}

} // namespace generic
} // namespace fakefactors
#endif /* GUARDCROWNFAKEFACTORS_CXX */
