#ifndef GUARDCROWNFAKEFACTORS_H
#define GUARDCROWNFAKEFACTORS_H

#include "ROOT/RDataFrame.hxx"
#include "correction.h"
#include "utility/CorrectionManager.hxx"
#include "utility/OnnxSessionManager.hxx"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

namespace fakefactors {
namespace generic {

// Maps the packed input column onto the inputs declared by a correction
struct ArgBuilder {
    enum class Kind { Value, Process, Syst };
    struct Slot {
        Kind kind;
        size_t index;
        correction::Variable::VarType type;
    };
    std::vector<Slot> slots;
    size_t syst_pos = 0;
    std::string process;

    ArgBuilder(const std::string &correction_name,
               const std::vector<correction::Variable> &inputs,
               const std::vector<std::string> &layout,
               const std::string &process);

    std::vector<correction::Variable::Type>
    build(const std::vector<float> &values, const std::string &syst) const;
};

struct NonClosureHandler {
    std::string prefix;
    const correction::CompoundCorrection *compound;
    size_t syst_pos;
    std::vector<std::string> variables;

    NonClosureHandler(const std::string &prefix,
                      const std::vector<std::string> &variables,
                      const correction::CompoundCorrection *c, size_t syst_pos);

    float evaluate(const std::string &systematic,
                   std::vector<correction::Variable::Type> args) const;
};

// Packs the columns into one std::vector<float>, casting each (int, short,
// double, ...) to float
ROOT::RDF::RNode build_inputs(ROOT::RDF::RNode df,
                              const std::string &outputname,
                              const std::vector<std::string> &input_columns);

// The terms are 0 (the fractions all 0) and not evaluated for events with guard
// < 0

// Evaluates the correction of a process (fake factor, DR->SR) from the
// correctionlib payloads
ROOT::RDF::RNode
correction_term(ROOT::RDF::RNode df,
                correctionManager::CorrectionManager &correctionManager,
                const std::string &outputname, const std::string &guard,
                const std::string &inputs_column,
                const std::vector<std::string> &layout,
                const std::string &process, const std::string &correction,
                const std::string &variation, const std::string &file);

// non_closure_prefix is the part of the non-closure variation names before the
// variable, e.g. "CMS_fake_t_QCD_non_closure_"
ROOT::RDF::RNode non_closure_term(
    ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correctionManager,
    const std::string &outputname, const std::string &guard,
    const std::string &inputs_column, const std::vector<std::string> &layout,
    const std::string &process, const std::string &compound_correction,
    const std::string &non_closure_prefix,
    const std::vector<std::string> &non_closure_variables,
    const std::string &variation, const std::string &file);

// The process fractions of a leg as one std::vector<float>, one entry per
// process
ROOT::RDF::RNode
fractions(ROOT::RDF::RNode df,
          correctionManager::CorrectionManager &correctionManager,
          const std::string &outputname, const std::string &guard,
          const std::string &inputs_column,
          const std::vector<std::string> &layout, const std::string &correction,
          const std::string &variation,
          const std::vector<std::string> &processes, const std::string &file);

// ONNX model with one output (fake factor, DR->SR) or one output per process
// (fractions)
ROOT::RDF::RNode
onnx_term(ROOT::RDF::RNode df, OnnxSessionManager &onnxSessionManager,
          const std::string &outputname, const std::string &guard,
          const std::string &inputs_column, const std::string &model_file_path);

ROOT::RDF::RNode
onnx_fractions(ROOT::RDF::RNode df, OnnxSessionManager &onnxSessionManager,
               const std::string &outputname, const std::string &guard,
               const std::string &inputs_column,
               const std::string &model_file_path, const size_t n_processes);

// Sum over the processes of max(fraction, 0) * prod(max(term, 0)); columns:
// guard, fractions, then n_terms terms per process
template <typename I> struct CombineHelper;

template <std::size_t... N> struct CombineHelper<std::index_sequence<N...>> {
    template <std::size_t> using Float = float;
    std::vector<size_t> n_terms;
    bool strict;

    float operator()(float guard, const std::vector<float> &fractions,
                     Float<N>... terms) const {
        const float values[] = {terms...};
        float sum = 0.0f;
        if (strict ? guard > 0.0f : guard >= 0.0f) {
            size_t k = 0;
            for (size_t p = 0; p < n_terms.size(); ++p) {
                float product = std::max(fractions.at(p), 0.0f);
                for (size_t t = 0; t < n_terms[p]; ++t)
                    product *= std::max(values[k++], 0.0f);
                sum += product;
            }
        }
        return sum;
    }
};

template <std::size_t N>
inline ROOT::RDF::RNode
combine(ROOT::RDF::RNode df, const std::string &outputname,
        const std::vector<std::string> &columns,
        const std::vector<size_t> &n_terms, const bool strict) {
    return df.Define(
        outputname, CombineHelper<std::make_index_sequence<N>>{n_terms, strict},
        columns);
}

} // namespace generic
} // namespace fakefactors
#endif /* GUARDCROWNFAKEFACTORS_H */
