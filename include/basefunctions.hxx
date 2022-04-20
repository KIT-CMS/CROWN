#ifndef GUARDBASEFUCTIONS_H
#define GUARDBASEFUCTIONS_H

#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "utility/RooFunctorThreadsafe.hxx"

enum Channel { MT = 0, ET = 1, TT = 2, EM = 3 };

namespace basefunctions {

// auto JSONFilter(auto &df, const std::string &json_path, const std::string
// &run,
//                 const std::string &luminosity, const std::string
//                 &filtername);
template <typename T>
ROOT::RDF::RNode rename(auto &df, const std::string &inputname,
                        const std::string &outputname);
template <typename T>
ROOT::RDF::RNode DefineQuantity(auto &df, const std::string &outputname,
                                T const &value);

template <class... Flags>
ROOT::RDF::RNode FilterFlagsAny(auto &df, const std::string &filtername,
                                const Flags &...flags);
template <class... Flags>
ROOT::RDF::RNode CombineFlagsAny(auto &df, const std::string &outputflag,
                                 const Flags &...flags);

template <typename T>
ROOT::RDF::RNode FilterIntSelection(auto &df, const std::string &quantity,
                                    const std::vector<T> &selection,
                                    const std::string &filtername);
/// Function to apply a maximal filter requirement to a quantity.
/// Returns true if the value is smaller than the given cut value
///
/// \param cut The cut value of the filter
///
/// \returns a lambda function to be used in RDF Define
inline auto FilterMax(const float &cut) {
    return [cut](const ROOT::RVec<float> &values) {
        ROOT::RVec<int> mask = values < cut;
        return mask;
    };
}

/// Function to apply a maximal filter requirement to a quantity.
/// Returns true if the absolute value is smaller than the given cut value
///
/// \param cut The cut value of the filter
///
/// \returns a lambda function to be used in RDF Define
inline auto FilterAbsMax(const float &cut) {
    return [cut](const ROOT::RVec<float> &values) {
        ROOT::RVec<int> mask = abs(values) < cut;
        return mask;
    };
}

/// Function to apply a minimal filter requirement to a quantity.
/// Returns true if the value is larger than the given cut value
///
/// \param cut The cut value of the filter
///
/// \returns a lambda function to be used in RDF Define
inline auto FilterMin(const float &cut) {
    // As in ROOT, for min we use >=
    return [cut](const ROOT::RVec<float> &values) {
        ROOT::RVec<int> mask = values >= cut;
        return mask;
    };
}

/// Function to apply a minimal filter requirement to an integer quantity.
/// Returns true if the value is larger than the given cut value
///
/// \param cut The cut value of the filter
///
/// \returns a lambda function to be used in RDF Define
inline auto FilterMinInt(const int &cut) {
    // As in ROOT, for min we use >=
    return [cut](const ROOT::RVec<int> &values) {
        ROOT::RVec<int> mask = values >= cut;
        return mask;
    };
}

/// Function to apply a minimal filter requirement to a quantity.
/// Returns true if the absolute value is larger than the given cut value
///
/// \param cut The cut value of the filter
///
/// \returns a lambda function to be used in RDF Define
inline auto FilterAbsMin(const float &cut) {
    return [cut](const ROOT::RVec<float> &values) {
        ROOT::RVec<int> mask = abs(values) >= cut;
        return mask;
    };
}

/// Function to combine two RVec Masks by multiplying the two RVec elementwise
///
/// \param mask_1 The first mask
/// \param mask_2 The second mask
///
/// \returns a lambda function which returns the multiplication of the two masks
inline auto MultiplyTwoMasks() {
    return [](const ROOT::RVec<Int_t> &mask_1,
              const ROOT::RVec<Int_t> &mask_2) { return mask_1 * mask_2; };
}

/// Function to filter the Tau ID in NanoAOD. The disciminator output is stored
/// as a bitmask. In order to filter based on a given WP, a bitwise comparison
/// between the given ID and the filtered index is performed. Return true if the
/// working point is passed by a tau.
///
/// \param index The bitmask index to be used for comparison
///
/// \returns a lambda function to be used in RDF Define
inline auto FilterID(const int &index) {
    return [index](const ROOT::RVec<UChar_t> &IDs) {
        ROOT::RVec<int> mask;
        for (auto const ID : IDs) {
            if (index > 0)
                mask.push_back(std::min(1, int(ID & 1 << index - 1)));
            else
                mask.push_back(int(1));
        }
        return mask;
    };
}
/// Function to filter the Jet ID in NanoAOD. Ths values are stored bitwise, so
/// the integer value has to be decoded into binary, and then the value of the
/// index bit has to be compared.
///
/// \param index The bitmask index to be used for comparison
///
/// \returns a lambda function to be used in RDF Define
inline auto FilterJetID(const int &index) {
    return [index](const ROOT::RVec<Int_t> &IDs) {
        ROOT::RVec<int> mask;
        for (auto const ID : IDs) {
            mask.push_back(std::min(1, (ID >> index) & 1));
        }
        return mask;
    };
}
/// Function to filter the Jet pileup ID in NanoAOD. This ID is applied on jets
/// below a given pt threshold. The jet pileup ID has 4 possible values:
/// 0==fail, 4==pass loose, 6==pass loose and medium, 7==pass loose, medium and
/// tight
///
/// \param PUindex The index to be used for comparison
/// \param PUptcut The pt threshold to be used for comparison
///
/// \returns a lambda function to be used in RDF Define
inline auto FilterJetPUID(const int &PUindex, const float &PUptcut) {
    return [PUindex, PUptcut](const ROOT::RVec<Int_t> &PUIDs,
                              const ROOT::RVec<float> &jet_pts) {
        ROOT::RVec<int> tmp_mask1 = PUIDs >= PUindex;
        ROOT::RVec<int> tmp_mask2 = jet_pts >= PUptcut;
        ROOT::RVec<int> mask = (tmp_mask1 + tmp_mask2) > 0;
        return mask;
    };
}
template <class... Inputs>
ROOT::RDF::RNode
evaluateWorkspaceFunction(auto &df, const std::string &outputname,
                          const std::shared_ptr<RooFunctorThreadsafe> &function,
                          const Inputs &...inputs);
template <typename T>
ROOT::RDF::RNode UnrollVectorQuantity(auto &df, const std::string &name,
                                      const std::vector<std::string> &names,
                                      const size_t &idx = 0);
} // namespace basefunctions

#endif /* GUARDBASEFUNCTIONS_H */
