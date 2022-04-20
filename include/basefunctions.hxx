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
auto FilterMax(const float &cut);
auto FilterAbsMax(const float &cut);
auto FilterMin(const float &cut);
auto FilterMinInt(const int &cut);
auto FilterAbsMin(const float &cut);
auto MultiplyTwoMasks();
auto FilterID(const int &index);
auto FilterJetID(const int &index);
auto FilterJetPUID(const int &PUindex, const float &PUptcut);
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
