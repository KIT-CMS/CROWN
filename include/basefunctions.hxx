#ifndef GUARDBASEFUCTIONS_H
#define GUARDBASEFUCTIONS_H

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
ROOT::RDF::RNode FilterMax(const float &cut);
ROOT::RDF::RNode FilterAbsMax(const float &cut);
ROOT::RDF::RNode FilterMin(const float &cut);
ROOT::RDF::RNode FilterMinInt(const int &cut);
ROOT::RDF::RNode FilterAbsMin(const float &cut);
ROOT::RDF::RNode MultiplyTwoMasks();
ROOT::RDF::RNode FilterID(const int &index);
ROOT::RDF::RNode FilterJetID(const int &index);
ROOT::RDF::RNode FilterJetPUID(const int &PUindex, const float &PUptcut);
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
