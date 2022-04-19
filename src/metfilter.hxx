/// The namespace that contains the metfilter function.

#include <string>

namespace metfilter {

/// Function to apply a metfilter to the dataframe
///
/// \param df the dataframe to add the quantity to
/// \param flagname Name of the Filterflag in NanoAOD
/// \param filtername Name of the Filter to be shown in the dataframe report
///
/// \returns a dataframe with the filter applied
auto ApplyMetFilter(auto &df, const std::string &flagname,
                    const std::string &filtername) {
    return df.Filter([](const bool flag) { return flag; }, {flagname},
                     filtername);
}

} // namespace metfilter