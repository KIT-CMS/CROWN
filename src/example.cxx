#ifndef GUARD_EXAMPLE_H
#define GUARD_EXAMPLE_H

#include "../include/defaults.hxx"
#include "../include/utility/CorrectionManager.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "utility/Logger.hxx"
#include "utility/RooFunctorThreadsafe.hxx"
#include "utility/utility.hxx"

namespace example {

// // Template for a Filter function
// ROOT::RDF::RNode <Filter func>(ROOT::RDF::RNode df, const std::string & <input>, <other inputs>) {
//     auto <lambda function> = [<static inputs>](<lambda inputs>) {};
//     auto df1 = df.Filter(<lambda function>, {<input>});
//     return df1;
// }

// // Template for a Define function
// ROOT::RDF::RNode <Define func>(ROOT::RDF::RNode df, const std::string & <input>, 
//                                const std::string & <output>, <other inputs>) {
//     auto <lambda function> = [<static inputs>](<inputs>) {};
//     auto df1 = df.Define(<output>, <lambda function>, {<input>});
//     return df1;
// }


} // namespace solution
#endif /* GUARD_EXAMPLE_H */
