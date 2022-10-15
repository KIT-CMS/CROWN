#ifndef GUARDMETFILTER_H
#define GUARDMETFILTER_H

namespace metfilter {

ROOT::RDF::RNode ApplyMetFilter(ROOT::RDF::RNode df,
                                const std::string &flagname,
                                const std::string &filtername);
} // namespace metfilter
#endif /* GUARDMETFILTER_H */