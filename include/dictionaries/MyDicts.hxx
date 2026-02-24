#ifndef MYDICTS_H
#define MYDICTS_H

#include <map>
#include <vector>
#include <string>

#ifdef __CLING__
// The '+' at the end tells ROOT to enable streamer/I/O support
#pragma link C++ class std::map<std::string, std::vector<std::string>>+;
#pragma link C++ class std::map<std::string, std::map<std::string, std::vector<std::string>>>+;
#endif

#endif