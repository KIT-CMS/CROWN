#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ class std::map < std::string, std::vector < std::string>> + ;
#pragma link C++ class std::map < std::string, std::map < std::string,         \
    std::vector < std::string>>> +                                             \
    ;
#endif