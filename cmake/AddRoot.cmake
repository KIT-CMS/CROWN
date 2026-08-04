find_package(ROOT 6.36 REQUIRED COMPONENTS ROOTVecOps ROOTDataFrame RooFit
                                           GenVector Minuit TMVA)

message(STATUS "")
message(STATUS "Found ROOT with following settings:")
message(STATUS "  Version: ${ROOT_VERSION}")
message(STATUS "  ROOT executable: ${ROOT_EXECUTABLE}")
message(STATUS "  Include directories: ${ROOT_INCLUDE_DIRS}")
message(STATUS "  Compiler flags: ${ROOT_CXX_FLAGS}")
message(STATUS "")

# Add ROOT flags to compile options, e.g. we have to use the same C++ standard
# Note that the flags from the build type, e.g. CMAKE_CXX_FLAGS_RELEASE, are
# automatically appended. You can check this during build time by enabling the
# verbose make output with "VERBOSE=1 make".
# Strip RPM hardening -specs= flags: they reference redhat-annobin-cc1/
# redhat-hardened-cc1, which point at an annobin.so relative to whichever
# compiler is invoked. ROOT_CXX_FLAGS carries these over from the host gcc
# ROOT was originally built with, but the CVMFS gcc release used here doesn't
# ship its own annobin.so, so cc1plus dies with "inaccessible plugin file".
string(REGEX REPLACE "-specs=[^ ]*" "" ROOT_CXX_FLAGS_FILTERED "${ROOT_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS_FILTERED}")

# Use -fconcepts with g++ to silence following warning: warning: use of 'auto'
# in parameter declaration only available with '-fconcepts
if(CMAKE_CXX_COMPILER_ID STREQUAL GNU)
  message(STATUS "Attach -fconcepts to the compiler flags to silence warnings.")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts")
endif()
# require the C++ standard 20 and don't want to fall back to lower
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 20)
