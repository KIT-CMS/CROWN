find_package(ROOT 6.36 REQUIRED COMPONENTS ROOTVecOps ROOTDataFrame RooFit
                                           GenVector Minuit TMVA)

message(STATUS "")
message(STATUS "Found ROOT with following settings:")
message(STATUS "  Version: ${ROOT_VERSION}")
message(STATUS "  ROOT executable: ${ROOT_EXECUTABLE}")
message(STATUS "  Include directories: ${ROOT_INCLUDE_DIRS}")
message(STATUS "  Compiler flags: ${ROOT_CXX_FLAGS}")
message(STATUS "")

# Strip -specs= flags from ROOT_CXX_FLAGS. These flags tell the compiler to load
# plugin shared libraries at paths relative to the compiler installation. Since
# ROOT_CXX_FLAGS reflects the toolchain ROOT was built with, these plugin paths
# may not exist or be incompatible when compiling with a different toolchain.
# This mismatch causes the compiler to fail when it cannot locate the referenced
# plugin files.
string(REGEX MATCHALL "-specs=[^ ]*" REMOVED_SPECS "${ROOT_CXX_FLAGS}")
string(REGEX REPLACE "-specs=[^ ]*" "" ROOT_CXX_FLAGS_FILTERED
                     "${ROOT_CXX_FLAGS}")
if(NOT ROOT_CXX_FLAGS STREQUAL ROOT_CXX_FLAGS_FILTERED)
  message(
    WARNING
      "Removed -specs= flags from ROOT_CXX_FLAGS to avoid plugin incompatibility: ${REMOVED_SPECS}"
  )
endif()

# Add ROOT flags to compile options, e.g. we have to use the same C++ standard
# Note that the flags from the build type, e.g. CMAKE_CXX_FLAGS_RELEASE, are
# automatically appended. You can check this during build time by enabling the
# verbose make output with "VERBOSE=1 make".
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
