cmake_minimum_required(VERSION 3.20)
# Unset some cached variables
unset(CMAKE_BUILD_TYPE CACHE)
unset(CMAKE_CXX_FLAGS_RELEASE CACHE)

# include cmake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
# parse args
include(ParseArguments)

# Set the default install directory to the build directory
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR} CACHE STRING "Set default install prefix to the build directory")
message(STATUS "Creating Project")
# Create the project
project(CROWN CXX)


# Print settings of the executable
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)
message(STATUS "")
message(STATUS "The executable is configured with following flags:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "  Base compiler flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "  Compiler flags from build type: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
message(STATUS "")

message(STATUS "Finding Packages")
# add the different dependencies
include(AddBaseDependencies)
include(AddRoot)
include(AddLogging)
include(AddCorrectionlib)
include(AddOnnxruntime)
include(ConfigurePython)

# installdir
if (NOT DEFINED INSTALLDIR)
    message(STATUS "No -DINSTALLDIR specified, using default: ${CMAKE_CURRENT_BINARY_DIR}/bin")
    set(INSTALLDIR ${CMAKE_CURRENT_BINARY_DIR}/bin)
endif()
set(GENERATE_CPP_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# crownlib setup
include(ConfigureCrownlib)

# run the code generation
include(CodeGeneration)

# run the build
include(Build)

# tests
# Include tests
enable_testing()
add_subdirectory(tests)