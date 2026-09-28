cmake_minimum_required(VERSION 3.20)
# Unset some cached variables
unset(CMAKE_BUILD_TYPE CACHE)
unset(CMAKE_CXX_FLAGS_RELEASE CACHE)

# include cmake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
# parse args
include(ParseArguments)

# Set the default install directory to the build directory
set(CMAKE_INSTALL_PREFIX
    ${CMAKE_BINARY_DIR}
    CACHE STRING "Set default install prefix to the build directory")

# Use the compilers from the active conda env (if any) unless the user
# explicitly set CMAKE_C_COMPILER / CMAKE_CXX_COMPILER themselves.
if(DEFINED ENV{CONDA_PREFIX})
  if(NOT CMAKE_C_COMPILER)
    file(GLOB CONDA_CC_CANDIDATES "$ENV{CONDA_PREFIX}/bin/*-cc")
    if(CONDA_CC_CANDIDATES)
      list(GET CONDA_CC_CANDIDATES 0 CONDA_CC)
      message(STATUS "Using conda env C compiler: ${CONDA_CC}")
      set(CMAKE_C_COMPILER
          "${CONDA_CC}"
          CACHE FILEPATH "C compiler from conda env")
    endif()
  endif()
  if(NOT CMAKE_CXX_COMPILER)
    file(GLOB CONDA_CXX_CANDIDATES "$ENV{CONDA_PREFIX}/bin/*-c++")
    if(CONDA_CXX_CANDIDATES)
      list(GET CONDA_CXX_CANDIDATES 0 CONDA_CXX)
      message(STATUS "Using conda env C++ compiler: ${CONDA_CXX}")
      set(CMAKE_CXX_COMPILER
          "${CONDA_CXX}"
          CACHE FILEPATH "C++ compiler from conda env")
    endif()
  endif()
endif()

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
message(
  STATUS
    "  Compiler flags from build type: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}"
)
message(STATUS "")

# Add conda library path early so all targets can find libstdc++
if(DEFINED ENV{CONDA_PREFIX})
  link_directories($ENV{CONDA_PREFIX}/lib)
  message(STATUS "Adding conda lib directory: $ENV{CONDA_PREFIX}/lib")
endif()

# installdir
if(NOT DEFINED INSTALLDIR)
  message(
    STATUS
      "No -DINSTALLDIR specified, using default: ${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  set(INSTALLDIR ${CMAKE_CURRENT_BINARY_DIR}/bin)
endif()
set(GENERATE_CPP_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

message(STATUS "Finding Packages")
# add the different dependencies
include(AddBaseDependencies)
include(AddRoot)
include(AddDicts)
include(AddLogging)
include(ConfigurePython)
include(AddCorrectionlib)
include(AddOnnxruntime)

# crownlib setup
include(ConfigureCrownlib)
if(BUILD_CROWNLIB_ONLY)
  # exit if only crownlib is being built
  return()
endif()
# run the code generation
include(CodeGeneration)

message(STATUS "Next is Build")
# run the build
include(Build)

# tests Include tests
enable_testing()
add_subdirectory(tests)
