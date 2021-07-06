cmake_minimum_required(VERSION 3.10)

# Check that config name is given with -DCONFIG, needed to generate the C++ code
if (NOT DEFINED ANALYSIS)
    message(FATAL_ERROR "Please specify the config to be used with -DANALYSIS=my_config_name")
endif()
message(STATUS "Set up analysis for config ${ANALYSIS}.")

if (NOT DEFINED CHANNELS)
    message(STATUS "No channel specificed, build the four main channels using -DCHANNELS=et;mt;tt;em")
    set(CHANNELS "et;mt;tt;em")
endif()

if (NOT DEFINED SHIFTS)
    message(STATUS "No shifts specificed, using -DSHIFTS=all")
    set(SHIFTS "all")
endif()
if (NOT DEFINED DEBUG)
    message(STATUS "No debug mode")
    set(DEBUG "false")
endif()

if (NOT DEFINED SAMPLES)
    message(FATAL_ERROR "Please specify the samples to be used with -DSAMPLES=samples")
endif()
message(STATUS "Set up analysis with --config ${ANALYSIS} --channels ${CHANNELS} --shifts ${SHIFTS} --samples ${SAMPLES} --debug ${DEBUG}")

# Define the default compiler flags for different build types, if different from the cmake defaults
set(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING "Set default compiler flags for build type Debug")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "Set default compiler flags for build type Release")

# The build type should be set so that the correct compiler flags are chosen
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel.")

# Set the default install directory to the build directory
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR} CACHE STRING "Set default install prefix to the build directory")

# Create the project
project(CROWN CXX)

# Find ROOT and print details
find_package(ROOT 6.24 REQUIRED COMPONENTS ROOTVecOps ROOTDataFrame RooFit)

message(STATUS "")
message(STATUS "Found ROOT with following settings:")
message(STATUS "  Version: ${ROOT_VERSION}")
message(STATUS "  ROOT executable: ${ROOT_EXECUTABLE}")
message(STATUS "  Include directories: ${ROOT_INCLUDE_DIRS}")
message(STATUS "  Compiler flags: ${ROOT_CXX_FLAGS}")
message(STATUS "")

# Add ROOT flags to compile options, e.g. we have to use the same C++ standard
# Note that the flags from the build type, e.g. CMAKE_CXX_FLAGS_RELEASE, are
# automatically appended. You can check this during build time by enabling
# the verbose make output with "VERBOSE=1 make".
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")

# Use -fconcepts with g++ to silence following warning:
# warning: use of 'auto' in parameter declaration only available with '-fconcepts
if (CMAKE_CXX_COMPILER_ID STREQUAL GNU)
    message(STATUS "Attach -fconcepts to the compiler flags to silence warnings.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts")
endif()

# Find the C++ standard from ROOT and set it as the standard of this project
# We require the C++ standard 17 or 20 and don't want to fall back to lower versions.
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if (${ROOT_CXX_FLAGS} MATCHES "\\-std\\=c\\+\\+17")
    message(STATUS "Set c++17 as the C++ standard.")
    set(CMAKE_CXX_STANDARD 17)
elseif (${ROOT_CXX_FLAGS} MATCHES "\\-std\\=c\\+\\+20")
    message(STATUS "Set c++20 as the C++ standard.")
    set(CMAKE_CXX_STANDARD 20)
else ()
    message(FATAL_ERROR "The standard c++17 or higher is required but not found in the ROOT flags: ${ROOT_CXX_FLAGS}")
endif()

# Build the logging library
include(ExternalProject)
ExternalProject_Add(spdlog
    PREFIX          spdlog
    GIT_REPOSITORY  https://github.com/gabime/spdlog.git
    GIT_TAG         v1.8.5
    CMAKE_ARGS      -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
                    -DCMAKE_BUILD_TYPE=Release
                    -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}
    LOG_DOWNLOAD 1 LOG_CONFIGURE 1 LOG_BUILD 1 LOG_INSTALL 1
)

# Make an imported target out of the build logging library
add_library(logging STATIC IMPORTED)
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/include") # required because the include dir must be existent for INTERFACE_INCLUDE_DIRECTORIES
include(GNUInstallDirs) # required to populate CMAKE_INSTALL_LIBDIR with lib or lib64 required for the destination of libspdlog.a
set_target_properties(logging PROPERTIES
    IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/libspdlog.a"
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/include")
add_dependencies(logging spdlog) # enforces to build spdlog before making the imported target

# Print settings of the executable
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)
message(STATUS "")
message(STATUS "The executable is configured with following flags:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "  Base compiler flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "  Compiler flags from build type: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
message(STATUS "")

# Find Python 3
find_package(Python 3.6 REQUIRED COMPONENTS Interpreter)

# detect virtualenv and set Pip args accordingly
if(DEFINED ENV{VIRTUAL_ENV} OR DEFINED ENV{CONDA_PREFIX})
  set(_pip_args)
else()
  set(_pip_args "--user")
endif()

function(find_python_package PYPINAME NAME MIN_VERSION)
    execute_process(COMMAND "${Python_EXECUTABLE}" "-c" "import ${NAME}; print(${NAME}.__version__)"
                    RESULT_VARIABLE PACKAGE_NOT_FOUND
                    OUTPUT_VARIABLE PACKAGE_VERSION
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(${PACKAGE_NOT_FOUND} EQUAL 1)
        execute_process(COMMAND ${Python_EXECUTABLE} -m pip install ${PYPINAME} ${_pip_args})
        execute_process(COMMAND "${Python_EXECUTABLE}" "-c" "import ${NAME}; print(${NAME}.__version__)"
                    RESULT_VARIABLE PACKAGE_NOT_FOUND
                    OUTPUT_VARIABLE PACKAGE_VERSION
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(${PACKAGE_NOT_FOUND} EQUAL 1)
            message(FATAL_ERROR "Failed to import ${PYPINAME} or get version.")
        endif()
    endif()
    if(PACKAGE_VERSION VERSION_LESS MIN_VERSION)
        message(FATAL_ERROR "The version of Python package ${PYPINAME} is too old (found ${PACKAGE_VERSION}, require at least ${MIN_VERSION}).")
    endif()
    message(STATUS "Found Python package ${PYPINAME} (require ${MIN_VERSION}, found ${PACKAGE_VERSION})")
endfunction()

find_python_package(GitPython git 3.1)

# Generate the C++ code
set(GENERATE_CPP_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(GENERATE_CPP_INPUT_TEMPLATE "${CMAKE_SOURCE_DIR}/code_generation/analysis_template.cxx")

message(STATUS "")
message(STATUS "Generate C++ code with following settings:")
message(STATUS "  Output directory: ${GENERATE_CPP_OUTPUT_DIRECTORY}")
message(STATUS "  Template: ${GENERATE_CPP_INPUT_TEMPLATE}")
message(STATUS "  Analysis: ${ANALYSIS}")
message(STATUS "  Channels: ${CHANNELS}")
message(STATUS "  Shifts: ${SHIFTS}")
message(STATUS "  Samples: ${SAMPLES}")
message(STATUS "")

file(MAKE_DIRECTORY ${GENERATE_CPP_OUTPUT_DIRECTORY})
execute_process(
    COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/generate.py --template ${GENERATE_CPP_INPUT_TEMPLATE} --output ${GENERATE_CPP_OUTPUT_DIRECTORY} --analysis ${ANALYSIS} --channels ${CHANNELS} --shifts ${SHIFTS} --samples ${SAMPLES} --debug ${DEBUG}
)

set(GENERATE_CPP_OUTPUT_FILELIST "${GENERATE_CPP_OUTPUT_DIRECTORY}/files.txt")
if(NOT EXISTS ${GENERATE_CPP_OUTPUT_FILELIST})
    message(FATAL_ERROR "List of generated C++ files in ${GENERATE_CPP_OUTPUT_FILELIST} does not exist.")
endif()

# Iterate over files from output filelist and add build and install targets
FILE(READ ${GENERATE_CPP_OUTPUT_FILELIST} FILELIST)
STRING(REGEX REPLACE "\n" ";" FILELIST ${FILELIST})
set(TARGET_NAMES "")
# copy all correction files into the install location
install(DIRECTORY data/ DESTINATION data)
foreach(FILENAME ${FILELIST})
    STRING(REGEX REPLACE ".cxx" "" TARGET_NAME ${FILENAME})
    list(APPEND TARGET_NAMES ${TARGET_NAME})
    set(FULL_PATH "${GENERATE_CPP_OUTPUT_DIRECTORY}/${FILENAME}")

    # Add build target
    message(STATUS "Add build target for file ${FILENAME}.")
    add_executable(${TARGET_NAME} ${FULL_PATH})
    # target_precompile_headers(${TARGET_NAME} PRIVATE jets.hxx lorentzvectors.hxx metfilter.hxx pairselection.hxx quantities.hxx)
    target_include_directories(${TARGET_NAME} PRIVATE ${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
    target_link_libraries(${TARGET_NAME} ROOT::ROOTVecOps ROOT::ROOTDataFrame ROOT::RooFit logging)
    # Add install target, basically just copying the executable around relative to CMAKE_INSTALL_PREFIX
    install(TARGETS ${TARGET_NAME} DESTINATION)
endforeach()

# Include tests
enable_testing()
add_subdirectory(tests)
