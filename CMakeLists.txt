cmake_minimum_required(VERSION 3.10)

# Define the default compiler flags for different build types, if different from the cmake defaults
set(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING "Set default compiler flags for build type Debug")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "Set default compiler flags for build type Release")

# The build type should be set so that the correct compiler flags are chosen
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel.")

# Create the project
project(ToBeNamedPrototype CXX)

# Find ROOT and print details
find_package(ROOT 6.24 REQUIRED COMPONENTS ROOTVecOps ROOTDataFrame)

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
set_target_properties(logging PROPERTIES
    IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/lib64/libspdlog.a"
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/include")
add_dependencies(logging spdlog) # enforces to build spdlog before making the imported target

# Register analysis executable target
add_executable(Analysis analysis.cxx)
target_link_libraries(Analysis ROOT::ROOTVecOps ROOT::ROOTDataFrame logging)

# Print settings of the executable
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)
message(STATUS "")
message(STATUS "The executable is configured with following flags:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "  Base compiler flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "  Compiler flags from build type: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
message(STATUS "")

# Add install target, basically just copying the executable around relative to CMAKE_INSTALL_PREFIX
install(TARGETS Analysis DESTINATION bin)
