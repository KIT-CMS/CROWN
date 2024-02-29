find_package(ROOT 6.26 REQUIRED COMPONENTS ROOTVecOps ROOTDataFrame RooFit
                                           GenVector)

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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")

# Use -fconcepts with g++ to silence following warning: warning: use of 'auto'
# in parameter declaration only available with '-fconcepts
if(CMAKE_CXX_COMPILER_ID STREQUAL GNU)
  message(STATUS "Attach -fconcepts to the compiler flags to silence warnings.")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts")
endif()

# Find the C++ standard from ROOT and set it as the standard of this project We
# require the C++ standard 17 or 20 and don't want to fall back to lower
# versions.
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(${ROOT_CXX_FLAGS} MATCHES "\\-std\\=c\\+\\+17")
  message(STATUS "Set c++17 as the C++ standard.")
  set(CMAKE_CXX_STANDARD 17)
elseif(${ROOT_CXX_FLAGS} MATCHES "\\-std\\=c\\+\\+20")
  message(STATUS "Set c++20 as the C++ standard.")
  set(CMAKE_CXX_STANDARD 20)
elseif(${ROOT_CXX_FLAGS} MATCHES "\\-std\\=c\\+\\+14")
  message(STATUS "c++14 found, setting c++17 as the C++ standard.")
  set(CMAKE_CXX_STANDARD 17)
else()
  message(
    FATAL_ERROR
      "The standard c++17 or higher is required but not found in the ROOT flags: ${ROOT_CXX_FLAGS}"
  )
endif()
