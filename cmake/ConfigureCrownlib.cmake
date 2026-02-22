# build a shared lib from all CROWN functions
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/include)

include_directories(${CMAKE_SOURCE_DIR}/analysis_configurations/${ANALYSIS}/cpp_addons/src)
include_directories(${CMAKE_SOURCE_DIR}/analysis_configurations/${ANALYSIS}/cpp_addons/include)

file(GLOB SOURCES_1 
  ${CMAKE_SOURCE_DIR}/src/*.cxx
  ${CMAKE_SOURCE_DIR}/src/*/*.cxx)

file(GLOB SOURCES_2 
  ${CMAKE_SOURCE_DIR}/analysis_configurations/${ANALYSIS}/cpp_addons/src/*.cxx
  ${CMAKE_SOURCE_DIR}/analysis_configurations/${ANALYSIS}/cpp_addons/src/*/*.cxx)

set(SOURCES ${SOURCES_1} ${SOURCES_2})

# Function to configure CROWNLIB to avoid code duplication in your if/else
macro(configure_crownlib_target)
    add_library(CROWNLIB SHARED ${SOURCES})
    target_include_directories(CROWNLIB PRIVATE ${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
    target_link_libraries(
        CROWNLIB
        "-Wl,--no-as-needed"
        MyDicts
        "-Wl,--as-needed"
        ROOT::Core
        ROOT::RIO
        ROOT::ROOTVecOps
        ROOT::ROOTDataFrame
        ROOT::RooFit
        ROOT::GenVector
        logging
        correctionlib
        nlohmann_json::nlohmann_json
        ${ONNX_RUNTIME_LIB_PATH}
    )
endmacro()

if(BUILD_CROWNLIB_ONLY)
    message(STATUS "Building only the CROWNLIB library")
    configure_crownlib_target()
    install(TARGETS CROWNLIB DESTINATION ${INSTALLDIR}/lib)
    return()
endif()

# Check for existing CROWNLIB
find_library(CROWNLIB_FOUND CROWNLIB HINTS ${INSTALLDIR}/lib ${CMAKE_CURRENT_BINARY_DIR})

if(NOT CROWNLIB_FOUND OR REBUILD_CROWN_LIB)
    message(STATUS "CROWNLIB not found or rebuild requested, building it")
    configure_crownlib_target()
    install(TARGETS CROWNLIB DESTINATION ${INSTALLDIR}/lib)
    set(CMAKE_BUILD_RPATH ${INSTALLDIR}/lib)
else()
    message(STATUS "Found CROWNLIB in ${CROWNLIB_FOUND}")
    install(FILES ${CROWNLIB_FOUND} DESTINATION ${INSTALLDIR}/lib)
    link_directories(${INSTALLDIR}/lib ${CMAKE_CURRENT_BINARY_DIR})
endif()