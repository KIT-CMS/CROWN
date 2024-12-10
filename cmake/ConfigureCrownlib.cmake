# build a shared lib from all CROWN functions
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/include)

include_directories(${CMAKE_SOURCE_DIR}/analysis_configurations/*/cpp_addons/src)
include_directories(${CMAKE_SOURCE_DIR}/analysis_configurations/*/cpp_addons/include)

file(GLOB SOURCES_1 ${CMAKE_SOURCE_DIR}/src/*.cxx)
file(GLOB SOURCES_2 ${CMAKE_SOURCE_DIR}/src/utility/*.cxx
     ${CMAKE_SOURCE_DIR}/src/RecoilCorrections/*.cxx
     ${CMAKE_SOURCE_DIR}/src/SVFit/*.cxx)

file(GLOB SOURCES_3 ${CMAKE_SOURCE_DIR}/analysis_configurations/*/cpp_addons/src/*.cxx)

set(SOURCES ${SOURCES_1} ${SOURCES_2} ${SOURCES_3})

if(BUILD_CROWNLIB_ONLY)
  message(STATUS "Building only the CROWNLIB library")
  add_library(CROWNLIB SHARED ${SOURCES})
  target_include_directories(CROWNLIB PRIVATE ${CMAKE_SOURCE_DIR}
                                              ${ROOT_INCLUDE_DIRS})
  target_link_libraries(
    CROWNLIB
    ROOT::ROOTVecOps
    ROOT::ROOTDataFrame
    ROOT::RooFit
    ROOT::GenVector
    logging
    correctionlib
    nlohmann_json::nlohmann_json
    ${ONNX_RUNTIME_LIB_PATH})
  install(TARGETS CROWNLIB DESTINATION ${INSTALLDIR}/lib)
  return()
endif()
# check if CROWNLIB is already installed
find_library(
  CROWNLIB_FOUND CROWNLIB HINTS ${INSTALLDIR}/lib ${CMAKE_CURRENT_BINARY_DIR}
                                ${CMAKE_CURRENT_BINARY_DIR}/lib)
if(NOT CROWNLIB_FOUND OR REBUILD_CROWN_LIB)
  message(STATUS "CROWNLIB not found, building it")
  # CROWNLIB not found, build it
  add_library(CROWNLIB SHARED ${SOURCES})
  target_include_directories(CROWNLIB PRIVATE ${CMAKE_SOURCE_DIR}
                                              ${ROOT_INCLUDE_DIRS})
  target_link_libraries(
    CROWNLIB
    ROOT::ROOTVecOps
    ROOT::ROOTDataFrame
    ROOT::RooFit
    ROOT::GenVector
    logging
    correctionlib
    nlohmann_json::nlohmann_json
    ${ONNX_RUNTIME_LIB_PATH})
  install(TARGETS CROWNLIB DESTINATION ${INSTALLDIR}/lib)
else()
  message(STATUS "Found CROWNLIB in ${CROWNLIB_FOUND}")
  install(FILES ${CROWNLIB_FOUND} DESTINATION ${INSTALLDIR}/lib)
  link_directories(${INSTALLDIR}/lib ${CMAKE_CURRENT_BINARY_DIR}
                   ${CMAKE_CURRENT_BINARY_DIR}/lib)
endif()
