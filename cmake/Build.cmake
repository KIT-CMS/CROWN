set(GENERATE_CPP_OUTPUT_FILELIST "${GENERATE_CPP_OUTPUT_DIRECTORY}/files.txt")
if(NOT EXISTS ${GENERATE_CPP_OUTPUT_FILELIST})
  message(
    FATAL_ERROR
      "List of generated C++ files in ${GENERATE_CPP_OUTPUT_FILELIST} does not exist."
  )
endif()

# Iterate over files from output filelist and add build and install targets
file(READ ${GENERATE_CPP_OUTPUT_FILELIST} FILELIST)
string(REGEX REPLACE "\n" ";" FILELIST ${FILELIST})
set(TARGET_NAMES "")
# copy all correction files into the install location
install(DIRECTORY data/ DESTINATION ${INSTALLDIR}/data)
if(PAYLOADS)
  install(
    DIRECTORY ${CMAKE_SOURCE_DIR}/analysis_configurations/${ANALYSIS}/payloads
    DESTINATION ${INSTALLDIR})
endif()

# also copy inish script needed for job tarball
install(FILES init.sh DESTINATION ${INSTALLDIR})
foreach(FILENAME ${FILELIST})
  # STRING(REGEX REPLACE ".cxx" "" TARGET_NAME ${FILENAME})
  cmake_path(GET FILENAME RELATIVE_PART RELATIVE_PATH)
  cmake_path(GET FILENAME FILENAME TARGET_FILENAMENAME)
  string(REGEX REPLACE ".cxx" "" TARGET_NAME ${TARGET_FILENAMENAME})
  string(REGEX REPLACE "/${TARGET_FILENAMENAME}" "" GENERATED_CODEBASE
                       ${RELATIVE_PATH})

  list(APPEND TARGET_NAMES ${TARGET_NAME})
  set(FULL_PATH "${GENERATE_CPP_OUTPUT_DIRECTORY}/${FILENAME}")

  # Add build target
  message(STATUS "Add build target for file ${FILENAME}.")

  # message(STATUS "FULL_PATH: ${FULL_PATH} / TARGET_NAME: ${TARGET_NAME}")
  # message(STATUS "Adding header files from
  # ${GENERATE_CPP_OUTPUT_DIRECTORY}/${GENERATED_CODEBASE}/include/*")
  file(
    GLOB GENERATED_HEADERS
    LIST_DIRECTORIES true
    "${GENERATE_CPP_OUTPUT_DIRECTORY}/${GENERATED_CODEBASE}/include/*")
  file(GLOB GENERATED_CXX_FILES
       "${GENERATE_CPP_OUTPUT_DIRECTORY}/${GENERATED_CODEBASE}/src/*/*.cxx")
  # message(STATUS "GENERATED_HEADERS ${GENERATED_HEADERS}")
  add_executable(${TARGET_NAME} ${FULL_PATH} ${GENERATED_CXX_FILES})
  # Adds a pre-build event to the Target copying the correctionlib.so file into
  # the /lib folder in the install directory
  target_include_directories(
    ${TARGET_NAME} PRIVATE ${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS}
                           $ORIGIN/lib/ lib/)
  target_link_libraries(
    ${TARGET_NAME}
    ROOT::ROOTVecOps
    ROOT::ROOTDataFrame
    ROOT::RooFit
    ROOT::GenVector
    logging
    correctionlib
    nlohmann_json::nlohmann_json
    CROWNLIB
    ${ONNX_RUNTIME_LIB_PATH})
  add_custom_command(
    TARGET ${TARGET_NAME}
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CORRECTION_LIB_PATH}"
            ${INSTALLDIR}/lib/libcorrectionlib.so)
  # Find shared libraries next to the executable in the /lib folder
  set_target_properties(
    ${TARGET_NAME} PROPERTIES BUILD_WITH_INSTALL_RPATH FALSE
                              LINK_FLAGS "-Wl,-rpath,$ORIGIN/lib")
  # Add install target, basically just copying the executable around relative to
  # CMAKE_INSTALL_PREFIX
  install(TARGETS ${TARGET_NAME} DESTINATION ${INSTALLDIR})
  install(
    CODE "execute_process(COMMAND ${CMAKE_SOURCE_DIR}/checks/get-diff.sh ${CMAKE_SOURCE_DIR} ${ANALYSIS} ${INSTALLDIR}/diff )"
  )

endforeach()
