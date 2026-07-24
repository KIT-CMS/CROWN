# Generate the C++ code
if(FRIENDS)
  set(GENERATE_CPP_INPUT_TEMPLATE
      "${CMAKE_SOURCE_DIR}/code_generation/analysis_template_friends.cxx")
else()
  set(GENERATE_CPP_INPUT_TEMPLATE
      "${CMAKE_SOURCE_DIR}/code_generation/analysis_template.cxx")
endif()
set(GENERATE_CPP_SUBSET_TEMPLATE
    "${CMAKE_SOURCE_DIR}/code_generation/subset_template.cxx")

message(STATUS "")
message(STATUS "Generate C++ code with following settings:")
message(STATUS "  Output directory: ${GENERATE_CPP_OUTPUT_DIRECTORY}")
message(STATUS "  Install directory: ${INSTALLDIR}")
message(STATUS "  Template: ${GENERATE_CPP_INPUT_TEMPLATE}")
message(STATUS "  Subset template: ${GENERATE_CPP_SUBSET_TEMPLATE}")
message(STATUS "  Analysis: ${ANALYSIS}")
message(STATUS "  Config: ${CONFIG}")
message(STATUS "  Channels: ${SCOPES}")
message(STATUS "  Shifts: ${SHIFTS}")
message(STATUS "  Samples: ${SAMPLES}")
message(STATUS "  Eras: ${ERAS}")
message(STATUS "")

file(MAKE_DIRECTORY ${GENERATE_CPP_OUTPUT_DIRECTORY})
if("${SAMPLES}" STREQUAL "all")
  execute_process(
    COMMAND
      ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/code_generation/expand_samples.py
      --analysis ${ANALYSIS} --config ${CONFIG}
    OUTPUT_VARIABLE SAMPLES_EXPANDED
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE ret
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  if(NOT ret EQUAL "0")
    message(FATAL_ERROR "SAMPLES=all expansion failed for config ${CONFIG}")
  endif()
  string(REPLACE "," ";" SAMPLES "${SAMPLES_EXPANDED}")
  message(STATUS "|> Expanded SAMPLES=all to: ${SAMPLES}")
endif()
# loop over all samples and eras and generate code for each one of them
foreach(ERA IN LISTS ERAS)
  foreach(SAMPLE IN LISTS SAMPLES)
    message("Generating code for sample ${SAMPLE} and era ${ERA}")
    message(
      "Running command: ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/generate.py --template ${GENERATE_CPP_INPUT_TEMPLATE} --subset-template ${GENERATE_CPP_SUBSET_TEMPLATE} --output ${GENERATE_CPP_OUTPUT_DIRECTORY} --analysis ${ANALYSIS} --config ${CONFIG} --scopes ${SCOPES} --shifts ${SHIFTS} --sample ${SAMPLE} --era ${ERA} --threads ${THREADS} --debug ${DEBUG_PARSED} --friends ${FRIENDS} --quantities-map ${QUANTITIESMAP}"
    )
    execute_process(
      COMMAND
        ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/generate.py --template
        ${GENERATE_CPP_INPUT_TEMPLATE} --subset-template
        ${GENERATE_CPP_SUBSET_TEMPLATE} --output
        ${GENERATE_CPP_OUTPUT_DIRECTORY} --analysis ${ANALYSIS} --config
        ${CONFIG} --scopes ${SCOPES} --shifts ${SHIFTS} --sample ${SAMPLE} --era
        ${ERA} --threads ${THREADS} --debug ${DEBUG_PARSED} --friends ${FRIENDS}
        --quantities-map ${QUANTITIESMAP}
      RESULT_VARIABLE ret
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    if(ret EQUAL "1")
      message(FATAL_ERROR "Code Generation Failed - Exiting !")
    endif()
  endforeach()
endforeach()
