# Executable to generate samples with dummy data
add_executable(generate_sample generate_sample.cxx)
target_link_libraries(generate_sample ROOT::ROOTVecOps ROOT::ROOTDataFrame logging)
target_include_directories(generate_sample PRIVATE ${CMAKE_SOURCE_DIR})

# Add test target to generate the testing sample mark this as the setup fixture
add_test(NAME generate_sample
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
         COMMAND generate_sample 100 sample.root)
set_tests_properties(generate_sample PROPERTIES FIXTURES_SETUP generate_input)

# Generate a test for each generated target and add the setup fixture as requirement
foreach(TARGET_NAME ${TARGET_NAMES})
    message(STATUS "Add test for target ${TARGET_NAME}")
    add_test(NAME ${TARGET_NAME}
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
             COMMAND ${TARGET_NAME} tests/sample.root output_${TARGET_NAME}.root)
    set_tests_properties(${TARGET_NAME} PROPERTIES FIXTURES_REQUIRED generate_input)
endforeach()
