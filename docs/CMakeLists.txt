cmake_minimum_required(VERSION 3.10)

# Create the project
project(Docs NONE)

# Search for dependencies
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
find_package(Doxygen 1.8 REQUIRED)
find_package(Sphinx 3.4 REQUIRED)
find_package(Python 3.6 REQUIRED COMPONENTS Interpreter)

execute_process(COMMAND "${Python_EXECUTABLE}" "-c" "import breathe; print(breathe.__version__)"
                RESULT_VARIABLE BREATHE_NOT_FOUND
                OUTPUT_VARIABLE BREATHE_VERSION
                OUTPUT_STRIP_TRAILING_WHITESPACE)
if(${BREATHE_NOT_FOUND} EQUAL 1)
    message(FATAL_ERROR "Failed to import breathe or get version.")
endif()
set(BREATHE_MIN_VERSION 4.29)
if(BREATHE_VERSION VERSION_LESS BREATHE_MIN_VERSION)
    message(FATAL_ERROR "The version of python-breathe is too old (found ${BREATHE_VERSION}, require at least ${BREATHE_MIN_VERSION}).")
endif()

# adapted from https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/

# Find all headers
file(GLOB_RECURSE HEADERS ${PROJECT_SOURCE_DIR}/../*.hxx)

# Write doxyfile
# Replace variables inside @@ with the current values
set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/..) # set in Doxyfile.in
set(DOXYGEN_OUTPUT_DIR ${PROJECT_BINARY_DIR}/docs/) # set in Doxyfile.in
set(DOXYFILE_IN ${PROJECT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${PROJECT_BINARY_DIR}/Doxyfile)
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY) # happens at configure time!

# Create target to generate doxygen docs
file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}) # Doxygen won't create this for us
set(DOXYGEN_INDEX_FILE ${DOXYGEN_BINARY_DIR}/docs/doxygen_index.html)
add_custom_target(Doxygen ALL
                   DEPENDS ${HEADERS} ${DOXYFILE_OUT} ${DOXYFILE_IN}
                   COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                   COMMENT "Generate doxygen index")

# Generate sphinx docs from doxygen xml
set(SPHINX_SOURCE ${PROJECT_SOURCE_DIR}/sphinx_source/)
set(SPHINX_OUTPUT_DIR ${PROJECT_BINARY_DIR}/docs)
set(SPHINX_INDEX_FILE ${PROJECT_SOURCE_DIR}/index.html)

add_custom_target(Sphinx ALL
                DEPENDS ${HEADERS}
                COMMAND ${SPHINX_EXECUTABLE} -b html
                # Tell Breathe where to find the Doxygen output
                -Dbreathe_projects.ToBeNamedPrototype=${DOXYGEN_OUTPUT_DIR}/xml
                ${SPHINX_SOURCE} ${SPHINX_OUTPUT_DIR}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMENT "Generating documentation with Sphinx")
add_dependencies(Sphinx Doxygen)
