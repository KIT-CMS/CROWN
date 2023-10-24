cmake_minimum_required(VERSION 3.10)

# Create the project
project(Docs NONE)

# Search for dependencies
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
find_package(Doxygen 1.8 REQUIRED)
find_package(Sphinx 3.4 REQUIRED)
find_package(Python 3.6 REQUIRED COMPONENTS Interpreter)

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
        message(STATUS "The version of Python package ${PYPINAME} is too old (found ${PACKAGE_VERSION}, require at least ${MIN_VERSION}). Starting update installation.")
        execute_process(COMMAND ${Python_EXECUTABLE} -m pip install --upgrade ${PYPINAME} ${_pip_args})
        execute_process(COMMAND "${Python_EXECUTABLE}" "-c" "import ${NAME}; print(${NAME}.__version__)"
                    RESULT_VARIABLE PACKAGE_NOT_FOUND
                    OUTPUT_VARIABLE PACKAGE_VERSION
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif()
    message(STATUS "Found Python package ${PYPINAME} (require ${MIN_VERSION}, found ${PACKAGE_VERSION})")
endfunction()

find_python_package(breathe breathe 4.35)
find_python_package(sphinx_rtd_theme sphinx_rtd_theme 1.2)

# adapted from https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/

# Find all headers
file(GLOB_RECURSE HEADERS
    ${PROJECT_SOURCE_DIR}/../include/*.hxx
    ${PROJECT_SOURCE_DIR}/../include/*/*.hxx)

# Write doxyfile
# Replace variables inside @@ with the current values
set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/../src ${PROJECT_SOURCE_DIR}/../include ${PROJECT_SOURCE_DIR}/../src/utility) # set in Doxyfile.in
list(JOIN DOXYGEN_INPUT_DIR "\" \"" DOXYGEN_INPUT_DIR)
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

# now run the update_namespace_docs.py script to add any missing namespace docs
add_custom_command(TARGET Doxygen POST_BUILD
                    COMMAND ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/update_namespace_docs.py --input ${PROJECT_BINARY_DIR} --output ${PROJECT_SOURCE_DIR}/sphinx_source/c_namespaces/
                   COMMENT "Updating namespace docs")
if(ret EQUAL "1")
    message( FATAL_ERROR "Namespace documentation update failed.")
endif()
# Generate sphinx docs from doxygen xml
set(SPHINX_SOURCE ${PROJECT_SOURCE_DIR}/sphinx_source/)
set(SPHINX_OUTPUT_DIR ${PROJECT_BINARY_DIR}/docs)
set(SPHINX_INDEX_FILE ${PROJECT_SOURCE_DIR}/index.html)

add_custom_target(Sphinx ALL
                DEPENDS ${HEADERS}
                COMMAND ${SPHINX_EXECUTABLE} -b html
                # Tell Breathe where to find the Doxygen output
                -Dbreathe_projects.CROWN=${DOXYGEN_OUTPUT_DIR}/xml
                ${SPHINX_SOURCE} ${SPHINX_OUTPUT_DIR}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMENT "Generating documentation with Sphinx")
add_dependencies(Sphinx Doxygen)
