find_package(Doxygen REQUIRED)
find_package(Sphinx REQUIRED)

# adapted from https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/

set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/)
set(DOXYGEN_OUTPUT_DIR ${PROJECT_SOURCE_DIR}/docs/doxygen)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)
set(DOXYFILE_IN ${PROJECT_SOURCE_DIR}/docs/Doxyfile.in)
set(DOXYFILE_OUT ${PROJECT_SOURCE_DIR}/Doxyfile)

set(SPHINX_SOURCE ${PROJECT_SOURCE_DIR}/docs/sphinx_source/)
set(SPHINX_BUILD ${PROJECT_SOURCE_DIR}/docs/sphinx/)


#Replace variables inside @@ with the current values
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}) #Doxygen won't create this for us
add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                   COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                   MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                   COMMENT "Generating docs")

add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})


add_custom_target(Sphinx ALL
                COMMAND ${SPHINX_EXECUTABLE} -b html
                # Tell Breathe where to find the Doxygen output
                -Dbreathe_projects.ToBeNamedPrototype=${DOXYGEN_OUTPUT_DIR}/xml
                ${SPHINX_SOURCE} ${SPHINX_BUILD}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMENT "Generating documentation with Sphinx")