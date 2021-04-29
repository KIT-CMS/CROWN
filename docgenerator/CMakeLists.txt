find_package(Doxygen REQUIRED)
find_package(Sphinx REQUIRED)

# adapted from https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/
# Find all headers
file(GLOB_RECURSE HEADERS ${PROJECT_SOURCE_DIR}/*.hxx)

set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/)
set(DOXYGEN_OUTPUT_DIR ${PROJECT_SOURCE_DIR}/docs/doxygen/)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)
set(DOXYFILE_IN ${PROJECT_SOURCE_DIR}/docgenerator/scripts/Doxyfile.in)
set(DOXYFILE_OUT ${PROJECT_SOURCE_DIR}/docgenerator/scripts/Doxyfile)

#Replace variables inside @@ with the current values
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}) #Doxygen won't create this for us
add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                   DEPENDS ${HEADERS}
                   COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                   MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                   COMMENT "Generating docs")

add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})

set(SPHINX_SOURCE ${PROJECT_SOURCE_DIR}/docgenerator/scripts/sphinx_source/)
set(SPHINX_BUILD ${PROJECT_SOURCE_DIR}/docs)
set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)

add_custom_target(Sphinx ALL
                DEPENDS ${HEADERS}
                COMMAND ${SPHINX_EXECUTABLE} -b html
                # Tell Breathe where to find the Doxygen output
                -Dbreathe_projects.ToBeNamedPrototype=${DOXYGEN_OUTPUT_DIR}/xml
                ${SPHINX_SOURCE} ${SPHINX_BUILD}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMENT "Generating documentation with Sphinx")

# # # Nice named target so we can run the job easily
# add_custom_target(Sphinx ALL DEPENDS ${SPHINX_INDEX_FILE})