# Find Python 3
find_package(Python 3.9 REQUIRED COMPONENTS Interpreter)

# detect virtualenv and set Pip args accordingly
if(DEFINED ENV{VIRTUAL_ENV} OR DEFINED ENV{CONDA_PREFIX})
  set(_pip_args)
else()
  set(_pip_args "--user")
endif()

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
        message(FATAL_ERROR "The version of Python package ${PYPINAME} is too old (found ${PACKAGE_VERSION}, require at least ${MIN_VERSION}).")
    endif()
    message(STATUS "Found Python package ${PYPINAME} (require ${MIN_VERSION}, found ${PACKAGE_VERSION})")
endfunction()