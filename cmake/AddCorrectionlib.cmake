function(install_correctionlib)
  execute_process(
    COMMAND "${Python_EXECUTABLE}" "-c"
            "import correctionlib; print(correctionlib.__version__)"
    RESULT_VARIABLE PACKAGE_NOT_FOUND
    OUTPUT_VARIABLE PACKAGE_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(${PACKAGE_NOT_FOUND} EQUAL 1)
    execute_process(
      COMMAND ${Python_EXECUTABLE} -m pip install ${_pip_args}
              git+https://github.com/cms-nanoAOD/correctionlib.git)
  endif()
  message(STATUS "Found correctionlib !")
endfunction()

# Adding correctionlib for scale factor evaluation for now the official pip
# package has some problem in the future "find_python_package(correctionlib
# correctionlib X.X)" should hopefully work
install_correctionlib()
message(STATUS "Setting up correctionlib ...")
execute_process(
  COMMAND correction config --cmake
  OUTPUT_VARIABLE CORRECTION_LIB_ARGS
  OUTPUT_STRIP_TRAILING_WHITESPACE)
string(REPLACE -Dcorrectionlib_DIR= "" CORRECTIONLIBPATH ${CORRECTION_LIB_ARGS})
# if correctionlib comes from cvmfs, change the correctionlibpath accordingly
if(${CORRECTIONLIBPATH} MATCHES "^/cvmfs/")
  message(STATUS "Setting up correctionlib from cvmfs ...")
  set(USING_CVMFS TRUE)
  find_package(correctionlib)
  find_library(CORRECTION_LIB_PATH correctionlib)
else()
  message(STATUS "Setting up correctionlib from local setup ...")
  set(USING_CVMFS FALSE)
  find_package(correctionlib REQUIRED PATHS ${CORRECTIONLIBPATH})
  set(CORRECTION_LIB_PATH "${CORRECTIONLIBPATH}/../lib/libcorrectionlib.so")
endif()
set(THREADS_PREFER_PTHREAD_FLAG ON)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
message(STATUS "Correctionlib library path: ${CORRECTION_LIB_PATH}")
