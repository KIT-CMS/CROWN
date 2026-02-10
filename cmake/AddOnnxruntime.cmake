# Find ONNXRuntime
if(DEFINED ENV{LCG_VERSION})
  string(REPLACE ":" ";" RUNTIME_PATH "$ENV{LD_LIBRARY_PATH}")
  message(STATUS "Found LCG stack, using it to find ONNXRuntime")

  find_library(
    ONNX_RUNTIME_LIB_PATH
    NAMES onnxruntime
    HINTS ${RUNTIME_PATH})

  if(ONNX_RUNTIME_LIB_PATH)
    get_filename_component(ONNX_RUNTIME_LIB_PATH ${ONNX_RUNTIME_LIB_PATH} REALPATH)
    get_filename_component(
      ONNX_RUNTIME_INCLUDE_PATH
      ${ONNX_RUNTIME_LIB_PATH}/../../include
      REALPATH)

    include_directories("${ONNX_RUNTIME_INCLUDE_PATH}/onnxruntime")
  endif()

  message(STATUS "ONNXRuntime library path: ${ONNX_RUNTIME_LIB_PATH}")

# Conda environment (pip-installed onnxruntime)
else()
  if(NOT DEFINED ENV{CONDA_PREFIX})
    message(FATAL_ERROR "Neither LCG stack nor Conda environment detected")
  endif()

  message(STATUS "Using Conda environment at $ENV{CONDA_PREFIX}")

  # Library is always in conda lib
  find_library(
    ONNX_RUNTIME_LIB_PATH
    NAMES onnxruntime
    HINTS $ENV{CONDA_PREFIX}/lib)

  if(NOT ONNX_RUNTIME_LIB_PATH)
    message(FATAL_ERROR "onnxruntime shared library not found in Conda env")
  endif()

  # Try Conda-style headers first
  set(ONNX_RUNTIME_INCLUDE_PATH
      $ENV{CONDA_PREFIX}/include)

  if(EXISTS ${ONNX_RUNTIME_INCLUDE_PATH}/onnxruntime/core/session/onnxruntime_cxx_api.h)
    message(STATUS "Found Conda-style ONNXRuntime headers")
  else()
    # Fallback: pip-installed headers in site-packages
    file(GLOB _ORT_PIP_INCLUDE
      $ENV{CONDA_PREFIX}/lib/python*/site-packages/onnxruntime/include)

    if(NOT _ORT_PIP_INCLUDE)
      message(FATAL_ERROR
        "onnxruntime headers not found (neither Conda nor pip layout)")
    endif()

    list(GET _ORT_PIP_INCLUDE 0 ONNX_RUNTIME_INCLUDE_PATH)
    message(STATUS "Found pip-installed ONNXRuntime headers")
  endif()

  include_directories(${ONNX_RUNTIME_INCLUDE_PATH})

  message(STATUS "ONNXRuntime include path: ${ONNX_RUNTIME_INCLUDE_PATH}")
  message(STATUS "ONNXRuntime library path: ${ONNX_RUNTIME_LIB_PATH}")
endif()
