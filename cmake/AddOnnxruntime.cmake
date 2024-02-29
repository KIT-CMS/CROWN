# Find ONNXRuntime
# first check if we have an LCG stack via LCG_VERSION environment variable
if (DEFINED ENV{LCG_VERSION})
string(REPLACE ":" ";" RUNTIME_PATH "$ENV{LD_LIBRARY_PATH}")
    message(STATUS "Found LCG stack, using it to find ONNXRuntime")
    find_library(ONNX_RUNTIME_LIB_PATH
        NAMES onnxruntime
        HINTS ${RUNTIME_PATH}
    )
    if (ONNX_RUNTIME_LIB_PATH)
        # get the real path of the library to find the include directory
        get_filename_component(ONNX_RUNTIME_LIB_PATH ${ONNX_RUNTIME_LIB_PATH} REALPATH)
        get_filename_component(ONNX_RUNTIME_INCLUDE_PATH ${ONNX_RUNTIME_LIB_PATH}/../../include REALPATH)
        message(STATUS "ONNXRuntime include path: ${ONNX_RUNTIME_INCLUDE_PATH}")
        include_directories("${ONNX_RUNTIME_INCLUDE_PATH}/core/session")
    endif()

    message(STATUS "ONNXRuntime library path: ${ONNX_RUNTIME_LIB_PATH}")
else()
    message(STATUS "No LCG stack found, not adding ONNXRuntime")
endif()