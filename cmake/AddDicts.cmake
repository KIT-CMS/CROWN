# Define Paths
set(CACHE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.cache")
set(PERSISTENT_LIB "${CACHE_DIR}/libMyDicts.so")
set(PERSISTENT_PCM "${CACHE_DIR}/libMyDicts_dict_rdict.pcm")
set(PERSISTENT_CC "${CACHE_DIR}/libMyDicts_dict.cc")
set(FINGERPRINT_FILE "${CACHE_DIR}/root_fingerprint.txt")

set(SRC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/include/dictionaries/MyDicts.hxx")
set(SRC_LINKDEF "${CMAKE_CURRENT_SOURCE_DIR}/include/dictionaries/LinkDef.hxx")

if(NOT EXISTS "${CACHE_DIR}")
  file(MAKE_DIRECTORY "${CACHE_DIR}")
endif()

# Fingerprint of the ROOT version used for the build. Relevant if the compiled files
# are used with a different ROOT version that can be incompatible.
execute_process(
  COMMAND root-config --prefix
  OUTPUT_VARIABLE ROOT_PREFIX
  OUTPUT_STRIP_TRAILING_WHITESPACE)
set(ROOT_FINGERPRINT "${ROOT_VERSION}|${ROOT_PREFIX}")

# Check if Rebuild is Needed
set(NEEDS_REBUILD FALSE)
if(NOT EXISTS "${PERSISTENT_LIB}"
   OR NOT EXISTS "${PERSISTENT_PCM}"
   OR NOT EXISTS "${PERSISTENT_CC}"
   OR NOT EXISTS "${FINGERPRINT_FILE}")
  set(NEEDS_REBUILD TRUE)
else()
  file(TIMESTAMP "${PERSISTENT_LIB}" BIN_TIME "%s")
  foreach(SRC_FILE "${SRC_HEADER}" "${SRC_LINKDEF}")
    file(TIMESTAMP "${SRC_FILE}" SRC_TIME "%s")
    if(${SRC_TIME} GREATER ${BIN_TIME})
      set(NEEDS_REBUILD TRUE)
    endif()
  endforeach()
  file(READ "${PERSISTENT_CC}" CACHED_CONTENT)
  string(FIND "${CACHED_CONTENT}" "${SRC_HEADER}" PATH_MATCH)
  if(${PATH_MATCH} EQUAL -1)
    message(
      STATUS
        "Cached file ${PERSISTENT_CC} includes a different header path than ${SRC_HEADER}"
    )
    set(NEEDS_REBUILD TRUE)
  endif()
  file(READ "${FINGERPRINT_FILE}" CACHED_FINGERPRINT)
  string(STRIP "${CACHED_FINGERPRINT}" CACHED_FINGERPRINT)
  if(NOT CACHED_FINGERPRINT STREQUAL ROOT_FINGERPRINT)
    message(
      STATUS
        "Cached ROOT dictionary was built against a different ROOT install "
        "(${CACHED_FINGERPRINT}) than the one currently in use (${ROOT_FINGERPRINT})"
    )
    set(NEEDS_REBUILD TRUE)
  endif()
endif()

# Build only if necessary
if(NEEDS_REBUILD)
  message(
    STATUS "Cache miss: Regenerating ROOT Dictionary and Bootstrap Lib...")

  # Remove the old fingerprint file to avoid confusion if the build fails
  # and the fingerprint is not updated
  file(REMOVE "${FINGERPRINT_FILE}")

  # Read in root-config flags
  execute_process(
    COMMAND root-config --cflags
    OUTPUT_VARIABLE ROOT_C_FLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
  execute_process(
    COMMAND root-config --libs
    OUTPUT_VARIABLE ROOT_L_FLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
  separate_arguments(C_FLAGS_LIST NATIVE_COMMAND "${ROOT_C_FLAGS}")
  separate_arguments(L_FLAGS_LIST NATIVE_COMMAND "${ROOT_L_FLAGS}")

  # Generate directly into the cache directory
  execute_process(
    COMMAND
      rootcling -f ${PERSISTENT_CC} -I${CMAKE_CURRENT_SOURCE_DIR}
      -I${CMAKE_CURRENT_SOURCE_DIR}/include -c ${SRC_HEADER} ${SRC_LINKDEF}
    WORKING_DIRECTORY ${CACHE_DIR} COMMAND_ERROR_IS_FATAL ANY)

  # Compile the bootstrap .so using the cached .cc
  execute_process(
    COMMAND
      g++ -shared -fPIC -o ${PERSISTENT_LIB} ${PERSISTENT_CC}
      -I${CMAKE_CURRENT_SOURCE_DIR}/include ${C_FLAGS_LIST} ${L_FLAGS_LIST}
      COMMAND_ERROR_IS_FATAL ANY)

  # Record which ROOT install these assets were built against
  file(WRITE "${FINGERPRINT_FILE}" "${ROOT_FINGERPRINT}")

  message(STATUS "Dictionary assets cached in .cache/")
else()
  message(STATUS "Using cached ROOT Dictionary assets")
endif()

# Deployment to Build Directory Copy the .so and .pcm so the CMake/Python
# process can load them
file(COPY "${PERSISTENT_LIB}" DESTINATION "${CMAKE_BINARY_DIR}")
file(COPY "${PERSISTENT_PCM}" DESTINATION "${CMAKE_BINARY_DIR}")

# Build MyDicts library
add_library(MyDicts SHARED ${PERSISTENT_CC})
target_include_directories(MyDicts PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(MyDicts PUBLIC ROOT::Core ROOT::RIO)

# Install dictionary and its PCM
install(TARGETS MyDicts DESTINATION ${INSTALLDIR}/lib)
install(FILES "${PERSISTENT_PCM}" DESTINATION ${INSTALLDIR}/lib)
