# Define Paths
set(CACHE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.cache")
set(PERSISTENT_LIB "${CACHE_DIR}/libMyDicts.so")
set(PERSISTENT_PCM "${CACHE_DIR}/libMyDicts_dict_rdict.pcm")
set(PERSISTENT_CC  "${CACHE_DIR}/libMyDicts_dict.cc")

set(SRC_HEADER  "${CMAKE_CURRENT_SOURCE_DIR}/include/dictionaries/MyDicts.hxx")
set(SRC_LINKDEF "${CMAKE_CURRENT_SOURCE_DIR}/include/dictionaries/LinkDef.hxx")

if(NOT EXISTS "${CACHE_DIR}")
    file(MAKE_DIRECTORY "${CACHE_DIR}")
endif()

# Check if Rebuild is Needed
set(NEEDS_REBUILD FALSE)
if(NOT EXISTS "${PERSISTENT_LIB}" OR NOT EXISTS "${PERSISTENT_PCM}" OR NOT EXISTS "${PERSISTENT_CC}")
    set(NEEDS_REBUILD TRUE)
else()
    file(TIMESTAMP "${SRC_HEADER}" SRC_TIME "%s")
    file(TIMESTAMP "${PERSISTENT_LIB}" BIN_TIME "%s")
    if(${SRC_TIME} GREATER ${BIN_TIME})
        set(NEEDS_REBUILD TRUE)
    endif()
endif()

# Build only if necessary
if(NEEDS_REBUILD)
    message(STATUS "Cache miss: Regenerating ROOT Dictionary and Bootstrap Lib...")

    # Read in root-config flags
    execute_process(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_C_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND root-config --libs OUTPUT_VARIABLE ROOT_L_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
    separate_arguments(C_FLAGS_LIST NATIVE_COMMAND "${ROOT_C_FLAGS}")
    separate_arguments(L_FLAGS_LIST NATIVE_COMMAND "${ROOT_L_FLAGS}")

    # Generate directly into the cache directory
    execute_process(
        COMMAND rootcling -f ${PERSISTENT_CC} 
                -I${CMAKE_CURRENT_SOURCE_DIR}
                -I${CMAKE_CURRENT_SOURCE_DIR}/include
                -c ${SRC_HEADER} ${SRC_LINKDEF}
        WORKING_DIRECTORY ${CACHE_DIR}
    )

    # Compile the bootstrap .so using the cached .cc
    execute_process(
        COMMAND g++ -shared -fPIC -o ${PERSISTENT_LIB} ${PERSISTENT_CC} 
                -I${CMAKE_CURRENT_SOURCE_DIR}/include ${C_FLAGS_LIST} ${L_FLAGS_LIST}
    )

    message(STATUS "Dictionary assets cached in .cache/")
else()
    message(STATUS "Using cached ROOT Dictionary assets")
endif()

# Deployment to Build Directory
# Copy the .so and .pcm so the CMake/Python process can load them
file(COPY "${PERSISTENT_LIB}" DESTINATION "${CMAKE_BINARY_DIR}")
file(COPY "${PERSISTENT_PCM}" DESTINATION "${CMAKE_BINARY_DIR}")

# Build MyDicts library
add_library(MyDicts SHARED ${PERSISTENT_CC})
target_include_directories(MyDicts PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(MyDicts PUBLIC ROOT::Core ROOT::RIO)

# Install dictionary and its PCM
install(TARGETS MyDicts DESTINATION ${INSTALLDIR}/lib)
install(FILES "${PERSISTENT_PCM}" DESTINATION ${INSTALLDIR}/lib)