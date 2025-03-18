# read debug and optimized flags from command line
option(BUILD_CROWNLIB_ONLY "Build only the CROWNLIB library" OFF)
set(REBUILD_CROWN_LIB "false") # used for non-production mode

if(NOT DEFINED DEBUG)
  message(
    STATUS
      "No Debug mode set, activate with -DDEBUG=true --> compile with debug symbols and run code generation with debug output"
  )
  set(DEBUG "false")
endif()

if(NOT DEFINED OPTIMIZED)
  message(
    STATUS
      "No Optimization not set, building with -DOPTIMIZED=true --> slower build times but faster runtimes"
  )
  set(OPTIMIZED "true")
endif()
# Convert args to lowercase
string(TOLOWER "${DEBUG}" DEBUG_PARSED)
string(TOLOWER "${OPTIMIZED}" OPTIMIZED_PARSED)
if(DEBUG_PARSED STREQUAL "true")
  message(STATUS "Debug mode")
  set(CMAKE_BUILD_TYPE
      "Debug"
      CACHE
        STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
  )
  set(CMAKE_CXX_FLAGS_DEBUG
      "-g"
      CACHE STRING "Set default compiler flags for build type Debug")
else()
  set(DEBUG_PARSED "false")
  if(OPTIMIZED_PARSED STREQUAL "true")
    message(STATUS "Optimized mode")
    set(CMAKE_BUILD_TYPE
        "Release"
        CACHE
          STRING
          "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
    )
    set(CMAKE_CXX_FLAGS_RELEASE
        "-O3 -DNDEBUG"
        CACHE STRING "Set default compiler flags for build type Release")
    find_program(CCACHE_FOUND ccache)
    if(CCACHE_FOUND)
      message(STATUS "ccache found at ${CCACHE_FOUND}, using it for compilation")
      set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_FOUND}")
    endif()
  else()
    message(STATUS "Unoptimized mode")
    set(CMAKE_BUILD_TYPE
        "Release"
        CACHE
          STRING
          "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
    )
    set(CMAKE_CXX_FLAGS_RELEASE
        "-DNDEBUG"
        CACHE STRING "Set default compiler flags for build type Release")
  endif()
endif()
# Only parse additional args if not building only the CROWNLIB library
if(NOT BUILD_CROWNLIB_ONLY)

  if(NOT DEFINED ANALYSIS)
    message(
      FATAL_ERROR
        "Please specify the Analysis to be used with -DANALYSIS=my_analysis_name"
    )
  endif()
  # if analysis is set, check the folder to find any potential payload files to
  # be used
  file(GLOB PAYLOADS
       ${CMAKE_SOURCE_DIR}/analysis_configurations/${ANALYSIS}/payloads/*)
  if(NOT PAYLOADS)
    message(
      STATUS
        "No payload files found in ${CMAKE_SOURCE_DIR}/analysis_configurations/${ANALYSIS}/payloads/ for analysis ${ANALYSIS}"
    )
  else()
    message(
      STATUS
        "Found payload files in ${CMAKE_SOURCE_DIR}/analysis_configurations/${ANALYSIS}/payloads/ for analysis ${ANALYSIS}"
    )
  endif()

  if(NOT DEFINED CONFIG)
    message(
      FATAL_ERROR
        "Please specify the config to be used with -DCONFIG=my_config_name")
  endif()

  if(NOT DEFINED SCOPES)
    message(
      FATAL_ERROR
        "No scope specificed, set the scopes via comma seperated list e.g. -DSCOPES=et,mt,tt,em"
    )
  endif()

  if(NOT DEFINED SHIFTS)
    message(
      STATUS
        "No shifts specificed, using -DSHIFTS=all. If you want to run nominal only, use -DSHIFTS=none"
    )
    set(SHIFTS "none")
  endif()

  if(NOT DEFINED QUANTITIESMAP)
    message(
      STATUS
        "No quantities map specified, none will be used. If you want to produce friends, you have to specify quantities maps for all friend files e.g. -DQUANTITIESMAP=quantities_map_1.json,quantities_map_2.json. The input can be a comma-separated list of JSON files and/or root files (for debugging purposes)."
    )
    set(FRIENDS "false")
    set(QUANTITIESMAP "none")
  else()
    set(FRIENDS "true")
  endif()

  if(NOT DEFINED SAMPLES)
    message(
      FATAL_ERROR "Please specify the samples to be used with -DSAMPLES=samples"
    )
  endif()

  if(NOT DEFINED ERAS)
    message(FATAL_ERROR "Please specify the eras to be used with -DERAS=eras")
  endif()

  if(NOT DEFINED PRODUCTION)
    message(
      STATUS
        "No production mode set --> will rebuild the CROWNLIB library if necessary"
    )
    set(REBUILD_CROWN_LIB "true")
  endif()
  if(NOT DEFINED THREADS)
    message(
      STATUS "No threads set, using single threaded mode with -DTHREADS=1")
    set(THREADS "1")
  endif()
  string(REPLACE "," ";" ERAS "${ERAS}")
  string(REPLACE "," ";" SAMPLES "${SAMPLES}")
  message(STATUS "---------------------------------------------")
  message(STATUS "|> Set up analysis for scopes ${SCOPES}.")
  message(STATUS "|> Set up analysis for ${ANALYSIS}.")
  message(STATUS "|> Set up analysis for config ${CONFIG}.")
  message(STATUS "|> Set up analysis for samples ${SAMPLES}.")
  message(STATUS "|> Set up analysis for eras ${ERAS}.")
  message(STATUS "|> Set up analysis for shifts ${SHIFTS}.")
  message(STATUS "|> Set up analysis with ${THREADS} threads.")
  message(STATUS "|> Set up analysis with debug mode : ${DEBUG_PARSED}.")
  message(
    STATUS "|> Set up analysis with optimization mode : ${OPTIMIZED_PARSED}.")
  message(STATUS "|> generator is set to ${CMAKE_GENERATOR}")
  # Define the default compiler flags for different build types, if different
  # from the cmake defaults The build type should be set so that the correct
  # compiler flags are chosen
  message(
    STATUS
      "|> Code generation arguments:  --analysis ${ANALYSIS} --config ${CONFIG} --scopes ${SCOPES} --shifts ${SHIFTS} --samples ${SAMPLES} --eras ${ERAS} --threads ${THREADS} --debug ${DEBUG_PARSED} --friends ${FRIENDS} --quantities_map ${QUANTITIESMAP}"
  )
  message(STATUS "---------------------------------------------")
else() # if building only the CROWNLIB library
  message(STATUS "Building only the CROWNLIB library")
  message(STATUS "No additional arguments parsed")
endif()
