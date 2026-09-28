message(STATUS "Including spdlog.")
include(GNUInstallDirs) # required to populate CMAKE_INSTALL_LIBDIR with lib or
                        # lib64 required for the destination of libspdlog.a. The
                        # value is passed on to the spdlog build below, so that
                        # both agree on where libspdlog.a ends up.

if(DEFINED SPDLOG_PREBUILT_LIB AND DEFINED SPDLOG_PREBUILT_INCLUDE
   AND EXISTS "${SPDLOG_PREBUILT_LIB}" AND EXISTS "${SPDLOG_PREBUILT_INCLUDE}")
  # Reuse a spdlog already built elsewhere (e.g. by a prior CROWNLIB-only build),
  # instead of fetching and compiling it again from scratch. Useful for build hosts
  # without (reliable) network access to GitHub.
  message(STATUS "Using prebuilt spdlog at ${SPDLOG_PREBUILT_LIB} (skipping fetch+build).")
  add_library(logging STATIC IMPORTED)
  set_target_properties(
    logging
    PROPERTIES IMPORTED_LOCATION "${SPDLOG_PREBUILT_LIB}"
               INTERFACE_INCLUDE_DIRECTORIES "${SPDLOG_PREBUILT_INCLUDE}")
else()
  # Build the logging library
  include(ExternalProject)
  ExternalProject_Add(
    spdlog
    PREFIX spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_SHALLOW 1
    GIT_TAG v1.14.1
    CMAKE_ARGS -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
               -DCMAKE_BUILD_TYPE=Release
               -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}
               -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
               -DCMAKE_CXX_FLAGS=-fpic
               -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    LOG_DOWNLOAD 1
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
    BUILD_BYPRODUCTS ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/libspdlog.a)

  message(STATUS "Configuring spdlog.")
  # Make an imported target out of the build logging library
  add_library(logging STATIC IMPORTED)
  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/include"
  )# required because the include dir must be existent for
   # INTERFACE_INCLUDE_DIRECTORIES
  set_target_properties(
    logging
    PROPERTIES IMPORTED_LOCATION
               "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/libspdlog.a"
               INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/include")
  add_dependencies(logging spdlog) # enforces to build spdlog before making the
                                   # imported target
endif()
