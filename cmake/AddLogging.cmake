message(STATUS "Including spdlog.")
# Build the logging library
include(ExternalProject)
ExternalProject_Add(spdlog
    PREFIX          spdlog
    GIT_REPOSITORY  https://github.com/gabime/spdlog.git
    GIT_SHALLOW     1
    GIT_TAG         v1.8.5
    CMAKE_ARGS      -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
                    -DCMAKE_BUILD_TYPE=Release
                    -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}
                    -DCMAKE_CXX_FLAGS=-fpic
    LOG_DOWNLOAD 1 LOG_CONFIGURE 1 LOG_BUILD 1 LOG_INSTALL 1
    BUILD_BYPRODUCTS ${CMAKE_INSTALL_PREFIX}/lib64/libspdlog.a
    BUILD_BYPRODUCTS ${CMAKE_INSTALL_PREFIX}/lib/libspdlog.a
)

message(STATUS "Configuring spdlog.")
# Make an imported target out of the build logging library
add_library(logging STATIC IMPORTED)
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/include") # required because the include dir must be existent for INTERFACE_INCLUDE_DIRECTORIES
include(GNUInstallDirs) # required to populate CMAKE_INSTALL_LIBDIR with lib or lib64 required for the destination of libspdlog.a
set_target_properties(logging PROPERTIES
    IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/libspdlog.a"
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/include")
add_dependencies(logging spdlog) # enforces to build spdlog before making the imported target