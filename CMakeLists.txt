CMAKE_MINIMUM_REQUIRED(VERSION 3.15)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")

# Set a default build type if none was specified
# ------------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

message(STATUS "Install directory: ${CMAKE_INSTALL_PREFIX}")
# About this project
# ------------------------------------------------------------------------------
project(moni)
SET(VERSION_MAJOR "0")
SET(VERSION_MINOR "1")
SET(VERSION_PATCH "0")
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

# Set environment
# ------------------------------------------------------------------------------

find_package(Git)
if(GIT_FOUND)
    message("git found: ${GIT_EXECUTABLE}")
else()
    message(WARNING "git not found. Cloning of submodules will not work.")
endif()



# Configure thirdparty
# ------------------------------------------------------------------------------
set(CMAKE_INSTALL_INCLUDEDIR "include") # This is an hack because include(GUIInstallDirs) doesn't work

add_subdirectory(thirdparty)


# Configure the compiler with the appropriate flags
# ------------------------------------------------------------------------------
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  # using Clang
  include(ConfigureCompilerClang)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  include(ConfigureCompilerGcc)
else ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	message(FATAL_ERROR "Only the compiler gcc and clang are supported")
endif()


add_subdirectory(src)
add_subdirectory(utils)

# Configure pipeline for build folder
set(USE_INSTALL_PATH False)
configure_file(${PROJECT_SOURCE_DIR}/pipeline/moni.in ${PROJECT_BINARY_DIR}/moni @ONLY)

# Configure pipeline for install folder
set(USE_INSTALL_PATH True)
configure_file(${PROJECT_SOURCE_DIR}/pipeline/moni.in ${PROJECT_BINARY_DIR}/moni.install @ONLY)


install(TARGETS ms rlebwt_ms_build extend_ksw2 compress_dictionary TYPE RUNTIME)
install(TARGETS SlpEncBuild pfp_thresholds pfp_thresholds64 TYPE RUNTIME)
install(PROGRAMS ${PROJECT_BINARY_DIR}/moni.install RENAME moni TYPE BIN)
# install(TARGETS ms rlebwt_ms_build extend_ksw2 DESTINATION bin)
# install(PROGRAMS ${PROJECT_SOURCE_DIR}/pipeline/moni DESTINATION bin)