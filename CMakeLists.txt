cmake_minimum_required(VERSION 3.18)

project(MuonTomography
    VERSION 0.1.0
    DESCRIPTION "Muon scattering tomography with 4 Micromegas chambers (X/Y strips)"
    LANGUAGES CXX)

# ---- Language standard (ACTS requires C++20) ----
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)  # for clangd / IDEs

# ---- Options ----
option(MTOMO_USE_ACTS            "Build against ACTS (external project)"                ON)
option(MTOMO_USE_ACTS_CLUSTERING "Use ACTS CCA clustering instead of the built-in one" OFF)
option(MTOMO_BUILD_TESTS         "Build unit tests"                                     OFF)

# Where to get ACTS:
#   system : use a LOCAL/pre-installed ACTS only (find_package); error if absent
#   fetch  : always download+build the pinned tag remotely (FetchContent)
#   auto   : try the local install first, and only if it is not found fall back
#            to fetching the pinned tag (that fallback builds ACTS from source)
set(MTOMO_ACTS_SOURCE "auto" CACHE STRING "How to obtain ACTS: auto | system | fetch")
set_property(CACHE MTOMO_ACTS_SOURCE PROPERTY STRINGS auto system fetch)

# For the LOCAL (system) install: optionally require a specific version.
#   -DMTOMO_ACTS_VERSION=46.8.0                         -> at least 46.8.0
#   -DMTOMO_ACTS_VERSION=46.8.0 -DMTOMO_ACTS_VERSION_EXACT=ON -> exactly 46.8.0
set(MTOMO_ACTS_VERSION "" CACHE STRING "Required version for a local ACTS (empty = any)")
option(MTOMO_ACTS_VERSION_EXACT "Require the local ACTS version to match exactly" OFF)

# For the REMOTE (fetch) path: which release tag to build.
set(MTOMO_ACTS_GIT_TAG "v46.7.0" CACHE STRING "ACTS git tag to fetch (fetch / auto-fallback)")

# ---- External dependencies ----
# ROOT: to read the .root data file from your professor.
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree TreePlayer Hist)

# A pre-installed ACTS is located with find_package. Point CMake at it with
#   -DCMAKE_PREFIX_PATH=<acts-install>        (the install prefix), or
#   -DActs_DIR=<acts-install>/lib/cmake/Acts  (the dir with ActsConfig.cmake).
# On CERN lxplus, sourcing an LCG view exports the prefix for you.
#
# Fetching a pinned tag builds ACTS from source (slow, needs Boost + Eigen, and
# possibly a newer CMake than 3.18). The fetched tree lands under build/_deps/.
function(mtomo_fetch_acts)
    include(FetchContent)
    set(ACTS_BUILD_EXAMPLES         ON CACHE BOOL "" FORCE)
    set(ACTS_BUILD_UNITTESTS        OFF CACHE BOOL "" FORCE)
    set(ACTS_BUILD_INTEGRATIONTESTS OFF CACHE BOOL "" FORCE)
    set(ACTS_BUILD_BENCHMARKS       OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(
        acts
        GIT_REPOSITORY https://github.com/acts-project/acts.git
        GIT_TAG        ${MTOMO_ACTS_GIT_TAG}
        GIT_SHALLOW    TRUE
    )
    FetchContent_MakeAvailable(acts)
endfunction()

if(MTOMO_USE_ACTS)
    if(NOT MTOMO_ACTS_SOURCE MATCHES "^(auto|system|fetch)$")
        message(FATAL_ERROR
            "MTOMO_ACTS_SOURCE must be 'auto', 'system' or 'fetch' (got '${MTOMO_ACTS_SOURCE}')")
    endif()

    set(_acts_exact "")
    if(MTOMO_ACTS_VERSION_EXACT)
        set(_acts_exact EXACT)
    endif()

    set(_acts_ready FALSE)

    # 1) Try the local install (for 'system' and 'auto').
    if(MTOMO_ACTS_SOURCE STREQUAL "system")
        # Must succeed -> REQUIRED (fails fast with the "ActsConfig.cmake not
        # found" message if the prefix is wrong).
        find_package(Acts ${MTOMO_ACTS_VERSION} ${_acts_exact} REQUIRED COMPONENTS Core)
        set(_acts_ready TRUE)
    elseif(MTOMO_ACTS_SOURCE STREQUAL "auto")
        # Probe quietly so we can fall back if it is missing.
        find_package(Acts ${MTOMO_ACTS_VERSION} ${_acts_exact} QUIET COMPONENTS Core)
        if(Acts_FOUND)
            set(_acts_ready TRUE)
        endif()
    endif()

    if(_acts_ready AND Acts_FOUND)
        message(STATUS "ACTS: using local install ${Acts_VERSION} (${Acts_DIR})")
    endif()

    # 2) Fetch the pinned tag (for 'fetch', or as the 'auto' fallback).
    if(NOT _acts_ready)
        if(MTOMO_ACTS_SOURCE STREQUAL "auto")
            message(STATUS "ACTS: no local install found -> fetching ${MTOMO_ACTS_GIT_TAG} "
                           "and building from source (can take tens of minutes). "
                           "Use -DMTOMO_ACTS_SOURCE=system to require a local install instead.")
        else()
            message(STATUS "ACTS: fetching ${MTOMO_ACTS_GIT_TAG} and building from source.")
        endif()
        mtomo_fetch_acts()
        set(_acts_ready TRUE)
    endif()
endif()

# ---- Targets ----
add_subdirectory(src)


if(MTOMO_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# ---- Configuration summary ----
message(STATUS "")
message(STATUS "== ${PROJECT_NAME} ${PROJECT_VERSION} ==")
message(STATUS "  C++ standard      : ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build type        : ${CMAKE_BUILD_TYPE}")
message(STATUS "  ROOT version      : ${ROOT_VERSION}")
message(STATUS "  Use ACTS          : ${MTOMO_USE_ACTS}")
if(MTOMO_USE_ACTS)
message(STATUS "  ACTS source mode  : ${MTOMO_ACTS_SOURCE}")
if(Acts_FOUND)
message(STATUS "  ACTS resolved     : local install ${Acts_VERSION} (${Acts_DIR})")
else()
message(STATUS "  ACTS resolved     : fetched tag ${MTOMO_ACTS_GIT_TAG} (built from source)")
endif()
endif()
message(STATUS "")
