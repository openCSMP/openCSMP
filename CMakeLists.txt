# SPDX-FileCopyrightText: © 2026 The openCSMP project
#
# SPDX-License-Identifier: LGPL-3.0-only

cmake_minimum_required(VERSION 3.29)

# ===============================================================================================
#  openCSMP — Library, Examples and Tests
# ===============================================================================================
#
#  Prerequisites:
#    - C++23 capable compiler (g++ 13.3 or later, Clang, Intel, NVC++)
#    - MPI
#    - git
#    - On Mac OSX with SAMG: you need Rosetta (Intel x86_64)
#                            Enter Rosetta shell first
#                            arch -x86_64 /bin/zsh
#
# cd /Users/skmatthai/OpenCSMP/open-csmp-2024/build
# rm -rf *
#
# cmake .. \
#    -DCMAKE_OSX_ARCHITECTURES=x86_64 \
#    -DCMAKE_BUILD_TYPE=Release \
#    -DSAMG_ENABLE=ON \
#    -DMPI_CXX_COMPILER=/usr/local/opt/open-mpi/bin/mpicxx \
#    -DMPI_C_COMPILER=/usr/local/opt/open-mpi/bin/mpicc \
#    -DMPI_Fortran_COMPILER=/usr/local/opt/open-mpi/bin/mpifort
#
#        (Note: this requires Intel MPI buid in Rosetta shell )
#
#  Usage:
#    cmake .. -DCMAKE_BUILD_TYPE=Release -DSAMG_ENABLE=ON
#    make -j$(nproc)
#
#    cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release -DSAMG_ENABLE=ON
#    cmake --build .
#
#  Safe commits (pre-push hook runs tests automatically):
#    .git/hooks/pre-push
# ===============================================================================================


# ===============================================================================================
# Build Options
# ===============================================================================================
option(BUILD_TESTING          "Enable building tests"                              ON )
option(BUILD_EXAMPLES         "Enable building examples"                           ON )
option(BUILD_SHARED_LIBS      "Build shared libraries (OFF = static)"              OFF)
option(CSMP_WITH_IMAGE_OUTPUT "Enable image output support via libturbojpeg"       OFF)
option(SAMG_ENABLE            "Enable Fraunhofer SAMG algebraic multigrid solver"  OFF )
option(CSMP_WITH_PETSC_SOLVER "Enable PETSc solver"                                ON )


# ===============================================================================================
# macOS / Rosetta 2 — must be set BEFORE project()
# use homebrew in the rosetta bash-3.2$ shell to install!
# ===============================================================================================
if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64" AND SAMG_ENABLE)
    set(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE STRING "Build architecture for macOS" FORCE)
    message(STATUS "Configuring for x86_64 (Rosetta 2) to support SAMG")
endif()


# ===============================================================================================
# Project Declaration
# ===============================================================================================
project(openCSMP
        VERSION     23
        DESCRIPTION "openCSMP library"
        LANGUAGES   CXX
)

set(CMAKE_CXX_STANDARD          23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)



# AddressSanitizer — toggle with -DCSMP_WITH_ASAN=ON when configuring.
# Applies globally (before any target is defined)
option(CSMP_WITH_ASAN "Build with AddressSanitizer" OFF)
if(CSMP_WITH_ASAN)
    add_compile_options(-fsanitize=address -g -O1 -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address)
endif()


# -----------------------------------------------------------------------
# External support libraries root — defined once to avoid path typos
#
# Moved up here so every dependency-discovery step that needs it (CGNS, SAMG) can rely on it being set.
# -----------------------------------------------------------------------
set(CSMP_EXTERNAL_LIBS "${CMAKE_CURRENT_SOURCE_DIR}/../external_support_libs"
        CACHE PATH "Path to external support libraries (samg, cgns, rapidxml etc.)")
message(STATUS "External support libs root: ${CSMP_EXTERNAL_LIBS}")
message(STATUS "  -> SAMG expected at: ${CSMP_EXTERNAL_LIBS}/samg/linux")
message(STATUS "  -> CGNS expected at:      ${CSMP_EXTERNAL_LIBS}/cgns")
message(STATUS "  -> RapidXML expected at:  ${CSMP_EXTERNAL_LIBS}/rapidxml")


# ===============================================================================================
# Dependencies (library discovery section) — Fetched at configure time
# ===============================================================================================

# 1. plf_colony (header-only, essential)
#    A fast, memory-efficient container for non-contiguous data.
# ------------------------------------------------------------------------------------------------
include(FetchContent)
FetchContent_Declare(
        plf_colony
        GIT_REPOSITORY https://github.com/mattreecebentley/plf_colony.git
    GIT_TAG        master
    GIT_SHALLOW    TRUE
    CONFIGURE_COMMAND ""
    BUILD_COMMAND     ""
    INSTALL_COMMAND   ""
)
FetchContent_Populate(plf_colony)
message(STATUS "plf_colony source: ${plf_colony_SOURCE_DIR}")


# 2. Eigen (header-only, essential)
#    Linear algebra library. Tests and docs are disabled to keep build times short.
# ------------------------------------------------------------------------------------------------
if ( APPLE )
    # -----------------------------------------------------------------------
    # macOS — Eigen via Homebrew
    # Homebrew Eigen does not ship CMake config files so we use the
    # include path directly rather than find_package.
    # -----------------------------------------------------------------------
    execute_process(
        COMMAND brew --prefix eigen
        OUTPUT_VARIABLE HOMEBREW_EIGEN_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )

    if ( HOMEBREW_EIGEN_PREFIX )
        set( EIGEN_INCLUDE_DIR "${HOMEBREW_EIGEN_PREFIX}/include/eigen3" )
        if ( EXISTS "${EIGEN_INCLUDE_DIR}/Eigen/Dense" )
            message( STATUS "Eigen found via Homebrew: ${EIGEN_INCLUDE_DIR}" )
        else()
            message( FATAL_ERROR
                "Eigen headers not found at ${EIGEN_INCLUDE_DIR}\n"
                "Install with: brew install eigen" )
        endif()
    else()
        message( FATAL_ERROR
            "Homebrew not found — cannot locate Eigen\n"
            "Install Homebrew and then: brew install eigen" )
    endif()

else()

set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_DOC     OFF CACHE BOOL "" FORCE)

# Temporarily override BUILD_TESTING so Eigen's CMakeLists.txt cannot
# register its tests even if it checks BUILD_TESTING directly
set(BUILD_TESTING OFF)

FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG        master
        GIT_SHALLOW    TRUE    # only fetch latest commit — much faster download
)
FetchContent_MakeAvailable(eigen)
message(STATUS "Eigen source: ${eigen_SOURCE_DIR}")
if(TARGET Eigen3::Eigen)
    message(STATUS "Eigen3::Eigen target found")
else()
    message(WARNING "Eigen3::Eigen CMake target not found — will fall back to manual include path")
endif()

# Restore BUILD_TESTING for openCSMP's own tests
set(BUILD_TESTING ON)

endif() # Eigen


# ===============================================================================================
# Dependencies — System Libraries (found via pkg-config / find_package)
# ===============================================================================================

# Find pkg-config once — used by PETSc and libturbojpeg
# install with: sudo apt install pkg-config
find_package(PkgConfig REQUIRED)

# 3. MPI (essential)
# ------------------------------------------------------------------------------------------------
if ( APPLE )
    if ( SAMG_ENABLE AND CMAKE_OSX_ARCHITECTURES STREQUAL "x86_64" )
        # SAMG requires x86_64 — use Rosetta Homebrew MPI at /usr/local/
        set( MPI_CXX_COMPILER     "/usr/local/opt/open-mpi/bin/mpicxx"  )
        set( MPI_C_COMPILER       "/usr/local/opt/open-mpi/bin/mpicc"   )
        set( MPI_Fortran_COMPILER "/usr/local/opt/open-mpi/bin/mpifort" )
        message( STATUS "Using x86_64 MPI for SAMG/Rosetta build" )
    else()
        # Native arm64 build — use arm64 Homebrew MPI at /opt/homebrew/
        execute_process(
            COMMAND brew --prefix open-mpi
            OUTPUT_VARIABLE HOMEBREW_MPI_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        set( MPI_CXX_COMPILER     "${HOMEBREW_MPI_PREFIX}/bin/mpicxx"  )
        set( MPI_C_COMPILER       "${HOMEBREW_MPI_PREFIX}/bin/mpicc"   )
        set( MPI_Fortran_COMPILER "${HOMEBREW_MPI_PREFIX}/bin/mpifort" )
        message( STATUS "Using arm64 MPI: ${HOMEBREW_MPI_PREFIX}" )
    endif()
endif()

find_package( MPI REQUIRED )
message( STATUS "MPI found: ${MPI_CXX_LIBRARIES}" )

# 4. PETSc (optional)
#    Install with: sudo apt install petsc-dev
# ------------------------------------------------------------------------------------------------
if(CSMP_WITH_PETSC_SOLVER)
    pkg_check_modules(PETSC REQUIRED PETSc)
    message(STATUS "PETSc found: ${PETSC_LIBRARIES}")
endif()


# 5. libjpeg (optional)
#    Install with: sudo apt install jpeg or using homebrew (OSX)
# ------------------------------------------------------------------------------------------------
if ( CSMP_WITH_IMAGE_OUTPUT )
    if ( APPLE AND CMAKE_OSX_ARCHITECTURES STREQUAL "x86_64" )
        # x86_64 Rosetta — jpeg installed via /usr/local/ Homebrew
        set( JPEG_INCLUDE_DIRS "/usr/local/opt/jpeg/include" )
        set( JPEG_LIBRARIES    "/usr/local/opt/jpeg/lib/libjpeg.dylib" )
        if ( NOT EXISTS "${JPEG_INCLUDE_DIRS}/jpeglib.h" )
            # Fall back to Cellar path if opt symlink not present
            set( JPEG_INCLUDE_DIRS "/usr/local/Cellar/jpeg/10/include" )
            set( JPEG_LIBRARIES    "/usr/local/Cellar/jpeg/10/lib/libjpeg.dylib" )
        endif()
        message( STATUS "JPEG x86_64: ${JPEG_INCLUDE_DIRS}" )

    elseif ( APPLE )
        # arm64 native build
        execute_process(
            COMMAND brew --prefix jpeg
            OUTPUT_VARIABLE HOMEBREW_JPEG_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        set( JPEG_INCLUDE_DIRS "${HOMEBREW_JPEG_PREFIX}/include" )
        set( JPEG_LIBRARIES    "${HOMEBREW_JPEG_PREFIX}/lib/libjpeg.dylib" )
        message( STATUS "JPEG arm64: ${JPEG_INCLUDE_DIRS}" )

    else()
        # Linux
        find_package( JPEG REQUIRED )
    endif()
endif()


# ===============================================================================================
# Dependencies — Third-Party Source Libraries (built as part of this project)
# ===============================================================================================

# 6. IAPS H2O EOS (optional, built from source if present)
#    Located in: thirdparty/iaps_h2o_eos/
# ------------------------------------------------------------------------------------------------
file(GLOB_RECURSE IAPS_SRC "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/iaps_h2o_eos/*.cpp")
list(FILTER IAPS_SRC EXCLUDE REGEX "/deprecated/")
if(IAPS_SRC)
    add_library(iaps_h2o_eos STATIC ${IAPS_SRC})
    target_include_directories(iaps_h2o_eos PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/iaps_h2o_eos
    )
    message(STATUS "IAPS H2O EOS found: ${IAPS_SRC}")
else()
    message(STATUS "IAPS H2O EOS not found — skipping")
endif()


# 7. CGNS (interface/header-only, from external_support_libs)
# -----------------------------------------------------------------------
# — install with: sudo apt install libcgns-dev
#
# CGNS — CFD General Notation System
#
#   macOS:  brew install cgns
#   Linux:  sudo apt install libcgns-dev
#   Fallback: headers from external_support_libs/cgns (no linking)
#
#   Discovery is attempted in order:
#     1. find_package  (works if CMake config files are present)
#     2. pkg-config    (works on most Linux systems)
#     3. Manual search (Homebrew path on macOS)
#     4. Direct system search
#     5. Fallback      (header-only from external_support_libs)
# -----------------------------------------------------------------------

set( CGNS_FOUND_USABLE FALSE )

# 1. Try find_package first
find_package( CGNS QUIET )
if ( CGNS_FOUND )
    message( STATUS "CGNS found via find_package: ${CGNS_LIBRARIES}" )
    set( CGNS_FOUND_USABLE TRUE )

# 2. Try pkg-config
elseif ( PKG_CONFIG_FOUND )
    pkg_check_modules( CGNS QUIET cgns )
    if ( CGNS_FOUND )
        message( STATUS "CGNS found via pkg-config: ${CGNS_LIBRARIES}" )
        set( CGNS_FOUND_USABLE TRUE )
    endif()
endif()

# 3. Try Homebrew manually on macOS
if ( APPLE AND NOT CGNS_FOUND_USABLE )
    execute_process(
        COMMAND brew --prefix cgns
        OUTPUT_VARIABLE HOMEBREW_CGNS_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if ( HOMEBREW_CGNS_PREFIX AND EXISTS "${HOMEBREW_CGNS_PREFIX}/include/cgnslib.h" )
        message( STATUS "CGNS found via Homebrew: ${HOMEBREW_CGNS_PREFIX}" )
        set( CGNS_INCLUDE_DIRS "${HOMEBREW_CGNS_PREFIX}/include" )
        set( CGNS_LIBRARIES    "${HOMEBREW_CGNS_PREFIX}/lib/libcgns.dylib" )
    set( CGNS_FOUND_USABLE TRUE )
endif()
endif()

# 4. Try a direct system search (covers apt-installed libcgns-dev on Ubuntu/Debian,
#    which ships no CMake config and no pkg-config file — just headers and a .so
#    dropped into standard system paths).
if ( NOT CGNS_FOUND_USABLE )
    find_path( CGNS_SYSTEM_INCLUDE_DIR NAMES cgnslib.h )
    find_library( CGNS_SYSTEM_LIBRARY NAMES cgns )
    if ( CGNS_SYSTEM_INCLUDE_DIR AND CGNS_SYSTEM_LIBRARY )
        message( STATUS "CGNS found via direct system search: ${CGNS_SYSTEM_LIBRARY}" )
        set( CGNS_INCLUDE_DIRS "${CGNS_SYSTEM_INCLUDE_DIR}" )
        set( CGNS_LIBRARIES    "${CGNS_SYSTEM_LIBRARY}" )
        set( CGNS_FOUND_USABLE TRUE )
    endif()
endif()

# 5. Fallback — use headers from external_support_libs (no linking)
if ( NOT CGNS_FOUND_USABLE )
    message( WARNING
        "CGNS not found via package manager — "
        "falling back to external_support_libs headers (no linking).\n"
        "  macOS:  brew install cgns\n"
        "  Linux:  sudo apt install libcgns-dev" )
    set( CGNS_INCLUDE_DIRS "${CSMP_EXTERNAL_LIBS}/cgns" )
    set( CGNS_LIBRARIES    "" )
endif()

# Create an interface target for consistent usage
add_library( cgns_iface INTERFACE )
target_include_directories( cgns_iface INTERFACE ${CGNS_INCLUDE_DIRS} )
if ( CGNS_LIBRARIES )
    target_link_libraries( cgns_iface INTERFACE ${CGNS_LIBRARIES} )
endif()


# 8. RapidXML (interface/header-only, from external_support_libs)
# ------------------------------------------------------------------------------------------------
# install via:  sudo apt install rapidxml-dev
# RapidXML (header-only, essential)
#   No apt package available — fetched automatically
FetchContent_Declare(
        rapidxml
        GIT_REPOSITORY https://github.com/discord/rapidxml.git
        GIT_TAG        master
)
FetchContent_MakeAvailable(rapidxml)
message(STATUS "RapidXML source: ${rapidxml_SOURCE_DIR}")


# ===============================================================================================
# Source and Header Discovery
# ===============================================================================================

# Collect all .cpp sources for the core library
file(GLOB_RECURSE LIBRARY_SOURCES CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/source/modules/**/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/modules/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/non-GIT-files/*.cpp"
)

# Collect all headers (used to derive include directories)
file(GLOB_RECURSE LIBRARY_HEADERS CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/source/modules/**/*.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/modules/**/*.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/**/*.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/**/*.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/non-GIT-files/*.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/non-GIT-files/*.hpp"
)

# Exclude deprecated files from both lists
list(FILTER LIBRARY_SOURCES EXCLUDE REGEX "/deprecated/")
list(FILTER LIBRARY_HEADERS EXCLUDE REGEX "/deprecated/")

# Exclude SAMG sources
# (SAMG sources are added back conditionally in the SAMG section below)
list(FILTER LIBRARY_SOURCES EXCLUDE REGEX ".*[Ss][Aa][Mm][Gg].*")

# All of hydrothermal_system_modelling (includes + modules) — ALWAYS discovered.
# What gets *compiled into a target* is controlled entirely by each consuming target's
# own option (BUILD_EXAMPLES for openCSMP-examples, BUILD_TESTING for CVFEM_tests_main).
#
#   cvfe_hydrothermal              — supporting library code, folded directly into openCSMP
#                                    itself (same treatment as the main source/ tree, inheriting
#                                    the project's flat-include model).
#   hydrothermal_system_simulator — supporting code for the examples and tests and future simulators.
#                                    NOT folded into openCSMP; consumed directly by openCSMP-examples
#                                    and CVFEM_tests_main.
set(CVFEM_HYDROTHERMAL_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/hydrothermal_system_modelling")

file(GLOB_RECURSE CVFEM_HYDROTHERMAL_SOURCES CONFIGURE_DEPENDS
        "${CVFEM_HYDROTHERMAL_ROOT}/modules/cvfe_hydrothermal/*.cpp"
)
file(GLOB_RECURSE CVFEM_HYDROTHERMAL_HEADERS CONFIGURE_DEPENDS
        "${CVFEM_HYDROTHERMAL_ROOT}/includes/cvfe_hydrothermal/*.h"
        "${CVFEM_HYDROTHERMAL_ROOT}/includes/cvfe_hydrothermal/*.hpp"
)
list(FILTER CVFEM_HYDROTHERMAL_SOURCES EXCLUDE REGEX "/deprecated/")
list(FILTER CVFEM_HYDROTHERMAL_HEADERS EXCLUDE REGEX "/deprecated/")

if(CVFEM_HYDROTHERMAL_SOURCES)
    list(APPEND LIBRARY_SOURCES ${CVFEM_HYDROTHERMAL_SOURCES})
    list(APPEND LIBRARY_HEADERS ${CVFEM_HYDROTHERMAL_HEADERS})
    message(STATUS "CVFEM Hydrothermal module (cvfe_hydrothermal) enabled — ${CVFEM_HYDROTHERMAL_ROOT}/modules/cvfe_hydrothermal")
else()
    message(WARNING "No sources found under ${CVFEM_HYDROTHERMAL_ROOT}/modules/cvfe_hydrothermal")
endif()

# Hydrothermal-System-Simulator (HSS) sources/headers — ALWAYS discovered, same as
# cvfe_hydrothermal above. NOT merged into LIBRARY_SOURCES/LIBRARY_HEADERS (i.e. not folded into
# the openCSMP library itself) — instead consumed directly by whichever targets need them
# (openCSMP-examples, CVFEM_tests_main), each gated by its own option (BUILD_EXAMPLES,
# BUILD_TESTING respectively). GLOB_RECURSE already walks subfolders under
# modules/hydrothermal_system_simulator and includes/hydrothermal_system_simulator, so nested
# subdirectories are picked up here too.
file(GLOB_RECURSE CVFEM_HSS_SOURCES CONFIGURE_DEPENDS
        "${CVFEM_HYDROTHERMAL_ROOT}/modules/hydrothermal_system_simulator/*.cpp"
)
file(GLOB_RECURSE CVFEM_HSS_HEADERS CONFIGURE_DEPENDS
        "${CVFEM_HYDROTHERMAL_ROOT}/includes/hydrothermal_system_simulator/*.h"
        "${CVFEM_HYDROTHERMAL_ROOT}/includes/hydrothermal_system_simulator/*.hpp"
)
list(FILTER CVFEM_HSS_SOURCES EXCLUDE REGEX "/deprecated/")
list(FILTER CVFEM_HSS_HEADERS EXCLUDE REGEX "/deprecated/")
if(NOT CVFEM_HSS_SOURCES)
    message(WARNING "No sources found under ${CVFEM_HYDROTHERMAL_ROOT}/modules/hydrothermal_system_simulator")
endif()

# This split is a safety net: CVFEM_HSS_LIB_SOURCES 
# Used for the future Hydrothermal Systems Simulator(HSS) driver source.
# When adding back main(), it will automatically land in CVFEM_HSS_MAIN_SOURCES.
set(CVFEM_HSS_MAIN_SOURCES "")
set(CVFEM_HSS_LIB_SOURCES  "")
foreach(_hss_src ${CVFEM_HSS_SOURCES})
    file(STRINGS "${_hss_src}" _hss_main_match REGEX "int[ \t]+main[ \t]*\\(")
    if(_hss_main_match)
        list(APPEND CVFEM_HSS_MAIN_SOURCES "${_hss_src}")
    else()
        list(APPEND CVFEM_HSS_LIB_SOURCES "${_hss_src}")
    endif()
endforeach()

if(CVFEM_HSS_MAIN_SOURCES)
    message(STATUS "CVFEM HSS driver source(s) with main(): ${CVFEM_HSS_MAIN_SOURCES}")
endif()

# Check that required non-GIT variable set headers exist
set(VARIABLE_INDEX_SETS
        ${CMAKE_CURRENT_SOURCE_DIR}/non-GIT-files/VariableSet_TracerTransfer.h
        ${CMAKE_CURRENT_SOURCE_DIR}/non-GIT-files/VariableSet_TwoPhaseFlow.h
        ${CMAKE_CURRENT_SOURCE_DIR}/non-GIT-files/VariableSet_CO2GeoSequestration.h
)
foreach(v_file ${VARIABLE_INDEX_SETS})
    if(NOT EXISTS ${v_file})
        message(WARNING "Required VariableSet file missing: ${v_file}")
    endif()
endforeach()


# ===============================================================================================
# Include Directory Extraction
#
#   Derives a flat list of unique directories from all discovered headers and sources.
#   This supports the project's flat-include model where headers are included by filename
#   only (e.g. #include "Node.h") rather than by relative path.
# ===============================================================================================
function(get_unique_directories result file_list)
    set(dir_list "")
    foreach(file_path ${file_list})
        if(EXISTS ${file_path})
            get_filename_component(dir_path ${file_path} DIRECTORY)
            list(APPEND dir_list ${dir_path})
        endif()
    endforeach()
    list(REMOVE_DUPLICATES dir_list)
    set(${result} ${dir_list} PARENT_SCOPE)
endfunction()

get_unique_directories(ALL_SUBDIRECTORIES "${LIBRARY_SOURCES};${LIBRARY_HEADERS}")

# Always include top-level roots for flat includes
list(APPEND ALL_SUBDIRECTORIES
        ${CMAKE_CURRENT_SOURCE_DIR}/source/includes
        ${CMAKE_CURRENT_SOURCE_DIR}/source/modules
        ${CMAKE_CURRENT_SOURCE_DIR}/source/non-GIT-files
        ${CMAKE_CURRENT_SOURCE_DIR}/non-GIT-files
)
list(REMOVE_DUPLICATES ALL_SUBDIRECTORIES)


# ===============================================================================================
# Core Library Target: openCSMP
# ===============================================================================================
add_library(openCSMP ${LIBRARY_SOURCES})

# Public include directories — propagated to all targets linking against openCSMP
target_include_directories(openCSMP PUBLIC
        $<BUILD_INTERFACE:${ALL_SUBDIRECTORIES}>
        $<BUILD_INTERFACE:${eigen_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include/openCSMP>
)

# Compile definitions
target_compile_definitions(openCSMP PUBLIC
        PYRAMID_TRIANGULAR_FACETS
)


# Core link dependencies
target_link_libraries(openCSMP PUBLIC
        MPI::MPI_CXX
)

# plf_colony — header-only, only plf_colony.h is needed
# plf_tools.h is an internal detail and must never be compiled directly
target_include_directories(openCSMP PUBLIC
        ${plf_colony_SOURCE_DIR}
)

# Eigen — linking differs between macOS (direct path) and Linux (FetchContent)
if ( APPLE )
    # Direct include path from Homebrew
    target_include_directories( openCSMP PUBLIC
        $<BUILD_INTERFACE:${EIGEN_INCLUDE_DIR}>
    )
    message( STATUS "Eigen include path: ${EIGEN_INCLUDE_DIR}" )
elseif ( TARGET Eigen3::Eigen )
    # FetchContent provided a proper CMake target
    target_link_libraries( openCSMP PUBLIC Eigen3::Eigen )
    message( STATUS "Eigen3::Eigen target found ✓" )
else()
    # FetchContent fallback — no CMake target
    target_include_directories( openCSMP PUBLIC
        $<BUILD_INTERFACE:${eigen_SOURCE_DIR}>
    )
endif()

# IAPS H2O EOS (optional)
if(TARGET iaps_h2o_eos)
    target_link_libraries(openCSMP PUBLIC iaps_h2o_eos)
endif()

# CGNS linking (discovery happened earlier)
if(CGNS_FOUND)
    target_link_libraries(openCSMP PUBLIC ${CGNS_LIBRARIES})
    target_include_directories(openCSMP PUBLIC ${CGNS_INCLUDE_DIRS})
else()
    target_link_libraries(openCSMP PUBLIC cgns)  # fallback interface library
endif()

# RapidXML — header-only, include path only needed
target_include_directories(openCSMP PUBLIC ${rapidxml_SOURCE_DIR})

# PETSc (optional)
if(CSMP_WITH_PETSC_SOLVER)
    target_include_directories(openCSMP PUBLIC ${PETSC_INCLUDE_DIRS})
    target_link_libraries(openCSMP PUBLIC ${PETSC_LIBRARIES})
    # PETSC_CFLAGS_OTHER contains raw -DFOO flags — must use compile_options not compile_definitions
    target_compile_options(openCSMP PRIVATE ${PETSC_CFLAGS_OTHER})
    target_compile_definitions(openCSMP PUBLIC CSMP_WITH_PETSC_SOLVER)
    endif()

# libjpeg (optional)
if(CSMP_WITH_IMAGE_OUTPUT)
    target_include_directories(openCSMP PRIVATE ${JPEG_INCLUDE_DIRS})
    target_link_libraries(openCSMP PRIVATE ${JPEG_LIBRARIES})
    target_compile_definitions(openCSMP PRIVATE CSMP_WITH_IMAGE_OUTPUT)
endif()

# ===============================================================================================
# SAMG Solver (commercial — Fraunhofer/SCAPOS)
#
#   Algebraic multigrid solver for systems of equations.
#   Requires prebuilt binaries in: ../external_support_libraries/samg/
#
#   Two modes:
#     SAMG_LEGACY=ON  (default) — original SAMG integration with SAMG_MULTIPLE_INSTANCES
#     SAMG_LEGACY=OFF           — modern SAMG integration, single instance
#
#   Toggle with:
#     cmake -DSAMG_ENABLE=ON ..
#     cmake -DSAMG_ENABLE=ON -DSAMG_LEGACY=OFF ..
# ===============================================================================================
option(SAMG_LEGACY "Use legacy SAMG integration (SAMG_MULTIPLE_INSTANCES)" OFF)

if(SAMG_ENABLE)
    message(STATUS "Enabling SAMG Solver...")

    # Locate prebuilt libraries first — no point proceeding without them
    set(SAMG_INTERFACE_REQUIRED FALSE)
    if(APPLE)
        set(SAMG_LIB_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../external_support_libs/samg/macosx")
        set(SAMG_LIB      "${SAMG_LIB_PATH}/libamg_mult.dylib")
        set(OMP_LIB       "${SAMG_LIB_PATH}/libiomp5.dylib")
    elseif(SAMG_LEGACY)
        find_library(SAMG_LIB NAMES amg_mult
            HINTS "${CSMP_EXTERNAL_LIBS}/samg/linux")
        find_library(OMP_LIB NAMES iomp5
            HINTS "${CSMP_EXTERNAL_LIBS}/samg/linux")
    else()
        find_library(SAMG_LIB NAMES samg
            HINTS "${CSMP_EXTERNAL_LIBS}/samg/linux")
        find_library(SAMG_INTERFACE_LIB NAMES samg_interface
            HINTS "${CSMP_EXTERNAL_LIBS}/samg/linux")
        find_library(OMP_LIB NAMES iomp5
            HINTS "${CSMP_EXTERNAL_LIBS}/samg/linux")
        set(SAMG_INTERFACE_REQUIRED TRUE)
    endif()

    if(SAMG_LIB AND OMP_LIB AND (SAMG_INTERFACE_LIB OR NOT SAMG_INTERFACE_REQUIRED))
        # Find SAMG source files
        file(GLOB_RECURSE SAMG_SRC
                "${CMAKE_CURRENT_SOURCE_DIR}/source/modules/integration/solver/SAMG_*.h"
                "${CMAKE_CURRENT_SOURCE_DIR}/source/modules/integration/solver/SAMG_*.cpp"
                "${CMAKE_CURRENT_SOURCE_DIR}/source/modules/integration/solver/SAMG_profiles/*.cpp"
        )

        # Common definitions for all SAMG builds
        target_compile_definitions(openCSMP PUBLIC
                CSMP_WITH_SAMG_SOLVER
                SAMG_UNIX_LINUX
                SAMG_LCASE_USCORE
        )

        # Legacy mode adds multiple-instance support and legacy API definitions
        if(SAMG_LEGACY)
            target_compile_definitions(openCSMP PUBLIC
                    SAMG_MULTIPLE_INSTANCES
                    SAMG_LEGACY
            )
            message(STATUS "SAMG mode: LEGACY (SAMG_MULTIPLE_INSTANCES enabled)")
        else()
            message(STATUS "SAMG mode: MODERN (single instance)")
        endif()

        # SAMG source files
        file(GLOB_RECURSE SAMG_SRC
                "${CMAKE_CURRENT_SOURCE_DIR}/source/modules/integration/solver/SAMG_*.cpp"
                "${CMAKE_CURRENT_SOURCE_DIR}/source/modules/integration/solver/SAMG_profiles/*.cpp"
        )
        if(SAMG_SRC)
            target_sources(openCSMP PRIVATE ${SAMG_SRC})
            target_include_directories(openCSMP PUBLIC
                    ${CSMP_EXTERNAL_LIBS}/samg # for samg.h
                    ${CMAKE_CURRENT_SOURCE_DIR}/source/modules/integration/solver/SAMG_profiles
            )
        else()
            message(WARNING "SAMG enabled but no SAMG source files found")
        endif()

        # SAMG_INTERFACE_LIB is empty/unset on APPLE (not searched for there — see above), so this
        # is harmless on macOS; CMake ignores the empty expansion.
        target_link_libraries(openCSMP PUBLIC ${SAMG_LIB} ${SAMG_INTERFACE_LIB} ${OMP_LIB})
        message(STATUS "SAMG library:   ${SAMG_LIB}")
        message(STATUS "SAMG interface library: ${SAMG_INTERFACE_LIB}")
        message(STATUS "OpenMP library: ${OMP_LIB}")

    else()
        message(WARNING "SAMG requested but libraries not found — disabling SAMG support")
        set(SAMG_ENABLE OFF)
    endif()
endif()


# ===============================================================================================
# Link-Time Optimisation (LTO / IPO)
# ===============================================================================================
include(CheckIPOSupported)
check_ipo_supported(RESULT IPO_SUPPORTED)
if(IPO_SUPPORTED)
    set_target_properties(openCSMP PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
    message(STATUS "Link-time optimisation enabled")
endif()


# ===============================================================================================
# Shared Include Directories for Examples and Tests
#
#   Both the examples and tests executables need the same broad set of include paths.
#   Defined once here as a list to avoid duplication.
# ===============================================================================================
set(CSMP_CONSUMER_INCLUDE_DIRS
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/analysis"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/discretisation"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/eos"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/integration"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/integration/cvfe_hydrothermal"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/integration/finite_elements"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/integration/finite_volumes"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/integration/finite_volumes/DES_transport_scheme"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/integration/finite_volumes/generic_transport_scheme"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/integration/finite_volumes/legacy_transport_scheme"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/integration/pde_integrators"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/integration/pde_operators"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/integration/pde_operators/analytical"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/integration/pde_operators/numerical"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/integration/solver"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/interfaces"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/interfaces/ANSYS"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/interfaces/CGNS"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/interfaces/Eclipse"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/interfaces/Fracman"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/interfaces/GeoModeller"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/interfaces/gOcad_SKUA"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/interfaces/JPEG"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/interfaces/Maple"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/interfaces/Matlab"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/interfaces/Rhinoceros"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/interfaces/Tecplot"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/interfaces/VTK"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/math"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/model"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/property_manipulation"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/regular_grids"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/run_control"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/thmc"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/thmc/fractures"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/thmc/geochemistry"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/thmc/mechanics"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/thmc/single_phase_flow"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/thmc/two_phase_flow"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/utilities"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/unit"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration"
        "${CMAKE_CURRENT_SOURCE_DIR}/examples/source/interrelations"
        "${CMAKE_CURRENT_SOURCE_DIR}/examples/source/visitors"
        ${plf_colony_SOURCE_DIR}
        ${eigen_SOURCE_DIR}
)
# Append optional include dirs only if the feature is enabled
if(CSMP_WITH_IMAGE_OUTPUT)
    list(APPEND CSMP_CONSUMER_INCLUDE_DIRS ${JPEG_INCLUDE_DIRS})
endif()

# Add Eigen include to consumer targets
if ( APPLE )
    list( APPEND CSMP_CONSUMER_INCLUDE_DIRS ${EIGEN_INCLUDE_DIR} )
else()
    list( APPEND CSMP_CONSUMER_INCLUDE_DIRS ${eigen_SOURCE_DIR} )
endif()



# ===============================================================================================
# Examples Executable: openCSMP-examples
# ===============================================================================================
if(BUILD_EXAMPLES)
    file(GLOB_RECURSE EXAMPLES_SRC
            "${CMAKE_CURRENT_SOURCE_DIR}/examples/source/*.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/examples/source/**/*.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/examples/source/*.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/examples/source/**/*.h"
    )
    # Include test infrastructure so examples can demonstrate test cases
    list(APPEND EXAMPLES_SRC
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/TestSuite.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/TestSuite.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/Test.cpp"
    ${CVFEM_HSS_LIB_SOURCES}
    )
    list(FILTER EXAMPLES_SRC EXCLUDE REGEX "/deprecated/")

    if(EXAMPLES_SRC)
        add_executable(openCSMP-examples ${EXAMPLES_SRC})

        target_include_directories(openCSMP-examples PRIVATE
                ${CSMP_CONSUMER_INCLUDE_DIRS}
                "${CMAKE_CURRENT_SOURCE_DIR}"
        )

        if(CVFEM_HSS_SOURCES OR CVFEM_HSS_HEADERS)
            get_unique_directories(CVFEM_HSS_INCLUDE_DIRS "${CVFEM_HSS_SOURCES};${CVFEM_HSS_HEADERS}")

            target_include_directories(openCSMP-examples PRIVATE
                    ${CVFEM_HSS_INCLUDE_DIRS}
                    "${CVFEM_HYDROTHERMAL_ROOT}/includes/hydrothermal_system_simulator"
            )
        endif()

        # Links against the core library (inherits all its PUBLIC dependencies)
        # jpeg is linked separately for any direct JPEG calls in examples
        target_link_libraries(openCSMP-examples PRIVATE openCSMP)
    else()
        message(STATUS "No example sources found — skipping openCSMP-examples")
    endif()
endif()


# ===============================================================================================
# Test Executable: openCSMP-tests
#
#   How to run:
#     ctest                                      — run all tests
#     ctest --output-on-failure                  — show output for failing tests
#     ctest --rerun-failed --output-on-failure   — rerun and report failures
#     ctest -R Suite_Fundamentals                — run a specific suite
#     ctest --output-junit test-dashboard.xml    — JUnit XML report
#
#   Test suites run in strict dependency order:
#     Fundamentals -> Interdependent1 -> Interdependent2 -> Interfaces -> Composite
# ===============================================================================================
if(BUILD_TESTING)
    include(CTest)
    enable_testing()

    file(GLOB_RECURSE TEST_FILES
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/unit/*.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/unit/*.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/*.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/*.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/Test.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/TestSuite.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/TestSuite.h"
    )

    list(FILTER TEST_FILES EXCLUDE REGEX "/deprecated/")

    add_executable(openCSMP-tests
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/unit_CTests_main.cpp
            ${TEST_FILES}
    )

    target_include_directories(openCSMP-tests PRIVATE
            ${CSMP_CONSUMER_INCLUDE_DIRS}
    )

    # Links against the core library (inherits all its PUBLIC dependencies)
    target_link_libraries(openCSMP-tests PRIVATE openCSMP)


    # ===============================================================================================
    # Test Executable: CVFEM_tests_main
    # ===============================================================================================
    # CVFEM_HSS_LIB_SOURCES is every HSS .cpp EXCEPT any that defines main() (auto-detected
    # above) — i.e. all of modules/hydrothermal_system_simulator, recursively, minus a driver
    # main() if one exists.
    add_executable(CVFEM_tests_main
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/CVFEM_tests_main.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/Test.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/TestSuite.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/CVFEM/CVFEM_1D_VVCase.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/CVFEM/CVFEM_2D_VVCase.cpp
            ${CVFEM_HSS_LIB_SOURCES}
    )

    target_include_directories(CVFEM_tests_main PRIVATE
            ${CSMP_CONSUMER_INCLUDE_DIRS}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/CVFEM
    )

    target_link_libraries(CVFEM_tests_main PRIVATE openCSMP)

    target_compile_definitions(CVFEM_tests_main PRIVATE
        CVFEM_BENCHMARK_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data-tests_CVFEM"
    )

    # Give the tests access to the HSS include tree — always available (see "Source and
    # Header Discovery" above).
    if(CVFEM_HSS_SOURCES OR CVFEM_HSS_HEADERS)
        get_unique_directories(CVFEM_HSS_INCLUDE_DIRS "${CVFEM_HSS_SOURCES};${CVFEM_HSS_HEADERS}")

        target_include_directories(CVFEM_tests_main PRIVATE
                ${CVFEM_HSS_INCLUDE_DIRS}
                "${CVFEM_HYDROTHERMAL_ROOT}/includes/hydrothermal_system_simulator"
                "${CMAKE_CURRENT_SOURCE_DIR}"
        )
    endif()

    add_test(NAME CVFEM_tests_main COMMAND CVFEM_tests_main)


    # -----------------------------------------------------------------------
    # Test Runtime Data
    #   Copies the tests/dataset directory into the build tree so tests can
    #   find their input files regardless of where ctest is invoked from.
    # -----------------------------------------------------------------------
    set(TEST_WORKING_DIR "${CMAKE_BINARY_DIR}/tests_runtime")

    add_custom_target(copy_test_data ALL
            COMMAND ${CMAKE_COMMAND} -E make_directory "${TEST_WORKING_DIR}"
            COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/dataset"
            "${TEST_WORKING_DIR}"
            COMMENT "Copying test dataset from tests/dataset to ${TEST_WORKING_DIR}"
            VERBATIM
    )
    add_dependencies(openCSMP-tests copy_test_data)
    if(TARGET CVFEM_tests_main)
        add_dependencies(CVFEM_tests_main copy_test_data)
    endif()

    # Also copy at configure time for IDE runs
    file(GLOB TEST_INPUT_FILES "${CMAKE_CURRENT_SOURCE_DIR}/tests/dataset/*")
    if(TEST_INPUT_FILES)
        file(COPY ${TEST_INPUT_FILES} DESTINATION ${TEST_WORKING_DIR})
    endif()

    # -----------------------------------------------------------------------
    # CTest Registration
    # -----------------------------------------------------------------------
    add_test(NAME openCSMP-tests COMMAND openCSMP-tests)
    set_tests_properties(openCSMP-tests PROPERTIES
            WORKING_DIRECTORY ${TEST_WORKING_DIR}
    )

    # Individual suites (run via --suite flag in the test executable)
    add_test(NAME Suite_Fundamentals    COMMAND openCSMP-tests --suite fundamentals)
    add_test(NAME Suite_Interdependent1 COMMAND openCSMP-tests --suite interdependent1)
    add_test(NAME Suite_Interdependent2 COMMAND openCSMP-tests --suite interdependent2)
    add_test(NAME Suite_Interfaces      COMMAND openCSMP-tests --suite interfaces)
    add_test(NAME Suite_Composite       COMMAND openCSMP-tests --suite composite)

    set_tests_properties(
                openCSMP-tests
                Suite_Fundamentals
                Suite_Interdependent1
                Suite_Interdependent2
                Suite_Interfaces
                Suite_Composite
                PROPERTIES
                WORKING_DIRECTORY  "${TEST_WORKING_DIR}"
    ENVIRONMENT       "DYLD_LIBRARY_PATH=${CSMP_EXTERNAL_LIBS}/samg/macosx"
    "LD_LIBRARY_PATH=${CSMP_EXTERNAL_LIBS}/samg/linux"
    "CSMP_DATA_ROOT=${CMAKE_CURRENT_SOURCE_DIR}/tests/dataset"
        )

    # Strict sequential dependency chain
    set_tests_properties(Suite_Interdependent1 PROPERTIES DEPENDS Suite_Fundamentals)
    set_tests_properties(Suite_Interdependent2 PROPERTIES DEPENDS Suite_Interdependent1)
    set_tests_properties(Suite_Interfaces      PROPERTIES DEPENDS Suite_Interdependent2)
    set_tests_properties(Suite_Composite       PROPERTIES DEPENDS Suite_Interfaces)

endif()


# ===============================================================================================
# Installation
# ===============================================================================================

# Core library
install(TARGETS openCSMP
        EXPORT  openCSMP-targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
)

# Optional executables (only installed if they were built)
if(TARGET openCSMP-examples)
    install(TARGETS openCSMP-examples RUNTIME DESTINATION bin)
endif()
if(TARGET openCSMP-tests)
    install(TARGETS openCSMP-tests RUNTIME DESTINATION bin)
endif()

# Public headers
install(DIRECTORY source/includes/      DESTINATION include/openCSMP)
install(DIRECTORY source/non-GIT-files/ DESTINATION include/openCSMP/non-GIT-files)

# Example inputs (copied to build tree at build time, installed separately)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/example_inputs")
    if(TARGET openCSMP-examples)
        add_custom_command(TARGET openCSMP-examples POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_directory
                "${CMAKE_CURRENT_SOURCE_DIR}/examples/example_inputs"
                "${CMAKE_BINARY_DIR}/example_inputs"
                COMMENT "Copying example_inputs to build directory"
                VERBATIM
        )
    endif()
    install(DIRECTORY examples/example_inputs/
            DESTINATION example_inputs
            PATTERN ".DS_Store"  EXCLUDE
            PATTERN "__MACOSX"   EXCLUDE
            PATTERN "Icon"       EXCLUDE
    )
endif()

# Continuous integration fixtures
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../continuous_integration_testing")
    install(DIRECTORY continuous_integration_testing/
            DESTINATION continuous_integration_testing
    )
endif()


# ===============================================================================================
# Configuration Summary
# ===============================================================================================
message(STATUS "------------------------------------------------------------")
message(STATUS "openCSMP configuration summary")
message(STATUS "------------------------------------------------------------")
message(STATUS "  C++ standard:       ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build type:         ${CMAKE_BUILD_TYPE}")
message(STATUS "  SAMG solver:        ${SAMG_ENABLE}")
message(STATUS "  PETSc solver:       ${CSMP_WITH_PETSC_SOLVER}")
message(STATUS "  Image output:       ${CSMP_WITH_IMAGE_OUTPUT}")
message(STATUS "  Build examples:     ${BUILD_EXAMPLES}")
message(STATUS "  Build tests:        ${BUILD_TESTING}")
message(STATUS "  LTO/IPO:            ${IPO_SUPPORTED}")
message(STATUS "  Source files:       ${LIBRARY_SOURCES}")
message(STATUS "  Include dirs:       ${ALL_SUBDIRECTORIES}")
message(STATUS "------------------------------------------------------------")
