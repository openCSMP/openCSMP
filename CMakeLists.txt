cmake_minimum_required(VERSION 3.29)

# ===============================================================================================
#  openCSMP — Library, Examples and Tests
# ===============================================================================================
#
# SPDX-FileCopyrightText: © 2026 The openCSMP project
# SPDX-License-Identifier: LGPL-3.0-only
#
#  Prerequisites (see README.md for the install commands per platform):
#    - C++23 capable compiler (g++ 13.3 or later, Clang, Intel, NVC++)
#    - CMake 3.29 or later, git, pkg-config
#    - MPI; PETSc, libjpeg and CGNS for the default configuration
#    - Eigen, plf_colony and RapidXML are downloaded at configure time
#      (an installed Eigen is used if one is found)
#    - On Mac OSX with SAMG: you need Rosetta (Intel x86_64)
#                            Enter Rosetta shell first
#                            arch -x86_64 /bin/zsh
#  Usage:
#    cd openCSMP
#    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
#    cmake --build build --parallel
#
#    with Ninja:   cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
#    with SAMG:    add -DSAMG_ENABLE=ON -DSAMG_ROOT=/path/to/samg
#
#  Running the examples (outside of the repository):
#    cmake --build build --target setup_examples     # copies the inputs to CSMP_EXAMPLES_WORKSPACE
#    build/openCSMP-examples                         # reads from / writes to that workspace
#
#  Tests:
#    ctest --test-dir build --output-on-failure
# ===============================================================================================


# ===============================================================================================
# Build Options
# ===============================================================================================
option(BUILD_TESTING          "Enable building tests"                              ON )
option(BUILD_EXAMPLES         "Enable building examples"                           ON )
option(BUILD_SHARED_LIBS      "Build shared libraries (OFF = static)"              OFF)
option(CSMP_WITH_IMAGE_OUTPUT "Enable image output support via libjpeg"            ON )
option(SAMG_ENABLE            "Enable Fraunhofer SAMG algebraic multigrid solver"  OFF)
option(CSMP_WITH_PETSC_SOLVER "Enable PETSc solver"                                ON )
option(OPENNURBS_ENABLE       "Enable OpenNURBS NURBS geometry kernel"             OFF)
option(BUILD_CVFEM_TESTS      "Enable building CVFEM verification & validation tests" OFF)


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
# enable testing (using Google's CTest framework)
include(CTest)

# ===============================================================================================
# C++ Standard — set globally as project default; also applied per-target via
# target_compile_features() after each target is defined.
# ===============================================================================================
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
# No directory outside of the repository is required. Open-source
# dependencies are either found on the system or downloaded into the
# build tree; only the commercial SAMG binaries need a location, which
# is given with -DSAMG_ROOT=... (see the SAMG section below).
# -----------------------------------------------------------------------


# ===============================================================================================
# Dependencies (library discovery section) — Fetched at configure time
# ===============================================================================================

include(FetchContent)

# 1. plf_colony (header-only, essential)
# ------------------------------------------------------------------------------------------------
FetchContent_Declare(
        plf_colony
        GIT_REPOSITORY https://github.com/mattreecebentley/plf_colony.git
        # pinned to a known-good commit for reproducible builds (was: master)
        GIT_TAG        b027591ab6280a131450a0b33e92e42676fc7219
        EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(plf_colony)
message(STATUS "plf_colony source: ${plf_colony_SOURCE_DIR}")


# 2. Eigen (header-only, essential)
#    Linear algebra library. An installed Eigen is used if one is found:
#      macOS:  brew install eigen           (located via `brew --prefix eigen`)
#      Linux:  sudo apt install libeigen3-dev
#    Otherwise Eigen 3.4.0 is downloaded at configure time (all platforms).
#    Use -DCSMP_FETCH_EIGEN=ON to force the download.
#    Result: CSMP_EIGEN_INCLUDE_DIR, the directory that contains Eigen/Dense.
# ------------------------------------------------------------------------------------------------
option(CSMP_FETCH_EIGEN "Download Eigen even if an installed copy is available" OFF)

if(CSMP_FETCH_EIGEN)
    unset(CSMP_EIGEN_INCLUDE_DIR CACHE)
else()
    set(_csmp_eigen_hints "")
    if(APPLE)
        # Homebrew: /opt/homebrew (arm64) or /usr/local (x86_64 / Rosetta)
        execute_process(
                COMMAND brew --prefix eigen
                OUTPUT_VARIABLE HOMEBREW_EIGEN_PREFIX
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
        )
        if(HOMEBREW_EIGEN_PREFIX)
            list(APPEND _csmp_eigen_hints "${HOMEBREW_EIGEN_PREFIX}/include")
        endif()
    endif()
    find_path(CSMP_EIGEN_INCLUDE_DIR
            NAMES Eigen/Dense
            HINTS ${_csmp_eigen_hints}
            PATH_SUFFIXES eigen3
            DOC   "Directory that contains Eigen/Dense"
    )
    unset(_csmp_eigen_hints)
endif()

if(CSMP_EIGEN_INCLUDE_DIR)
    message(STATUS "Eigen found (installed): ${CSMP_EIGEN_INCLUDE_DIR}")
else()
    message(STATUS "Eigen not found on this system — downloading Eigen 3.4.0")

    set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_DOC     OFF CACHE BOOL "" FORCE)

    # Temporarily override BUILD_TESTING so Eigen's CMakeLists.txt cannot
    # register its tests even if it checks BUILD_TESTING directly.
    set(_saved_build_testing ${BUILD_TESTING})
    set(BUILD_TESTING OFF)

    FetchContent_Declare(
            eigen
            GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
            GIT_TAG        3.4.0
            GIT_SHALLOW    TRUE
            EXCLUDE_FROM_ALL
    )
    FetchContent_MakeAvailable(eigen)

    # Restore BUILD_TESTING for openCSMP's own tests
    set(BUILD_TESTING ${_saved_build_testing})
    unset(_saved_build_testing)

    # not cached: the download is re-evaluated (cheaply) at every configure
    set(CSMP_EIGEN_INCLUDE_DIR "${eigen_SOURCE_DIR}")
    message(STATUS "Eigen downloaded: ${CSMP_EIGEN_INCLUDE_DIR}")
endif()

if(NOT EXISTS "${CSMP_EIGEN_INCLUDE_DIR}/Eigen/Dense")
    message(FATAL_ERROR
            "Eigen headers not found at ${CSMP_EIGEN_INCLUDE_DIR}\n"
            "Install Eigen (macOS: brew install eigen, Linux: sudo apt install libeigen3-dev)\n"
            "or configure with -DCSMP_FETCH_EIGEN=ON")
endif() # Eigen



# ===============================================================================================
# Dependencies — System Libraries (found via pkg-config / find_package)
# ===============================================================================================

# Find pkg-config once — used by PETSc and CGNS
# Install with: sudo apt install pkg-config
find_package(PkgConfig REQUIRED)


# 3. MPI (essential)
# ------------------------------------------------------------------------------------------------
if(APPLE)
    if(SAMG_ENABLE AND CMAKE_OSX_ARCHITECTURES STREQUAL "x86_64")
        # SAMG requires x86_64 — use Rosetta Homebrew MPI at /usr/local/
        set(MPI_CXX_COMPILER     "/usr/local/opt/open-mpi/bin/mpicxx" )
        set(MPI_C_COMPILER       "/usr/local/opt/open-mpi/bin/mpicc"  )
        set(MPI_Fortran_COMPILER "/usr/local/opt/open-mpi/bin/mpifort")
        message(STATUS "Using x86_64 MPI for SAMG/Rosetta build")
    else()
        # Native arm64 build — use arm64 Homebrew MPI at /opt/homebrew/
        execute_process(
                COMMAND brew --prefix open-mpi
                OUTPUT_VARIABLE HOMEBREW_MPI_PREFIX
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
        )
        if(HOMEBREW_MPI_PREFIX AND EXISTS "${HOMEBREW_MPI_PREFIX}/bin/mpicxx")
            set(MPI_CXX_COMPILER     "${HOMEBREW_MPI_PREFIX}/bin/mpicxx" )
            set(MPI_C_COMPILER       "${HOMEBREW_MPI_PREFIX}/bin/mpicc"  )
            set(MPI_Fortran_COMPILER "${HOMEBREW_MPI_PREFIX}/bin/mpifort")
            message(STATUS "Using arm64 MPI: ${HOMEBREW_MPI_PREFIX}")
        else()
            # no Homebrew Open MPI: let find_package(MPI) search the standard locations
            message(STATUS "Homebrew Open MPI not found (brew install open-mpi) — searching for another MPI")
        endif()
    endif()
endif()

find_package(MPI QUIET)
if(NOT MPI_CXX_FOUND)
    message(FATAL_ERROR
            "MPI (C++) not found.\n"
            "  macOS:  brew install open-mpi\n"
            "  Linux:  sudo apt install libopenmpi-dev openmpi-bin")
endif()
message(STATUS "MPI found: ${MPI_CXX_LIBRARIES}")


# 4. PETSc (optional, default ON)
#    Linux:  sudo apt install petsc-dev
#    macOS:  brew install petsc
# ------------------------------------------------------------------------------------------------
if(CSMP_WITH_PETSC_SOLVER)
    if(APPLE)
        if(SAMG_ENABLE AND CMAKE_OSX_ARCHITECTURES STREQUAL "x86_64")
            # x86_64 Rosetta — Homebrew installed at /usr/local
            set(ENV{PKG_CONFIG_PATH} "/usr/local/opt/petsc/lib/pkgconfig:/usr/local/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
            message(STATUS "Searching for PETSc (x86_64 / Rosetta)...")
        else()
            # Native arm64 — Homebrew installed at /opt/homebrew
            execute_process(
                    COMMAND brew --prefix petsc
                    OUTPUT_VARIABLE HOMEBREW_PETSC_PREFIX
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET
            )
            if(HOMEBREW_PETSC_PREFIX)
                set(ENV{PKG_CONFIG_PATH} "${HOMEBREW_PETSC_PREFIX}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
            else()
                set(ENV{PKG_CONFIG_PATH} "/opt/homebrew/opt/petsc/lib/pkgconfig:/opt/homebrew/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
            endif()
            message(STATUS "Searching for PETSc (arm64): ${HOMEBREW_PETSC_PREFIX}")
        endif()
    endif()

    # IMPORTED_TARGET generates PkgConfig::PETSC target with correct -I, -L, and -l flags attached
    pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)
    message(STATUS "PETSc include dirs: ${PETSC_INCLUDE_DIRS}")
    message(STATUS "PETSc libraries:    ${PETSC_LIBRARIES}")
endif()


# 5. libjpeg (optional)
#    Linux:  sudo apt install libjpeg-dev
#    macOS:  brew install jpeg
# ------------------------------------------------------------------------------------------------
if(CSMP_WITH_IMAGE_OUTPUT)
    if(APPLE AND CMAKE_OSX_ARCHITECTURES STREQUAL "x86_64")
        # x86_64 Rosetta — jpeg installed via /usr/local/ Homebrew
        set(JPEG_INCLUDE_DIRS "/usr/local/opt/jpeg/include")
        set(JPEG_LIBRARIES    "/usr/local/opt/jpeg/lib/libjpeg.dylib")
        if(NOT EXISTS "${JPEG_INCLUDE_DIRS}/jpeglib.h")
            message(FATAL_ERROR
                    "JPEG headers not found at ${JPEG_INCLUDE_DIRS}\n"
                    "Install with: brew install jpeg  (in the Rosetta shell)")
        endif()
        message(STATUS "JPEG x86_64: ${JPEG_INCLUDE_DIRS}")

    elseif(APPLE)
        # arm64 native build
        execute_process(
                COMMAND brew --prefix jpeg
                OUTPUT_VARIABLE HOMEBREW_JPEG_PREFIX
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
        )
        if(HOMEBREW_JPEG_PREFIX AND EXISTS "${HOMEBREW_JPEG_PREFIX}/include/jpeglib.h")
            set(JPEG_INCLUDE_DIRS "${HOMEBREW_JPEG_PREFIX}/include")
            set(JPEG_LIBRARIES    "${HOMEBREW_JPEG_PREFIX}/lib/libjpeg.dylib")
        else()
            # e.g. jpeg-turbo, MacPorts or conda: let CMake search
            find_package(JPEG QUIET)
            if(NOT JPEG_FOUND)
                message(FATAL_ERROR
                        "libjpeg not found.\n"
                        "Install with: brew install jpeg\n"
                        "or configure with -DCSMP_WITH_IMAGE_OUTPUT=OFF")
            endif()
        endif()
        message(STATUS "JPEG arm64: ${JPEG_INCLUDE_DIRS}")

    else()
        # Linux
        find_package(JPEG REQUIRED)
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


# 7. CGNS — CFD General Notation System
#    Linux:  sudo apt install libcgns-dev
#    macOS:  brew install cgns
#    Fallback: header bundled with openCSMP (no linking)
#
#   Discovery order:
#     7.1 find_package  (CMake config files)
#     7.2 pkg-config
#     7.3 Homebrew manual search (macOS)
#     7.4 Direct system search   (covers apt-installed libcgns-dev)
#     7.5 Fallback header-only   (source/includes/interfaces/CGNS)
# -----------------------------------------------------------------------
set(CGNS_FOUND_USABLE FALSE)

# 7.1 Try find_package first.
#     The CGNS config package defines only imported targets (CGNS::cgns_shared,
#     CGNS::cgns_static) and no CGNS_LIBRARIES / CGNS_INCLUDE_DIRS. Its targets
#     also name an "hdf5-shared" dependency that the package does not load, so
#     linking the target directly can fail. Take the shared library path and
#     include directory from the target instead; libcgns itself records its
#     HDF5 dependency. Without a usable shared target, fall through to 7.2+.
find_package(CGNS QUIET)
if(CGNS_FOUND AND TARGET CGNS::cgns_shared)
    get_target_property(_cgns_configs CGNS::cgns_shared IMPORTED_CONFIGURATIONS)
    set(_cgns_lib "")
    foreach(_cfg IN LISTS _cgns_configs)
        get_target_property(_cgns_lib CGNS::cgns_shared IMPORTED_LOCATION_${_cfg})
        if(_cgns_lib)
            break()
        endif()
    endforeach()
    if(NOT _cgns_lib)
        get_target_property(_cgns_lib CGNS::cgns_shared IMPORTED_LOCATION)
    endif()
    get_target_property(_cgns_inc CGNS::cgns_shared INTERFACE_INCLUDE_DIRECTORIES)
    if(_cgns_lib AND _cgns_inc)
        list(REMOVE_DUPLICATES _cgns_inc)
        set(CGNS_LIBRARIES    "${_cgns_lib}")
        set(CGNS_INCLUDE_DIRS "${_cgns_inc}")
        message(STATUS "CGNS found via find_package: ${CGNS_LIBRARIES}")
        set(CGNS_FOUND_USABLE TRUE)
    endif()
endif()

# 7.2 Try pkg-config
if(NOT CGNS_FOUND_USABLE AND PKG_CONFIG_FOUND)
    pkg_check_modules(CGNS QUIET cgns)
    if(CGNS_FOUND)
        message(STATUS "CGNS found via pkg-config: ${CGNS_LIBRARIES}")
        set(CGNS_FOUND_USABLE TRUE)
    endif()
endif()

# 7.3 Try Homebrew manually on macOS
if(APPLE AND NOT CGNS_FOUND_USABLE)
    execute_process(
            COMMAND brew --prefix cgns
            OUTPUT_VARIABLE HOMEBREW_CGNS_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
    )
    if(HOMEBREW_CGNS_PREFIX AND EXISTS "${HOMEBREW_CGNS_PREFIX}/include/cgnslib.h")
        message(STATUS "CGNS found via Homebrew: ${HOMEBREW_CGNS_PREFIX}")
        set(CGNS_INCLUDE_DIRS "${HOMEBREW_CGNS_PREFIX}/include")
        set(CGNS_LIBRARIES    "${HOMEBREW_CGNS_PREFIX}/lib/libcgns.dylib")
        set(CGNS_FOUND_USABLE TRUE)
    endif()
endif()

# 7.4 Direct system search (covers apt-installed libcgns-dev on Ubuntu/Debian)
if(NOT CGNS_FOUND_USABLE)
    find_path(CGNS_SYSTEM_INCLUDE_DIR NAMES cgnslib.h)
    find_library(CGNS_SYSTEM_LIBRARY  NAMES cgns)
    if(CGNS_SYSTEM_INCLUDE_DIR AND CGNS_SYSTEM_LIBRARY)
        message(STATUS "CGNS found via direct system search: ${CGNS_SYSTEM_LIBRARY}")
        set(CGNS_INCLUDE_DIRS "${CGNS_SYSTEM_INCLUDE_DIR}")
        set(CGNS_LIBRARIES    "${CGNS_SYSTEM_LIBRARY}")
        set(CGNS_FOUND_USABLE TRUE)
    endif()
endif()

# 7.5 Fallback — use the cgnslib.h bundled with openCSMP (no linking)
if(NOT CGNS_FOUND_USABLE)
    message(WARNING
            "CGNS library not found — using the bundled header only.\n"
            "openCSMP itself will build, but executables that use the CGNS interface\n"
            "will fail to link with undefined cg_* symbols. To avoid this install CGNS:\n"
            "  macOS:  brew install cgns\n"
            "  Linux:  sudo apt install libcgns-dev")
    set(CGNS_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/interfaces/CGNS")
    set(CGNS_LIBRARIES    "")
endif()

# 7.6 Interface target for consistent usage across all consumers
add_library(cgns_iface INTERFACE)
target_include_directories(cgns_iface INTERFACE ${CGNS_INCLUDE_DIRS})
if(CGNS_LIBRARIES)
    target_link_libraries(cgns_iface INTERFACE ${CGNS_LIBRARIES})
endif()


# 8. RapidXML (header-only, essential)
#    No apt package available — fetched automatically.
# ------------------------------------------------------------------------------------------------
# 8. rapidxml
FetchContent_Declare(
        rapidxml
        GIT_REPOSITORY https://github.com/discord/rapidxml.git
        # pinned to a known-good commit for reproducible builds (was: master)
        GIT_TAG        2ae4b2888165a393dfb6382168825fddf00c27b9
        EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(rapidxml)
message(STATUS ">>> rapidxml_SOURCE_DIR = ${rapidxml_SOURCE_DIR}")
message(STATUS ">>> rapidxml.hpp exists = ${rapidxml_SOURCE_DIR}/rapidxml.hpp")


# 9. OpenNURBS (McNeel & Associates — NURBS geometry kernel)
#    Cloned automatically into <build>/_deps/opennurbs-src if not present.
#    A CMakeLists.txt wrapper is generated automatically.
#
#    System prerequisites:
#      Linux:  sudo apt install zlib1g-dev
#      macOS:  brew install zlib  (usually already present via Xcode)
# ------------------------------------------------------------------------------------------------
#    Use -DOPENNURBS_ROOT=/path/to/opennurbs to use an existing checkout instead.
set(OPENNURBS_ROOT "" CACHE PATH "Existing OpenNURBS checkout (empty: clone into the build tree)")
set(OPENNURBS_GIT_BRANCH "8.x" CACHE STRING "Branch of https://github.com/mcneel/opennurbs to clone")

get_filename_component(_csmp_former_opennurbs_dir
        "${CMAKE_CURRENT_SOURCE_DIR}/../external_support_libs/opennurbs" ABSOLUTE)
if(OPENNURBS_ROOT)
    set(OPENNURBS_SOURCE_DIR "${OPENNURBS_ROOT}")
elseif(OPENNURBS_ENABLE AND EXISTS "${_csmp_former_opennurbs_dir}/.git")
    # backward compatibility: checkout next to the repository, as used previously
    set(OPENNURBS_SOURCE_DIR "${_csmp_former_opennurbs_dir}")
    message(STATUS "OPENNURBS_ROOT not set — using the former default location: ${OPENNURBS_SOURCE_DIR}")
else()
    set(OPENNURBS_SOURCE_DIR "${CMAKE_BINARY_DIR}/_deps/opennurbs-src")
endif()
unset(_csmp_former_opennurbs_dir)
set(OPENNURBS_AVAILABLE  FALSE)

if(OPENNURBS_ENABLE)
    # Step 1 — Make sure the parent directory in the build tree exists
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/_deps")

    # Step 2 — Clone OpenNURBS if not already present
    if(OPENNURBS_ROOT AND NOT EXISTS "${OPENNURBS_SOURCE_DIR}")
        message(FATAL_ERROR "OPENNURBS_ROOT does not exist: ${OPENNURBS_ROOT}")
    endif()
    if(NOT OPENNURBS_ROOT AND NOT EXISTS "${OPENNURBS_SOURCE_DIR}/.git")
        message(STATUS "Cloning OpenNURBS into ${OPENNURBS_SOURCE_DIR} ...")
        execute_process(
                COMMAND git clone
                --depth 1
                --branch ${OPENNURBS_GIT_BRANCH}
                https://github.com/mcneel/opennurbs.git
                "${OPENNURBS_SOURCE_DIR}"
                RESULT_VARIABLE _opennurbs_clone_result
                OUTPUT_VARIABLE _opennurbs_clone_output
                ERROR_VARIABLE  _opennurbs_clone_error
        )
        if(NOT _opennurbs_clone_result EQUAL 0)
            message(FATAL_ERROR
                    "Failed to clone OpenNURBS:\n${_opennurbs_clone_error}\n"
                    "Check your network connection and that git is installed.")
        endif()
        message(STATUS "OpenNURBS cloned successfully")
    else()
        message(STATUS "OpenNURBS already cloned: ${OPENNURBS_SOURCE_DIR}")
    endif()

    # Step 3 — Generate a CMakeLists.txt wrapper if not already present
    #          (never overwrites an existing one so you can customise it)
    if(NOT EXISTS "${OPENNURBS_SOURCE_DIR}/CMakeLists.txt")
        message(STATUS "Generating OpenNURBS CMakeLists.txt wrapper...")
        file(WRITE "${OPENNURBS_SOURCE_DIR}/CMakeLists.txt" [=[
cmake_minimum_required(VERSION 3.29)
project(opennurbs CXX)

file(GLOB_RECURSE OPENNURBS_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

# Exclude example and test files shipped with the SDK
list(FILTER OPENNURBS_SRC EXCLUDE REGEX "[Ee]xample")
list(FILTER OPENNURBS_SRC EXCLUDE REGEX "[Tt]est")

add_library(opennurbs STATIC ${OPENNURBS_SRC})

target_include_directories(opennurbs PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/opennurbs>
)

target_compile_features(opennurbs PUBLIC cxx_std_23)
set_target_properties(opennurbs PROPERTIES
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS        OFF
)

    if(OPENNURBS_AVAILABLE AND APPLE)
        target_link_libraries(opennurbs PUBLIC
                "-framework CoreFoundation"
                "-framework Foundation"
        )
    endif()

find_package(ZLIB REQUIRED)
target_link_libraries(opennurbs PUBLIC ZLIB::ZLIB)
]=])
        message(STATUS "OpenNURBS CMakeLists.txt wrapper generated")
    endif()

    # Step 4 — Add as a subdirectory now that the wrapper exists
    add_subdirectory(
            "${OPENNURBS_SOURCE_DIR}"
            "${CMAKE_BINARY_DIR}/opennurbs"
    )

	if(NOT TARGET opennurbs)
		message(FATAL_ERROR
				"OpenNURBS target was not created by add_subdirectory.\n"
				"Check the generated CMakeLists.txt at: ${OPENNURBS_SOURCE_DIR}/CMakeLists.txt")
	endif()

    # macOS framework dependencies — inside the guard where the target exists
    if(APPLE)
        target_link_libraries(opennurbs PUBLIC
                "-framework CoreFoundation"
                "-framework Foundation"
        )
    endif()

    set(OPENNURBS_AVAILABLE TRUE)
    message(STATUS "OpenNURBS available: ${OPENNURBS_SOURCE_DIR}")

else()
    message(STATUS "OpenNURBS disabled")
endif()


# ===============================================================================================
# Source and Header Discovery
# ===============================================================================================

# Collect all .cpp sources for the core library
file(GLOB_RECURSE LIBRARY_SOURCES CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/source/modules/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/non-GIT-files/*.cpp"
)

# Collect all headers (used to derive include directories)
file(GLOB_RECURSE LIBRARY_HEADERS CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/source/modules/*.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/modules/*.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/*.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/*.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/non-GIT-files/*.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/non-GIT-files/*.hpp"
)

# Exclude deprecated files from both lists
list(FILTER LIBRARY_SOURCES EXCLUDE REGEX "/deprecated/")
list(FILTER LIBRARY_HEADERS EXCLUDE REGEX "/deprecated/")

# Exclude SAMG sources
# (SAMG sources are added back conditionally in the SAMG section below)
list(FILTER LIBRARY_SOURCES EXCLUDE REGEX ".*[Ss][Aa][Mm][Gg].*")

# Hydrothermal_system_modelling (includes + modules) — ALWAYS discovered.
# What gets *compiled into a target* is controlled entirely by each consuming
# target's own option (BUILD_EXAMPLES for openCSMP-examples, BUILD_TESTING for
# CVFEM_tests_main).
#
#   cvfe_hydrothermal             — supporting library code, folded directly into
#                                   openCSMP itself (same treatment as the main
#                                   source/ tree, inheriting the project's
#                                   flat-include model).
#   hydrothermal_system_simulator — supporting code for examples, tests and future
#                                   simulators. NOT folded into openCSMP; consumed
#                                   directly by openCSMP-examples and CVFEM_tests_main.
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

# Hydrothermal-System-Simulator (HSS) sources/headers — ALWAYS discovered.
# NOT merged into LIBRARY_SOURCES/LIBRARY_HEADERS (not folded into the openCSMP
# library itself) — consumed directly by whichever targets need them
# (openCSMP-examples, CVFEM_tests_main), each gated by its own option.
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

# Split HSS sources into driver (contains main()) and library sources.
# When a main() is added back it will automatically land in CVFEM_HSS_MAIN_SOURCES.
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
            cmake_path(GET file_path PARENT_PATH dir_path)
            list(APPEND dir_list ${dir_path})
        endif()
    endforeach()
    list(REMOVE_DUPLICATES dir_list)
    set(${result} ${dir_list} PARENT_SCOPE)
endfunction()

get_unique_directories(ALL_SUBDIRECTORIES "${LIBRARY_SOURCES};${LIBRARY_HEADERS}")

if(CVFEM_HSS_SOURCES OR CVFEM_HSS_HEADERS)
    get_unique_directories(CVFEM_HSS_INCLUDE_DIRS "${CVFEM_HSS_SOURCES};${CVFEM_HSS_HEADERS}")
    message(STATUS "CVFEM HSS include dirs: ${CVFEM_HSS_INCLUDE_DIRS}")
endif()

# Always include top-level roots for flat includes
list(APPEND ALL_SUBDIRECTORIES
        ${CMAKE_CURRENT_SOURCE_DIR}/source/includes
        ${CMAKE_CURRENT_SOURCE_DIR}/source/modules
        ${CMAKE_CURRENT_SOURCE_DIR}/source/non-GIT-files
        ${CMAKE_CURRENT_SOURCE_DIR}/non-GIT-files
)

# FetchContent directories are not on disk at GLOB time so they are
# not picked up by get_unique_directories() — append them explicitly.
list(APPEND ALL_SUBDIRECTORIES
        ${plf_colony_SOURCE_DIR}
        ${rapidxml_SOURCE_DIR}
)

list(REMOVE_DUPLICATES ALL_SUBDIRECTORIES)



# ===============================================================================================
# Core Library Target: openCSMP
# ===============================================================================================
add_library(openCSMP ${LIBRARY_SOURCES})

# C++ standard applied per-target (global defaults set above are the fallback)
target_compile_features(openCSMP PUBLIC cxx_std_23)
set_target_properties(openCSMP PROPERTIES
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS        OFF
)

# Public include directories — propagated to all targets linking against openCSMP
target_include_directories(openCSMP PUBLIC
        $<BUILD_INTERFACE:${ALL_SUBDIRECTORIES}>
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
target_include_directories(openCSMP PUBLIC
        ${plf_colony_SOURCE_DIR}
        ${rapidxml_SOURCE_DIR}
)

# Eigen — header-only; installed or downloaded, see dependency section 2
target_include_directories(openCSMP PUBLIC
        $<BUILD_INTERFACE:${CSMP_EIGEN_INCLUDE_DIR}>
)
message(STATUS "Eigen include path: ${CSMP_EIGEN_INCLUDE_DIR}")

# IAPS H2O EOS (optional)
if(TARGET iaps_h2o_eos)
    target_link_libraries(openCSMP PUBLIC iaps_h2o_eos)
endif()

# CGNS — cgns_iface carries include dirs and link libraries via INTERFACE
target_link_libraries(openCSMP PUBLIC cgns_iface)

# RapidXML — header-only
target_include_directories(openCSMP PUBLIC ${rapidxml_SOURCE_DIR})

# OpenNURBS (optional)
if(OPENNURBS_AVAILABLE)
    target_link_libraries(openCSMP PUBLIC opennurbs)
    target_include_directories(openCSMP PUBLIC
            $<BUILD_INTERFACE:${OPENNURBS_SOURCE_DIR}>
    )
endif()



# PETSc (optional)
if(CSMP_WITH_PETSC_SOLVER)
    target_link_libraries(openCSMP PUBLIC PkgConfig::PETSC)
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
#   Requires the licensed, prebuilt binaries. Tell CMake where they are with
#     -DSAMG_ROOT=/path/to/samg
#   where that directory contains linux/ and/or macosx/ with the libraries
#   (and the SAMG headers). For backward compatibility the former location
#   ../external_support_libs/samg is still used if SAMG_ROOT is not given.
#
#   Two modes:
#     SAMG_LEGACY=OFF (default) — modern SAMG integration, single instance
#     SAMG_LEGACY=ON            — original SAMG integration with SAMG_MULTIPLE_INSTANCES
#
#   Toggle with:
#     cmake -DSAMG_ENABLE=ON -DSAMG_ROOT=/path/to/samg ..
#     cmake -DSAMG_ENABLE=ON -DSAMG_ROOT=/path/to/samg -DSAMG_LEGACY=ON ..
# ===============================================================================================
option(SAMG_LEGACY "Use legacy SAMG integration (SAMG_MULTIPLE_INSTANCES)" OFF)
set(SAMG_ROOT "" CACHE PATH "Directory with the licensed SAMG binaries (contains linux/ and/or macosx/)")

if(SAMG_ENABLE)
    message(STATUS "Enabling SAMG Solver...")

    # Where are the SAMG binaries?
    if(NOT SAMG_ROOT)
        get_filename_component(_csmp_former_samg_dir
                "${CMAKE_CURRENT_SOURCE_DIR}/../external_support_libs/samg" ABSOLUTE)
        if(EXISTS "${_csmp_former_samg_dir}")
            set(SAMG_ROOT "${_csmp_former_samg_dir}")
            message(STATUS "SAMG_ROOT not set — using the former default location: ${SAMG_ROOT}")
        endif()
        unset(_csmp_former_samg_dir)
    endif()
    message(STATUS "SAMG root: ${SAMG_ROOT}")

    # Locate prebuilt libraries — no point proceeding without them
    set(SAMG_INTERFACE_REQUIRED FALSE)
    if(APPLE)
        set(SAMG_LIB_PATH "${SAMG_ROOT}/macosx")
        set(SAMG_LIB      "${SAMG_LIB_PATH}/libamg_mult.dylib")
        set(OMP_LIB       "${SAMG_LIB_PATH}/libiomp5.dylib")
        # these are fixed paths, not search results: make sure the files are there
        if(NOT EXISTS "${SAMG_LIB}" OR NOT EXISTS "${OMP_LIB}")
            set(SAMG_LIB "")
            set(OMP_LIB  "")
        endif()
    elseif(SAMG_LEGACY)
        find_library(SAMG_LIB NAMES amg_mult
                HINTS "${SAMG_ROOT}/linux")
        find_library(OMP_LIB NAMES iomp5
                HINTS "${SAMG_ROOT}/linux")
    else()
        find_library(SAMG_LIB NAMES samg
                HINTS "${SAMG_ROOT}/linux")
        find_library(SAMG_INTERFACE_LIB NAMES samg_interface
                HINTS "${SAMG_ROOT}/linux")
        find_library(OMP_LIB NAMES iomp5
                HINTS "${SAMG_ROOT}/linux")
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

        if(SAMG_LEGACY)
            target_compile_definitions(openCSMP PUBLIC
                    SAMG_MULTIPLE_INSTANCES
                    SAMG_LEGACY
            )
            message(STATUS "SAMG mode: LEGACY (SAMG_MULTIPLE_INSTANCES enabled)")
        else()
            message(STATUS "SAMG mode: MODERN (single instance)")
        endif()

        if(SAMG_SRC)
            target_sources(openCSMP PRIVATE ${SAMG_SRC})
            target_include_directories(openCSMP PUBLIC
                    ${SAMG_ROOT}
                    ${CMAKE_CURRENT_SOURCE_DIR}/source/modules/integration/solver/SAMG_profiles
            )
        else()
            message(WARNING "SAMG enabled but no SAMG source files found")
        endif()

        # SAMG_INTERFACE_LIB is empty/unset on APPLE — CMake ignores the empty expansion
        target_link_libraries(openCSMP PUBLIC ${SAMG_LIB} ${SAMG_INTERFACE_LIB} ${OMP_LIB})
        message(STATUS "SAMG library:           ${SAMG_LIB}")
        message(STATUS "SAMG interface library: ${SAMG_INTERFACE_LIB}")
        message(STATUS "OpenMP library:         ${OMP_LIB}")

    else()
        message(WARNING "SAMG requested but its libraries were not found under SAMG_ROOT "
                "('${SAMG_ROOT}') — disabling SAMG support.\n"
                "Configure with -DSAMG_ROOT=/path/to/samg (the directory that contains linux/ or macosx/)")
        set(SAMG_ENABLE OFF)
    endif()
endif()

if(APPLE AND SAMG_ENABLE)
    if(NOT CMAKE_OSX_ARCHITECTURES STREQUAL "x86_64")
        message(FATAL_ERROR
                "SAMG requires an x86_64 build on macOS. "
                "Configure with -DCMAKE_OSX_ARCHITECTURES=x86_64")
    endif()
endif()

# ===============================================================================================
# Link-Time Optimisation (LTO / IPO)
# ===============================================================================================
option(CSMP_WITH_LTO "Enable link-time optimisation" OFF)

if(CSMP_WITH_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT IPO_SUPPORTED)

    if(IPO_SUPPORTED)
        set_target_properties(openCSMP PROPERTIES
                INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
endif()


# ===============================================================================================
# Shared Include Directories for Examples and Tests
#
#   Both the examples and tests executables need the same broad set of include paths.
#   Defined once here as a list to avoid duplication.
#   Note: targets that link openCSMP PRIVATE already inherit its PUBLIC include dirs
#   transitively. This list covers direct #include references in example/test sources
#   that are not routed through the library headers.
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
        ${rapidxml_SOURCE_DIR}
)

if(OPENNURBS_AVAILABLE)
    list(APPEND CSMP_CONSUMER_INCLUDE_DIRS ${OPENNURBS_SOURCE_DIR})
endif()

# Append optional include dirs only if the feature is enabled
if(CSMP_WITH_IMAGE_OUTPUT)
    list(APPEND CSMP_CONSUMER_INCLUDE_DIRS ${JPEG_INCLUDE_DIRS})
endif()

# Eigen include path for consumer targets
list(APPEND CSMP_CONSUMER_INCLUDE_DIRS ${CSMP_EIGEN_INCLUDE_DIR})


# ===============================================================================================
# Examples Workspace
#
#   The examples are run in a workspace directory OUTSIDE of the repository and of the build
#   tree, so that simulation output is neither committed by accident nor lost when the build
#   directory is deleted:
#
#     <workspace>/example_inputs/                  inputs, copied from examples/example_inputs
#     <workspace>/example_outputs/<example name>/  results, e.g. for visualisation with ParaView
#
#   Create / update it with (never overwrites files that already exist there):
#     cmake --build <build dir> --target setup_examples
#
#   Default location: ~/openCSMP/examples; change with -DCSMP_EXAMPLES_WORKSPACE=/some/dir
#   At run time the executable uses, in this order: the option --workspace <dir>, the
#   environment variable CSMP_EXAMPLES_DIR, the current directory if it is an
#   example_inputs/ directory, and finally the location configured here.
# ===============================================================================================
if(DEFINED ENV{HOME})
    file(TO_CMAKE_PATH "$ENV{HOME}" _csmp_home)
elseif(DEFINED ENV{USERPROFILE})
    file(TO_CMAKE_PATH "$ENV{USERPROFILE}" _csmp_home)
else()
    set(_csmp_home "${CMAKE_BINARY_DIR}")
endif()
set(CSMP_EXAMPLES_WORKSPACE "${_csmp_home}/openCSMP/examples"
        CACHE PATH "Directory outside of the repository in which the examples are run")
unset(_csmp_home)
message(STATUS "Examples workspace: ${CSMP_EXAMPLES_WORKSPACE}")

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/example_inputs")
    add_custom_target(setup_examples
            COMMAND ${CMAKE_COMMAND}
            "-DCSMP_EXAMPLE_INPUTS=${CMAKE_CURRENT_SOURCE_DIR}/examples/example_inputs"
            "-DCSMP_EXAMPLES_WORKSPACE=${CSMP_EXAMPLES_WORKSPACE}"
            "-DCSMP_REPOSITORY=${CMAKE_CURRENT_SOURCE_DIR}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/CsmpSetupExamples.cmake"
            COMMENT "Setting up the examples workspace in ${CSMP_EXAMPLES_WORKSPACE}"
            VERBATIM
            USES_TERMINAL
    )
endif()


# ===============================================================================================
# Examples Executable: openCSMP-examples
# ===============================================================================================
if(BUILD_EXAMPLES)
    file(GLOB_RECURSE EXAMPLES_SRC
            "${CMAKE_CURRENT_SOURCE_DIR}/examples/source/*.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/examples/source/*.h"
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

        target_compile_features(openCSMP-examples PRIVATE cxx_std_23)
        set_target_properties(openCSMP-examples PROPERTIES
                CXX_STANDARD_REQUIRED ON
                CXX_EXTENSIONS        OFF
        )

        target_include_directories(openCSMP-examples PRIVATE
                ${CSMP_CONSUMER_INCLUDE_DIRS}
                "${CMAKE_CURRENT_SOURCE_DIR}"
        )

        if(CVFEM_HSS_SOURCES OR CVFEM_HSS_HEADERS)
            target_include_directories(openCSMP-examples PRIVATE
                    ${CVFEM_HSS_INCLUDE_DIRS}
                    "${CVFEM_HYDROTHERMAL_ROOT}/includes/hydrothermal_system_simulator"
            )
        endif()

        # Links against the core library (inherits all its PUBLIC dependencies)
        target_link_libraries(openCSMP-examples PRIVATE openCSMP)

        # Default examples workspace, see "Examples Workspace" below
        target_compile_definitions(openCSMP-examples PRIVATE
                CSMP_EXAMPLES_WORKSPACE_DEFAULT="${CSMP_EXAMPLES_WORKSPACE}"
        )
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
#     Suite_Refactoring -> Suite_Fundamentals -> Suite_Interdependent1
#       -> Suite_Interdependent2 -> Suite_Interfaces -> Suite_Composite
# ===============================================================================================
if(BUILD_TESTING)

# Individual test infrastructure files (not covered by the recursive globs)
set(TEST_INFRASTRUCTURE_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/Test.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/Test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/TestSuite.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/TestSuite.cpp"
)

# Recursive glob for unit and integration test sources
file(GLOB_RECURSE TEST_UNIT_FILES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/unit/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/unit/*.h"
)
file(GLOB_RECURSE TEST_INTEGRATION_FILES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/*.h"
)

set(TEST_FILES
    ${TEST_INFRASTRUCTURE_FILES}
    ${TEST_UNIT_FILES}
    ${TEST_INTEGRATION_FILES}
)

list(FILTER TEST_FILES EXCLUDE REGEX "/deprecated/")

# Exclude the main entry point from the glob — it is listed explicitly below
list(FILTER TEST_FILES EXCLUDE REGEX "unit_CTests_main\\.cpp")

add_executable(openCSMP-tests
"${CMAKE_CURRENT_SOURCE_DIR}/tests/unit_CTests_main.cpp"
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/Test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/TestSuite.cpp
        ${TEST_FILES}
)

target_compile_features(openCSMP-tests PRIVATE cxx_std_23)
set_target_properties(openCSMP-tests PROPERTIES
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS        OFF
)

target_include_directories(openCSMP-tests PRIVATE
    ${CSMP_CONSUMER_INCLUDE_DIRS}
    "${CMAKE_CURRENT_SOURCE_DIR}/tests"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/unit"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration"
)

target_link_libraries(openCSMP-tests PRIVATE openCSMP)

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

# Also copy at configure time for IDE runs (XCode, CLion etc.)
file(GLOB TEST_INPUT_FILES "${CMAKE_CURRENT_SOURCE_DIR}/tests/dataset/*")
if(TEST_INPUT_FILES)
    file(COPY ${TEST_INPUT_FILES} DESTINATION ${TEST_WORKING_DIR})
endif()

# -----------------------------------------------------------------------
# CTest Registration
# -----------------------------------------------------------------------
add_test(NAME Suite_Fundamentals    COMMAND openCSMP-tests --suite fundamentals)
add_test(NAME Suite_Interdependent1 COMMAND openCSMP-tests --suite interdependent1)
add_test(NAME Suite_Interdependent2 COMMAND openCSMP-tests --suite interdependent2)
add_test(NAME Suite_Interfaces      COMMAND openCSMP-tests --suite interfaces)
add_test(NAME Suite_Composite       COMMAND openCSMP-tests --suite composite)
add_test(NAME Suite_Refactoring     COMMAND openCSMP-tests --suite refactored)

set(_test_environment
"CSMP_DATA_ROOT=${CMAKE_CURRENT_SOURCE_DIR}/tests/dataset"
)
if(SAMG_ENABLE)
    # the SAMG shared libraries must be found at run time
    list(APPEND _test_environment
        "DYLD_LIBRARY_PATH=${SAMG_ROOT}/macosx"
        "LD_LIBRARY_PATH=${SAMG_ROOT}/linux"
    )
endif()

set_tests_properties(
Suite_Fundamentals
Suite_Interdependent1
Suite_Interdependent2
Suite_Interfaces
Suite_Composite
Suite_Refactoring
PROPERTIES
WORKING_DIRECTORY "${TEST_WORKING_DIR}"
ENVIRONMENT       "${_test_environment}"
)

# Strict sequential dependency chain — Suite_Refactoring always runs first
set_tests_properties(Suite_Fundamentals    PROPERTIES DEPENDS Suite_Refactoring)
set_tests_properties(Suite_Interdependent1 PROPERTIES DEPENDS Suite_Fundamentals)
set_tests_properties(Suite_Interdependent2 PROPERTIES DEPENDS Suite_Interdependent1)
set_tests_properties(Suite_Interfaces      PROPERTIES DEPENDS Suite_Interdependent2)
set_tests_properties(Suite_Composite       PROPERTIES DEPENDS Suite_Interfaces)

endif() # BUILD_TESTING


# ===============================================================================================
# CVFEM Verification & Validation Test Executable: CVFEM_tests_main
#
#   Stand-alone executable — built and registered independently of BUILD_TESTING.
#   Runs hydrothermal CVFEM V&V cases against benchmark data.
#
#   Enable with:
#     cmake .. -DBUILD_CVFEM_TESTS=ON
#
#   How to run:
#     ctest -R CVFEM                             — run all CVFEM tests
#     ctest -R CVFEM --output-on-failure         — show output for failures
#     ./CVFEM_tests_main                         — run directly
#
#   Benchmark data expected at:
#     tests/data-tests_CVFEM/
# ===============================================================================================
if(BUILD_CVFEM_TESTS)

    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/data-tests_CVFEM")
        message(WARNING
                "BUILD_CVFEM_TESTS is ON but benchmark data directory not found:\n"
                "  ${CMAKE_CURRENT_SOURCE_DIR}/tests/data-tests_CVFEM\n"
                "Tests will build but may fail at runtime.")
    endif()

    add_executable(CVFEM_tests_main
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/CVFEM_tests_main.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/Test.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/TestSuite.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/CVFEM/CVFEM_1D_VVCase.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/CVFEM/CVFEM_2D_VVCase.cpp
            ${CVFEM_HSS_LIB_SOURCES}
    )

    target_compile_features(CVFEM_tests_main PRIVATE cxx_std_23)
    set_target_properties(CVFEM_tests_main PROPERTIES
            CXX_STANDARD_REQUIRED ON
            CXX_EXTENSIONS        OFF
    )

    target_include_directories(CVFEM_tests_main PRIVATE
            ${CSMP_CONSUMER_INCLUDE_DIRS}
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/CVFEM"
    )

    target_link_libraries(CVFEM_tests_main PRIVATE openCSMP)

    target_compile_definitions(CVFEM_tests_main PRIVATE
            CVFEM_BENCHMARK_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data-tests_CVFEM"
    )

    if(CVFEM_HSS_SOURCES OR CVFEM_HSS_HEADERS)
        target_include_directories(CVFEM_tests_main PRIVATE
                ${CVFEM_HSS_INCLUDE_DIRS}
                "${CVFEM_HYDROTHERMAL_ROOT}/includes/hydrothermal_system_simulator"
                "${CMAKE_CURRENT_SOURCE_DIR}"
        )
    endif()

    # -----------------------------------------------------------------------
    # CVFEM Test Runtime Data
    #   Separate working directory from openCSMP-tests to avoid interference
    #   with the main test suite.
    # -----------------------------------------------------------------------
    set(CVFEM_TEST_WORKING_DIR "${CMAKE_BINARY_DIR}/cvfem_tests_runtime")

    add_custom_target(copy_cvfem_test_data ALL
            COMMAND ${CMAKE_COMMAND} -E make_directory "${CVFEM_TEST_WORKING_DIR}"
            COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/data-tests_CVFEM"
            "${CVFEM_TEST_WORKING_DIR}"
            COMMENT "Copying CVFEM benchmark data to ${CVFEM_TEST_WORKING_DIR}"
            VERBATIM
    )
    add_dependencies(CVFEM_tests_main copy_cvfem_test_data)

    # Also copy at configure time for IDE runs
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/data-tests_CVFEM")
        file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/tests/data-tests_CVFEM/"
                DESTINATION "${CVFEM_TEST_WORKING_DIR}"
        )
    endif()

    # -----------------------------------------------------------------------
    # CTest Registration
    # -----------------------------------------------------------------------
    add_test(NAME CVFEM_1D_VVCase
            COMMAND CVFEM_tests_main --suite cvfem_1d)
    add_test(NAME CVFEM_2D_VVCase
            COMMAND CVFEM_tests_main --suite cvfem_2d)

    set(_cvfem_environment
            "CSMP_DATA_ROOT=${CMAKE_CURRENT_SOURCE_DIR}/tests/data-tests_CVFEM"
    )
    if(SAMG_ENABLE)
        list(APPEND _cvfem_environment
                "DYLD_LIBRARY_PATH=${SAMG_ROOT}/macosx"
                "LD_LIBRARY_PATH=${SAMG_ROOT}/linux"
        )
    endif()

    set_tests_properties(
            CVFEM_1D_VVCase
            CVFEM_2D_VVCase
            PROPERTIES
            WORKING_DIRECTORY "${CVFEM_TEST_WORKING_DIR}"
            ENVIRONMENT       "${_cvfem_environment}"
    )

    # 2D case depends on 1D passing
    set_tests_properties(CVFEM_2D_VVCase PROPERTIES DEPENDS CVFEM_1D_VVCase)

    # Installation (optional)
    install(TARGETS CVFEM_tests_main RUNTIME DESTINATION bin OPTIONAL)

endif() # BUILD_CVFEM_TESTS


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
install(DIRECTORY non-GIT-files/        DESTINATION include/openCSMP/non-GIT-files)

# Example inputs — installed with the executables; to run the examples, copy them into a
# workspace with the 'setup_examples' target (see "Examples Workspace" above)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/example_inputs")
    install(DIRECTORY examples/example_inputs/
            DESTINATION example_inputs
            PATTERN ".DS_Store" EXCLUDE
            PATTERN "__MACOSX"  EXCLUDE
            PATTERN "Icon"      EXCLUDE
    )
endif()

# Continuous integration fixtures
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/continuous_integration_testing")
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
if(SAMG_ENABLE)
    message(STATUS "  SAMG root:          ${SAMG_ROOT}")
endif()
message(STATUS "  PETSc solver:       ${CSMP_WITH_PETSC_SOLVER}")
message(STATUS "  Image output:       ${CSMP_WITH_IMAGE_OUTPUT}")
message(STATUS "  Eigen:              ${CSMP_EIGEN_INCLUDE_DIR}")
message(STATUS "  Build examples:     ${BUILD_EXAMPLES}")
message(STATUS "  Examples workspace: ${CSMP_EXAMPLES_WORKSPACE}")
message(STATUS "  Build tests:        ${BUILD_TESTING}")
message(STATUS "  Build CVFEM tests:  ${BUILD_CVFEM_TESTS}")
message(STATUS "  LTO/IPO:            ${IPO_SUPPORTED}")
list(LENGTH LIBRARY_SOURCES _src_count)
message(STATUS "  Source files:       ${_src_count} files")
message(STATUS "------------------------------------------------------------")

