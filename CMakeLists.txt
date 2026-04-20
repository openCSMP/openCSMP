cmake_minimum_required(VERSION 3.29)

# ===============================================================================================
#  Setup of open-csmp (library, examples and tests)
# ===============================================================================================
# prerequisites
# - C++23 capabable compiler (g++ 13.3 or later, CLang, Intel, NVC++)
# - MPI 
# - git
# - Rosetta on OSX if you want to use SAMG (Intel (x86_64))
#     cmake -DCMAKE_OSX_ARCHITECTURES=x86_64  (note that this will also require Intel MPI!)
# ===============================================================================================
# usage alternatives
#
#   cmake .. -DCMAKE_BUILD_TYPE="Release" -Dsamg=ON
#   make -j$(sysctl -n hw.ncpu)
#
#   cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release -Dsamg=ON
#   cmake --build .
#
# ===============================================================================================

# toggle Fraunhofer's SAMG Solver and its configuration ON if you have access to it!
option(SAMG_ENABLE "Enable SAMG Solver" OFF)
option(BUILD_TESTING "Enable building tests" ON)
option(BUILD_EXAMPLES "Enable building examples" ON)

if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64" AND SAMG_ENABLE)
    set(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE STRING "Build architecture for OSX" FORCE)
    message(STATUS "Configuring for x86_64 architecture to use SAMG with Rosetta 2")
endif()

project(open-csmp
        VERSION 23
        DESCRIPTION "Open-CSMP++ library"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# OFF sets default to 'lib' for shared libs as opposed to .so, .dll etc.
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)


# ===============================================================================================
# Dependencies: MPI
# ===============================================================================================
find_package(MPI REQUIRED)
message(STATUS "MPI found: ${MPI_CXX_LIBRARIES}")


# ===============================================================================================
# Source/Headers discovery (recursive)
# ===============================================================================================
# Collect .cpp sources (compile these). Collect headers for determining include dirs.
file(GLOB_RECURSE LIBRARY_SOURCES CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/source/modules/**/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/modules/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/non-GIT-files/*.cpp"
)

file(GLOB_RECURSE LIBRARY_HEADERS CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/source/modules/**/*.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/modules/**/*.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/**/*.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/source/includes/**/*.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/non-GIT-files/*.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/non-GIT-files/*.hpp"
)

# Always include plf_colony header explicitly (you had this earlier)
list(APPEND LIBRARY_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/plf_colony/plf_colony.h")

# Combine for some operations (but only .cpp go into the lib target initially)
set(LIBRARY_SOURCE_FILES ${LIBRARY_SOURCES} ${LIBRARY_HEADERS})

# Variable Sets (Warning: WILL BE DEPRECATED)
set(VARIABLE_INDEX_SETS
        ${CMAKE_CURRENT_SOURCE_DIR}/non-GIT-files/VariableSet_TracerTransfer.h
        ${CMAKE_CURRENT_SOURCE_DIR}/non-GIT-files/VariableSet_TwoPhaseFlow.h
        ${CMAKE_CURRENT_SOURCE_DIR}/non-GIT-files/VariableSet_CO2GeoSequestration.h
)
# Check if they exist to prevent obscure compile errors
foreach(v_file ${VARIABLE_INDEX_SETS})
    if(NOT EXISTS ${v_file})
        message(WARNING "External VariableSet file missing: ${v_file}")
    endif()
endforeach()


# ===============================================================================================
# Extract unique directories from headers/sources for flat includes
# ===============================================================================================
function(get_directories result file_list)
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

# Remove deprecated files (if present)
list(FILTER LIBRARY_SOURCE_FILES EXCLUDE REGEX "/deprecated/")

# Extract all directories that contain headers/sources
get_directories(ALL_SUBDIRECTORIES "${LIBRARY_SOURCE_FILES}")

# Always add top-level include roots that should be searched (flat includes)
list(APPEND ALL_SUBDIRECTORIES
        ${CMAKE_CURRENT_SOURCE_DIR}/source/includes
        ${CMAKE_CURRENT_SOURCE_DIR}/source/modules
        ${CMAKE_CURRENT_SOURCE_DIR}/source/non-GIT-files
        ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/plf_colony
        ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/jpeg
)
list(REMOVE_DUPLICATES ALL_SUBDIRECTORIES)


# ===============================================================================================
# Filter out SAMG source files by default (they are added back conditionally later)
# ===============================================================================================
if(NOT SAMG_ENABLE)
    list(FILTER LIBRARY_SOURCE_FILES EXCLUDE REGEX ".*SAMG.*")
    list(FILTER LIBRARY_SOURCE_FILES EXCLUDE REGEX ".*samg.*")
endif()


# ===============================================================================================
# Core Library Definition
# ===============================================================================================
add_library(open-csmp-core ${LIBRARY_SOURCE_FILES})

# Export include dirs for build & install (flat include model)
target_include_directories(open-csmp-core PUBLIC
        $<BUILD_INTERFACE:${ALL_SUBDIRECTORIES}>
        $<INSTALL_INTERFACE:include/open-csmp>
)


# ===============================================================================================
# Third Party Libraries (Built as part of the project)
# ===============================================================================================

# 1. IAPS H2O EOS (if sources exist)
# ----------------------------------
file(GLOB_RECURSE IAPS_SRC "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/iaps_h2o_eos/*.cpp")
list(FILTER IAPS_SRC EXCLUDE REGEX "/deprecated/")
if(IAPS_SRC)
    add_library(iaps_h2o_eos STATIC ${IAPS_SRC})
    target_include_directories(iaps_h2o_eos PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/iaps_h2o_eos)
    target_link_libraries(open-csmp-core PUBLIC iaps_h2o_eos)
else()
    message(STATUS "IAPS H2O EOS not found - skipping.")
endif()

# 2. JPEG (libjpeg-turbo would be better but does not lnk with Rosetta) 
# ---------------------------------------------------------------------
file(GLOB JPEG_SOURCE_FILES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/jpeg/*.cpp"
    RECURSE
)
list(FILTER JPEG_SOURCE_FILES EXCLUDE REGEX "thirdparty/jpeg/deprecated/.*")
add_library(jpeg STATIC ${JPEG_SOURCE_FILES})
target_include_directories(jpeg PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/jpeg")

# 3. Eigen (header-only)
# ---------------------------------------------------------------------
set(EIGEN_PATH "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/eigen")
if(CMAKE_GENERATOR MATCHES "Xcode")
    set(EIGEN_PATH "${EIGEN_PATH}/eigen")
endif()
add_library(eigen_lib INTERFACE)
target_include_directories(eigen_lib INTERFACE ${EIGEN_PATH})
target_link_libraries(open-csmp-core PUBLIC eigen_lib)

# 4. CGNS & RapidXML (Interface only)
# ---------------------------------------------------------------------
add_library(cgns INTERFACE)
target_include_directories(cgns INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/cgns)
target_link_libraries(open-csmp-core PUBLIC cgns)

add_library(rapidxml INTERFACE)
target_include_directories(rapidxml INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/rapidxml)
target_link_libraries(open-csmp-core PUBLIC rapidxml)

# Link Core to basic dependencies
target_link_libraries(open-csmp-core PUBLIC MPI::MPI_CXX)


# ===============================================================================================
# Feature: SAMG Solver
# ===============================================================================================
if(SAMG_ENABLE)
    message(STATUS "Enabling SAMG Solver...")
    target_compile_definitions(open-csmp-core PUBLIC
            CSMP_WITH_SAMG_SOLVER
            SAMG_UNIX_LINUX
            SAMG_LCASE_USCORE
            SAMG_MULTIPLE_INSTANCES
            LEGACY_SAMG
    )

    file(GLOB_RECURSE SAMG_SRC
            "${CMAKE_CURRENT_SOURCE_DIR}/source/modules/integration/solver/SAMG_*.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/source/modules/integration/solver/SAMG_profiles/*.cpp"
    )
    if(SAMG_SRC)
        target_sources(open-csmp-core PRIVATE ${SAMG_SRC})

        target_include_directories(open-csmp-core PUBLIC
                ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/samg
                ${CMAKE_CURRENT_SOURCE_DIR}/source/modules/integration/solver/SAMG_profiles
        )
    endif()

    if(APPLE)
        set(SAMG_LIB_PATH "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/samg/macosx")
        set(SAMG_LIB "${SAMG_LIB_PATH}/libamg_mult.dylib")
        set(OMP_LIB "${SAMG_LIB_PATH}/libiomp5.dylib")
    else()
        find_library(SAMG_LIB NAMES amg_mult HINTS "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/samg/linux")
        find_library(OMP_LIB NAMES iomp5 HINTS "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/samg/linux")
    endif()

    if(SAMG_LIB AND OMP_LIB)
        target_link_libraries(open-csmp-core PUBLIC ${SAMG_LIB} ${OMP_LIB})
    else()
        message(WARNING "SAMG requested but SAMG_LIB or OMP_LIB is missing. Disabling SAMG support.")
        set(WITH_SAMG OFF)
    endif()
endif()


# ===============================================================================================
# Optimization
# ===============================================================================================
include(CheckIPOSupported)
check_ipo_supported(RESULT IPO_SUPPORTED)
if(IPO_SUPPORTED)
    set_target_properties(open-csmp-core PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()


# ===============================================================================================
# Examples Executable
# ===============================================================================================
if(BUILD_EXAMPLES)
    # 1. Collect Example Sources
    file(GLOB_RECURSE EXAMPLES_SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/examples/source/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/examples/source/**/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/examples/source/*.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/examples/source/**/*.h"
    )

    # 2. Including CSMP API Library Headers + Tests because they are demo'ed in examples
    list(APPEND EXAMPLES_SRC 
        ${LIBRARY_HEADERS}
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/TestSuite.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/TestSuite.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/Test.cpp"
    )

    # 3. Exclude deprecated directory
    list(FILTER EXAMPLES_SRC EXCLUDE REGEX "${CMAKE_CURRENT_SOURCE_DIR}/examples/source/deprecated/.*")

    # 4. Create Executable
    if(EXAMPLES_SRC)
        add_executable(open-csmp-examples ${EXAMPLES_SRC})
        
        target_include_directories(open-csmp-examples PRIVATE
          "${CMAKE_CURRENT_SOURCE_DIR}"
          "${CMAKE_CURRENT_SOURCE_DIR}/examples/source/interrelations"
          "${CMAKE_CURRENT_SOURCE_DIR}/examples/source/visitors"
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
          "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/jpeg"
        )
        # This links libjpeg-turbo or standard libjpeg
        target_link_libraries(open-csmp-examples PRIVATE open-csmp-core jpeg )
    else()
        message(STATUS "No example sources found.")
    endif()
endif()


# ===============================================================================================
# Testing (unit and integration tests)
# ===============================================================================================
if(BUILD_TESTING)
    include(CTest)
    enable_testing()

    file(GLOB_RECURSE TEST_FILES
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/unit/*.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/*.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/unit/*.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/*.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/Test.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/TestSuite.cpp"
    )
    list(FILTER TEST_FILES EXCLUDE REGEX "/deprecated/")

    add_executable(open-csmp-tests ${CMAKE_CURRENT_SOURCE_DIR}/tests/unit_CTests_main.cpp ${TEST_FILES})
    target_include_directories(open-csmp-tests PRIVATE
          "${CMAKE_CURRENT_SOURCE_DIR}/tests"
          "${CMAKE_CURRENT_SOURCE_DIR}/tests/unit"
          "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration"
          "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/plf_colony"
          "${CMAKE_CURRENT_SOURCE_DIR}"
          "${CMAKE_CURRENT_SOURCE_DIR}/examples/source/interrelations"
          "${CMAKE_CURRENT_SOURCE_DIR}/examples/source/visitors"
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
          "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/jpeg"
    )
    target_link_libraries( open-csmp-tests PRIVATE open-csmp-core )

    # Directory to place test runtime data
    set(TEST_WORKING_DIR "${CMAKE_BINARY_DIR}/tests_runtime")

    # Create the directory at build time
    add_custom_target(copy_test_data ALL
            COMMAND ${CMAKE_COMMAND} -E make_directory "${TEST_WORKING_DIR}"
            COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/dataset"
            "${TEST_WORKING_DIR}"
            COMMENT "Copying test dataset to ${TEST_WORKING_DIR}"
    )
    add_dependencies(open-csmp-tests copy_test_data)

    file(GLOB TEST_INPUT_FILES "${CMAKE_CURRENT_SOURCE_DIR}/tests/dataset/*")
    if(TEST_INPUT_FILES)
        file(COPY ${TEST_INPUT_FILES} DESTINATION ${TEST_WORKING_DIR})
    endif()

    add_test(NAME open-csmp-tests COMMAND open-csmp-tests)
    set_tests_properties(open-csmp-tests PROPERTIES WORKING_DIRECTORY ${TEST_WORKING_DIR})


# Enable testing
enable_testing()

# Define the individual test suites
    add_test(NAME Suite_Fundamentals     COMMAND open-csmp-tests --suite fundamentals)
    add_test(NAME Suite_Interdependent1  COMMAND open-csmp-tests --suite interdependent1)
    add_test(NAME Suite_Interdependent2  COMMAND open-csmp-tests --suite interdependent2)
    add_test(NAME Suite_Interfaces       COMMAND open-csmp-tests --suite interfaces)
    add_test(NAME Suite_Composite        COMMAND open-csmp-tests --suite composite)

# Enforce the strict execution sequence
    set_tests_properties(Suite_Interdependent1 PROPERTIES DEPENDS Suite_Fundamentals)
    set_tests_properties(Suite_Interdependent2 PROPERTIES DEPENDS Suite_Interdependent1)
    set_tests_properties(Suite_Interfaces      PROPERTIES DEPENDS Suite_Interdependent2)
    set_tests_properties(Suite_Composite       PROPERTIES DEPENDS Suite_Interfaces)

# NOTE - how to perform the testing:
#    compile your code and run ctest in your build directory.
#
#    To run normally: ctest
#
#    To see the detailed output of a failing test: ctest --output-on-failure
#
#    To run a specific suite manually: ctest -R Suite_Fundamentals
#
#    CTest will give a pass/fail summary broken down by suite, 
#    while the standard test output will print test-by-test Passed/Failed counts 
#    generated by your custom TestSuite objects.
#
#    get nice dashboard with JUnit XML format using option: ctest --output-junit test-dashboard.xml

endif()


# ===============================================================================================
# Installation
# ===============================================================================================
install(TARGETS open-csmp-core
        EXPORT open-csmp-targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
)

if(TARGET open-csmp-examples)
    install(TARGETS open-csmp-examples RUNTIME DESTINATION bin)
endif()

if(TARGET open-csmp-tests)
    install(TARGETS open-csmp-tests RUNTIME DESTINATION bin)
endif()

# Install headers (keep install layout simple: include/open-csmp/*)
install(DIRECTORY source/includes/ DESTINATION include/open-csmp)
install(DIRECTORY source/non-GIT-files/ DESTINATION include/open-csmp/non-GIT-files)

# Install example inputs and CI fixtures (if present)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/example_inputs")
    if(TARGET open-csmp-examples)
        add_custom_command(TARGET open-csmp-examples POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_directory
                "${CMAKE_CURRENT_SOURCE_DIR}/examples/example_inputs"
                "${CMAKE_BINARY_DIR}/example_inputs"
                || ${CMAKE_COMMAND} -E echo "WARNING: Could not update example_inputs – probably locked by another process"
                COMMENT "Copying example_inputs (non-fatal)"
                VERBATIM
        )
    endif()

    # Keep the install rule unchanged
    install(DIRECTORY examples/example_inputs/
            DESTINATION example_inputs
            PATTERN ".DS_Store" EXCLUDE
            PATTERN "__MACOSX" EXCLUDE
            PATTERN "Icon" EXCLUDE)
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/continuous_integration_testing")
install(DIRECTORY continuous_integration_testing/ DESTINATION continuous_integration_testing)
endif()


# ===============================================================================================
# Diagnostics / Summary of what was found
# ===============================================================================================
message(STATUS "Found ${LIBRARY_SOURCES} source files.")
message(STATUS "Header dirs added to include path: ${ALL_SUBDIRECTORIES}")
message(STATUS "SAMG Enabled: ${SAMG_ENABLE}")
message(STATUS "Build examples: ${BUILD_EXAMPLES}")
message(STATUS "Build tests: ${BUILD_TESTING}")
