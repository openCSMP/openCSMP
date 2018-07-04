####################################################
## CMAKE Makefile builder
## CSMP++
## Created by: Julian E. Mindel
## Creation Date: 19.12.2013
## Update: 02.04.2014
## For support please email: julian.mindel@gmail.com
## If you are not sure of what you are doing, please ask! 
## you can type:  
## ccmake . (unix)
## cmake . (unix)
## cmake-gui (linux , windows)
## IMPORTANT NOTE: This cmake lists file will generate makefiles for all projects
## In unix/linux, using parallel building via make -j will likely cause linking errors
## since make (independently of cmake), attempts to build more targets in parallel
## and does not wait for a dependency to be built.
## This does seem to work WITHIN one subdirectory though, however in our case, each
## application or support library has its own, so using subdirectories is quite useful.
## The only way found so far to get around this is to "make -j targetname" (one by one)
## Or simply run make -j several times untill all targets are built so the dependencies are found.
#####################################################   

cmake_minimum_required(VERSION 2.8)
project(csmpALL)

SET(BUILD_UNIT_TESTS ON CACHE BOOL OFF)
SET(BUILD_TEST_CASES ON CACHE BOOL OFF)
SET(BUILD_VV_CASES ON CACHE BOOL OFF)

SET(BUILD_EXAMPLES ON CACHE BOOL OFF)

SET(BUILD_FECFV_SIMULATOR2D ON CACHE BOOL OFF)
SET(BUILD_RESERVOIR_SIMULATOR ON CACHE BOOL OFF)
SET(BUILD_DFMUPSCALING ON CACHE BOOL OFF)
SET(BUILD_EQPERMTENSOR OFF CACHE BOOL ON)
SET(BUILD_BOUSSINESQ_GEOTHERMAL_SIMULATOR ON CACHE BOOL OFF)
SET(BUILD_GEOTHERMAL_SIMULATOR ON CACHE BOOL OFF)
SET(BUILD_GEOCCS_SIMULATOR OFF CACHE BOOL OFF)
SET(BUILD_HYDROGEOCHEM_SIMULATOR ON CACHE BOOL OFF) ## necessary for the HydroGeoChemicalSimulator. Needs GEMIPM2K, below.

SET(BUILD_GEMIPM2K ON CACHE BOOL OFF) ## necessary for the HydroGeoChemicalSimulator.
SET(BUILD_MESCHACH OFF CACHE BOOL OFF)

IF(BUILD_GEMIPM2K)
   SET(CSMP_WITH_GEM2K 1)
ENDIF(BUILD_GEMIPM2K)

IF(BUILD_MESCHACH)
   SET(CSMP_WITH_MESCHACH 1)
ENDIF(BUILD_MESCHACH)

SET(BUILD_EOS_H2O_LIBS ON CACHE BOOL OFF)

SET(CMAKE_DEBUG_POSTFIX "")
IF(${CMAKE_BUILD_TYPE} MATCHES "Debug")
   SET(CMAKE_DEBUG_POSTFIX "-debug")
ENDIF(${CMAKE_BUILD_TYPE} MATCHES "Debug")

SET(OBJ_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/../lib/obj${CMAKE_DEBUG_POSTFIX}")

# Build CSMP library (default)
ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/../trunk/CMake ${OBJ_OUTPUT_PATH}/csmp_core)

IF(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/../support_libraries/gemipm2k AND BUILD_GEMIPM2K)
    message(STATUS "Will configure the gemipm2k")
    ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/../support_libraries/gemipm2k/CMake ${OBJ_OUTPUT_PATH}/gemipm2k)
ENDIF(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/../support_libraries/gemipm2k AND BUILD_GEMIPM2K)

# Build related applications.
IF(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/../tests/unit_tests AND BUILD_UNIT_TESTS)
    message(STATUS "Will configure the unit_tests")
    ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/../tests/unit_tests/CMake ${OBJ_OUTPUT_PATH}/unit_tests)
ENDIF()

IF(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/../tests/test_cases AND BUILD_TEST_CASES)
    message(STATUS "Will configure the test_cases")
    ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/../tests/test_cases/CMake ${OBJ_OUTPUT_PATH}/test_cases)
ENDIF()

IF(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/../tests/vv_cases AND BUILD_VV_CASES)
    message(STATUS "Will configure the vv_cases")
    ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/../tests/vv_cases/CMake ${OBJ_OUTPUT_PATH}/vv_cases)
ENDIF()

IF(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/../applications/examples AND BUILD_EXAMPLES)
    message(STATUS "Will configure the Examples")
    ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/../applications/examples/CMake ${OBJ_OUTPUT_PATH}/examples)
ENDIF()

IF(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/../applications/reservoir_simulator AND BUILD_RESERVOIR_SIMULATOR)
    message(STATUS "Will configure the reservoir_simulator")
    ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/../applications/reservoir_simulator/CMake ${OBJ_OUTPUT_PATH}/reservoir_simulator)
ENDIF()

IF(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/../applications/FECFV_Simulator AND BUILD_FECFV_SIMULATOR2D)
    message(STATUS "Will configure the FECFV_Simulator")
    ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/../applications/FECFV_Simulator/CMake ${OBJ_OUTPUT_PATH}/FECFV_Simulator)
ENDIF()

IF(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/../applications/BoussinesqGeothermalSimulator AND BUILD_BOUSSINESQ_GEOTHERMAL_SIMULATOR)
    message(STATUS "Will configure the BoussinesqGeothermalSimulator")
    ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/../applications/BoussinesqGeothermalSimulator/CMake ${OBJ_OUTPUT_PATH}/BoussinesqGeothermalSimulator)
ENDIF()

IF(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/../applications/GTS AND BUILD_GEOTHERMAL_SIMULATOR)
    message(STATUS "Will configure the GTS")
    ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/../applications/GTS/CMake ${OBJ_OUTPUT_PATH}/GTS)
ENDIF()

IF(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/../applications/GeologicalCarbonCaptureAndStorageSimulator AND BUILD_GEOCCS_SIMULATOR)
    message(STATUS "Will configure the GEOCCS Simulator")
    ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/../applications/GeologicalCarbonCaptureAndStorageSimulator/CMake ${OBJ_OUTPUT_PATH}/GeologicalCarbonCaptureAndStorageSimulator)
ENDIF()

IF(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/../applications/EqPermTensor AND BUILD_EQPERMTENSOR)
    message(STATUS "Will configure the EqPermTensor applications")
    ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/../applications/EqPermTensor/CMake ${OBJ_OUTPUT_PATH}/EqPermTensor)
ENDIF()

IF(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/../applications/DFM_upscaling AND BUILD_DFMUPSCALING)
    message(STATUS "Will configure the DFM_upscaling applications")
    ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/../applications/DFM_upscaling/CMake ${OBJ_OUTPUT_PATH}/DFM_upscaling)
ENDIF()

IF(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/../applications/HydroGeoChemicalSimulator AND BUILD_HYDROGEOCHEM_SIMULATOR)
    message(STATUS "Will configure the HydroGeoChemicalSimulator applications")
    ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/../applications/HydroGeoChemicalSimulator/CMake ${OBJ_OUTPUT_PATH}/HydroGeoChemicalSimulator)
ENDIF()
