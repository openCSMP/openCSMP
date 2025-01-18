cmake_minimum_required(VERSION 3.10)

project(open-csmp)
set(PACKAGE_NAME "csmp-api")
set(PACKAGE_VERSION "UoM1-alpha")
set(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")
set(PACKAGE_TARNAME "${PACKAGE_NAME}")


if(DEFINED ENV{csmp_api_SOURCE_DIR})
    message("Environment Variable csmp_api_SOURCE_DIR is defined to $ENV{csmp_api_SOURCE_DIR}")
else()
    message(WARNING "Environment Variable csmp_api_SOURCE_DIR NOT defined: This will likely lead to compilation issues. Please define variable to location of open-csmp")
endif()

set(CMAKE_MODULE_PATH $ENV{csmp_api_SOURCE_DIR}/cmake)

message("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
message("csmp_api_SOURCE_DIR: $ENV{csmp_api_SOURCE_DIR}")

include(CTest)
include(CsmpPlatform)
include(CsmpDependencies)

add_subdirectory(source)
add_subdirectory(thirdparty)
#JC: skip temperally since some examples related to pde_integrator_uom need to be updated
#QS: should be fixed now
add_subdirectory(examples)
add_subdirectory(tests)

FILE(COPY cmake/CsmpTestResultTransform.cmake DESTINATION ${CMAKE_BINARY_DIR})
FILE(WRITE "${CMAKE_BINARY_DIR}/CMakeLists.txt" "")
FILE(COPY ctest_runner.cmake DESTINATION ${CMAKE_BINARY_DIR})
FILE(COPY CTestConfig.cmake DESTINATION ${CMAKE_BINARY_DIR})
CONFIGURE_FILE(ctest_runner.vars.in ${CMAKE_BINARY_DIR}/ctest_runner.vars)
