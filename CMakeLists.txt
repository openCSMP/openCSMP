project(csmp_api)

cmake_minimum_required(VERSION 3.1)
set(PACKAGE_NAME "csmp-api")
set(PACKAGE_VERSION "UoM1-alpha")
set(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")
set(PACKAGE_TARNAME "${PACKAGE_NAME}")

set(CMAKE_MODULE_PATH ${csmp_api_SOURCE_DIR}/cmake)

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
