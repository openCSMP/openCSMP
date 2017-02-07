project(csmp_api)

cmake_minimum_required(VERSION 3.1)
set(PACKAGE_NAME "csmp-api")
set(PACKAGE_VERSION "0.0.0a")
set(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")
set(PACKAGE_TARNAME "${PACKAGE_NAME}")

set(CMAKE_MODULE_PATH ${csmp_api_SOURCE_DIR}/cmake)

include(CTest)
include(CsmpPlatform)
include(CsmpDependencies)

add_subdirectory(source)



FILE(COPY cmake/CsmpTestResultTransform.cmake DESTINATION ${CMAKE_BINARY_DIR})
CONFIGURE_FILE(ctest_runner.cmake.in ${CMAKE_BINARY_DIR}/ctest_runner.cmake)
add_subdirectory(tests)

