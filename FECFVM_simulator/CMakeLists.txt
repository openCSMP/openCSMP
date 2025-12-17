project(csmp_FECFVM_simulator)

cmake_minimum_required(VERSION 3.1)
set(PACKAGE_NAME "csmp-FECFVM-simulator")
set(PACKAGE_VERSION "UoM1-alpha")
set(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")
set(PACKAGE_TARNAME "${PACKAGE_NAME}")

set(CMAKE_MODULE_PATH ${csmp_FECFVM_simulator_SOURCE_DIR}/cmake)

include(CsmpPlatform)
include(CsmpDependencies)

add_subdirectory(modules)
