include("ctest_runner.vars")

set(CTEST_PROJECT_NAME "CSMP-API")
set(CTEST_CMAKE_GENERATOR "${CMAKE_GENERATOR}")
set(CTEST_SOURCE_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CTEST_BINARY_DIRECTORY "${CMAKE_BINARY_DIR}")

set(to_xunit_xsl "${CMAKE_SOURCE_DIR}/tests/ctest2xunit.xsl")
set(xunit_xml "${CMAKE_BINARY_DIR}/xunit.xml")

include("CsmpTestResultTransform.cmake")

ctest_start("Test")
ctest_configure(BUILD "${CMAKE_BINARY_DIR}")
ctest_build(BUILD "${CMAKE_BINARY_DIR}")
ctest_test(BUILD "${CMAKE_BINARY_DIR}")

CTEST_GET_RESULTDIR(result_dir "${CMAKE_BINARY_DIR}")
message(STATUS "result_dir is ${result_dir}")

# ctest_generate_resultindex(index_xml "${result_dir}")
file(GLOB input_xml "${result_dir}/Test.xml")

ctest_xsl_transform("${to_xunit_xsl}" "${input_xml}" "${xunit_xml}")

