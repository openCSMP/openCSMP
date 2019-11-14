include("ctest_runner.vars")

set(CTEST_PROJECT_NAME "CSMP-API")
set(CTEST_CMAKE_GENERATOR "${CMAKE_GENERATOR}")
set(CTEST_SOURCE_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CTEST_BINARY_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE 1000000)
set(CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE 1000000)

set(to_xunit_xsl "${CMAKE_SOURCE_DIR}/tests/ctest2junit.xsl")
set(xunit_xml "${CMAKE_BINARY_DIR}/xunit.xml")

include("${CMAKE_SOURCE_DIR}/cmake/CsmpTestResultTransform.cmake")

ctest_start("Test")
ctest_configure(BUILD "${CMAKE_BINARY_DIR}")
ctest_build(BUILD "${CMAKE_BINARY_DIR}")
ctest_test(BUILD "${CMAKE_BINARY_DIR}")

CTEST_GET_RESULTDIR(result_dir "${CMAKE_BINARY_DIR}")
message("The result is ${result_dir}")
CTEST_GENERATE_RESULTINDEX(input_xml "${result_dir}")
file(GLOB input_xml "${result_dir}/Test.xml")
message("The xunit is ${CMAKE_BINARY_DIR}/xunit.xml")
ctest_xsl_transform("${to_xunit_xsl}" "${input_xml}" "${xunit_xml}")

