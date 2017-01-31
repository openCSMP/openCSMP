MACRO(SUBDIRLIST result curdirs)
  SET(dirlist "")
  SET(curdirlist ${curdirs} ${ARGN})
  FOREACH(curdir ${curdirlist})
    FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
    FOREACH(child ${children})
      IF(IS_DIRECTORY ${curdir}/${child})
	  LIST(APPEND dirlist ${curdir}/${child})
      ENDIF()
    ENDFOREACH()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

SUBDIRLIST(BASE_INCLUDES_DIR
	${CMAKE_SOURCE_DIR}/source/includes
)

SUBDIRLIST(BASE_INCLUDES_SUBDIR ${BASE_INCLUDES_DIR})
SUBDIRLIST(BASE_INCLUDES_SUBSUBDIR ${BASE_INCLUDES_SUBDIR})

SUBDIRLIST(MODULES_SOURCE_DIR
	${CMAKE_SOURCE_DIR}/source/modules
)

INCLUDE_DIRECTORIES(
	${BASE_INCLUDES_DIR}
	${BASE_INCLUDES_SUBDIR}
	${BASE_INCLUDES_SUBSUBDIR}
	${MODULES_SOURCE_DIR}
	"${CMAKE_SOURCE_DIR}/lib/catch/include"
	"${CMAKE_BINARY_DIR}/source/generated_includes"
	"${CMAKE_CURRENT_SOURCE_DIR}/unit"
)


ADD_DEFINITIONS(-DRUNNING_UNDER_CATCH)

FILE(GLOB TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/unit/*.cpp)
ADD_EXECUTABLE(csmptest unit_tests_catch.cpp ${TEST_SOURCES})
TARGET_LINK_LIBRARIES(csmptest csmpcore ${SAMG_LIBRARIES})
ADD_TEST(
	NAME csmptest
	COMMAND csmptest -r xml
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/data/fixtures
)

