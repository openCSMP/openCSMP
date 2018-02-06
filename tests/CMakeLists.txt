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
	"${CMAKE_SOURCE_DIR}/thirdparty/catch/include"
	"${CMAKE_SOURCE_DIR}/thirdparty/iaps_h2o_eos"
	"${CMAKE_BINARY_DIR}/source/generated_includes"
	"${CMAKE_CURRENT_SOURCE_DIR}/unit"
	"${CMAKE_CURRENT_SOURCE_DIR}/integration"
)

if(CMAKE_HOST_WIN32)
	LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/bin/windows)
	LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/build/windows)
endif()

FILE(GLOB TEST_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/unit/*.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/integration/*.cpp
)

ADD_DEFINITIONS(-DRUNNING_UNDER_CATCH)

if(BUILD_TESTING)
    ADD_EXECUTABLE(csmptest unit_tests_catch.cpp ${TEST_SOURCES})
    TARGET_LINK_LIBRARIES(csmptest csmpmeschach csmpcore csmpiaps ${SAMG_LIBRARIES})
    if(CMAKE_HOST_WIN32)
		ENABLE_TESTING()
        ADD_TEST(
        NAME csmptest
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/data/fixtures
	    COMMAND ${CMAKE_BINARY_DIR}/tests/Release/csmptest.exe -r junit		
        )
    else()
        ADD_TEST(
	    NAME csmptest
	    COMMAND csmptest -r junit
	    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/data/fixtures
        )
    endif()

    ADD_EXECUTABLE(csmptest_main unit_tests_main.cpp ${TEST_SOURCES})
    TARGET_LINK_LIBRARIES(csmptest_main csmpmeschach csmpcore csmpiaps ${SAMG_LIBRARIES})
    if(CMAKE_HOST_WIN32)
	ENABLE_TESTING()
        ADD_TEST(
        NAME csmptest_main
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/data/fixtures
	    COMMAND ${CMAKE_BINARY_DIR}/tests/Release/csmptest_main.exe
        )
    else()
        ADD_TEST(
	    NAME csmptest_main
	    COMMAND csmptest_main
	    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/data/fixtures
        )
    endif()

endif(BUILD_TESTING)

ADD_EXECUTABLE(csmptest_development unit_tests_development.cpp ${TEST_SOURCES})
TARGET_LINK_LIBRARIES(csmptest_development csmpmeschach csmpcore csmpiaps ${SAMG_LIBRARIES})

