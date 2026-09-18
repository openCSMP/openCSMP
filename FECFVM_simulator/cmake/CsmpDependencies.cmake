# SPDX-FileCopyrightText: © 2026 The openCSMP project
#
# SPDX-License-Identifier: LGPL-3.0-only

# CSMP directories

set(CSMP_SOURCE_DIRECTORY_INIT "${csmp-FECFVM-simulator_SOURCE_DIR}/../csmp-api-library")

set(CSMP_SOURCE_DIRECTORY "${CSMP_SOURCE_DIRECTORY_INIT}" CACHE PATH "Path of CSMP source directory")

if (CMAKE_HOST_WIN32)
    set(CSMP_BUILD_DIRECTORY_GUESS "${CSMP_SOURCE_DIRECTORY}/../csmp-api-library-build")
else()
    set(CSMP_BUILD_DIRECTORY_GUESS "${CSMP_SOURCE_DIRECTORY}/thebuild")
endif()
set(CSMP_BUILD_DIRECTORY_INIT "${CSMP_BUILD_DIRECTORY_GUESS}")

# Search likely locations for a directory
set(D)
foreach(D
    ${CSMP_BUILD_DIRECTORY_GUESS}
    ${CSMP_SOURCE_DIRECTORY}/thebuild
    ${CSMP_SOURCE_DIRECTORY}/../csmp-api-library-build
    )
    if (EXISTS "$D")
        set(CSMP_BUILD_DIRECTORY_INIT ${D})
	break()
    endif()
endforeach()

set(CSMP_BUILD_DIRECTORY ${CSMP_BUILD_DIRECTORY_INIT} CACHE PATH "Path of CSMP build directory")

MESSAGE(STATUS "Assuming CSMP_SOURCE_DIRECTORY is ${CSMP_SOURCE_DIRECTORY}")
MESSAGE(STATUS "Assuming CSMP_BUILD_DIRECTORY is ${CSMP_BUILD_DIRECTORY}")

# SAMG

set(SAMG_LIBRARIES)
set(I)
foreach(I ${PLATFORM_EXTERNAL_LIBRARIES})
    set(L "NOTFOUND")
    find_library(L
	"${I}"
	PATHS ${CSMP_SOURCE_DIRECTORY}/${SAMG_SEARCH_SUFFIX}
    )
    if(${L} EQUAL "NOTFOUND")
	MESSAGE(FATAL_ERROR "Cannot find library ${I} in directory ${CSMP_SOURCE_DIRECTORY}/${SAMG_SEARCH_SUFFIX}")
    endif()
    list(APPEND SAMG_LIBRARIES "${L}")
endforeach()


# CSMP

# Find the generated includes directory
set(CSMP_GENERATED_INCLUDES_DIR "${CSMP_BUILD_DIRECTORY}/source/generated_includes")
if(NOT EXISTS "${CSMP_GENERATED_INCLUDES_DIR}/CSMP_number_types.h")
    MESSAGE(FATAL_ERROR "Cannot find CSMP_number_types.h in directory ${CSMP_GENERATED_INCLUDES_DIR}. Was CSMP_BUILD_DIRECTORY set correctly?")
endif()

if(NOT EXISTS "${CSMP_SOURCE_DIRECTORY}/source/includes/model/Node.h")
    MESSAGE(FATAL_ERROR "Cannot find Node.h in directory ${CSMP_SOURCE_DIRECTORY}. Was CSMP_SOURCE_DIRECTORY set correctly?")
else()

endif()

set(CSMP_LIBRARIES)
set(I)

if (CMAKE_BUILD_TYPE STREQUAL "Release")
	foreach(I csmpcore csmpiaps csmpjpeg)
		set(L "NOTFOUND")
		find_library(L
		"${I}"
		PATHS
		${CSMP_BUILD_DIRECTORY}/source
		${CSMP_BUILD_DIRECTORY}/source/Release
		${CSMP_BUILD_DIRECTORY}/source/Debug
		${CSMP_BUILD_DIRECTORY}/thirdparty
		${CSMP_BUILD_DIRECTORY}/thirdparty/Release
		${CSMP_BUILD_DIRECTORY}/thirdparty/Debug
		)
		if(${L} EQUAL "NOTFOUND")
			MESSAGE(FATAL_ERROR "Cannot find library ${I} in directory ${CSMP_BUILD_DIRECTORY}")
		endif()
		list(APPEND CSMP_LIBRARIES "${L}")
	endforeach()
else()
	foreach(I csmpcore_debug csmpiaps_debug csmpjpeg_debug)
		set(L "NOTFOUND")
		find_library(L
		"${I}"
		PATHS
		${CSMP_BUILD_DIRECTORY}/source
		${CSMP_BUILD_DIRECTORY}/source/Release
		${CSMP_BUILD_DIRECTORY}/source/Debug
		${CSMP_BUILD_DIRECTORY}/thirdparty
		${CSMP_BUILD_DIRECTORY}/thirdparty/Release
		${CSMP_BUILD_DIRECTORY}/thirdparty/Debug
		)
		if(${L} EQUAL "NOTFOUND")
			MESSAGE(FATAL_ERROR "Cannot find library ${I} in directory ${CSMP_BUILD_DIRECTORY}")
		endif()
		list(APPEND CSMP_LIBRARIES "${L}")
	endforeach()
endif()

MESSAGE(STATUS "Using SAMG libraries: ${SAMG_LIBRARIES}")
MESSAGE(STATUS "Using CSMP libraries: ${CSMP_LIBRARIES}")
MESSAGE(STATUS "Using CSMP generated includes directory: ${CSMP_GENERATED_INCLUDES_DIR}")

mark_as_advanced(CSMP_INCLUDE_DIR CSMP_LIBRARIES)

