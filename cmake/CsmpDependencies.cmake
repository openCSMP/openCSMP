# SAMG

if (SAMG_INCLUDE_DIR AND SAMG_LIBRARIES)
    # in cache already
    set(SAMG_FIND_QUIETLY TRUE)
endif (SAMG_INCLUDE_DIR AND SAMG_LIBRARIES)

find_package(PkgConfig)

set(SAMG_DEFINITIONS ${PC_SAMG_CFLAGS_OTHER})

find_path(SAMG_INCLUDE_DIR NAMES samg.h
    PATHS
    ${PLATFORM_INCLUDES}
    ${CMAKE_SOURCE_DIR}/lib/samg
)

find_library(SAMG_LIBRARIES
    NAMES amg_mult ifcoremt ifportmt imf intlc iomp5 irc svml
    PATHS
    ${PLATFORM_LIBS}
)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SAMG DEFAULT_MSG SAMG_INCLUDE_DIR SAMG_LIBRARIES)
mark_as_advanced(SAMG_INCLUDE_DIR SAMG_LIBRARIES)
include_directories(SAMG_INCLUDE_DIR)


