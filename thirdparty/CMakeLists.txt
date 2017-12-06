

# iaps_h2o_eos
FILE(GLOB IAPS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/iaps_h2o_eos/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/iaps_h2o_eos/*.h )
ADD_LIBRARY(csmpiaps ${IAPS_SRC})
TARGET_INCLUDE_DIRECTORIES(csmpiaps PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/iaps_h2o_eos)
install(TARGETS csmpiaps ARCHIVE DESTINATION lib)

# jpeg
FILE(GLOB JPEG_SRC ${CMAKE_CURRENT_SOURCE_DIR}/jpeg/*.cpp  ${CMAKE_CURRENT_SOURCE_DIR}/jpeg/*.h )
ADD_LIBRARY(csmpjpeg ${JPEG_SRC})
TARGET_INCLUDE_DIRECTORIES(csmpjpeg PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/jpeg)
install(TARGETS csmpjpeg ARCHIVE DESTINATION lib)

# meschach
FILE(GLOB MESCHACH_SRC ${CMAKE_CURRENT_SOURCE_DIR}/meschach/source/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/meschach/source/*.h )
ADD_LIBRARY(csmpmeschach ${MESCHACH_SRC})
TARGET_INCLUDE_DIRECTORIES(csmpjpeg PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/meschach/source)
install(TARGETS csmpmeschach ARCHIVE DESTINATION lib)


# Header files
FILE(GLOB PUBLIC_INCLUDE_FILES
	RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
	"${CMAKE_CURRENT_SOURCE_DIR}/*/*.h"
	)
foreach ( I ${PUBLIC_INCLUDE_FILES} )
    if (IS_DIRECTORY "${I}")
    elseif("${I}" MATCHES ".h$")
        get_filename_component( dir ${I} DIRECTORY )
	install( FILES "${CMAKE_CURRENT_SOURCE_DIR}/${I}"
		DESTINATION "include/thirdparty/${dir}")
    else()
	message(WARNING "Excluding file ${I} from installation")
    endif()
endforeach()


