

# iaps_h2o_eos
FILE(GLOB IAPS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/iaps_h2o_eos/*.cpp )
ADD_LIBRARY(csmpiaps ${IAPS_SRC})
TARGET_INCLUDE_DIRECTORIES(csmpiaps PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/iaps_h2o_eos)
install(TARGETS csmpiaps ARCHIVE DESTINATION lib)

# jpeg
FILE(GLOB JPEG_SRC ${CMAKE_CURRENT_SOURCE_DIR}/jpeg/*.cpp )
ADD_LIBRARY(csmpjpeg ${JPEG_SRC})
TARGET_INCLUDE_DIRECTORIES(csmpjpeg PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/jpeg)
install(TARGETS csmpjpeg ARCHIVE DESTINATION lib)


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


