# Install examples

file(GLOB_RECURSE EXAMPLES_FILES_AND_DIRS
	RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
	"${CMAKE_CURRENT_SOURCE_DIR}/*")
foreach (I ${EXAMPLES_FILES_AND_DIRS})
    if (IS_DIRECTORY "${I}")
	message(WARNING "Examining examples directory ${I}")
    elseif("${I}" MATCHES "CMakeLists.txt")
	message(WARNING "Excluding file ${I} from installation")
    else()
	message(WARNING "Examining examples file ${I}")
        get_filename_component( dir ${I} DIRECTORY )
	install( FILES "${CMAKE_CURRENT_SOURCE_DIR}/${I}"
		DESTINATION "examples/${dir}")
    endif()
endforeach()

