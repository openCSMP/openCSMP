# Install examples

file(GLOB_RECURSE EXAMPLES_FILES_AND_DIRS
	RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
	"${CMAKE_CURRENT_SOURCE_DIR}/*")
foreach (I ${EXAMPLES_FILES_AND_DIRS})
    if (IS_DIRECTORY "${I}")
	message(STATUS "Examining examples directory ${I}")
    elseif("${I}" MATCHES "CMakeLists.txt")
	message(WARNING "Excluding file ${I} from installation")
    else()
        get_filename_component( dir ${I} DIRECTORY )
	install( FILES "${CMAKE_CURRENT_SOURCE_DIR}/${I}"
		DESTINATION "examples/${dir}")
    endif()
endforeach()


SET(INCLUDES_DIR "${CMAKE_SOURCE_DIR}/source/includes")
SET(THIRDPARTY_DIR "${CMAKE_SOURCE_DIR}/thirdparty")
SET(GEN_INCLUDES_DIR "${CMAKE_BINARY_DIR}/source/generated_includes")

FILE(GLOB EXAMPLES_SRC
	${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/source/*.h
	)

FILE(GLOB EXAMPLES_EXCLUSIONS
	${CMAKE_CURRENT_SOURCE_DIR}/source/*.old.cpp
)
LIST(REMOVE_ITEM EXAMPLES_SRC ${EXAMPLES_EXCLUSIONS} )


ADD_EXECUTABLE(examples
	${EXAMPLES_SRC}
)

# XXX This should not require anything in MODULES_DIR
TARGET_INCLUDE_DIRECTORIES(examples PRIVATE
	${INCLUDES_DIR}/run_control
	${INCLUDES_DIR}/utilities
	${INCLUDES_DIR}/math
	${INCLUDES_DIR}/model
	${INCLUDES_DIR}/analysis
	${INCLUDES_DIR}/discretisation
	${INCLUDES_DIR}/property_manipulation
	${INCLUDES_DIR}/regular_grids
	${INCLUDES_DIR}/eos
	${INCLUDES_DIR}/interfaces
	${INCLUDES_DIR}/interfaces/ANSYS
	${INCLUDES_DIR}/interfaces/Eclipse
	${INCLUDES_DIR}/interfaces/VTK
	${INCLUDES_DIR}/interfaces/JPEG
	${INCLUDES_DIR}/interfaces/Matlab
	${INCLUDES_DIR}/interfaces/Maple
	${INCLUDES_DIR}/interfaces/Rhinoceros
	${INCLUDES_DIR}/integration/finite_elements
	${INCLUDES_DIR}/integration/finite_volumes
	${INCLUDES_DIR}/integration/finite_volumes/legacy_transport_scheme
	${INCLUDES_DIR}/integration/solver
	${INCLUDES_DIR}/integration/pde_integrators
	${INCLUDES_DIR}/integration/pde_operators
	${INCLUDES_DIR}/integration/pde_operators/numerical
	${INCLUDES_DIR}/integration/pde_operators/analytical
	${INCLUDES_DIR}/thmc
	${INCLUDES_DIR}/thmc/mechanics
	${INCLUDES_DIR}/thmc/single_phase_flow
	${INCLUDES_DIR}/thmc/two_phase_flow
	${THIRDPARTY_DIR}/catch/include
	${THIRDPARTY_DIR}/iaps_h2o_eos
	${GEN_INCLUDES_DIR}
)

TARGET_LINK_LIBRARIES(examples csmpmeschach csmpcore csmpiaps csmpjpeg)

