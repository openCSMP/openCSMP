
SET(INCLUDES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/includes")
SET(MODULES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/modules")
SET(GEN_INCLUDES_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated_includes")
FILE(MAKE_DIRECTORY ${GEN_INCLUDES_DIR})

# Preprocessors and generators
SET_SOURCE_FILES_PROPERTIES(CSMP_number_types_main PROPERTIES GENERATED TRUE)
ADD_EXECUTABLE(CSMP_number_types_main
	${CMAKE_CURRENT_SOURCE_DIR}/generate_number_types.cpp
)

SET_SOURCE_FILES_PROPERTIES(CSMP_variable_set_processor PROPERTIES GENERATED TRUE)
ADD_EXECUTABLE(CSMP_variable_set_processor
	variable_set_processor.cpp csv.h
)

# Generated includes
ADD_CUSTOM_TARGET(GeneratedIncludes)
SET_SOURCE_FILES_PROPERTIES(GeneratedIncludes PROPERTIES SYMBOLIC TRUE)

# CSMP number types header

SET_SOURCE_FILES_PROPERTIES(CSMP_number_types.h PROPERTIES GENERATED TRUE)
# SET_SOURCE_FILES_PROPERTIES(GenerateNumberTypes PROPERTIES SYMBOLIC TRUE)
ADD_CUSTOM_TARGET(GenerateNumberTypes DEPENDS ${GEN_INCLUDES_DIR}/CSMP_number_types.h)
ADD_CUSTOM_COMMAND(
    OUTPUT ${GEN_INCLUDES_DIR}/CSMP_number_types.h
    COMMAND $<TARGET_FILE:CSMP_number_types_main>
    MAIN_DEPENDENCY CSMP_number_types_main
    WORKING_DIRECTORY ${GEN_INCLUDES_DIR}
    COMMENT "Generating CSMP_number_types.h"
)
ADD_DEPENDENCIES(GeneratedIncludes GenerateNumberTypes)
INSTALL(FILES "${GEN_INCLUDES_DIR}/CSMP_number_types.h" DESTINATION include)


# Variables
macro(GenerateVariables structname source_csv target_header)
    # SET_SOURCE_FILES_PROPERTIES(Generate${structname} PROPERTIES SYMBOLIC TRUE)
    SET_SOURCE_FILES_PROPERTIES(${target_header} PROPERTIES GENERATED TRUE)
    ADD_CUSTOM_TARGET(Generate${structname} DEPENDS ${GEN_INCLUDES_DIR}/${target_header})
    ADD_CUSTOM_COMMAND(
        OUTPUT ${GEN_INCLUDES_DIR}/${target_header}
        COMMAND $<TARGET_FILE:CSMP_variable_set_processor>
        ARGS ${structname} ${source_csv} ${target_header}
	      DEPENDS CSMP_variable_set_processor ${source_csv}
        WORKING_DIRECTORY ${GEN_INCLUDES_DIR}
        COMMENT "Generating ${target_header}"
    )
    ADD_DEPENDENCIES(GeneratedIncludes Generate${structname})
    INSTALL(FILES "${GEN_INCLUDES_DIR}/${target_header}" DESTINATION include)
endmacro()

GenerateVariables(VariableSet_TracerTransfer ${CMAKE_CURRENT_SOURCE_DIR}/includes/integration/finite_volumes/generic_transport_scheme/VariableSet_TracerTransfer.csv VariableSet_TracerTransfer.h)
GenerateVariables(Variables_TwoPhaseFlow ${CMAKE_CURRENT_SOURCE_DIR}/includes/integration/finite_volumes/generic_transport_scheme/VariableSet_TwoPhaseFlow.csv Variables_TwoPhaseFlow.h)
GenerateVariables(VariableSet_CO2GeoSequestration ${CMAKE_CURRENT_SOURCE_DIR}/includes/integration/finite_volumes/generic_transport_scheme/VariableSet_CO2GeoSequestration.csv VariableSet_CO2GeoSequestration.h)



# CSMP core library

file(GLOB_RECURSE CSMP_CORE_SOURCE_FILES "*.cpp")
file(GLOB_RECURSE CSMP_CORE_SOURCE_EXCLUSIONS_OLD "*.old.cpp")
file(GLOB_RECURSE CSMP_CORE_SOURCE_EXCLUSIONS_NEW "*.new.cpp")
FILE(GLOB CSMP_CORE_SOURCE_EXCLUSIONS
	${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
	${MODULES_DIR}/interfaces/CGNS/CGNS_*.cpp
	${MODULES_DIR}/interfaces/CGNS/CGNS_*.cpp
	${MODULES_DIR}/model/SplitBoundaryInterface.cpp
	${MODULES_DIR}/run_control/ExampleSuiteMainWindow.cpp
	#${MODULES_DIR}/integration/finite_volumes/generic_transport_scheme/source/*.cpp
)
LIST(REMOVE_ITEM CSMP_CORE_SOURCE_FILES ${CSMP_CORE_SOURCE_EXCLUSIONS} ${CSMP_CORE_SOURCE_EXCLUSIONS_OLD} ${CSMP_CORE_SOURCE_EXCLUSIONS_NEW})

file(GLOB_RECURSE CSMP_CORE_HEADER_FILES "*.h")

macro(GroupSources curdir)
  file(GLOB children RELATIVE ${PROJECT_SOURCE_DIR}/source/${curdir} ${PROJECT_SOURCE_DIR}/source/${curdir}/*)
  foreach(child ${children})
    if(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/source/${curdir}/${child})
      GroupSources(${curdir}/${child})
    else()
      if (CMAKE_HOST_APPLE)
	string(REGEX REPLACE "^modules/" "Source Files/" groupname ${curdir})
        string(REPLACE "/" "\\\\" groupname ${groupname})
      elseif(CMAKE_HOST_WIN32)
	string(REGEX REPLACE "^modules/" "Source/" groupname ${curdir})
        string(REPLACE "/" "\\" groupname ${groupname})
      endif()

      source_group("${groupname}" FILES ${PROJECT_SOURCE_DIR}/source/${curdir}/${child})
    endif()
  endforeach()
endmacro()

GroupSources(modules)

macro(GroupHeaders curdir)
  file(GLOB children RELATIVE ${PROJECT_SOURCE_DIR}/source/${curdir} ${PROJECT_SOURCE_DIR}/source/${curdir}/*)
  foreach(child ${children})
    if(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/source/${curdir}/${child})
      GroupHeaders(${curdir}/${child})
    else()
      if (CMAKE_HOST_APPLE)
	string(REGEX REPLACE "^includes/" "Header Files/" groupname ${curdir})
        string(REPLACE "/" "\\\\" groupname ${groupname})
      elseif(CMAKE_HOST_WIN32)
	string(REGEX REPLACE "^includes/" "Include/" groupname ${curdir})
        string(REPLACE "/" "\\" groupname ${groupname})
      endif()

      source_group("${groupname}" FILES ${PROJECT_SOURCE_DIR}/source/${curdir}/${child})
    endif()
  endforeach()
endmacro()

GroupHeaders(includes)

ADD_LIBRARY(csmpcore ${CSMP_CORE_SOURCE_FILES} ${CSMP_CORE_HEADER_FILES})

ADD_DEPENDENCIES(csmpcore GeneratedIncludes)

TARGET_LINK_LIBRARIES(csmpcore ${SAMG_LIBRARIES})

# XXX We shouldn't depend on the examples directory.
TARGET_INCLUDE_DIRECTORIES(csmpcore PRIVATE
    ${GEN_INCLUDES_DIR}
    ${INCLUDES_DIR}/analysis
    ${INCLUDES_DIR}/thmc/fractures
    ${INCLUDES_DIR}/thmc/geochemistry
    ${INCLUDES_DIR}/thmc/mechanics
    ${INCLUDES_DIR}/thmc/single_phase_flow
    ${INCLUDES_DIR}/thmc/two_phase_flow
    ${INCLUDES_DIR}/interfaces/ANSYS
    ${INCLUDES_DIR}/interfaces/CGNS
    ${INCLUDES_DIR}/interfaces/Eclipse
    ${INCLUDES_DIR}/interfaces/Eclipse_UoM
    ${INCLUDES_DIR}/interfaces/Fracman
    ${INCLUDES_DIR}/interfaces/GeoModeller
    ${INCLUDES_DIR}/interfaces/gOcad_SKUA
    ${INCLUDES_DIR}/interfaces/JPEG
    ${INCLUDES_DIR}/interfaces/Maple
    ${INCLUDES_DIR}/interfaces/Matlab
    ${INCLUDES_DIR}/interfaces/Rhinoceros
    ${INCLUDES_DIR}/interfaces/Tecplot
    ${INCLUDES_DIR}/interfaces/VTK
    ${INCLUDES_DIR}/integration/finite_elements
    ${INCLUDES_DIR}/integration/finite_volumes
    ${INCLUDES_DIR}/integration/finite_volumes/generic_transport_scheme
    ${INCLUDES_DIR}/integration/finite_volumes/legacy_transport_scheme
    ${INCLUDES_DIR}/integration/finite_volumes/DES_transport_scheme
    ${INCLUDES_DIR}/integration/pde_integrators
    ${INCLUDES_DIR}/integration/pde_operators
    ${INCLUDES_DIR}/integration/pde_operators/analytical
    ${INCLUDES_DIR}/integration/pde_operators/numerical
    ${INCLUDES_DIR}/integration/solver
    ${INCLUDES_DIR}/discretisation
    ${INCLUDES_DIR}/eos
    ${INCLUDES_DIR}/integration
    ${INCLUDES_DIR}/interfaces
    ${INCLUDES_DIR}/math
    ${INCLUDES_DIR}/model
    ${INCLUDES_DIR}/property_manipulation
    ${INCLUDES_DIR}/regular_grids
    ${INCLUDES_DIR}/run_control
    ${INCLUDES_DIR}/thmc
    ${INCLUDES_DIR}/utilities
    ${MODULES_DIR}/analysis
    ${MODULES_DIR}/eos/h2o
    ${MODULES_DIR}/eos/h2o_nacl
    ${MODULES_DIR}/eos/h2o_nacl_lookup
    ${MODULES_DIR}/integration/cvfe_hydrothermal
    ${MODULES_DIR}/integration/finite_elements
    ${MODULES_DIR}/integration/finite_volumes
    ${MODULES_DIR}/integration/finite_volumes/generic_transport_sch
    ${MODULES_DIR}/integration/finite_volumes/generic_transport_sch
    ${MODULES_DIR}/integration/finite_volumes/legacy_transport_sche
    ${MODULES_DIR}/integration/pde_integrators
    ${MODULES_DIR}/integration/pde_operators
    ${MODULES_DIR}/integration/pde_operators/analytical_integration
    ${MODULES_DIR}/integration/pde_operators/fe_numerical_integrati
    ${MODULES_DIR}/integration/pde_operators/fv_integration
    ${MODULES_DIR}/integration/solver
    ${MODULES_DIR}/integration/solver/SAMG_profiles
    ${MODULES_DIR}/interfaces/ANSYS
    ${MODULES_DIR}/interfaces/CGNS
    ${MODULES_DIR}/interfaces/Eclipse
    ${MODULES_DIR}/interfaces/Eclipse_UoM
    ${MODULES_DIR}/interfaces/Fracman
    ${MODULES_DIR}/interfaces/GeoModeller
    ${MODULES_DIR}/interfaces/gOcad_SKUA
    ${MODULES_DIR}/interfaces/JPEG
    ${MODULES_DIR}/interfaces/Maple
    ${MODULES_DIR}/interfaces/Matlab
    ${MODULES_DIR}/interfaces/Rhinoceros
    ${MODULES_DIR}/interfaces/Tecplot
    ${MODULES_DIR}/interfaces/VTK
    ${MODULES_DIR}/property_manipulation/interrelations
    ${MODULES_DIR}/thmc/fractures
    ${MODULES_DIR}/thmc/geochemistry
    ${MODULES_DIR}/thmc/mechanics
    ${MODULES_DIR}/thmc/single_phase_flow
    ${MODULES_DIR}/thmc/two_phase_flow
    ${MODULES_DIR}/discretisation
    ${MODULES_DIR}/eos
    ${MODULES_DIR}/integration
    ${MODULES_DIR}/interfaces
    ${MODULES_DIR}/math
    ${MODULES_DIR}/model
    ${MODULES_DIR}/property_manipulation
    ${MODULES_DIR}/regular_grids
    ${MODULES_DIR}/run_control
    ${MODULES_DIR}/thmc
    ${MODULES_DIR}/utilities
    ${CMAKE_SOURCE_DIR}/lib/catch/include
    ${CMAKE_SOURCE_DIR}/thirdparty/iaps_h2o_eos
    ${CMAKE_SOURCE_DIR}/thirdparty/meschach/source
    ${CMAKE_SOURCE_DIR}/thirdparty/jpeg
    ${CMAKE_SOURCE_DIR}/thirdparty/samg
    ${CMAKE_SOURCE_DIR}/examples/source
)


# Install targets

# Header files
FILE(GLOB_RECURSE PUBLIC_INCLUDE_FILES
	RELATIVE "${INCLUDES_DIR}"
	"${INCLUDES_DIR}/*"
	)
foreach ( I ${PUBLIC_INCLUDE_FILES} )
    if (IS_DIRECTORY "${I}")
	message(STATUS "Examining include directory ${I}")
    elseif("${I}" MATCHES ".h$")
        get_filename_component( dir ${I} DIRECTORY )
	install( FILES "${INCLUDES_DIR}/${I}"
		DESTINATION "include/${dir}")
    else()
	message(WARNING "Excluding file ${I} from installation")
    endif()
endforeach()


# Libraries

install(TARGETS csmpcore ARCHIVE DESTINATION lib)

