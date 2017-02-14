
SET(INCLUDES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/includes")
SET(MODULES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/modules")
SET(GEN_INCLUDES_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated_includes")

FILE(MAKE_DIRECTORY ${GEN_INCLUDES_DIR})

ADD_EXECUTABLE(CSMP_number_types_main
	generate_number_types.cpp
)

ADD_CUSTOM_COMMAND(
    OUTPUT "${GEN_INCLUDES_DIR}/CSMP_number_types.h"
    COMMAND $<TARGET_FILE:CSMP_number_types_main>
    MAIN_DEPENDENCY CSMP_number_types_main
    WORKING_DIRECTORY "${GEN_INCLUDES_DIR}"
)

ADD_CUSTOM_TARGET(number_types_header DEPENDS "${GEN_INCLUDES_DIR}/CSMP_number_types.h")


FILE(GLOB ANALYSIS_SRC ${MODULES_DIR}/analysis/*.cpp )
FILE(GLOB DISCRETISATION_SRC ${MODULES_DIR}/discretisation/*.cpp )
FILE(GLOB EOS_SRC ${MODULES_DIR}/eos/*.cpp )
FILE(GLOB EOS_H2O_SRC ${MODULES_DIR}/eos/h2o/*.cpp )
FILE(GLOB EOS_H2O_NACL_SRC ${MODULES_DIR}/eos/h2o_nacl/*.cpp )
FILE(GLOB EOS_H2O_NACL_LOOKUP_SRC ${MODULES_DIR}/eos/h2o_nacl_lookup/*.cpp )
FILE(GLOB INTEGRATION_CVFE_HYDROTHERMAL_SRC ${MODULES_DIR}/integration/cvfe_hydrothermal/*.cpp )
FILE(GLOB INTEGRATION_FINITE_ELEMENTS_SRC ${MODULES_DIR}/integration/finite_elements/*.cpp )
FILE(GLOB INTEGRATION_FINITE_VOLUMES_SRC ${MODULES_DIR}/integration/finite_volumes/*.cpp )
FILE(GLOB INTEGRATION_FINITE_VOLUMES_GENERIC_SRC ${MODULES_DIR}/integration/finite_volumes/generic_transport_scheme/source/*.cpp )
FILE(GLOB INTEGRATION_FINITE_VOLUMES_LEGACY_SRC ${MODULES_DIR}/integration/finite_volumes/legacy_transport_scheme/*.cpp )
FILE(GLOB INTEGRATION_PDE_INTEGRATORS_SRC ${MODULES_DIR}/integration/pde_integrators/*.cpp )
FILE(GLOB INTEGRATION_PDE_OPERATORS_SRC ${MODULES_DIR}/integration/pde_operators/*.cpp )
FILE(GLOB INTEGRATION_PDE_OPERATORS_ANALYTIC_SRC ${MODULES_DIR}/integration/pde_operators/analytical_integration/*.cpp )
FILE(GLOB INTEGRATION_PDE_OPERATORS_FE_SRC ${MODULES_DIR}/integration/pde_operators/fe_numerical_integration/*.cpp )
FILE(GLOB INTEGRATION_PDE_OPERATORS_FV_SRC ${MODULES_DIR}/integration/pde_operators/fv_integration/*.cpp )
FILE(GLOB INTEGRATION_SOLVER_SRC ${MODULES_DIR}/integration/solver/*.cpp )
FILE(GLOB INTEGRATION_SOLVER_SAMG_SRC ${MODULES_DIR}/integration/solver/SAMG_profiles/*.cpp )
FILE(GLOB INTERFACES_SRC ${MODULES_DIR}/interfaces/*.cpp )
FILE(GLOB INTERFACES_ANSYS_SRC ${MODULES_DIR}/interfaces/ANSYS/*.cpp )
FILE(GLOB CGNS_SRC ${MODULES_DIR}/interfaces/CGNS/*.cpp )
    FILE(GLOB CGNS_SRC_EXCLUSIONS ${MODULES_DIR}/interfaces/CGNS/CGNS_*.cpp )
    LIST(REMOVE_ITEM CGNS_SRC ${CGNS_SRC_EXCLUSIONS})
FILE(GLOB ECLIPSE_SRC ${MODULES_DIR}/interfaces/Eclipse/*.cpp )
FILE(GLOB FRACMAN_SRC ${MODULES_DIR}/interfaces/Fracman/*.cpp )
FILE(GLOB GEOMODELLER_SRC ${MODULES_DIR}/interfaces/GeoModeller/*.cpp )
FILE(GLOB GOCAD_SKUA_SRC ${MODULES_DIR}/interfaces/gOcad_SKUA/*.cpp )
FILE(GLOB JPEG_SRC ${MODULES_DIR}/interfaces/JPEG/*.cpp )
FILE(GLOB MAPLE_SRC ${MODULES_DIR}/interfaces/Maple/*.cpp )
FILE(GLOB MATLAB_SRC ${MODULES_DIR}/interfaces/Matlab/*.cpp )
FILE(GLOB RHINOCEROS_SRC ${MODULES_DIR}/interfaces/Rhinoceros/*.cpp )
FILE(GLOB TECPLOT_SRC ${MODULES_DIR}/interfaces/Tecplot/*.cpp )
FILE(GLOB VTK_SRC ${MODULES_DIR}/interfaces/VTK/*.cpp )
FILE(GLOB MATH_SRC ${MODULES_DIR}/math/*.cpp )
FILE(GLOB MODEL_SRC ${MODULES_DIR}/model/*.cpp )
FILE(GLOB PROPERTY_MANIPULATION_SRC ${MODULES_DIR}/property_manipulation/*.cpp )
FILE(GLOB PROPERTY_MANIPULATION_INTER_SRC ${MODULES_DIR}/property_manipulation/interrelations/*.cpp )
FILE(GLOB REGULAR_GRIDS_SRC ${MODULES_DIR}/regular_grids/*.cpp )
FILE(GLOB RUN_CONTROL_SRC ${MODULES_DIR}/run_control/*.cpp : $(MODULE_BASE)/run_control/ExampleSuiteMainWindow.cpp )
    FILE(GLOB RUN_CONTROL_SRC_EXCLUSIONS ${MODULES_DIR}/run_control/ExampleSuiteMainWindow.cpp )
    LIST(REMOVE_ITEM RUN_CONTROL_SRC ${RUN_CONTROL_SRC_EXCLUSIONS})
FILE(GLOB THMC_FRACTURES_SRC ${MODULES_DIR}/thmc/fractures/*.cpp )
FILE(GLOB THMC_GEOCHEMISTRY_SRC ${MODULES_DIR}/thmc/geochemistry/*.cpp )
FILE(GLOB THMC_MECHANICS_SRC ${MODULES_DIR}/thmc/mechanics/*.cpp )
FILE(GLOB THMC_SINGLE_PHASE_FLOW_SRC ${MODULES_DIR}/thmc/single_phase_flow/*.cpp )
FILE(GLOB THMC_TWO_PHASE_FLOW_SRC ${MODULES_DIR}/thmc/two_phase_flow/*.cpp )
FILE(GLOB UTILITIES_SRC ${MODULES_DIR}/utilities/*.cpp )

ADD_LIBRARY(csmpcore
	${ANALYSIS_SRC}
	${DISCRETISATION_SRC}
	${EOS_SRC}
	${EOS_H2O_SRC}
	${EOS_H2O_NACL_SRC}
	${EOS_H2O_NACL_LOOKUP_SRC}
	${INTEGRATION_CVFE_HYDROTHERMAL_SRC}
	${INTEGRATION_FINITE_ELEMENTS_SRC}
	${INTEGRATION_FINITE_VOLUMES_SRC}
	${INTEGRATION_FINITE_VOLUMES_GENERIC_SRC}
	${INTEGRATION_FINITE_VOLUMES_LEGACY_SRC}
	${INTEGRATION_PDE_INTEGRATORS_SRC}
	${INTEGRATION_PDE_OPERATORS_SRC}
	${INTEGRATION_PDE_OPERATORS_ANALYTIC_SRC}
	${INTEGRATION_PDE_OPERATORS_FE_SRC}
	${INTEGRATION_PDE_OPERATORS_FV_SRC}
	${INTEGRATION_SOLVER_SRC}
	${INTEGRATION_SOLVER_SAMG_SRC}
	${INTERFACES_SRC}
	${INTERFACES_ANSYS_SRC}
	${CGNS_SRC}
	${ECLIPSE_SRC}
	${FRACMAN_SRC}
	${GEOMODELLER_SRC}
	${GOCAD_SKUA_SRC}
	${JPEG_SRC}
	${MAPLE_SRC}
	${MATLAB_SRC}
	${RHINOCEROS_SRC}
	${TECPLOT_SRC}
	${VTK_SRC}
	${MATH_SRC}
	${MODEL_SRC}
	${PROPERTY_MANIPULATION_SRC}
	${PROPERTY_MANIPULATION_INTER_SRC}
	${REGULAR_GRIDS_SRC}
	${RUN_CONTROL_SRC}
	${THMC_FRACTURES_SRC}
	${THMC_GEOCHEMISTRY_SRC}
	${THMC_MECHANICS_SRC}
	${THMC_SINGLE_PHASE_FLOW_SRC}
	${THMC_TWO_PHASE_FLOW_SRC}
	${UTILITIES_SRC}
)

ADD_DEPENDENCIES(csmpcore number_types_header)

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
    ${CMAKE_SOURCE_DIR}/lib/jpeg
    ${CMAKE_SOURCE_DIR}/lib/samg
    ${CMAKE_SOURCE_DIR}/lib/iaps_h2o_eos
    ${CMAKE_SOURCE_DIR}/examples/source
    ${CMAKE_SOURCE_DIR}/tests
)


# Install targets

# Header files
FILE(GLOB_RECURSE PUBLIC_INCLUDE_FILES
	RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/includes"
	"${CMAKE_CURRENT_SOURCE_DIR}/includes/*"
	)
foreach ( I ${PUBLIC_INCLUDE_FILES} )
    if (IS_DIRECTORY "${I}")
	message(STATUS "Examining include directory ${I}")
    elseif("${I}" MATCHES ".h$")
        get_filename_component( dir ${I} DIRECTORY )
	install( FILES "${CMAKE_CURRENT_SOURCE_DIR}/includes/${I}"
		DESTINATION "include/${dir}")
    else()
	message(WARNING "Excluding file ${I} from installation")
    endif()
endforeach()

INSTALL(FILES "${GEN_INCLUDES_DIR}/CSMP_number_types.h" DESTINATION include)

# Libraries

install(TARGETS csmpcore ARCHIVE DESTINATION lib)

