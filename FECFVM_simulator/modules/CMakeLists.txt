
SET(INCLUDES_DIR "${PROJECT_SOURCE_DIR}/includes")
SET(MODULES_DIR "${PROJECT_SOURCE_DIR}/modules")

file(GLOB_RECURSE CSMP_FECFVM_SIMULATOR_SOURCE_FILES "${MODULES_DIR}/*.cpp")
file(GLOB_RECURSE CSMP_FECFVM_SIMULATOR_HEADER_FILES "${INCLUDES_DIR}/*.h")


macro(GroupSources curdir)
  file(GLOB children RELATIVE ${PROJECT_SOURCE_DIR}/modules/${curdir} ${PROJECT_SOURCE_DIR}/modules/${curdir}/*)
  foreach(child ${children})
    if(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/modules/${curdir}/${child})
      GroupSources(${curdir}/${child})
    else()
      if (CMAKE_HOST_APPLE)
	string(REGEX REPLACE "^modules/" "Source Files/" groupname ${curdir})
        string(REPLACE "/" "\\\\" groupname ${groupname})
      elseif(CMAKE_HOST_WIN32)
	string(REGEX REPLACE "^modules/" "Source/" groupname ${curdir})
        string(REPLACE "/" "\\" groupname ${groupname})
      endif()

      source_group("${groupname}" FILES ${PROJECT_SOURCE_DIR}/modules/${curdir}/${child})
    endif()
  endforeach()
endmacro()

GroupSources(modules)

macro(GroupHeaders curdir)
  file(GLOB children RELATIVE ${PROJECT_SOURCE_DIR}/includes/${curdir} ${PROJECT_SOURCE_DIR}/includes/${curdir}/*)
  foreach(child ${children})
    if(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/includes/${curdir}/${child})
      GroupHeaders(${curdir}/${child})
    else()
      if (CMAKE_HOST_APPLE)
	string(REGEX REPLACE "^includes/" "Header Files/" groupname ${curdir})
        string(REPLACE "/" "\\\\" groupname ${groupname})
      elseif(CMAKE_HOST_WIN32)
	string(REGEX REPLACE "^includes/" "Include/" groupname ${curdir})
        string(REPLACE "/" "\\" groupname ${groupname})
      endif()

      source_group("${groupname}" FILES ${PROJECT_SOURCE_DIR}/includes/${curdir}/${child})
    endif()
  endforeach()
endmacro()

GroupHeaders(includes)

ADD_EXECUTABLE(csmp_FECFVM_simulator
    ${CSMP_FECFVM_SIMULATOR_SOURCE_FILES}
    ${CSMP_FECFVM_SIMULATOR_HEADER_FILES}
    ${PROJECT_SOURCE_DIR}/IMPES_main.cpp
)

TARGET_LINK_LIBRARIES(csmp_FECFVM_simulator ${CSMP_LIBRARIES} ${SAMG_LIBRARIES})

set(CSMP_INCLUDE_DIR "${CSMP_SOURCE_DIRECTORY}/source/includes")
TARGET_INCLUDE_DIRECTORIES(csmp_FECFVM_simulator PRIVATE    
    ${INCLUDES_DIR}
    ${CSMP_GENERATED_INCLUDES_DIR}
    ${CSMP_SOURCE_DIRECTORY}/thirdparty/iaps_h2o_eos
    ${CSMP_SOURCE_DIRECTORY}/thirdparty/samg
    ${CSMP_SOURCE_DIRECTORY}/tests/unit
    ${CSMP_SOURCE_DIRECTORY}/tests/integration
    ${CSMP_INCLUDE_DIR}/analysis
    ${CSMP_INCLUDE_DIR}/thmc/fractures
    ${CSMP_INCLUDE_DIR}/thmc/geochemistry
    ${CSMP_INCLUDE_DIR}/thmc/mechanics
    ${CSMP_INCLUDE_DIR}/thmc/single_phase_flow
    ${CSMP_INCLUDE_DIR}/thmc/two_phase_flow
    ${CSMP_INCLUDE_DIR}/interfaces/ANSYS
    ${CSMP_INCLUDE_DIR}/interfaces/CGNS
    ${CSMP_INCLUDE_DIR}/interfaces/Eclipse    
    ${CSMP_INCLUDE_DIR}/interfaces/Fracman
    ${CSMP_INCLUDE_DIR}/interfaces/GeoModeller
    ${CSMP_INCLUDE_DIR}/interfaces/gOcad_SKUA
    ${CSMP_INCLUDE_DIR}/interfaces/JPEG
    ${CSMP_INCLUDE_DIR}/interfaces/Maple
    ${CSMP_INCLUDE_DIR}/interfaces/Matlab
    ${CSMP_INCLUDE_DIR}/interfaces/Rhinoceros
    ${CSMP_INCLUDE_DIR}/interfaces/Tecplot
    ${CSMP_INCLUDE_DIR}/interfaces/VTK
    ${CSMP_INCLUDE_DIR}/integration/finite_elements
    ${CSMP_INCLUDE_DIR}/integration/finite_volumes
    ${CSMP_INCLUDE_DIR}/integration/finite_volumes/generic_transport_scheme
    ${CSMP_INCLUDE_DIR}/integration/finite_volumes/legacy_transport_scheme
    ${CSMP_INCLUDE_DIR}/integration/finite_volumes/DES_transport_scheme
    ${CSMP_INCLUDE_DIR}/integration/pde_integrators
    ${CSMP_INCLUDE_DIR}/integration/pde_operators
    ${CSMP_INCLUDE_DIR}/integration/pde_operators/analytical
    ${CSMP_INCLUDE_DIR}/integration/pde_operators/numerical
    ${CSMP_INCLUDE_DIR}/integration/solver
    ${CSMP_INCLUDE_DIR}/discretisation
    ${CSMP_INCLUDE_DIR}/eos
    ${CSMP_INCLUDE_DIR}/integration
    ${CSMP_INCLUDE_DIR}/interfaces
    ${CSMP_INCLUDE_DIR}/math
    ${CSMP_INCLUDE_DIR}/model
    ${CSMP_INCLUDE_DIR}/property_manipulation
    ${CSMP_INCLUDE_DIR}/regular_grids
    ${CSMP_INCLUDE_DIR}/run_control
    ${CSMP_INCLUDE_DIR}/thmc
    ${CSMP_INCLUDE_DIR}/utilities
)


# Install targets


# Libraries

#install(TARGETS csmp_FECFVM_simulator ARCHIVE DESTINATION bin)

