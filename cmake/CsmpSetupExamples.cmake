# SPDX-FileCopyrightText: © 2026 The openCSMP project
#
# SPDX-License-Identifier: LGPL-3.0-only

# Sets up the directory in which the openCSMP examples are run (the "examples workspace").
#
# Run through the build target
#     cmake --build <build dir> --target setup_examples
# or directly
#     cmake -DCSMP_EXAMPLE_INPUTS=<repository>/examples/example_inputs
#           -DCSMP_EXAMPLES_WORKSPACE=<workspace> [-DCSMP_REPOSITORY=<repository>]
#           -P cmake/CsmpSetupExamples.cmake
#
# Copies the example inputs to <workspace>/example_inputs. Files that already exist there
# are left untouched, so that inputs edited by the user and models created with the
# Interfaces example survive a repeated setup.

if(NOT CSMP_EXAMPLE_INPUTS OR NOT CSMP_EXAMPLES_WORKSPACE)
    message(FATAL_ERROR "CsmpSetupExamples.cmake: CSMP_EXAMPLE_INPUTS and CSMP_EXAMPLES_WORKSPACE must be set")
endif()

# a leading ~ is not expanded when the path comes from a CMake cache entry
if(CSMP_EXAMPLES_WORKSPACE MATCHES "^~" AND DEFINED ENV{HOME})
    string(REGEX REPLACE "^~" "$ENV{HOME}" CSMP_EXAMPLES_WORKSPACE "${CSMP_EXAMPLES_WORKSPACE}")
endif()
get_filename_component(_inputs    "${CSMP_EXAMPLE_INPUTS}"     ABSOLUTE)
get_filename_component(_workspace "${CSMP_EXAMPLES_WORKSPACE}" ABSOLUTE)

if(NOT IS_DIRECTORY "${_inputs}")
    message(FATAL_ERROR "Example inputs not found: ${_inputs}")
endif()

# simulation output does not belong into the repository
if(CSMP_REPOSITORY)
    get_filename_component(_repository "${CSMP_REPOSITORY}" ABSOLUTE)
    string(FIND "${_workspace}/" "${_repository}/" _position)
    if(_position EQUAL 0)
        message(FATAL_ERROR
                "The examples workspace must be located outside of the repository.\n"
                "  workspace:  ${_workspace}\n"
                "  repository: ${_repository}\n"
                "Configure with -DCSMP_EXAMPLES_WORKSPACE=<directory outside of the repository>")
    endif()
endif()

set(_destination "${_workspace}/example_inputs")
file(MAKE_DIRECTORY "${_destination}")
file(MAKE_DIRECTORY "${_workspace}/example_outputs")

file(GLOB_RECURSE _files RELATIVE "${_inputs}" "${_inputs}/*")
set(_copied 0)
set(_kept   0)
foreach(_file IN LISTS _files)
    if(_file MATCHES "(^|/)\\.DS_Store$" OR _file MATCHES "(^|/)__MACOSX/")
        continue()
    endif()
    if(EXISTS "${_destination}/${_file}")
        math(EXPR _kept "${_kept} + 1")
    else()
        get_filename_component(_subdirectory "${_file}" DIRECTORY)
        file(COPY "${_inputs}/${_file}" DESTINATION "${_destination}/${_subdirectory}")
        math(EXPR _copied "${_copied} + 1")
    endif()
endforeach()

message("")
message("openCSMP examples workspace: ${_workspace}")
message("  inputs:  ${_destination}")
message("           (${_copied} files copied, ${_kept} existing files left untouched)")
message("  outputs: ${_workspace}/example_outputs/<example name>/")
message("")
message("Run the examples with:   openCSMP-examples --workspace \"${_workspace}\"")
message("In an IDE, alternatively set the working directory to the inputs directory above.")
message("")
