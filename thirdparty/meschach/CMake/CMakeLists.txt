####################################################
## CMAKE Makefile builder
## CSMP++
## Created by: Julian Mindel
## Creation Date: 17.07.2013
## type:
## ccmake . (unix)
## cmake-gui (linux , windows)
#####################################################
cmake_minimum_required(VERSION 2.8)
project(meschach)

INCLUDE_DIRECTORIES(../source)

FILE(GLOB SOURCE_FILES "../source/*.cpp")
FILE(GLOB HEADER_FILES "../source/*.h")
FILE(GLOB HPP_FILES "../source/*.hpp")
SET(SOURCE ${SOURCE} ${SOURCE_FILES} ${HEADER_FILES} ${HPP_FILES})

IF("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
   IF(WIN32)
      SET(CMAKE_EXECUTABLE_SUFFIX "-debug.exe") #for executables in windows.
   ELSE(WIN32)
      SET(CMAKE_EXECUTABLE_SUFFIX "-debug") #for executables in linux/unix
   ENDIF(WIN32)

   SET(CMAKE_DEBUG_POSTFIX "-debug") # For libraries
ENDIF("${CMAKE_BUILD_TYPE}" MATCHES "Debug")

IF(${CMAKE_GENERATOR} MATCHES "Visual Studio")
   SET(BUILD_SHARED_LIBS OFF CACHE BOOL OFF)
ELSE()
   SET(BUILD_SHARED_LIBS ON CACHE BOOL OFF)
ENDIF(${CMAKE_GENERATOR} MATCHES "Visual Studio")

IF(BUILD_SHARED_LIBS)
  message("Building MESCHACH Shared Libraries")
ENDIF()

SET(LIBRARY_OUTPUT_PATH "${CMAKE_CURRENT_LIST_DIR}/../../../lib")
SET(EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_LIST_DIR}/../../../bin")

IF("${CMAKE_C_COMPILER}" MATCHES "clang")
  SET(COMPILER_NAME "clang++")
ELSEIF("${CMAKE_C_COMPILER}" MATCHES "gcc")
  SET(COMPILER_NAME "g++")
ELSEIF("${CMAKE_C_COMPILER}" MATCHES "icc")
  ADD_DEFINITIONS(-D__INTEL_COMPILER)
  SET(COMPILER_NAME "icpc")
ELSEIF("${CMAKE_C_COMPILER}" MATCHES "cl.exe")
  SET(COMPILER_NAME "VS2012")
ENDIF("${CMAKE_C_COMPILER}" MATCHES "clang")

IF(CMAKE_GENERATOR MATCHES "Visual Studio 11")
  SET(COMPILER_NAME "VS2012")
ENDIF(CMAKE_GENERATOR MATCHES "Visual Studio 11")

IF(CMAKE_GENERATOR MATCHES "Visual Studio 12")
  SET(COMPILER_NAME "VS2013")
ENDIF(CMAKE_GENERATOR MATCHES "Visual Studio 12")

ADD_LIBRARY("${PROJECT_NAME}-${COMPILER_NAME}" ${SOURCE})

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY LINK_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "meschach link directory : '${dir}'")
endforeach(dir ${dirs})

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "meschach include directory : '${dir}'")
endforeach(dir ${dirs})
