#///////////////////////////////////////////////////////////////////////////
#//
#//    Copyright 2010
#//
#//    This file is part of rootpwa
#//
#//    rootpwa is free software: you can redistribute it and/or modify
#//    it under the terms of the GNU General Public License as published by
#//    the Free Software Foundation, either version 3 of the License, or
#//    (at your option) any later version.
#//
#//    rootpwa is distributed in the hope that it will be useful,
#//    but WITHOUT ANY WARRANTY; without even the implied warranty of
#//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#//    GNU General Public License for more details.
#//
#//    You should have received a copy of the GNU General Public License
#//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
#//
#///////////////////////////////////////////////////////////////////////////
#//-------------------------------------------------------------------------
#// File and Version Information:
#// $Rev::                             $: revision of last commit
#// $Author::                          $: author of last commit
#// $Date::                            $: date of last commit
#//
#// Description:
#//      main build file for ROOTPWA
#//
#//
#// Author List:
#//      Sebastian Neubert    TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


# check if cmake has the required version
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)


# set verbosity
set(CMAKE_VERBOSE_MAKEFILE 0)  # if set to 1 compile and link commands are displayed during build
# the same effect can be achieved by calling 'make VERBOSE=1'


# define project
project(ROOTPWA)


# load some common cmake macros
# set path, where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ is checked
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmakeModules")
message(STATUS "Using cmake module path '${CMAKE_MODULE_PATH}'")
include(CommonMacros)


# force out-of-source builds.
enforce_out_of_source_build()


# warn user if system is not UNIX
if(NOT UNIX)
  message(FATAL_ERROR "This is an unsupported system.")
endif()
message(STATUS "Detected host system '${CMAKE_HOST_SYSTEM_NAME}' version '${CMAKE_HOST_SYSTEM_VERSION}' architecture '${CMAKE_HOST_SYSTEM_PROCESSOR}'")
message(STATUS "Compiling for system '${CMAKE_SYSTEM_NAME}' version '${CMAKE_SYSTEM_VERSION}' architecture '${CMAKE_SYSTEM_PROCESSOR}'")


# define build types
# set a default build type for single-configuration CMake generators, if no build type is set.
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type was specified. Setting build type to 'Release'.")
  set(CMAKE_BUILD_TYPE Release)
endif()
# common compiler flags
#set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -pedantic")
#set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror")
set(CMAKE_CXX_FLAGS "-Wall -Werror")
# flags for specific build types
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_LDFLAGS_DEBUG "-g")
# report global build settings
foreach(_LANG "C" "CXX")
  message(STATUS "Using ${_LANG} compiler '${CMAKE_${_LANG}_COMPILER}'")
  message(STATUS "Using ${_LANG} compiler flags '${CMAKE_${_LANG}_FLAGS}'")
endforeach()
foreach(_BUILD_TYPE "DEBUG" "MINSIZEREL" "RELEASE" "RELWITHDEBINFO")
  message(STATUS "Using CXX compiler flags '${CMAKE_CXX_FLAGS_${_BUILD_TYPE}}' for build type ${_BUILD_TYPE}")
endforeach()
message(STATUS "Build type is '${CMAKE_BUILD_TYPE}'")


# redirect output files
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
message(STATUS "Using library output path '${LIBRARY_OUTPUT_PATH}'")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
message(STATUS "Using executable output path '${EXECUTABLE_OUTPUT_PATH}'")


# make CMAKE_SOURCE_DIR accessible via predefined macro CMAKE_SOURCE_DIR
if(CMAKE_SOURCE_DIR)
  add_definitions(-D'CMAKE_SOURCE_DIR=\"${CMAKE_SOURCE_DIR}\"')
else()
  add_definitions(-D'CMAKE_SOURCE_DIR=\"\"')
endif()


# setup Boost
# environment variable $BOOST_ROOT is expected to point to non-standard locations
# set(Boost_DEBUG 1)
find_package(Boost 1.42.0 REQUIRED)
if(Boost_FOUND)
  message(STATUS "Using Boost include directory '${Boost_INCLUDE_DIR}'")
  message(STATUS "Using Boost library directory '${Boost_LIBRARY_DIRS}'")
else()
  message(SEND_ERROR "Could not find Boost installation. Is environment variable BOOST_ROOT=${BOOST_ROOT} set correctly?")
endif()


# setup libconfig
# environment variable $LIBCONFIG is expected to point to libconfig location
set(LIBCONFIG $ENV{LIBCONFIG})
if(LIBCONFIG)
  message(STATUS "Using libconfig in directory '${LIBCONFIG}'")
  set(LIBCONFIG_INCLUDE_DIR "${LIBCONFIG}/include")
  set(LIBCONFIG_LIBRARY_DIR "${LIBCONFIG}/lib")
  set(_LIBCONFIG_LIB_NAMES "config++")
  find_library(LIBCONFIG_LIBS
    NAMES ${_LIBCONFIG_LIB_NAMES}
    PATHS ${LIBCONFIG_LIBRARY_DIR}
    NO_DEFAULT_PATH)
  if(NOT LIBCONFIG_LIBS)
    message(SEND_ERROR "Cannot find libconfig library '${_LIBCONFIG_LIB_NAMES}' in '${LIBCONFIG_LIBRARY_DIR}'")
  else()
    message(STATUS "Using libconfig library '${LIBCONFIG_LIBS}'")
  endif()
else()
  message(SEND_ERROR "Environment variable LIBCONFIG='${LIBCONFIG}' is not set correctly")
endif()


# setup ROOT includes and libraries
find_package(ROOT 5.26 REQUIRED Minuit2)


# setup CUDA
set(USE_CUDA FALSE)
#set(CUDA_SHADER_MODELS "13" "20")
set(CUDA_SHADER_MODELS "13")
if(NOT DEFINED ENV{NVSDKCUDA_ROOT})
  # use default SDK path, if not defined by user
  set(ENV{NVSDKCUDA_ROOT} "$ENV{HOME}/NVIDIA_GPU_Computing_SDK/C")
endif()
find_package(CUDA 3.0)
if(NOT CUDA_FOUND)
  message(STATUS "Cannot find CUDA installation. No CUDA components will be generated.")
  set(USE_CUDA FALSE)
else()
  message(STATUS "Found CUDA version ${CUDA_VERSION}")
  message(STATUS "Using CUDA toolkit root directory '${CUDA_TOOLKIT_ROOT_DIR}'")
  message(STATUS "Using CUDA include directories '${CUDA_INCLUDE_DIRS}'")
  message(STATUS "Using CUDA libraries '${CUDA_LIBRARIES}'")
  if(NOT CUDA_SDK_ROOT_DIR)
    message(SEND_ERROR "Could not find CUDA SDK directory in '$ENV{NVSDKCUDA_ROOT}'. Is environment variable $NVSDKCUDA_ROOT set correctly.")
  else()
		# find cuda utility library
		find_library(CUDA_CUTIL_LIBRARY
			NAMES cutil_${CMAKE_SYSTEM_PROCESSOR}
			PATHS ${CUDA_SDK_ROOT_DIR}/lib
			NO_DEFAULT_PATH
			)
		if(NOT CUDA_CUTIL_LIBRARY)
			message(SEND_ERROR "Could not find CUDA utility library 'cutil_${CMAKE_SYSTEM_PROCESSOR}' in '${CUDA_SDK_ROOT_DIR}/lib'")
		else()
			set(USE_CUDA TRUE)
			message(STATUS "Using CUDA utility library '${CUDA_CUTIL_LIBRARY}'")
			message(STATUS "Using CUDA SDK root directory '${CUDA_SDK_ROOT_DIR}'")
			message(STATUS "Enabling compilation of CUDA components")
			# define for which shader models code is generated
			message(STATUS "Generating binary files for shader model(s): ${CUDA_SHADER_MODELS}")
			foreach(_SHADER_MODEL ${CUDA_SHADER_MODELS})
				list(APPEND CUDA_NVCC_FLAGS
					"--generate-code=arch=compute_${_SHADER_MODEL},code=\\\"sm_${_SHADER_MODEL},compute_${_SHADER_MODEL}\\\"")
			endforeach()
			#list(APPEND CUDA_NVCC_FLAGS "--ptxas-options=-v") # prints kernel's resource usage
			# set some global switches
			set(CUDA_BUILD_CUBIN     OFF)
			set(CUDA_BUILD_EMULATION OFF)
			set(CUDA_VERBOSE_BUILD   OFF)
			# enable CUDA code in other parts of the project
			add_definitions(-DUSE_CUDA)
		endif()
  endif()
endif()


# make SVN version string accessible via predefined macro SVN_VERSION
find_package(Subversion)
if(Subversion_FOUND)
  # unfortunately CMAKE only parses 'svn info'
  find_program(SVNVERSION_EXECUTABLE
    svnversion
		)
  if(NOT SVNVERSION_EXECUTABLE)
    message(STATUS "Could not find subversion command 'svnversion'. Repository version unknown.")
  else()
    execute_process(
      COMMAND ${SVNVERSION_EXECUTABLE} "${CMAKE_SOURCE_DIR}"
      OUTPUT_VARIABLE SVN_VERSION
      RESULT_VARIABLE _SVNVERSION_RETURN
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT ${_SVNVERSION_RETURN})
      message(STATUS "Subversion repository revision is '${SVN_VERSION}'")
    else()
      message(STATUS "Error running 'svnversion'. Repository version unknown.")
      set(SVN_VERSION "")
    endif()
  endif()
else()
  message(STATUS "Could not find subversion installation. Repository version unknown.")
endif()
if(SVN_VERSION)
  add_definitions(-D'SVN_VERSION=\"${SVN_VERSION}\"')
else()
  add_definitions(-D'SVN_VERSION=\"\"')
endif()


# setup doxygen
find_package(Doxygen)
if(NOT DOXYGEN_FOUND)
  message(WARNING "Cannot find Doxygen. No documentation will be generated.")
else()
  set(DOXYGEN_TARGET  "doxygen")
  set(DOXYGEN_DOC_DIR "${CMAKE_SOURCE_DIR}/html-doc")
  set(DOXYGEN_CONF    "${CMAKE_SOURCE_DIR}/rootpwaDoxyfile.conf")
  message(STATUS "Run 'make ${DOXYGEN_TARGET}' to create Doxygen documentation files in '${DOXYGEN_DOC_DIR}'")
  add_custom_target(${DOXYGEN_TARGET}
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_CONF}
    DEPENDS ${DOXYGEN_CONF}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		)
endif()


enable_testing()


# build subprojects

# PWA2000
add_subdirectory(pwa2000)

# Utilities
set(UTILITIES_INCLUDE_DIR
  ${CMAKE_SOURCE_DIR}/utilities
	)
add_subdirectory(utilities)

# CUDA stuff
if(USE_CUDA)
	set(CUDA_INCLUDE_DIR
		${CMAKE_SOURCE_DIR}/cuda
		)
	set(CUDA_LIBS
		cuda
		)
	add_subdirectory(cuda)
else()
	set(CUDA_INCLUDE_DIR)
	set(CUDA_LIBS)
endif()

# RootPwa
set(ROOTPWA_INCLUDE_DIR
  ${CMAKE_SOURCE_DIR}/src
	)
set(ROOTPWA_LIBS
  RootPwa
	)
add_subdirectory(src)

# RootPwaAmp
# add_subdirectory(amplitude)

# RootPwaTools
add_subdirectory(tools)

# RootPwaGen
add_subdirectory(generators)

# RootPwaKey
add_subdirectory(keygen)

# Fortran
add_subdirectory(fortran)



# file(GLOB _CC_LIST ${CMAKE_CURRENT_SOURCE_DIR} "amplitude/*.cc")
# #message("foo ${_CC_LIST}")
# switch_file_extension("${_CC_LIST}" ".cc" ".h" _H_LIST)
# #message("bar ${_H_LIST}")
# root_generate_dictionary("testDict.cc" "${ROOT_INCLUDE_DIR}" "${_H_LIST}" "testLinkDef.h")
