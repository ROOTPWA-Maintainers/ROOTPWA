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
include(FeatureSummary)


# force out-of-source builds.
enforce_out_of_source_build()


# warn user if system is not UNIX
if(NOT UNIX)
  message(FATAL_ERROR "This is an unsupported system.")
endif()
# requires Cmake V2.8.5
#include(ProcessorCount)
#ProcessorCount(NMB_CPU_CORES)
# message(STATUS "Detected host system '${CMAKE_HOST_SYSTEM_NAME}' version "
# 	"'${CMAKE_HOST_SYSTEM_VERSION}', architecture '${CMAKE_HOST_SYSTEM_PROCESSOR}', "
# 	"number of CPU cores ${NMB_CPU_CORES}")
message(STATUS "Detected host system '${CMAKE_HOST_SYSTEM_NAME}' version "
	"'${CMAKE_HOST_SYSTEM_VERSION}', architecture '${CMAKE_HOST_SYSTEM_PROCESSOR}'")
message(STATUS "Compiling for system '${CMAKE_SYSTEM_NAME}' version '${CMAKE_SYSTEM_VERSION}' "
	"architecture '${CMAKE_SYSTEM_PROCESSOR}'")


# define build types
# predefined build types: "DEBUG" "RELEASE" "RELWITHDEBINFO" "MINSIZEREL"
# set a default build type for single-configuration CMake generators, if no build type is set.
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type was specified. Setting build type to 'RELEASE'.")
  set(CMAKE_BUILD_TYPE RELEASE)
endif()
# common compiler flags
#set(CMAKE_CXX_FLAGS "-Wall -Woverloaded-virtual -Werror -Wextra -pedantic")
#set(CMAKE_CXX_FLAGS "-Wall -Woverloaded-virtual -Werror -Wextra")
set(CMAKE_CXX_FLAGS "-Wall -Woverloaded-virtual -Werror")
# flags for specific build types
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
# `gcc -march=native -Q --help=target` shows the relevant compiler switches
# consider adding -DNODEBUG (disables all asserts) and -ffast-math
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g ${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_LDFLAGS_DEBUG "-g")
set(CMAKE_CXX_LDFLAGS_RELWITHDEBINFO "-g")
# report global build settings
foreach(_LANG "C" "CXX")
  message(STATUS "Using ${_LANG} compiler '${CMAKE_${_LANG}_COMPILER}'")
  message(STATUS "Using ${_LANG} compiler flags '${CMAKE_${_LANG}_FLAGS}'")
endforeach()
message(STATUS "Build type is '${CMAKE_BUILD_TYPE}'")
message(STATUS "Using CXX compiler flags '${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}' "
	"for build type ${CMAKE_BUILD_TYPE}")
message(STATUS "Using linker flags '${CMAKE_CXX_LDFLAGS_${CMAKE_BUILD_TYPE}}' "
	"for build type ${CMAKE_BUILD_TYPE}")
#include(CMakePrintSystemInformation)


# redirect output files
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
message(STATUS "Using library output path '${LIBRARY_OUTPUT_PATH}'")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
message(STATUS "Using executable output path '${EXECUTABLE_OUTPUT_PATH}'")


# setup Boost
# environment variable $BOOST_ROOT is expected to point to non-standard locations
# set(Boost_DEBUG 1)
find_package(Boost 1.43.0 REQUIRED)
set_feature_info(Boost "C++ Template Library" "http://www.boost.org/" "required")
if(Boost_FOUND)
	set(Boost_LIBRARY_VERSION "${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
  message(STATUS "Using Boost ${Boost_LIBRARY_VERSION} include directories '${Boost_INCLUDE_DIRS}'")
  message(STATUS "Using Boost ${Boost_LIBRARY_VERSION} library directories '${Boost_LIBRARY_DIRS}'")
else()
  message(FATAL_ERROR "Could not find Boost installation. "
		"Is environment variable BOOST_ROOT=${BOOST_ROOT} set correctly? Please read INSTALL.")
endif()
# setup uprooted Boost libraries
set(uprootedBoost_INCLUDE_DIRS
  # local copy of Boost Graph Library from V1.43.0; overides graph library in Boost supplied by user
  ${CMAKE_SOURCE_DIR}/uprootedBoost/boost_1_43_0
  )
message(STATUS "using uprooted Boost include directories '${uprootedBoost_INCLUDE_DIRS}'")
set(origBoost_INCLUDE_DIRS
	${Boost_INCLUDE_DIRS}
	)
# use uprooted Boost libraries by default
set(Boost_INCLUDE_DIRS
	${uprootedBoost_INCLUDE_DIRS}
	${origBoost_INCLUDE_DIRS}
	)


# setup libconfig
# environment variable $LIBCONFIG is expected to point to libconfig location
find_package(Libconfig 1.4 REQUIRED)
set_feature_info(Libconfig "simple C++ library for processing structured configuration files"
	"http://www.hyperrealm.com/libconfig/" "required")


# setup ROOT includes and libraries
find_package(ROOT 5.26 REQUIRED Minuit2 MathMore TreePlayer)
set_feature_info(ROOT "C++ HEP analysis framework" "http://root.cern.ch/" "required")
# some code switches that depend on ROOT version
# std::complex is not supported as TTree leaf in ROOT versions below 5.27.06
if(NOT ROOT_VERSION VERSION_LESS "5.27.06")
	set(USE_STD_COMPLEX_TREE_LEAFS TRUE)
	add_definitions(-DUSE_STD_COMPLEX_TREE_LEAFS)
endif()
# TFitResult produces a name clash for ROOT versions from 5.25.0 on
if(ROOT_VERSION VERSION_LESS "5.25.00")
	set(USE_TFITRESULT TRUE)
	add_definitions(-DUSE_TFITRESULT)
endif()


# setup CUDA
set(CUDA_SHADER_MODELS "13" "20")
if(NOT DEFINED ENV{NVSDKCUDA_ROOT})
  # use default SDK path, if not defined by user
  set(ENV{NVSDKCUDA_ROOT} "$ENV{HOME}/NVIDIA_GPU_Computing_SDK/C")
endif()
find_package(CUDA 4.0)
set_feature_info(CUDA "nvidia GPGPU framework" "http://developer.nvidia.com/category/zone/cuda-zone"
	"optional")
set(USE_CUDA FALSE)
if(NOT CUDA_FOUND)
  set(USE_CUDA FALSE)
  message(STATUS "Cannot find CUDA installation. No CUDA components will be generated.")
else()
  message(STATUS "Found CUDA version ${CUDA_VERSION}")
  message(STATUS "Using CUDA toolkit root directory '${CUDA_TOOLKIT_ROOT_DIR}'")
  message(STATUS "Using CUDA include directories '${CUDA_INCLUDE_DIRS}'")
  message(STATUS "Using CUDA libraries '${CUDA_LIBRARIES}'")
  if(NOT CUDA_SDK_ROOT_DIR)
    set(USE_CUDA FALSE)
    message(STATUS "Found CUDA intallation, but could not find CUDA SDK in '$ENV{NVSDKCUDA_ROOT}'. "
			"Check that $NVSDKCUDA_ROOT is set correctly. No CUDA components will be generated.")
  else()
    # find cuda utility library
    set(CUDA_CUTIL_LIB_SUFFIX "i386")
    if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
      set(CUDA_CUTIL_LIB_SUFFIX ${CMAKE_SYSTEM_PROCESSOR})
    endif()
    find_library(CUDA_CUTIL_LIBRARY
      NAMES cutil_${CUDA_CUTIL_LIB_SUFFIX}
      PATHS ${CUDA_SDK_ROOT_DIR}/lib
      NO_DEFAULT_PATH
      )
    if(NOT CUDA_CUTIL_LIBRARY)
      set(USE_CUDA FALSE)
      message(STATUS "Found CUDA intallation, but could not find CUDA utility library "
				"'cutil_${CUDA_CUTIL_LIB_SUFFIX}' in '${CUDA_SDK_ROOT_DIR}/lib'. "
				"No CUDA components will be generated.")
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
			set(CUDA_LIB_DIRS "${CUDA_TOOLKIT_ROOT_DIR}, ${CUDA_SDK_ROOT_DIR}")
    endif()
  endif()
endif()
if(CUDA_FOUND AND NOT USE_CUDA)
	disable_feature(CUDA)
endif()


# setup MPI
find_package(MPI)
set_feature_info(MPI "standardized and portable message-passing system for parallel computing"
	"http://www.boost.org/doc/libs/1_47_0/doc/html/mpi.html" "optional")
set(USE_MPI FALSE)
if(NOT MPI_FOUND)
  set(USE_MPI FALSE)
  message(STATUS "Cannot find MPI installation. No MPI components will be generated.")
else()
  set(USE_MPI TRUE)
  message(STATUS "Using MPI include directories '${MPI_INCLUDE_PATH}'")
  message(STATUS "Using MPI libraries '${MPI_LIBRARIES}'")
  # setup Boost.MPI library
  set(Boost_MPI_LIBRARY_DIR "${BOOST_ROOT}/stage/lib")
  set(_Boost_MPI_LIB_NAMES "boost_mpi" "boost_serialization")
  foreach(_Boost_MPI_LIB_NAME ${_Boost_MPI_LIB_NAMES})
    find_library(_Boost_MPI_LIB_${_Boost_MPI_LIB_NAME}
      NAMES ${_Boost_MPI_LIB_NAME}
      PATHS ${Boost_MPI_LIBRARY_DIR}
      NO_DEFAULT_PATH)
    if(NOT _Boost_MPI_LIB_${_Boost_MPI_LIB_NAME})
      set(USE_MPI FALSE)
      message(STATUS "Cannot find Boost.MPI library '${_Boost_MPI_LIB_NAME}' in "
				"'${Boost_MPI_LIBRARY_DIR}'. Please read INSTALL on how to build Boost.MPI.")
    else()
      list(APPEND Boost_MPI_LIBS ${_Boost_MPI_LIB_${_Boost_MPI_LIB_NAME}})
    endif()
  endforeach()
  if(USE_MPI)
    message(STATUS "Using Boost.MPI libraries '${Boost_MPI_LIBS}'")
    message(STATUS "Enabling compilation of MPI components")
    add_definitions(-DUSE_MPI)
  else()
		disable_feature(MPI)
    message(STATUS "Found MPI installation, but Boost.MPI is not built correctly (consult INSTALL). "
			"No MPI components will be generated.")
  endif()
endif()


# make SVN version string accessible via predefined macro SVN_VERSION
find_package(Subversion)
set_feature_info(Subversion "software versioning and revision control system"
	"http://subversion.apache.org/" "needed only to determine repository revision")
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
    if(NOT _SVNVERSION_RETURN)
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
set_feature_info(Doxygen "documentation generator" "http://www.doxygen.org/" "optional")
if(NOT DOXYGEN_FOUND)
  message(WARNING "Cannot find Doxygen. No documentation will be generated.")
else()
  set(DOXYGEN_TARGET  "doxygen")
  set(DOXYGEN_DOC_DIR "${CMAKE_SOURCE_DIR}/html-doc")
  set(DOXYGEN_CONF    "${CMAKE_SOURCE_DIR}/rootpwaDoxyfile.conf")
  message(STATUS "Run 'make ${DOXYGEN_TARGET}' to create Doxygen documentation files "
		"in '${DOXYGEN_DOC_DIR}'")
  add_custom_target(${DOXYGEN_TARGET}
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_CONF}
    DEPENDS ${DOXYGEN_CONF}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
endif()


# make some environment variables accessible via predefined macro variables
execute_process(COMMAND hostname
  OUTPUT_VARIABLE HOSTNAME
  RESULT_VARIABLE _HOSTNAME_RETURN
  OUTPUT_STRIP_TRAILING_WHITESPACE)  # HOSTNAME is a shell builtin, _not_ an environment variable
if(_HOSTNAME_RETURN)
  set(HOSTNAME "")
endif()
foreach(_CMAKEVAR "CMAKE_HOST_SYSTEM_NAME" "CMAKE_HOST_SYSTEM_PROCESSOR" "CMAKE_HOST_SYSTEM_VERSION"
		"CMAKE_SOURCE_DIR" "CMAKE_BUILD_TYPE" "HOSTNAME" "ROOTSYS"
		"Boost_LIBRARY_VERSION" "Boost_INCLUDE_DIRS" "Boost_MPI_LIBRARY_DIR"
		"Libconfig_VERSION" "Libconfig_DIR" "CUDA_VERSION" "CUDA_LIB_DIRS")
  if(_CMAKEVAR)
    #message(STATUS "${_CMAKEVAR}=${${_CMAKEVAR}}")
    add_definitions(-D'${_CMAKEVAR}=\"${${_CMAKEVAR}}\"')
  else()
    add_definitions(-D'${_CMAKEVAR}=\"\"')
  endif()
endforeach()
foreach(_ENVVAR "USER")
  set(_ENVVARVAL $ENV{${_ENVVAR}})
  if(_ENVVARVAL)
    #message(STATUS "${_ENVVAR}=${_ENVVARVAL}")
    add_definitions(-D'${_ENVVAR}=\"${_ENVVARVAL}\"')
  else()
    add_definitions(-D'${_ENVVAR}=\"\"')
  endif()
endforeach()


enable_testing()


# build subprojects

set(UTILITIES_INCLUDE_DIR
  ${CMAKE_SOURCE_DIR}/utilities
  )

# PWA2000
add_subdirectory(pwa2000)

# RootPwaUtil
set(UTILITIES_LIBS
  RootPwaUtil
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

# RootPwaAmp
set(AMPLITUDE_INCLUDE_DIR
  ${CMAKE_SOURCE_DIR}/amplitude
  )
set(AMPLITUDE_LIBS
  RootPwaAmp
  )
add_subdirectory(amplitude)

# RootPwa
set(ROOTPWA_INCLUDE_DIR
  ${CMAKE_SOURCE_DIR}/src
  )
set(ROOTPWA_LIBS
  RootPwa
  )
add_subdirectory(src)

# compile rootscripts into a library since they are linked by the gui
set(SCRIPTS_INCLUDE_DIR
  ${CMAKE_SOURCE_DIR}/src/rootscripts
  )
add_subdirectory(src/rootscripts)
set(ROOTPWASCRIPTS_LIBS
  RootPwaScripts
  )

# RootPwaTools
set(TOOLS_INCLUDE_DIR
  ${CMAKE_SOURCE_DIR}/tools
  )
add_subdirectory(tools)
set(ROOTPWATOOLS_LIBS
  RootPwaTools
  )

# RootPwaGen
set(GEN_INCLUDE_DIR
  ${CMAKE_SOURCE_DIR}/generators
  )
add_subdirectory(generators)
set(ROOTPWAGEN_LIBS
  RootPwaGen
  )

# RootPwaKey
add_subdirectory(keygen)

# Fortran
#add_subdirectory(fortran)

# Relativic Amplitudes
add_subdirectory(relampl)

# gui
add_subdirectory(gui/src)

# massdep
#add_subdirectory(massdep)

# MPI
if(USE_MPI)
  add_subdirectory(mpi)
endif()


print_enabled_features()
print_disabled_features()
#feature_summary(WHAT ALL FATAL_ON_MISSING_REQUIRED_PACKAGES) requires V2.8.5+


message(STATUS "Run 'make' to build the project.")
message(STATUS "On multi-core machines 'make -j#', where # is the number of parallel jobs, "
	"can speedup compilation considerably.")
