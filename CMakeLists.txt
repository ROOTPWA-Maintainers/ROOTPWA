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


option(DEBUG_OUTPUT "en/disable debug output" OFF)


# define project
message(STATUS "")
message(STATUS ">>> Setting up project 'ROOTPWA'.")
project(ROOTPWA)


# load some common cmake macros
# set path, where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ is checked
message(STATUS "")
message(STATUS ">>> Setting up Cmake modules.")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmakeModules")
message(STATUS "Using cmake module path '${CMAKE_MODULE_PATH}'.")
include(CommonMacros)
include(FeatureSummary)


# force out-of-source builds.
enforce_out_of_source_build()


# warn user if system is not UNIX
message(STATUS "")
message(STATUS ">>> Setting up system environment.")
if(NOT UNIX)
	message(FATAL_ERROR "This is an unsupported system.")
endif()
if(CMAKE_VERSION VERSION_GREATER 2.8.4)
	include(ProcessorCount)  # requires Cmake V2.8.5+
	ProcessorCount(NMB_CPU_CORES)
	message(STATUS "Detected host system '${CMAKE_HOST_SYSTEM_NAME}' version "
		"'${CMAKE_HOST_SYSTEM_VERSION}', architecture '${CMAKE_HOST_SYSTEM_PROCESSOR}', "
		"${NMB_CPU_CORES} CPU core(s).")
else()
	message(STATUS "Detected host system '${CMAKE_HOST_SYSTEM_NAME}' version "
		"'${CMAKE_HOST_SYSTEM_VERSION}', architecture '${CMAKE_HOST_SYSTEM_PROCESSOR}'.")
endif()
message(STATUS "Compiling for system '${CMAKE_SYSTEM_NAME}' version '${CMAKE_SYSTEM_VERSION}', "
	"architecture '${CMAKE_SYSTEM_PROCESSOR}'.")


# define build types
# predefined build types: "DEBUG" "RELEASE" "RELWITHDEBINFO" "MINSIZEREL"
# set a default build type for single-configuration CMake generators, if no build type is set.
message(STATUS "")
message(STATUS ">>> Setting up compiler environment.")
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
if(NOT CMAKE_COMPILER_IS_GNUCXX)
	# use gcc instead of llvm compiler on MacOS 10.7
	if(${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Darwin" AND ${CMAKE_HOST_SYSTEM_VERSION} STREQUAL "11.3.0")
		message(STATUS "Detected MacOS 10.7: using gcc instead of llvm compiler.")
		set(CMAKE_C_COMPILER gcc-4.2)
		set(CMAKE_CXX_COMPILER g++-4.2)
		set(CMAKE_COMPILER_IS_GNUCXX TRUE)
	endif()
	# use gcc instead of llvm compiler on MacOS 10.8
	if(${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Darwin" AND ${CMAKE_HOST_SYSTEM_VERSION} STREQUAL "12.0.0")
		message(STATUS "Detected MacOS 10.8: using gcc instead of llvm compiler.")
		set(CMAKE_C_COMPILER gcc)
		set(CMAKE_CXX_COMPILER g++)
		set(CMAKE_COMPILER_IS_GNUCXX TRUE)
	endif()
endif()
if(CMAKE_COMPILER_IS_GNUCXX)
	execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
	# workaround for bug concerning -march=native switch in gcc 4.2.1 shipped with MacOS 10.7
	# see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=33144
	if(GCC_VERSION VERSION_EQUAL "4.2.1")
		message(STATUS "Detected gcc version 4.2.1; -march=native switch is disabled because of gcc bug.")
		set(CMAKE_CXX_FLAGS_RELEASE "-O3")
	endif()
endif()
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g ${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_LDFLAGS_DEBUG "-g")
set(CMAKE_CXX_LDFLAGS_RELWITHDEBINFO "-g")
# report global build settings
foreach(_LANG "C" "CXX")
	message(STATUS "Using ${_LANG} compiler '${CMAKE_${_LANG}_COMPILER}'.")
	message(STATUS "Using ${_LANG} compiler flags '${CMAKE_${_LANG}_FLAGS}'.")
endforeach()
unset(_LANG)
message(STATUS "Build type is '${CMAKE_BUILD_TYPE}'.")
message(STATUS "Using CXX compiler flags '${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}' "
	"for build type ${CMAKE_BUILD_TYPE}.")
message(STATUS "Using linker flags '${CMAKE_CXX_LDFLAGS_${CMAKE_BUILD_TYPE}}' "
	"for build type ${CMAKE_BUILD_TYPE}.")
#include(CMakePrintSystemInformation)


# redirect output files
message(STATUS "")
message(STATUS ">>> Setting up output paths.")
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
message(STATUS "Using library output path '${LIBRARY_OUTPUT_PATH}'.")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
message(STATUS "Using executable output path '${EXECUTABLE_OUTPUT_PATH}'.")

# setup Boost
# environment variable $BOOST_ROOT is expected to point to non-standard locations
message(STATUS "")
message(STATUS ">>> Setting up Boost library.")
set(_BOOST_COMPONENTS "mpi" "serialization" "python")
# set(Boost_DEBUG 1)
set(Boost_USE_STATIC_LIBS    OFF)
set(Boost_USE_MULTITHREADED  OFF)
set(Boost_USE_STATIC_RUNTIME OFF)
if(     ("$ENV{BOOST_ROOT}"       STREQUAL "")
		AND ("$ENV{BOOSTROOT}"        STREQUAL "")
		AND ("$ENV{Boost_DIR}"        STREQUAL "")
		AND ("$ENV{BOOST_INCLUDEDIR}" STREQUAL "")
		AND ("$ENV{BOOST_LIBRARYDIR}" STREQUAL ""))
	set(Boost_NO_SYSTEM_PATHS OFF)
else()
	set(Boost_NO_SYSTEM_PATHS ON)
endif()
# this is a somewhat ugly hack
# the problem is that components cannot be defined as optional while
# at the same time the library is required. the third find_package
# line is needed in case the components are not found, because
# Boost_FOUND is set to FALSE.
find_package(Boost 1.50.0 REQUIRED)
find_package(Boost 1.50.0 QUIET COMPONENTS ${_BOOST_COMPONENTS})
foreach(_BOOST_COMPONENT ${_BOOST_COMPONENTS})
	string(TOUPPER ${_BOOST_COMPONENT} _BOOST_COMPONENT)
	if(Boost_${_BOOST_COMPONENT}_FOUND)
		message(STATUS "    Found Boost component ${_BOOST_COMPONENT} at "
			"'${Boost_${_BOOST_COMPONENT}_LIBRARY}'.")
	endif()
endforeach()
unset(_BOOST_COMPONENT)
unset(_BOOST_COMPONENTS)
find_package(Boost 1.50.0 REQUIRED)
set_feature_info(Boost
	"C++ Template Library"
	"http://www.boost.org/"
	"required"
	)
if(Boost_FOUND)
	set(Boost_LIBRARY_VERSION "${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
	message(STATUS "Using Boost ${Boost_LIBRARY_VERSION} include directories '${Boost_INCLUDE_DIRS}'.")
	message(STATUS "Using Boost ${Boost_LIBRARY_VERSION} library directories '${Boost_LIBRARY_DIRS}'.")
else()
	message(FATAL_ERROR "Could not find Boost installation. "
		"Is environment variable BOOST_ROOT=${BOOST_ROOT} set correctly? Please read INSTALL.")
endif()


# setup libconfig
# environment variable $LIBCONFIG is expected to point to libconfig location
message(STATUS "")
message(STATUS ">>> Setting up libconfig library.")
find_package(Libconfig 1.4 REQUIRED)
set_feature_info(Libconfig
	"Simple C++ library for processing structured configuration files"
	"http://www.hyperrealm.com/libconfig/"
	"required"
	)


# setup ROOT includes and libraries
message(STATUS "")
message(STATUS ">>> Setting up ROOT.")
find_package(ROOT 5.26 REQUIRED Minuit2 MathMore TreePlayer)
set_feature_info(ROOT
	"C++ HEP analysis framework"
	"http://root.cern.ch/"
	"required"
	)
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
message(STATUS "")
message(STATUS ">>> Setting up CUDA.")
option(USE_CUDA "En/Disables compilation of CUDA components" OFF)
set(CUDA_SHADER_MODELS "13" "20")
if(NOT DEFINED ENV{NVSDKCUDA_ROOT})
	# use default SDK path, if not defined by user
	set(ENV{NVSDKCUDA_ROOT} "$ENV{HOME}/NVIDIA_GPU_Computing_SDK/C")
endif()
find_package(CUDA 4.0)
set_feature_info(CUDA
	"nvidia GPGPU framework"
	"http://developer.nvidia.com/category/zone/cuda-zone"
	"optional"
	)
if(NOT CUDA_FOUND)
	set(USE_CUDA OFF)
	message(STATUS "Cannot find CUDA installation. No CUDA components will be generated.")
else()
	message(STATUS "Found CUDA version ${CUDA_VERSION}")
	message(STATUS "Using CUDA toolkit root directory '${CUDA_TOOLKIT_ROOT_DIR}'.")
	message(STATUS "Using CUDA include directories '${CUDA_INCLUDE_DIRS}'.")
	message(STATUS "Using CUDA libraries '${CUDA_LIBRARIES}'.")
	if(NOT CUDA_SDK_ROOT_DIR)
		set(USE_CUDA OFF)
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
			set(USE_CUDA OFF)
			message(STATUS "Found CUDA intallation, but could not find CUDA utility library "
				"'cutil_${CUDA_CUTIL_LIB_SUFFIX}' in '${CUDA_SDK_ROOT_DIR}/lib'. "
				"No CUDA components will be generated.")
		else()
			set(USE_CUDA ON)
			message(STATUS "Using CUDA utility library '${CUDA_CUTIL_LIBRARY}'.")
			message(STATUS "Using CUDA SDK root directory '${CUDA_SDK_ROOT_DIR}'.")
			message(STATUS "Enabling compilation of CUDA components.")
			# define for which shader models code is generated
			message(STATUS "Generating binary files for shader model(s): ${CUDA_SHADER_MODELS}")
			foreach(_SHADER_MODEL ${CUDA_SHADER_MODELS})
				list(APPEND CUDA_NVCC_FLAGS
					"--generate-code=arch=compute_${_SHADER_MODEL},code=\\\"sm_${_SHADER_MODEL},compute_${_SHADER_MODEL}\\\"")
			endforeach()
			unset(_SHADER_MODEL)
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
message(STATUS "")
message(STATUS ">>> Setting up MPI.")
option(USE_MPI "En/Disables compilation of MPI components" OFF)
find_package(MPI)
set_feature_info(MPI
	"Standardized and portable message-passing system for parallel computing"
	"http://www.boost.org/doc/libs/1_50_0/doc/html/mpi.html"
	"optional"
	)
if(NOT MPI_CXX_FOUND)
	set(USE_MPI OFF)
	message(STATUS "Cannot find MPI installation. No MPI components will be generated.")
else()
	set(USE_MPI ON)
	message(STATUS "Using MPI compiler wrapper '${MPI_CXX_COMPILER}'.")
	message(STATUS "Using CXX compiler flags '${MPI_CXX_COMPILE_FLAGS}' for MPI code.")
	message(STATUS "Using CXX linker flags '${MPI_CXX_LINK_FLAGS}' for MPI code.")
	message(STATUS "Using MPI include directories '${MPI_CXX_INCLUDE_PATH}'.")
	message(STATUS "Use '${MPIEXEC}' executable to run MPI programs.")
	#message(STATUS "Using MPI libraries '${MPI_CXX_LIBRARIES}'.")
	# check that Boost.MPI library was found
	set(_BOOST_MPI_COMPONENTS "mpi" "serialization")
	foreach(_BOOST_MPI_COMPONENT ${_BOOST_MPI_COMPONENTS})
		string(TOUPPER ${_BOOST_MPI_COMPONENT} _BOOST_MPI_COMPONENT)
		if(NOT Boost_${_BOOST_MPI_COMPONENT}_FOUND)
			set(USE_MPI OFF)
			message(STATUS "Cannot find Boost.MPI library '${_BOOST_MPI_COMPONENT}' in "
				"'${Boost_LIBRARY_DIRS}'. Please read INSTALL on how to build Boost.MPI.")
		else()
			list(APPEND Boost_MPI_LIBS ${Boost_${_BOOST_MPI_COMPONENT}_LIBRARY})
		endif()
	endforeach()
	unset(_BOOST_MPI_COMPONENTS)
	if(USE_MPI)
		message(STATUS "Using Boost.MPI libraries '${Boost_MPI_LIBS}'.")
		message(STATUS "Enabling compilation of MPI components.")
		add_definitions(-DUSE_MPI)
	else()
		disable_feature(MPI)
		message(STATUS "Found MPI installation, but Boost.MPI is not built correctly (consult INSTALL). "
			"No MPI components will be generated.")
	endif()
endif()


# setup Python
message(STATUS "")
message(STATUS ">>> Setting up Python.")
option(USE_PYTHON "En/Disables Python bindings" OFF)
find_package(Python 2.7)
set_feature_info(Python
	"General-purpose, high-level programming language"
	"http://www.python.org"
	"optional"
	)
if(NOT PYTHONLIBS_FOUND)
	set(USE_PYTHON OFF)
	message(STATUS "Cannot find Python libraries. No Python bindings will be generated.")
else()
	set(USE_PYTHON ON)
	# check that Boost.Python library was found
	set(_BOOST_PYTHON_COMPONENTS "python")
	foreach(_BOOST_PYTHON_COMPONENT ${_BOOST_PYTHON_COMPONENTS})
		string(TOUPPER ${_BOOST_PYTHON_COMPONENT} _BOOST_PYTHON_COMPONENT)
		if(NOT Boost_${_BOOST_PYTHON_COMPONENT}_FOUND)
			set(USE_PYTHON OFF)
			message(STATUS "Cannot find Boost.Python library in "
				"'${Boost_LIBRARY_DIRS}'. Please read INSTALL on how to build Boost.Python.")
		else()
			list(APPEND Boost_PYTHON_LIBS ${Boost_${_BOOST_PYTHON_COMPONENT}_LIBRARY})
		endif()
	endforeach()
	unset(_BOOST_PYTHON_COMPONENTS)
	if(USE_PYTHON)
		message(STATUS "Using Boost.Python libraries '${Boost_PYTHON_LIBS}'.")
		message(STATUS "Enabling compilation of Python components.")
		add_definitions(-DUSE_PYTHON)
	else()
		disable_feature(Python)
		message(STATUS "Found Python installation, but Boost.Python is not built correctly (consult INSTALL). "
			"No Python bindings will be generated.")
	endif()
endif()


# make SVN version string accessible via predefined macro SVN_VERSION
message(STATUS "")
message(STATUS ">>> Setting up Subversion.")
find_package(Subversion)
set_feature_info(Subversion
	"Software versioning and revision control system"
	"http://subversion.apache.org/"
	"needed only to determine repository revision"
	)
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
			message(STATUS "Subversion repository revision is '${SVN_VERSION}'.")
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
message(STATUS "")
message(STATUS ">>> Setting up Doxygen.")
find_package(Doxygen)
set_feature_info(Doxygen
	"Documentation generator"
	"http://www.doxygen.org/"
	"optional"
	)
if(NOT DOXYGEN_FOUND)
	message(WARNING "Cannot find Doxygen. Documentation cannot be generated.")
else()
	set(DOXYGEN_TARGET  "doxygen")
	set(DOXYGEN_DOC_DIR "${CMAKE_SOURCE_DIR}/html-doc")
	set(DOXYGEN_CONF    "${CMAKE_SOURCE_DIR}/rootpwaDoxyfile.conf")
	message(STATUS "Run 'make ${DOXYGEN_TARGET}' to create Doxygen documentation files "
		"in '${DOXYGEN_DOC_DIR}'.")
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
	OUTPUT_STRIP_TRAILING_WHITESPACE)
if(_HOSTNAME_RETURN)
	set(HOSTNAME "")
endif()
unset(_HOSTNAME_RETURN)
if(DEBUG_OUTPUT)
	message(STATUS "Adding definitions for internal variables:")
endif()
foreach(_CMAKEVAR
		"CMAKE_HOST_SYSTEM_NAME"
		"CMAKE_HOST_SYSTEM_PROCESSOR"
		"CMAKE_HOST_SYSTEM_VERSION"
		"NMB_CPU_CORES"
		"HOSTNAME"
		"CMAKE_SOURCE_DIR"
		"CMAKE_BUILD_TYPE"
		"Boost_LIBRARY_VERSION"
		"Boost_INCLUDE_DIRS"
		"Libconfig_VERSION"
		"Libconfig_DIR"
		"ROOTSYS"
		"CUDA_VERSION"
		"CUDA_LIB_DIRS"
		"PYTHONLIBS_VERSION_STRING"
		"PYTHON_INCLUDE_DIRS"
		)
	if(_CMAKEVAR)
		if(DEBUG_OUTPUT)
			message(STATUS "        ${_CMAKEVAR} = ${${_CMAKEVAR}}")
		endif()
		add_definitions(-D'${_CMAKEVAR}=\"${${_CMAKEVAR}}\"')
	else()
		add_definitions(-D'${_CMAKEVAR}=\"\"')
	endif()
endforeach()
unset(_CMAKEVAR)
if(DEBUG_OUTPUT)
	message(STATUS "Adding definitions for environment variables:")
endif()
foreach(_ENVVAR
		"USER"
		)
	set(_ENVVARVAL $ENV{${_ENVVAR}})
	if(_ENVVARVAL)
		if(DEBUG_OUTPUT)
			message(STATUS "        ${_ENVVAR} = ${_ENVVARVAL}")
		endif()
		add_definitions(-D'${_ENVVAR}=\"${_ENVVARVAL}\"')
	else()
		add_definitions(-D'${_ENVVAR}=\"\"')
	endif()
	unset(_ENVVARVAL)
endforeach()
unset(_ENVVAR)


enable_testing()


# build subprojects

set(UTILITIES_INCLUDE_DIR
	${CMAKE_SOURCE_DIR}/utilities
	)

# PWA2000
add_subdirectory(pwa2000)

# Utilities
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

# Python bindings
if(USE_PYTHON)
	set(PYINTERFACE_LIBS
		RootPwaPy
		)
	add_subdirectory(pyInterface)
endif()


# print feature summary
message(STATUS "")
message(STATUS ">>> Feature summary:")
#if(CMAKE_VERSION VERSION_GREATER 2.8.4)
#	feature_summary(WHAT ALL FATAL_ON_MISSING_REQUIRED_PACKAGES)  # requires Cmake V2.8.5+
#else()
print_enabled_features()
print_disabled_features()
#endif()


message(STATUS "")
message(STATUS ">>> Done.")
message(STATUS "Run 'make' to build the project.")
message(STATUS "On multi-core machines 'make -j#', where # is the number of parallel jobs, "
	"can speedup compilation considerably.")
