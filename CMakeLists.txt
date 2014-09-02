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
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)


# set verbosity
set(CMAKE_VERBOSE_MAKEFILE 0)  # if set to 1 compile and link commands are displayed during build
                               # the same effect can be achieved by calling 'make VERBOSE=1'


option(DEBUG_OUTPUT "en/disable debug output" OFF)


# define project
message(STATUS "")
message(STATUS ">>> Setting up project 'ROOTPWA'.")
project(ROOTPWA)


# set path, where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ is checked
message(STATUS "")
message(STATUS ">>> Setting up Cmake modules.")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmakeModules")
message(STATUS "Using cmake module path '${CMAKE_MODULE_PATH}'.")
# load cmake macros
include(CommonMacros)
include(FeatureSummary)
include(ProcessorCount)
include(CheckCXXSourceCompiles)


# force out-of-source builds.
enforce_out_of_source_build()


# warn user if system is not UNIX
message(STATUS "")
message(STATUS ">>> Setting up system environment.")
if(NOT UNIX)
	message(FATAL_ERROR "This is an unsupported system.")
endif()
ProcessorCount(NMB_CPU_CORES)
message(STATUS "Detected host system '${CMAKE_HOST_SYSTEM_NAME}' version "
	"'${CMAKE_HOST_SYSTEM_VERSION}', architecture '${CMAKE_HOST_SYSTEM_PROCESSOR}', "
	"${NMB_CPU_CORES} CPU core(s).")
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


# set common compiler flags
# consider flags -Wextra -pedantic -Wsuggest-attribute=pure -Wsuggest-attribute=const -Wsuggest-attribute=noreturn
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Woverloaded-virtual -Werror")
# set additional flags for specific build types
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
# -march=native switch was removed, because is caused trouble in heterogeneous environments
# `gcc -march=native -Q --help=target` shows the relevant compiler switches
# consider adding -DNODEBUG (disables all asserts) and -ffast-math
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
# code should now compile with clang so forcing of gcc should not be necessary anymore
# if(NOT CMAKE_COMPILER_IS_GNUCXX AND ${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Darwin")
# 	# use gcc instead of llvm compiler on MacOS 10.7 (Darwin 11.x.y)
# 	if(${CMAKE_HOST_SYSTEM_VERSION} VERSION_GREATER "10.999.999" AND ${CMAKE_HOST_SYSTEM_VERSION} VERSION_LESS "12.0.0")
# 		message(STATUS "Detected MacOS 10.7: using gcc instead of llvm compiler.")
# 		set(CMAKE_C_COMPILER gcc-4.2)
# 		set(CMAKE_CXX_COMPILER g++-4.2)
# 		set(CMAKE_COMPILER_IS_GNUCXX TRUE)
# 	endif()
# 	# use gcc instead of llvm compiler on MacOS 10.8 (Darwin 12.x.y)
# 	if(${CMAKE_HOST_SYSTEM_VERSION} VERSION_GREATER "11.999.999" AND ${CMAKE_HOST_SYSTEM_VERSION} VERSION_LESS "13.0.0")
# 		message(STATUS "Detected MacOS 10.8: using gcc instead of llvm compiler.")
# 		set(CMAKE_C_COMPILER gcc)
# 		set(CMAKE_CXX_COMPILER g++)
# 		set(CMAKE_COMPILER_IS_GNUCXX TRUE)
# 	endif()
# endif()
# -march=native switch was disabled so no need anymore to handle the gcc 4.2.1 bug
# if(CMAKE_COMPILER_IS_GNUCXX)
# 	execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
# 	# workaround for bug concerning -march=native switch in gcc 4.2.1 shipped with MacOS 10.7
# 	# see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=33144
# 	if(GCC_VERSION VERSION_EQUAL "4.2.1")
# 		message(STATUS "Detected gcc version 4.2.1; -march=native switch is disabled because of gcc bug.")
# 		set(CMAKE_CXX_FLAGS_RELEASE "-O3")
# 	endif()
# endif()
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -g")
# report global build settings
foreach(_LANG "C" "CXX")
	message(STATUS "${_LANG} compiler ID is '${CMAKE_${_LANG}_COMPILER_ID} V${CMAKE_${_LANG}_COMPILER_VERSION}'.")
	message(STATUS "Using ${_LANG} compiler '${CMAKE_${_LANG}_COMPILER}'.")
	message(STATUS "Using general ${_LANG} compiler flags '${CMAKE_${_LANG}_FLAGS}'.")
endforeach()
unset(_LANG)
message(STATUS "Build type is '${CMAKE_BUILD_TYPE}'.")
message(STATUS "Using additional CXX compiler flags '${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}' "
	"for build type ${CMAKE_BUILD_TYPE}.")
message(STATUS "Using additional linker flags '${CMAKE_CXX_LDFLAGS_${CMAKE_BUILD_TYPE}}' "
	"for build type ${CMAKE_BUILD_TYPE}.")
if(DEBUG_OUTPUT)
	include(CMakePrintSystemInformation)
endif()


# detect some compiler/STL features
# test whether std::complex provides real and imag accessors
check_cxx_source_compiles("#include <complex>
int main() { std::complex<double> a; a.real() = 0; }" COMPILER_PROVIDES_STD_COMPLEX_ACCESSORS)
if(COMPILER_PROVIDES_STD_COMPLEX_ACCESSORS)
	add_definitions(-DCOMPILER_PROVIDES_STD_COMPLEX_ACCESSORS)
endif()


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
set(_BOOST_COMPONENTS "mpi" "serialization" "python" "timer" "system")
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
set_package_properties(Boost
	PROPERTIES
	DESCRIPTION "C++ Template Library"
	URL "http://www.boost.org/"
	TYPE REQUIRED
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
set_package_properties(Libconfig
	PROPERTIES
	DESCRIPTION "Simple C++ library for processing structured configuration files"
	URL "http://www.hyperrealm.com/libconfig/"
	TYPE REQUIRED
	)


# setup ROOT includes and libraries
message(STATUS "")
message(STATUS ">>> Setting up ROOT.")
set(_ROOT_COMPONENTS "Minuit2" "MathMore" "TreePlayer")
find_package(ROOT 5.26 REQUIRED ${_ROOT_COMPONENTS})
set_package_properties(ROOT
	PROPERTIES
	DESCRIPTION "C++ HEP analysis framework"
	URL "http://root.cern.ch/"
	TYPE REQUIRED
	)
# some code switches that depend on ROOT version
# std::complex is not supported as TTree leaf in ROOT versions below 5.27.06
if(NOT ROOT_VERSION VERSION_LESS "5.27.06")
	set(USE_STD_COMPLEX_TREE_LEAFS TRUE)
	add_definitions(-DUSE_STD_COMPLEX_TREE_LEAFS)
endif()
add_feature_info(Complex_tree_leaf USE_STD_COMPLEX_TREE_LEAFS "The Complex_tree_leaf feature provides ROOT data formats for amplitudes and amplitude integrals.")


# setup CUDA
message(STATUS "")
message(STATUS ">>> Setting up CUDA.")
option(USE_CUDA "En/Disables compilation of CUDA components" OFF)
set(CUDA_SHADER_MODELS "20")
find_package(CUDA 5.5)
set_package_properties(CUDA
        PROPERTIES
        DESCRIPTION "nvidia GPGPU framework"
	URL "http://developer.nvidia.com/category/zone/cuda-zone"
	PURPOSE "Used to accelerate calculations, if system is equipped with suitable GPU"
	TYPE OPTIONAL
	)
if(NOT CUDA_FOUND)
	set(USE_CUDA OFF)
	message(STATUS "Cannot find CUDA installation. No CUDA components will be generated.")
else()
	message(STATUS "Found CUDA version ${CUDA_VERSION}")
        message(STATUS "Using CUDA toolkit root directory '${CUDA_TOOLKIT_ROOT_DIR}'.")
        message(STATUS "Using CUDA include directories '${CUDA_INCLUDE_DIRS}'.")
        message(STATUS "Using CUDA libraries '${CUDA_LIBRARIES}'.")

        if(NOT DEFINED ENV{CUDA_SAMPLES_ROOT_DIR})
                # use default path for CUDA samples, if not defined by user
                set(ENV{CUDA_SAMPLES_ROOT_DIR} "$ENV{HOME}/NVIDIA_CUDA-${CUDA_VERSION}_Samples")
                endif()
        # find include directory
        find_file(CUDA_SAMPLES_INC_DIR
                        NAMES inc
            PATHS $ENV{CUDA_SAMPLES_ROOT_DIR}/common
            NO_DEFAULT_PATH)
        if(NOT CUDA_SAMPLES_INC_DIR)
                        set(USE_CUDA OFF)
                message(STATUS "Found CUDA intallation, but could not find include directory in CUDA samples directory in '$ENV{CUDA_SAMPLES_ROOT_DIR}'. "
                        "Check that environment variable $CUDA_SAMPLES_ROOT_DIR is set correctly. No CUDA components will be generated.")
                else()
                        set(USE_CUDA ON)
                message(STATUS "Using CUDA samples include directory in '${CUDA_SAMPLES_INC_DIR}'.")
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
        endif()
endif()
add_feature_info(CUDA_acceleration USE_CUDA "The CUDA_acceleration feature provides GPU acceleration.")


# setup MPI
message(STATUS "")
message(STATUS ">>> Setting up MPI.")
option(USE_MPI "En/Disables compilation of MPI components" OFF)
find_package(MPI)
set_package_properties(MPI
	PROPERTIES
	DESCRIPTION "Standardized and portable message-passing system for parallel computing"
	URL "http://www.boost.org/doc/libs/1_50_0/doc/html/mpi.html"
	PURPOSE "At the moment just used for tests"
	TYPE OPTIONAL
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
		message(STATUS "Found MPI installation, but Boost.MPI is not built correctly (consult INSTALL). "
			"No MPI components will be generated.")
	endif()
endif()
add_feature_info(MPI_parallization USE_MPI "The MPI_parallization feature will allow parallelized execution of some tools.")


# setup Python
message(STATUS "")
message(STATUS ">>> Setting up Python.")
option(USE_PYTHON "En/Disables Python bindings" OFF)
if(ROOT_VERSION VERSION_LESS "5.34.04")
	message(STATUS "Python bindings require ROOT version 5.34.04 or newer "
		"(found ${ROOT_VERSION}).")
	message(STATUS "No Python bindings will be generated")
	set(USE_PYTHON OFF)
else()
	find_package(PYTHON 2.6)
	set_package_properties(PYTHON
		PROPERTIES
		DESCRIPTION "General-purpose, high-level programming language"
		URL "http://www.python.org"
		PURPOSE "Provides scripting interface to ROOTPWA classes"
		TYPE OPTIONAL
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
			message(STATUS "Setting up ROOT again now with PyROOT library.")
			find_package(ROOT 5.34.04 REQUIRED ${_ROOT_COMPONENTS} PyROOT)
			message(STATUS "Using Boost.Python libraries '${Boost_PYTHON_LIBS}'.")
			message(STATUS "Enabling compilation of Python components.")
			add_definitions(-DUSE_PYTHON)
			set(PYTHON_LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/pyLib")
			file(MAKE_DIRECTORY ${PYTHON_LIBRARY_OUTPUT_PATH})
			file(MAKE_DIRECTORY ${PYTHON_LIBRARY_OUTPUT_PATH}/pyRootPwa)
			message(STATUS "Using python bindings output path '${PYTHON_LIBRARY_OUTPUT_PATH}'.")
		else()
			message(STATUS "Found Python installation, but Boost.Python is not built correctly (consult INSTALL). "
				"No Python bindings will be generated.")
		endif()
	endif()
endif()
unset(_ROOT_COMPONENTS)
add_feature_info(Python_scripting USE_PYTHON "The Python_scripting feature provides a scripting interface for some classes.")


# make git hash accessible via predefined macro GIT_HASH
message(STATUS "")
message(STATUS ">>> Setting up Git.")
find_package(Git)
set_package_properties(Git
	PROPERTIES
	DESCRIPTION "Free and open source distributed version control system"
	URL "http://git-scm.com/"
	PURPOSE "Used to determine repository hash that is put into executables"
	TYPE RECOMMENDED
	)
if(GIT_FOUND)
	execute_process(
		COMMAND ${GIT_EXECUTABLE} log --pretty="%H" -n1 -- ${CMAKE_SOURCE_DIR}
		OUTPUT_VARIABLE GIT_HASH
		RESULT_VARIABLE _GIT_LOG_RETURN
		OUTPUT_STRIP_TRAILING_WHITESPACE
		)
	if(NOT _GIT_LOG_RETURN)
		message(STATUS "Git repository hash is '${GIT_HASH}'.")
	else()
		message(STATUS "Error running 'git'. Repository hash unknown.")
	endif()
endif()


# setup doxygen
message(STATUS "")
message(STATUS ">>> Setting up Doxygen.")
find_package(Doxygen)
set_package_properties(Doxygen
	PROPERTIES
	DESCRIPTION "Documentation generator"
	URL "http://www.doxygen.org/"
  PURPOSE "Generates HTML source-code documentation"
	TYPE RECOMMENDED
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


enable_testing()


# set include and library paths
set(UTILITIES_INCLUDE_DIR
	${CMAKE_SOURCE_DIR}/utilities
	)
set(UTILITIES_LIBS
	RootPwaUtil
	)

set(STORAGEFORMATS_INCLUDE_DIR
	${CMAKE_SOURCE_DIR}/storageFormats
	)
set(STORAGEFORMATS_LIBS
	RootPwaMeta
	)

set(TOOLS_INCLUDE_DIR
	${CMAKE_SOURCE_DIR}/tools
	)
set(ROOTPWATOOLS_LIBS
	RootPwaTools
	)

set(AMPLITUDE_INCLUDE_DIR
	${CMAKE_SOURCE_DIR}/amplitude
	)
set(AMPLITUDE_LIBS
	RootPwaAmp
	)

set(ROOTPWA_INCLUDE_DIR
	${CMAKE_SOURCE_DIR}/src
	)
set(ROOTPWA_LIBS
	RootPwa
	)

set(GEN_INCLUDE_DIR
	${CMAKE_SOURCE_DIR}/generators
	)
set(ROOTPWA_NBPSGEN_LIBS
	RootPwaNBPSGen
	)
set(ROOTPWAGEN_LIBS
	RootPwaGen
	)

if(USE_PYTHON)
	set(PYINTERFACE_LIBS
		RootPwaPy
		)
endif()

if(USE_CUDA)
        set(CUDACODE_INCLUDE_DIR
                ${CMAKE_SOURCE_DIR}/cuda
                )
        set(CUDACODE_LIBS
                RootPwaCuda
                )
else()
        set(CUDACODE_INCLUDE_DIR)
        set(CUDACODE_LIBS)
endif()


# build subprojects

# Utilities
add_subdirectory(utilities)

# Metadata
add_subdirectory(storageFormats)

# CUDA stuff
if(USE_CUDA)
	add_subdirectory(cuda)
endif()

# RootPwaAmp
add_subdirectory(amplitude)

# RootPwa
add_subdirectory(src)

# RootPwaTools
add_subdirectory(tools)

# RootPwaGen
add_subdirectory(generators)

# Fortran
#add_subdirectory(fortran)

# Relativic Amplitudes
add_subdirectory(relativisticAmpCorrections)

# massdep
#add_subdirectory(massdep)

# MPI
if(USE_MPI)
	add_subdirectory(mpi)
endif()

# Python bindings
if(USE_PYTHON)
	add_subdirectory(pyInterface)
endif()


# print package and feature summary
message(STATUS "")
message(STATUS ">>> Package summary:")
feature_summary(WHAT PACKAGES_FOUND
	FATAL_ON_MISSING_REQUIRED_PACKAGES
	DESCRIPTION "Found packages:")
feature_summary(WHAT PACKAGES_NOT_FOUND
	FATAL_ON_MISSING_REQUIRED_PACKAGES
	DESCRIPTION "Packages not found:")
message(STATUS "")
message(STATUS ">>> Feature summary:")
feature_summary(WHAT ENABLED_FEATURES
	FATAL_ON_MISSING_REQUIRED_PACKAGES
	DESCRIPTION "Enabled features:")
feature_summary(WHAT DISABLED_FEATURES
	FATAL_ON_MISSING_REQUIRED_PACKAGES
	DESCRIPTION "Disabled features:")


message(STATUS "")
message(STATUS ">>> Done.")
message(STATUS "Run 'make' to build the project.")
if(${NMB_CPU_CORES} GREATER 1)
	message(STATUS "On this machine using 'make -j${NMB_CPU_CORES}' can speedup compilation considerably.")
endif()
