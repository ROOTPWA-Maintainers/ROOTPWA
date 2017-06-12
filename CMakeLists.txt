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
message(STATUS ">>> This is CMake version ${CMAKE_VERSION}.")
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
cmake_policy(VERSION 3.1)


# set verbosity
set(CMAKE_VERBOSE_MAKEFILE 0)  # if set to 1 compile and link commands are displayed during build
                               # the same effect can be achieved by calling 'make VERBOSE=1'


option(DEBUG_OUTPUT "en/disable debug output" OFF)


# define project
message(STATUS "")
message(STATUS ">>> Setting up project 'ROOTPWA'.")
project(ROOTPWA)


# require a C++11 compatible compiler
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# set path, where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ is checked
message(STATUS "")
message(STATUS ">>> Setting up Cmake modules.")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmakeModules")
message(STATUS "Using cmake module path '${CMAKE_MODULE_PATH}'.")
# load cmake macros
include(CommonMacros)
include(FeatureSummary)
include(ProcessorCount)


# On Apples the default filename extension of a shared library is ".dylib",
# unfortunately ROOT libraries end with ".so" so Cmake's mechanism of
# determining whether to add a rpath breaks. Also add ".so" to the filetypes
# to be considered as a shared library.
if(APPLE AND "${CMAKE_SHARED_LIBRARY_SUFFIX}" STREQUAL ".dylib")
	list(APPEND CMAKE_EXTRA_SHARED_LIBRARY_SUFFIXES ".so")
endif()


# force out-of-source builds.
enforce_out_of_source_build()


# abort, if system is not UNIX
message(STATUS "")
message(STATUS ">>> Setting up system environment.")
if(NOT UNIX)
	message(FATAL_ERROR " !!! This is an unsupported system.")
endif()
ProcessorCount(NMB_CPU_CORES)
message(STATUS "Detected host system '${CMAKE_HOST_SYSTEM_NAME}' version "
	"'${CMAKE_HOST_SYSTEM_VERSION}', architecture '${CMAKE_HOST_SYSTEM_PROCESSOR}', "
	"${NMB_CPU_CORES} CPU core(s).")
message(STATUS "Compiling for system '${CMAKE_SYSTEM_NAME}' version '${CMAKE_SYSTEM_VERSION}', "
	"architecture '${CMAKE_SYSTEM_PROCESSOR}'.")


# warn user about outdated or unsupported compiler
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
    message(WARNING " ??? Your GCC version ${CMAKE_CXX_COMPILER_VERSION} is too low. GCC version 4.8 or higher is required. Proceed at your own risk.")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.4)
    message(WARNING " ??? Your Clang version ${CMAKE_CXX_COMPILER_VERSION} is too low. "
			"Clang version 3.4 or higher is required. Proceed at your own risk.")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0)
    message(WARNING " ??? Your AppleClang version ${CMAKE_CXX_COMPILER_VERSION} is too low. "
			"AppleClang version 8.0 or higher is required. Proceed at your own risk.")
  endif()
else()
  message(WARNING " ??? '${CMAKE_CXX_COMPILER_ID} V${CMAKE_CXX_COMPILER_VERSION}' "
		"is an unsupported compiler. Compilation has only been tested with GCC, Clang and AppleClang. "
		"Proceed at your own risk.")
endif()


# define build types
# predefined build types: "DEBUG" "RELEASE" "RELWITHDEBINFO" "MINSIZEREL"
# set a default build type for single-configuration CMake generators, if no build type is set.
message(STATUS "")
message(STATUS ">>> Setting up compiler environment.")
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type was specified. Setting build type to 'RELWITHDEBINFO'.")
	set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
endif()


# set common compiler flags
# consider flags -pedantic -Wsuggest-attribute=pure -Wsuggest-attribute=noreturn -Wsuggest-attribute=const
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Woverloaded-virtual -Werror")
# options switchable via command line
# enable option by running cmake .. -D <VAR>=<any value>
if(ARCH_NATIVE)
	# !Note! may cause trouble in heterogeneous environments
	# `gcc -march=native -Q --help=target` shows the compiler switches set by this option on your machine
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()
# set additional flags for specific build types
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
# consider adding -DNODEBUG (disables all asserts) and -ffast-math
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
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
set(Boost_USE_MULTITHREADED  ON)
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
find_package(Boost 1.63.0 REQUIRED)
find_package(Boost 1.63.0 QUIET COMPONENTS ${_BOOST_COMPONENTS})
foreach(_BOOST_COMPONENT ${_BOOST_COMPONENTS})
	string(TOUPPER ${_BOOST_COMPONENT} _BOOST_COMPONENT)
	if(Boost_${_BOOST_COMPONENT}_FOUND)
		message(STATUS "    Found Boost component ${_BOOST_COMPONENT} at "
			"'${Boost_${_BOOST_COMPONENT}_LIBRARY}'.")
	endif()
endforeach()
unset(_BOOST_COMPONENT)
unset(_BOOST_COMPONENTS)
find_package(Boost 1.63.0 REQUIRED)
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
	message(FATAL_ERROR " !!! Could not find Boost installation. "
		"Is environment variable BOOST_ROOT=${BOOST_ROOT} set correctly? Please consult README.md.")
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


# setup yaml-cpp
# environment variable $YAML_CPP is expected to point to yaml-cpp location
message(STATUS "")
message(STATUS ">>> Setting up yaml-cpp library.")
find_package(YamlCpp 0.5 REQUIRED)
set_package_properties(YamlCpp
	PROPERTIES
	DESCRIPTION "Simple C++ library for processing structured configuration files"
	URL "https://github.com/jbeder/yaml-cpp"
	TYPE REQUIRED
	)


# setup ROOT includes and libraries
message(STATUS "")
message(STATUS ">>> Setting up ROOT.")
set(_ROOT_COMPONENTS "Minuit2" "MathMore" "TreePlayer" "PyROOT")
find_package(ROOT 5.34.04 REQUIRED ${_ROOT_COMPONENTS})
set_package_properties(ROOT
	PROPERTIES
	DESCRIPTION "C++ HEP analysis framework"
	URL "http://root.cern.ch/"
	TYPE REQUIRED
	)


# setup Python
message(STATUS "")
message(STATUS ">>> Setting up Python.")
find_package(PYTHON 2.6 REQUIRED)
set_package_properties(PYTHON
	PROPERTIES
	DESCRIPTION "General-purpose, high-level programming language"
	URL "http://www.python.org"
	PURPOSE "Provides scripting interface to ROOTPWA classes"
	TYPE REQUIRED
	)
# check that Boost.Python library was found
if(NOT Boost_PYTHON_FOUND)
	message(FATAL_ERROR " !!! Cannot find Boost.Python library in "
		"'${Boost_LIBRARY_DIRS}'. Please consult README.md on how to build Boost.Python.")
endif()
set(Boost_PYTHON_LIBS "${Boost_PYTHON_LIBRARY}")
message(STATUS "Using Boost.Python libraries '${Boost_PYTHON_LIBS}'.")
set(PYTHON_LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/pyLib")
message(STATUS "Using python bindings output path '${PYTHON_LIBRARY_OUTPUT_PATH}'.")


# setup numpy
message(STATUS "")
message(STATUS ">>> Setting up numpy.")
find_package(NumPy 1.8 REQUIRED)
set_package_properties(NumPy
	PROPERTIES
	DESCRIPTION "numerical calculations in Python"
	URL "http://www.numpy.org"
	PURPOSE "Used for handling of multi-dimensional arrays in Python."
	TYPE REQUIRED
	)
add_definitions(-DNPY_NO_DEPRECATED_API=NPY_1_8_API_VERSION)


# setup pylint (optional)
message(STATUS "")
message(STATUS ">>> Setting up pylint.")
find_package(PYLINT 1.5.5)
set_package_properties(PYLINT
        PROPERTIES
        DESCRIPTION "code analysis for Python"
        URL "https://www.pylint.org"
        PURPOSE "Used to assess quality of Python code."
        TYPE OPTIONAL
        )
if(NOT PYLINT_FOUND)
        message(STATUS "Cannot find pylint executable. Quality of python code will no be assessed.")
endif()


# setup CUDA (optional)
message(STATUS "")
message(STATUS ">>> Setting up CUDA.")
option(USE_CUDA "En/Disables compilation of CUDA components" OFF)
set(CUDA_SHADER_MODELS "20")
find_package(CUDA 8.0)
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

	# find CUDA sample directory
	if(NOT DEFINED ENV{CUDA_SAMPLES_ROOT_DIR})
		# use default path for CUDA samples, if not defined by user
		set(ENV{CUDA_SAMPLES_ROOT_DIR} "${CUDA_TOOLKIT_ROOT_DIR}/samples")
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
if (NOT USE_CUDA)
	set(DO_NOT_COMPILE_RPWA_CUDA TRUE)
endif()


# setup MPI (optional)
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
				"'${Boost_LIBRARY_DIRS}'. Please consult README.md on how to build Boost.MPI.")
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
		message(STATUS "Found MPI installation, but Boost.MPI is not built correctly. "
			"Please consult README.md on how to build Boost.MPI. No MPI components will be generated.")
	endif()
endif()
add_feature_info(MPI_parallization USE_MPI "The MPI_parallization feature will allow parallelized execution of some tools.")
if (NOT USE_MPI)
	set(DO_NOT_COMPILE_RPWA_MPI TRUE)
endif()


# setup BAT (optional)
message(STATUS "")
message(STATUS ">>> Setting up BAT.")
option(USE_BAT "En/Disables usage of BAT" OFF)
find_package(BAT 1.0)
set_package_properties(BAT
	PROPERTIES
	DESCRIPTION "Bayesian Analysis Toolkit"
	URL "http://mpp.mpg.de/bat/"
	PURPOSE "Used in the mass-dependent fitter"
	)
if(NOT BAT_FOUND)
	set(USE_BAT OFF)
	message(STATUS "Cannot find BAT installation. BAT will not be used.")
else()
	set(USE_BAT ON)
	add_definitions(-DUSE_BAT)
endif()
add_feature_info(Bayesian_Analysis_Toolkit USE_BAT "Allow usage of the Bayesian Analysis Toolkit.")


# setup nlopt (optional)
message(STATUS "")
message(STATUS ">>> Setting up NLopt.")
find_package(NLopt)
set_package_properties(NLopt
	PROPERTIES
	DESCRIPTION "Library for nonlinear optimization"
	URL "http://ab-initio.mit.edu/wiki/index.php/NLopt"
	TYPE OPTIONAL
	)
if(NOT NLopt_FOUND)
	set(USE_NLOPT OFF)
	message(STATUS "Cannot find NLopt installation. No NLopt components will be generated.")
else()
	set(USE_NLOPT ON)
	add_definitions(-DUSE_NLOPT)
endif()
add_feature_info(NLopt USE_NLOPT "NLopt is a library for nonlinear optimization.")


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
		message(Warning " ??? Running 'git' did not succeed. Repository hash unknown.")
	endif()
endif()


# setup doxygen (optional)
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
	message(WARNING " ??? Cannot find Doxygen. Documentation cannot be generated.")
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


# define include and library paths for project-internal libraries
#
# the following scheme is used:
#
# * for each source folder a shared library is generated
#
# * the name 'libRootPwa<Folder>.so' of the libary is derived from the
#   directory name, e.g. the code in the 'foo' directory will be
#   compiled into 'libRootPwaFoo.so'.
#
# * for each source folder two CMake variables are defined:
#   RPWA_<FOLDER>_INCLUDE_DIR that points to the source directory and
#   RPWA_<FOLDER>_LIB that holds the library name
#
# * compilation of each directory can be switched off by defining a
#   variable with name DO_NOT_COMPILE_RPWA_<FOLDER> and setting it to
#   a non-empty value
#
# * each folder may contain one linkdef.h file, which will be used to
#   generate ROOT dictionaries (dict.cc)

message(STATUS "")
message(STATUS ">>> Defining include and library paths.")
set(_SRC_DIRS
	cuda
	decayAmplitude
	fortran
	generators
	highLevelInterface
	nBodyPhaseSpace
	partialWaveFit
	particleData
	pyInterface
	relativisticAmpCorrections
	resonanceFit
	storageFormats
	utilities
	)
# loop over all source directories
foreach(_SRC_DIR ${_SRC_DIRS})
	var_name_from_dir_name(${_SRC_DIR} _VAR_NAME)
	if(NOT DO_NOT_COMPILE_${_VAR_NAME})
		lib_name_from_dir_name(${_SRC_DIR} _LIB_NAME)
		set(${_VAR_NAME}_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/${_SRC_DIR}")
		set(${_VAR_NAME}_LIB "${_LIB_NAME}")
		message(STATUS ">>> Directory '${_SRC_DIR}' will be compiled into library '${_LIB_NAME}', "
			"setting up variables '${_VAR_NAME}_INCLUDE_DIR' and '${_VAR_NAME}_LIB'.")
		unset(_LIB_NAME)
	endif()
	unset(_VAR_NAME)
endforeach()


# build project-internal libraries
foreach(_SRC_DIR ${_SRC_DIRS})
	var_name_from_dir_name(${_SRC_DIR} _VAR_NAME)
	if(NOT DO_NOT_COMPILE_${_VAR_NAME})
		add_subdirectory("${_SRC_DIR}")
	endif()
	unset(_VAR_NAME)
endforeach()
unset(_SRC_DIRS)

# build tests
add_subdirectory(test)


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
