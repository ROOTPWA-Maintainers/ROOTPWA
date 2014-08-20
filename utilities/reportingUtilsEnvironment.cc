
#include <iostream>

#include "reportingUtils.hpp"
#include "reportingUtilsEnvironment.h"


using namespace rpwa;


	//////////////////////////////////////////////////////////////////////////////
	// functions to print version and compilation info

	// check macro variables set by Makefile
#ifndef CMAKE_HOST_SYSTEM_NAME
#define CMAKE_HOST_SYSTEM_NAME "undefined"
#endif
#ifndef CMAKE_HOST_SYSTEM_PROCESSOR
#define CMAKE_HOST_SYSTEM_PROCESSOR "undefined"
#endif
#ifndef CMAKE_HOST_SYSTEM_VERSION
#define CMAKE_HOST_SYSTEM_VERSION "undefined"
#endif
#ifndef NMB_CPU_CORES
#define NMB_CPU_CORES "undefined"
#endif
#ifndef HOSTNAME
#define HOSTNAME "undefined"
#endif
#ifndef USER
#define USER "undefined"
#endif
#ifndef CMAKE_SOURCE_DIR
#define CMAKE_SOURCE_DIR "undefined"
#endif
#ifndef CMAKE_BUILD_TYPE
#define CMAKE_BUILD_TYPE "undefined"
#endif
#ifndef GIT_HASH
#define GIT_HASH "undefined"
#endif
#ifndef Boost_LIBRARY_VERSION
#define Boost_LIBRARY_VERSION "undefined"
#endif
#ifndef Boost_INCLUDE_DIRS
#define Boost_INCLUDE_DIRS "undefined"
#endif
#ifndef Libconfig_VERSION
#define Libconfig_VERSION "undefined"
#endif
#ifndef Libconfig_DIR
#define Libconfig_DIR "undefined"
#endif
#ifndef ROOTSYS
#define ROOTSYS "undefined"
#endif
#ifndef CUDA_VERSION
#define CUDA_VERSION "undefined"
#endif
#ifndef CUDA_LIB_DIRS
#define CUDA_LIB_DIRS "undefined"
#endif
#ifndef PYTHONLIBS_VERSION_STRING
#define PYTHONLIBS_VERSION_STRING "undefined"
#endif
#ifndef PYTHON_INCLUDE_DIRS
#define PYTHON_INCLUDE_DIRS "undefined"
#endif


void rpwa::printGitHash()
{
	printInfo << "git repository hash at compile time was "
	          << "'" << GIT_HASH << "'" << std::endl;
}


void rpwa::printCompilerInfo()
{
	printInfo << "this executable was compiled" << std::endl
	          << "    on a '" << CMAKE_HOST_SYSTEM_NAME << " " << CMAKE_HOST_SYSTEM_VERSION << " "
	          << CMAKE_HOST_SYSTEM_PROCESSOR << "' machine" << std::endl
	          << "    in directory '" << USER << "@" << HOSTNAME ":"
	          << CMAKE_SOURCE_DIR << "'" << std::endl
	          << "    on " << __DATE__ << " " << __TIME__
	          << " by compiler version " << __VERSION__ << ","
	          << " build type is '" << CMAKE_BUILD_TYPE << "'" << std::endl;
}


void rpwa::printLibraryInfo()
{
	printInfo << "this project was linked against" << std::endl
	          << "    Boost version " << Boost_LIBRARY_VERSION << " in '" << Boost_INCLUDE_DIRS << "'"  << std::endl
	          << "    libConfig version " << Libconfig_VERSION << " in '" << Libconfig_DIR << "'" << std::endl
	          << "    ROOT version " << ROOT_RELEASE << " in '" << ROOTSYS << "'"  << std::endl;
#ifdef USE_CUDA
	std::cout << "    CUDA version " << CUDA_VERSION << " in '" << CUDA_LIB_DIRS << "'" << std::endl;
#endif
#ifdef USE_MPI
	std::cout << "    MPI libraries; using Boost.MPI" << std::endl;
#endif
#ifdef USE_PYTHON
	std::cout << "    Python version " << PYTHONLIBS_VERSION_STRING << " in '" << PYTHON_INCLUDE_DIRS << "'; using Boost.Python" << std::endl;
#endif
}
