
#include <iostream>

#include "environment.h"
#include "reportingUtils.hpp"
#include "reportingUtilsEnvironment.h"


using namespace rpwa;


	//////////////////////////////////////////////////////////////////////////////
	// functions to print version and compilation info

void rpwa::printGitHash()
{
	printInfo << "git repository hash at compile time was "
	          << "'" << gitHash() << "'" << std::endl;
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
	std::cout << "    Python version " << PYTHONLIBS_VERSION_STRING << " in '" << PYTHON_INCLUDE_DIRS << "'; using Boost.Python" << std::endl;
}


std::string rpwa::gitHash()
{
	return GIT_HASH;
}
