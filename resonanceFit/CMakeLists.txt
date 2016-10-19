#///////////////////////////////////////////////////////////////////////////
#//
#//    Copyright 2010-2012 Sebastian Neubert (TUM)
#//    Copyright 2014 Sebastian Uhl (TUM)
#//
#//    This file is part of ROOTPWA
#//
#//    ROOTPWA is free software: you can redistribute it and/or modify
#//    it under the terms of the GNU General Public License as published by
#//    the Free Software Foundation, either version 3 of the License, or
#//    (at your option) any later version.
#//
#//    ROOTPWA is distributed in the hope that it will be useful,
#//    but WITHOUT ANY WARRANTY; without even the implied warranty of
#//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#//    GNU General Public License for more details.
#//
#//    You should have received a copy of the GNU General Public License
#//    along with ROOTPWA.  If not, see <http://www.gnu.org/licenses/>.
#//
#///////////////////////////////////////////////////////////////////////////
#//-------------------------------------------------------------------------
#//
#// Description:
#//      build file for RootPwaResonanceFit library
#//
#//-------------------------------------------------------------------------


message_setup_this_dir()


# set include directories
include_directories(
	${RPWA_DECAYAMPLITUDE_INCLUDE_DIR}
	${RPWA_PARTIALWAVEFIT_INCLUDE_DIR}
	${RPWA_PARTICLEDATA_INCLUDE_DIR}
	${RPWA_UTILITIES_INCLUDE_DIR}
	SYSTEM
	${Boost_INCLUDE_DIRS}
	${YamlCpp_INCLUDE_DIR}
	${ROOT_INCLUDE_DIR}
	)


# source files that are compiled into library
set(SOURCES
	massDepFit.cc
	cache.cc
	components.cc
	fsmd.cc
	function.cc
	minimizerRoot.cc
	model.cc
	parameters.cc
	)


# library
lib_name_for_current_dir(THIS_LIB)
make_shared_library(
	"${THIS_LIB}"
	"${SOURCES}"
	"${YamlCpp_LIBS}"
	"${ROOT_LIBS}"
	"${RPWA_PARTICLEDATA_LIB}"
	"${RPWA_UTILITIES_LIB}"
	"${RPWA_PARTIALWAVEFIT_LIB}"
	"${RPWA_DECAYAMPLITUDE_LIB}"
	)


make_executable(resonanceFit resonanceFitExe.cc ${THIS_LIB})
