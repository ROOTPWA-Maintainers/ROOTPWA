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
#//      build file for RootPwa library
#//
#//
#// Author List:
#//      Sebastian Neubert    TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


message_setup_this_dir()


# set include directories
set(INCLUDE_DIRECTORIES
	${CMAKE_CURRENT_SOURCE_DIR}
	${RPWA_DECAYAMPLITUDE_INCLUDE_DIR}
	${PARTIALWAVEFIT_INCLUDE_DIR}
	${RPWA_PARTICLEDATA_INCLUDE_DIR}
	${RPWA_UTILITIES_INCLUDE_DIR}
	SYSTEM
	${Boost_INCLUDE_DIRS}
	${Libconfig_INCLUDE_DIR}
	${ROOT_INCLUDE_DIR}
	)
include_directories(${INCLUDE_DIRECTORIES})


# source files that are compiled into library
set(SOURCES
	massDepFit.cc
	massDepFitComponents.cc
	massDepFitLikeli.cc
	massDepFitModel.cc
	)


# library
set(THIS_LIB "RootPwaResonanceFit")
make_shared_library(
	"${THIS_LIB}"
	"${SOURCES}"
	"${Libconfig_LIBS}"
	"${ROOT_LIBS}"
	"${RPWA_PARTICLEDATA_LIB}"
	"${RPWA_UTILITIES_LIB}"
	"${PARTIALWAVEFIT_LIBS}"
	"${RPWA_DECAYAMPLITUDE_LIB}"
	)


make_executable(pwaMassFit pwaMassFit.cc ${THIS_LIB})
