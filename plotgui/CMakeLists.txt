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
#//      build file for plotgui 
#//
#//
#// Author List:
#//      Stephan Schmeing          TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------



# set include directories
set(INCLUDE_DIRECTORIES
	${uprootedBoost_INCLUDE_DIRS}
	${Boost_INCLUDE_DIRS}
	${ROOT_INCLUDE_DIR}
	${ROOTPWA_INCLUDE_DIR}
	${UTILITIES_INCLUDE_DIR}
	${AMPLITUDE_INCLUDE_DIR}
	)
include_directories(${INCLUDE_DIRECTORIES})


# source files that are compiled into executable
set(SOURCES
	CompassPwaFileBase.cc
	CompassPwaFileFitResults.cc
	CompassPwaFileNormIntegrals.cc
	CompassPwaFilePhaseSpaceIntegrals.cc
	CompassPwaFileObject.cc
	CompassPwaFileLoader.cc
	)

# library
set(THIS_LIB "RootPwaPlotGui")
make_shared_library("${THIS_LIB}" "${SOURCES}"
	"${BOOST_LIBS}"
	"${ROOT_LIBS}"
	"${ROOTPWA_LIBS}"
	)

# executables
make_executable(TestCompassPwaFileparsers     main_TestCompassPwaFileparsers.cc   ${THIS_LIB}  )
make_executable(CompassPwaFilesToRoot     main_CompassPwaFilesToRootTree.cc   ${THIS_LIB}  )
