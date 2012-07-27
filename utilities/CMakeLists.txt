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
#//      build file for utilities
#//
#//
#// Author List:
#//      Boris Grube          TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


# set include directories
set(INCLUDE_DIRECTORIES
	${Boost_INCLUDE_DIRS}
	${ROOT_INCLUDE_DIR}
	${PWA2000_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	)
include_directories(${INCLUDE_DIRECTORIES})


# executables
make_executable(testSumAccumulators testSumAccumulators.cc)
make_executable(testLibppFunctions  testLibppFunctions.cc "${ROOT_LIBS};${PWA2000_LIBS}")
make_executable(testSpinUtils       testSpinUtils.cc)
