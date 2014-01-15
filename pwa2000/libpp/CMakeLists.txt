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
#//      build file for pp library
#//
#//
#// Author List:
#//      Sebastian Neubert    TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


message(STATUS "")
message(STATUS ">>> Setting up 'libpp' directory.")


set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -O2 -Wall -Woverloaded-virtual -Werror")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -Wall -Woverloaded-virtual -Werror")


# set include directories
set(INCLUDE_DIRECTORIES
	${CMAKE_CURRENT_SOURCE_DIR}
	)
include_directories(${INCLUDE_DIRECTORIES})


# source files that are compiled into library
set(SOURCES
	event.cc
	massDep.cc
	particle.cc
	pputil.cc
	lorentz.cc
	matrix.cc
	particleData.cc
	Vec.cc
	)


# library
set(THIS_LIB "pp")
make_shared_library(
	"${THIS_LIB}"
	"${SOURCES}"
	)
