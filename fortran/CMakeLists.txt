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
#//      build file for Fortran library
#//
#// If CMake does not find your fortran compiler you may have to set
#// the FC environment variable to your gfortran executable
#//
#// !!! important !!!
#//
#// Beware that ccmake seems to have problems parsing comment lines in
#// fixed-format sources, so it may pick up 'use' or 'module' keywords
#// from the comments that generate stray dependencies and leading to
#// make errors. One (ugly) workaround is to edit the comments in the
#// Fortran sources.
#//
#// As far as I can tell, the problem is not yet solved. See
#// http://public.kitware.com/Bug/view.php?id=8785
#//
#// Author List:
#//      Boris Grube          TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


message_setup_this_dir()


# enable Fortran support
enable_language(Fortran)


# add compiler options
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-automatic -ffixed-line-length-none")
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE} -g")
message(STATUS "Using Fortran compiler '${CMAKE_Fortran_COMPILER}'")
message(STATUS "Using Fortran compiler flags '${CMAKE_Fortran_FLAGS}'")
message(STATUS "Using additional Fortran compiler flags '${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}' "
	"for build type ${CMAKE_BUILD_TYPE}.")


# source files that are compiled into library
set(SOURCES
	bw_example.f  # example file; not used in any other place
	)


# library
lib_name_for_current_dir(THIS_LIB)
make_shared_library(
	"${THIS_LIB}"
	"${SOURCES}"
	)
