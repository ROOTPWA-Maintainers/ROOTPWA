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
#//      build file for Fortran library 
#//
#// !!! important !!!
#//
#// Beware that ccmake seems to have problems parsing comment lines in
#// fixed-format sources, so it may pick up 'use' or 'module' keywords
#// from the comments that generate stray dependencies and leading to
#// make errors. One (ugly) workaround is to edit comments in the
#// Fortran sources.
#//
#// As far as I can tell the problem is not yet solved. See
#// http://public.kitware.com/Bug/view.php?id=8785
#//
#//
#// Author List:
#//      Boris Grube          TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


# enable Fortran support
enable_language(Fortran)


# add compiler options
#set(CMAKE_Fortran_FLAGS "-O2 -m32 -fno-automatic -ffixed-line-length-none -finit-local-zero")
#set(CMAKE_Fortran_FLAGS "-O2 -m32 -fno-automatic -ffixed-line-length-none")
set(CMAKE_Fortran_FLAGS "-O2 -fno-automatic -ffixed-line-length-none")


# source files that are compiled into library
set(SOURCES
  eps_k1.f
  eps_k1p.f
  eps_m.f
  eps_mx.f
)
   

# library
set(THIS_LIB "Fortran")
make_shared_library("${THIS_LIB}" "${SOURCES}")
