#///////////////////////////////////////////////////////////////////////////
#//
#//    Copyright 2009 Sebastian Neubert
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


# !!! important !!!
#
# Beware that ccmake seems to have problems parsing comment lines in
# fixed-format sources, so it may pick up 'use' or 'module' keywords
# from the comments that generate stray dependencies and leading to
# make errors. One (ugly) workaround is to edit comments in the
# Fortran sources.
#
# As far as I can tell the problem is not yet solved. see
# http://public.kitware.com/Bug/view.php?id=8785


PROJECT(rootpwa)


ENABLE_TESTING()


# enable Fortran support
ENABLE_LANGUAGE(Fortran)


# add compiler options
#SET(CMAKE_Fortran_FLAGS "-O2 -m32 -fno-automatic -ffixed-line-length-none -finit-local-zero")
SET(CMAKE_Fortran_FLAGS "-O2 -m32 -fno-automatic -ffixed-line-length-none")


# SET(INCLUDE_DIRECTORIES
#   $ENV{ROOTSYS}/include
#   $ENV{PWA2000}/include
#   $ENV{ROOTPWA}/src
#   $ENV{LIBCONFIG}/include
#   ${CMAKE_SOURCE_DIR}/src
#   ${CMAKE_SOURCE_DIR}/generators
# )

# include_directories(${INCLUDE_DIRECTORIES})

SET(SOURCES
  eps_k1.f
  eps_k1p.f
  eps_m.f
  eps_mx.f
)
   

# define ROOT libraries
# SET(ROOT_LIB $ENV{ROOTSYS}/lib)
# SET(ROOT_BIN $ENV{ROOTSYS}/bin)
# SET(ROOT_LIBS
#   ${ROOT_LIB}/libCore.so
#   ${ROOT_LIB}/libCint.so
#   ${ROOT_LIB}/libRIO.so
#   ${ROOT_LIB}/libNet.so
#   ${ROOT_LIB}/libHist.so
#   ${ROOT_LIB}/libGraf.so
#   ${ROOT_LIB}/libGraf3d.so
#   ${ROOT_LIB}/libGpad.so
#   ${ROOT_LIB}/libTree.so
#   ${ROOT_LIB}/libThread.so
#   ${ROOT_LIB}/libRint.so
#   ${ROOT_LIB}/libPostscript.so
#   ${ROOT_LIB}/libMatrix.so
#   ${ROOT_LIB}/libPhysics.so
#   ${ROOT_LIB}/libMathCore.so
#   ${ROOT_LIB}/libMinuit.so
#   libdl.so
#   libm.so
# )

# define PWA libraries
# SET(PWA_LIB $ENV{PWA2000}/lib)
# SET(PWA_LIBS
#   ${PWA_LIB}/libpp.a
#   ${PWA_LIB}/libpwa.a
# )


#-------------------------------------------------
#------- The library ---- ------------------------
#-------------------------------------------------

SET(SOURCES ${SOURCES})

ADD_LIBRARY(fortran SHARED ${SOURCES})

# SET_TARGET_PROPERTIES(fortran PROPERTIES LINKER_LANGUAGE Fortran)

# TARGET_LINK_LIBRARIES(fortran ${ROOT_LIBS})
# TARGET_LINK_LIBRARIES(fortran ${PWA_LIBS})
# TARGET_LINK_LIBRARIES(fortran ${ROOTPWA_LIBS})


#-------------------------------------------------
#------- Install options -------------------------
#-------------------------------------------------

INSTALL(FILES libfortran.so DESTINATION lib)
# INSTALL(FILES DESTINATION bin PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE)
