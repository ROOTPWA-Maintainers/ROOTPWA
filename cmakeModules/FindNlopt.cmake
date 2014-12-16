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
#//      cmake module for finding NLopt installation
#//      NLopt installation location is defined by environment variable $LIBCONFIG
#//
#//      following variables are defined:
#//      Nlopt_DIR              - NLopt installation directory
#//      Nlopt_INCLUDE_DIR      - NLopt header directory
#//      Nlopt_LIBRARY_DIR      - NLopt library directory
#//      Nlopt_LIBS             - NLopt library files
#//
#//      Example usage:
#//          find_package(Nlopt 1.4 REQUIRED)
#//
#//
#//-------------------------------------------------------------------------


set(Nlopt_FOUND        FALSE)
set(Nlopt_ERROR_REASON "")
set(Nlopt_DEFINITIONS  "")
set(Nlopt_LIBS)


set(Nlopt_DIR $ENV{NLOPT})
if(NOT Nlopt_DIR)

	set(Nlopt_FOUND TRUE)

	set(_Nlopt_LIB_NAMES "nlopt")
	find_library(Nlopt_LIBS
		NAMES ${_Nlopt_LIB_NAMES})
	if(NOT Nlopt_LIBS)
		set(Nlopt_FOUND FALSE)
		set(Nlopt_ERROR_REASON "${Nlopt_ERROR_REASON} Cannot find NLopt library '${_Nlopt_LIB_NAMES}'.")
	else()
		get_filename_component(Nlopt_DIR ${Nlopt_LIBS} PATH)
	endif()
	unset(_Nlopt_LIB_NAMES)

	set(_Nlopt_HEADER_FILE_NAME "NLopt.h")
	find_file(_Nlopt_HEADER_FILE
		NAMES ${_Nlopt_HEADER_FILE_NAME})
	if(NOT _Nlopt_HEADER_FILE)
		set(Nlopt_FOUND FALSE)
		set(Nlopt_ERROR_REASON "${Nlopt_ERROR_REASON} Cannot find NLopt header file '${_Nlopt_HEADER_FILE_NAME}'.")
	endif()
	unset(_Nlopt_HEADER_FILE_NAME)
	unset(_Nlopt_HEADER_FILE)

	if(NOT Nlopt_FOUND)
		set(Nlopt_ERROR_REASON "${Nlopt_ERROR_REASON} Nlopt not found in system directories (and environment variable NLOPT is not set).")
	endif()



else()

	set(Nlopt_FOUND TRUE)

	set(Nlopt_INCLUDE_DIR "${Nlopt_DIR}/include")
	if(NOT EXISTS "${Nlopt_INCLUDE_DIR}")
		set(Nlopt_FOUND FALSE)
		set(Nlopt_ERROR_REASON "${Nlopt_ERROR_REASON} Directory '${Nlopt_INCLUDE_DIR}' does not exist.")
	endif()

	set(Nlopt_LIBRARY_DIR "${Nlopt_DIR}/lib")
	if(NOT EXISTS "${Nlopt_LIBRARY_DIR}")
		set(Nlopt_FOUND FALSE)
		set(Nlopt_ERROR_REASON "${Nlopt_ERROR_REASON} Directory '${Nlopt_LIBRARY_DIR}' does not exist.")
	endif()

	set(_Nlopt_LIB_NAMES "nlopt_cxx")
	find_library(Nlopt_LIBS
		NAMES ${_Nlopt_LIB_NAMES}
		PATHS ${Nlopt_LIBRARY_DIR}
		NO_DEFAULT_PATH)
	if(NOT Nlopt_LIBS)
		set(Nlopt_FOUND FALSE)
		set(Nlopt_ERROR_REASON "${Nlopt_ERROR_REASON} Cannot find NLopt library '${_Nlopt_LIB_NAMES}' in '${Nlopt_LIBRARY_DIR}'.")
	endif()
	unset(_Nlopt_LIB_NAMES)

	set(_Nlopt_HEADER_FILE_NAME "nlopt.hpp")
	find_file(_Nlopt_HEADER_FILE
		NAMES ${_Nlopt_HEADER_FILE_NAME}
		PATHS ${Nlopt_INCLUDE_DIR}
		NO_DEFAULT_PATH)
	if(NOT _Nlopt_HEADER_FILE)
		set(Nlopt_FOUND FALSE)
		set(Nlopt_ERROR_REASON "${Nlopt_ERROR_REASON} Cannot find NLopt header file '${_Nlopt_HEADER_FILE_NAME}' in '${Nlopt_INCLUDE_DIR}'.")
	endif()
	unset(_Nlopt_HEADER_FILE_NAME)
	unset(_Nlopt_HEADER_FILE)

endif()


# make variables changeable
mark_as_advanced(
	Nlopt_INCLUDE_DIR
	Nlopt_LIBRARY_DIR
	Nlopt_LIBS
	Nlopt_DEFINITIONS
	)


# report result
if(Nlopt_FOUND)
	message(STATUS "Found NLopt in '${Nlopt_DIR}'.")
	message(STATUS "Using NLopt include directory '${Nlopt_INCLUDE_DIR}'.")
	message(STATUS "Using NLopt library '${Nlopt_LIBS}'.")
else()
	if(Nlopt_FIND_REQUIRED)
		message(FATAL_ERROR "Unable to find requested NLopt installation:${Nlopt_ERROR_REASON}")
	else()
		if(NOT Nlopt_FIND_QUIETLY)
			message(STATUS "Nlopt was not found:${Nlopt_ERROR_REASON}")
		endif()
	endif()
endif()
