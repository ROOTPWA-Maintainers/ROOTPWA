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
#//      cmake module for finding Libconfig installation
#//      Libconfig installation location is defined by environment variable $LIBCONFIG
#//
#//      following variables are defined:
#//      Libconfig_VERSION          - Libconfig version
#//      Libconfig_MAJOR_VERSION    - Libconfig major version
#//      Libconfig_MINOR_VERSION    - Libconfig minor version
#//      Libconfig_SUBMINOR_VERSION - Libconfig patch level
#//      Libconfig_DIR              - Libconfig installation directory
#//      Libconfig_INCLUDE_DIR      - Libconfig header directory
#//      Libconfig_LIBRARY_DIR      - Libconfig library directory
#//      Libconfig_LIBS             - Libconfig library files
#//
#//      Example usage:
#//          find_package(Libconfig 1.4 REQUIRED)
#//
#//
#//-------------------------------------------------------------------------


set(Libconfig_FOUND        TRUE)
set(Libconfig_ERROR_REASON "")

set(Libconfig_VERSION          NOTFOUND)
set(Libconfig_MAJOR_VERSION    NOTFOUND)
set(Libconfig_MINOR_VERSION    NOTFOUND)
set(Libconfig_SUBMINOR_VERSION NOTFOUND)
set(Libconfig_DIR              NOTFOUND)
set(Libconfig_INCLUDE_DIR      NOTFOUND)
set(Libconfig_LIBRARY_DIR      NOTFOUND)
set(Libconfig_LIBS             NOTFOUND)


# try to get the environment variable pointing to the Libconfig installation
# directory
set(Libconfig_DIR $ENV{LIBCONFIG})


# find the library
set(_Libconfig_LIBRARY_NAMES "config++")
if(Libconfig_DIR)
	# search only in Libconfig_DIR
	find_library(Libconfig_LIBS
		NAMES ${_Libconfig_LIBRARY_NAMES}
		PATHS "${Libconfig_DIR}/lib"
		NO_DEFAULT_PATH)
else()
	# search system-wide
	find_library(Libconfig_LIBS
		NAMES ${_Libconfig_LIBRARY_NAMES})
endif()
if(Libconfig_LIBS)
	get_filename_component(Libconfig_LIBRARY_DIR ${Libconfig_LIBS} DIRECTORY)
else()
	set(Libconfig_FOUND FALSE)
	if(Libconfig_DIR)
		set(Libconfig_ERROR_REASON "${Libconfig_ERROR_REASON} Did not find Libconfig library '${_Libconfig_LIBRARY_NAMES}' in directory '${Libconfig_DIR}/lib'.")
	else()
		set(Libconfig_ERROR_REASON "${Libconfig_ERROR_REASON} Did not find Libconfig library '${_Libconfig_LIBRARY_NAMES}' in any standard library directory.")
	endif()
endif()
unset(_Libconfig_LIBRARY_NAMES)


# find the include directory
set(_Libconfig_HEADER_FILE_NAME "libconfig.h++")
if(Libconfig_DIR)
	# search only in Libconfig_DIR
	find_path(Libconfig_INCLUDE_DIR
		NAMES ${_Libconfig_HEADER_FILE_NAME}
		PATHS "${Libconfig_DIR}/include"
		NO_DEFAULT_PATH)
else()
	# search system-wide
	find_path(Libconfig_INCLUDE_DIR
		NAMES ${_Libconfig_HEADER_FILE_NAME})
endif()
if(NOT Libconfig_INCLUDE_DIR)
	set(Libconfig_FOUND FALSE)
	if(Libconfig_DIR)
		set(Libconfig_ERROR_REASON "${Libconfig_ERROR_REASON} Did not find Libconfig include file '${_Libconfig_HEADER_FILE_NAME}' in directory '${Libconfig_DIR}/include'.")
	else()
		set(Libconfig_ERROR_REASON "${Libconfig_ERROR_REASON} Did not find Libconfig include file '${_Libconfig_HEADER_FILE_NAME}' in any standard include directory.")
	endif()
endif()


# get version from header file
if(Libconfig_DIR)
	# search only in Libconfig_DIR
	find_file(_Libconfig_HEADER_FILE
		NAMES ${_Libconfig_HEADER_FILE_NAME}
		PATHS "${Libconfig_DIR}/include"
		NO_DEFAULT_PATH)
else()
	# search system-wide
	find_file(_Libconfig_HEADER_FILE
		NAMES ${_Libconfig_HEADER_FILE_NAME})
endif()
if(_Libconfig_HEADER_FILE)
	# parse version string
	file(STRINGS ${_Libconfig_HEADER_FILE} _Libconfig_VERSION_LINES
		REGEX "LIBCONFIGXX_VER_(MAJOR|MINOR|REVISION)")
	list(LENGTH _Libconfig_VERSION_LINES _NMB_Libconfig_VERSION_LINES)
	if(NOT _NMB_Libconfig_VERSION_LINES EQUAL 3)
		set(Libconfig_FOUND FALSE)
		set(Libconfig_ERROR_REASON "${Libconfig_ERROR_REASON} Cannot determine Libconfig version: file '${_Libconfig_HEADER_FILE}' contains ${_NMB_Libconfig_VERSION_LINES} instead of 3 version lines.")
	else()
		string(REGEX REPLACE
			"[A-Za-z0-9_;# \t]*#define[ \t]+LIBCONFIGXX_VER_MAJOR[ \t]+([0-9]+)[A-Za-z0-9_;# \t]*"
			"\\1"	Libconfig_MAJOR_VERSION "${_Libconfig_VERSION_LINES}")
		string(REGEX REPLACE
			"[A-Za-z0-9_;# \t]*#define[ \t]+LIBCONFIGXX_VER_MINOR[ \t]+([0-9]+)[A-Za-z0-9_;# \t]*"
			"\\1" Libconfig_MINOR_VERSION "${_Libconfig_VERSION_LINES}")
		string(REGEX REPLACE
			"[A-Za-z0-9_;# \t]*#define[ \t]+LIBCONFIGXX_VER_REVISION[ \t]+([0-9]+)[A-Za-z0-9_;# \t]*"
			"\\1"	Libconfig_SUBMINOR_VERSION "${_Libconfig_VERSION_LINES}")
	endif()
	set(Libconfig_VERSION
		"${Libconfig_MAJOR_VERSION}.${Libconfig_MINOR_VERSION}.${Libconfig_SUBMINOR_VERSION}")
	unset(_Libconfig_VERSION_LINES)
	unset(_NMB_Libconfig_VERSION_LINES)
endif()
unset(_Libconfig_HEADER_FILE)
unset(_Libconfig_HEADER_FILE_NAME)


if(Libconfig_ERROR_REASON AND NOT Libconfig_FIND_QUIETLY)
	message(STATUS "Problems while finding the requested Libconfig installation:${Libconfig_ERROR_REASON}")
endif()
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Libconfig
	FOUND_VAR Libconfig_FOUND
	REQUIRED_VARS Libconfig_VERSION Libconfig_INCLUDE_DIR Libconfig_LIBRARY_DIR Libconfig_LIBS
	VERSION_VAR Libconfig_VERSION)
# additional reporting
if(Libconfig_FOUND AND NOT Libconfig_FIND_QUIETLY)
	message(STATUS "Using Libconfig include directory '${Libconfig_INCLUDE_DIR}'.")
	message(STATUS "Using Libconfig library '${Libconfig_LIBS}'.")
endif()


# hide variables from normal GUI
mark_as_advanced(
	Libconfig_VERSION
	Libconfig_MAJOR_VERSION
	Libconfig_MINOR_VERSION
	Libconfig_SUBMINOR_VERSION
	Libconfig_DIR
	Libconfig_INCLUDE_DIR
	Libconfig_LIBRARY_DIR
	Libconfig_LIBS
	)


if(NOT Libconfig_FOUND)
	unset(Libconfig_VERSION)
	unset(Libconfig_MAJOR_VERSION)
	unset(Libconfig_MINOR_VERSION)
	unset(Libconfig_SUBMINOR_VERSION)
	unset(Libconfig_DIR)
	unset(Libconfig_INCLUDE_DIR)
	unset(Libconfig_LIBRARY_DIR)
	unset(Libconfig_LIBS)
endif()
