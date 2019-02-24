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
#//      NLopt installation location is defined by environment variable $NLOPT
#//
#//      following variables are defined:
#//      NLopt_VERSION     - NLopt version
#//      NLopt_DIR         - NLopt installation directory
#//      NLopt_INCLUDE_DIR - NLopt header directory
#//      NLopt_LIBRARY_DIR - NLopt library directory
#//      NLopt_LIBS        - NLopt library files
#//
#//      Example usage:
#//          find_package(NLopt 1.4 REQUIRED)
#//
#//
#//-------------------------------------------------------------------------


set(NLopt_FOUND        TRUE)
set(NLopt_ERROR_REASON "")

set(NLopt_VERSION     NOTFOUND)
set(NLopt_DIR         NOTFOUND)
set(NLopt_INCLUDE_DIR NOTFOUND)
set(NLopt_LIBRARY_DIR NOTFOUND)
set(NLopt_LIBS        NOTFOUND)


# try to get the environment variable pointing to the NLopt installation
# directory
set(NLopt_DIR $ENV{NLOPT})

# find the library
set(_NLopt_LIBRARY_NAMES "nlopt" "nlopt_cxx")
if(NLopt_DIR)
	# search only in NLopt_DIR
	find_library(NLopt_LIBS
		NAMES ${_NLopt_LIBRARY_NAMES}
		PATHS "${NLopt_DIR}/lib"
		      "${NLopt_DIR}/build"
		NO_DEFAULT_PATH)
else()
	# search system-wide
	find_library(NLopt_LIBS
		NAMES ${_NLopt_LIBRARY_NAMES})
endif()
if(NLopt_LIBS)
	get_filename_component(NLopt_LIBRARY_DIR ${NLopt_LIBS} DIRECTORY)
else()
	set(NLopt_FOUND FALSE)
	if(NLopt_DIR)
		set(NLopt_ERROR_REASON "${NLopt_ERROR_REASON} Did not find NLopt library '${_NLopt_LIBRARY_NAMES}' in directories '${NLopt_DIR}/lib' or '${NLopt_DIR}/build'.")
	else()
		set(NLopt_ERROR_REASON "${NLopt_ERROR_REASON} Did not find NLopt library '${_NLopt_LIBRARY_NAMES}' in any standard library directory.")
	endif()
endif()
unset(_NLopt_LIBRARY_NAMES)


# find the include directory
set(_NLopt_HEADER_FILE_NAME "nlopt.h")
if(NLopt_DIR)
	# search only in NLopt_DIR
	find_path(NLopt_INCLUDE_DIR
		NAMES ${_NLopt_HEADER_FILE_NAME}
		PATHS "${NLopt_DIR}/include"
		      "${NLopt_DIR}/build/src/api"
		NO_DEFAULT_PATH)
else()
	# search system-wide
	find_path(NLopt_INCLUDE_DIR
		NAMES ${_NLopt_HEADER_FILE_NAME})
endif()
if(NOT NLopt_INCLUDE_DIR)
	set(NLopt_FOUND FALSE)
	if(NLopt_DIR)
		set(NLopt_ERROR_REASON "${NLopt_ERROR_REASON} Did not find NLopt include file '${_NLopt_HEADER_FILE_NAME}' in directories '${NLopt_DIR}/include' or '${NLopt_DIR}/build/src/api'.")
	else()
		set(NLopt_ERROR_REASON "${NLopt_ERROR_REASON} Did not find NLopt include file '${_NLopt_HEADER_FILE_NAME}' in any standard include directory.")
	endif()
endif()
unset(_NLopt_HEADER_FILE_NAME)


# try to get the version from the pkgconfig file
set(_NLopt_PKG_CONFIG_FILE_NAME "nlopt.pc")
if(NLopt_DIR)
	# search only in NLopt_DIR
	find_file(_NLopt_PC_FILE
		NAMES ${_NLopt_PKG_CONFIG_FILE_NAME}
		PATHS "${NLopt_DIR}/lib/pkgconfig"
		      "${NLopt_DIR}/build"
		NO_DEFAULT_PATH)
else()
	# search system-wide
	find_file(_NLopt_PC_FILE
		NAMES ${_NLopt_PKG_CONFIG_FILE_NAME}
		PATHS "${NLopt_LIBRARY_DIR}/pkgconfig"
		      "${NLopt_LIBRARY_DIR}")
endif()
if(_NLopt_PC_FILE)
	parse_version_from_pkg_config_file("${_NLopt_PC_FILE}" NLopt_VERSION)
endif()
unset(_NLopt_PKG_CONFIG_FILE_NAME)
unset(_NLopt_PC_FILE)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NLopt
	FOUND_VAR NLopt_FOUND
	REQUIRED_VARS NLopt_VERSION NLopt_INCLUDE_DIR NLopt_LIBRARY_DIR NLopt_LIBS
	VERSION_VAR NLopt_VERSION)
# additional reporting
if(NLopt_FOUND)
	message(STATUS "Using NLopt include directory '${NLopt_INCLUDE_DIR}'.")
	message(STATUS "Using NLopt library '${NLopt_LIBS}'.")
else()
	message(STATUS "Unable to find requested NLopt installation:${NLopt_ERROR_REASON}")
endif()


# hide variables from normal GUI
mark_as_advanced(
	NLopt_VERSION
	NLopt_DIR
	NLopt_INCLUDE_DIR
	NLopt_LIBRARY_DIR
	NLopt_LIBS
	)


if(NOT NLopt_FOUND)
	unset(NLopt_VERSION)
	unset(NLopt_DIR)
	unset(NLopt_INCLUDE_DIR)
	unset(NLopt_LIBRARY_DIR)
	unset(NLopt_LIBS)
endif()
