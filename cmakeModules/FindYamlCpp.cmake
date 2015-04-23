#///////////////////////////////////////////////////////////////////////////
#//
#//    Copyright 2014
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
#//      cmake module for finding yaml-cpp installation
#//      yaml-cpp installation location is defined by environment variable
#//      $YAML_CPP
#//
#//      following variables are defined:
#//      YamlCpp_DIR              - yaml-cpp installation directory
#//      YamlCpp_INCLUDE_DIR      - yaml-cpp header directory
#//      YamlCpp_LIBS             - yaml-cpp library files
#//      YamlCpp_VERSION          - yaml-cpp version
#//
#//      Example usage:
#//          find_package(YamlCpp 0.5 REQUIRED)
#//
#//
#//-------------------------------------------------------------------------


set(YamlCpp_FOUND        True)
set(YamlCpp_ERROR_REASON "")

set(YamlCpp_DIR)
set(YamlCpp_INCLUDE_DIR)
set(YamlCpp_LIBS)
set(YamlCpp_VERSION)


# try to get the environment variable pointing to the yaml-cpp installation
# directory
set(YamlCpp_DIR $ENV{YAML_CPP})


# find the include directory
set(_YamlCpp_HEADER_FILE_NAME "yaml-cpp/yaml.h")
if(YamlCpp_DIR)
	find_path(YamlCpp_INCLUDE_DIR
		NAMES ${_YamlCpp_HEADER_FILE_NAME}
		PATHS ${YamlCpp_DIR}/include
		NO_DEFAULT_PATH
		)
else()
	find_path(YamlCpp_INCLUDE_DIR
		NAMES ${_YamlCpp_HEADER_FILE_NAME}
		)
endif()

if(NOT YamlCpp_INCLUDE_DIR)
	set(YamlCpp_FOUND False)
	if(YamlCpp_DIR)
		set(YamlCpp_ERROR_REASON "${YamlCpp_ERROR_REASON} Did not find yaml-cpp include file '${_YamlCpp_HEADER_FILE_NAME}' in directory '${YamlCpp_DIR}/include'.")
	else()
		set(YamlCpp_ERROR_REASON "${YamlCpp_ERROR_REASON} Did not find yaml-cpp include file '${_YamlCpp_HEADER_FILE_NAME}' in any standard include directory.")
	endif()
endif()
unset(_YamlCpp_HEADER_FILE_NAME)


# find the library
set(_YamlCpp_LIBRARY_NAME yaml-cpp)
if(YamlCpp_DIR)
	find_library(YamlCpp_LIBS
		NAMES ${_YamlCpp_LIBRARY_NAME}
		PATHS ${YamlCpp_DIR}/lib
		      ${YamlCpp_DIR}/build
		NO_DEFAULT_PATH
		)
else()
	find_library(YamlCpp_LIBS
		NAMES ${_YamlCpp_LIBRARY_NAME}
		)
endif()

if(NOT YamlCpp_LIBS)
	set(YamlCpp_FOUND False)
	if(YamlCpp_DIR)
		set(YamlCpp_ERROR_REASON "${YamlCpp_ERROR_REASON} Did not find yaml-cpp library '${_YamlCpp_LIBRARY_NAME}' in directory '${YamlCpp_DIR}/lib' or '${YamlCpp_DIR}/build'.")
	else()
		set(YamlCpp_ERROR_REASON "${YamlCpp_ERROR_REASON} Did not find yaml-cpp library '${_YamlCpp_LIBRARY_NAME}' in any standard library directory.")
	endif()
endif()
unset(_YamlCpp_LIBRARY_NAME)


# check the version
if(NOT YamlCpp_VERSION)
	if(YamlCpp_FIND_VERSION_EXACT)
		set(YamlCpp_FOUND False)
		set(YamlCpp_ERROR_REASON "${YamlCpp_ERROR_REASON} Exact version of yaml-cpp required, but the version could not be extracted.")
	else()
		set(YamlCpp_VERSION "unknown")
		message(STATUS "Could not extract version for yaml-cpp, let us hope it still works.")
	endif()
else()
	if(YamlCpp_FIND_VERSION_EXACT)
		if(NOT YamlCpp_VERSION VERSION_EQUAL YamlCpp_FIND_VERSION)
			set(YamlCpp_FOUND False)
			set(YamlCpp_ERROR_REASON "${YamlCpp_ERROR_REASON} Did not find version ${YamlCpp_FIND_VERSION} of yaml-cpp.")
		endif()
	else()
		if(YamlCpp_VERSION VERSION_LESS YamlCpp_FIND_VERSION)
			set(YamlCpp_FOUND False)
			set(YamlCpp_ERROR_REASON "${YamlCpp_ERROR_REASON} Did not find version ${YamlCpp_FIND_VERSION} or later of yaml-cpp.")
		endif()
	endif()
endif()


# report result
if(YamlCpp_FOUND)
	message(STATUS "Found yaml-cpp version ${YamlCpp_VERSION} in '${YamlCpp_DIR}'.")
	message(STATUS "Using yaml-cpp include directory '${YamlCpp_INCLUDE_DIR}'.")
	message(STATUS "Using yaml-cpp library '${YamlCpp_LIBS}'.")
else()
	unset(YamlCpp_DIR)
	unset(YamlCpp_INCLUDE_DIR)
	unset(YamlCpp_LIBS)
	unset(YamlCpp_VERSION)

	if(YamlCpp_FIND_REQUIRED)
		message(FATAL_ERROR "Unable to find requested yaml-cpp installation:${YamlCpp_ERROR_REASON}")
	else()
		if(NOT YamlCpp_FIND_QUIETLY)
			if(YamlCpp_FIND_VERSION_EXACT)
				message(STATUS "yaml-cpp version ${YamlCpp_FIND_VERSION} was not found:${YamlCpp_ERROR_REASON}")
			else()
				message(STATUS "yaml-cpp version ${YamlCpp_FIND_VERSION} or later was not found:${YamlCpp_ERROR_REASON}")
			endif()
		endif()
	endif()
endif()


# make variables changeable
mark_as_advanced(
	YamlCpp_DIR
	YamlCpp_INCLUDE_DIR
	YamlCpp_LIBS
	)
