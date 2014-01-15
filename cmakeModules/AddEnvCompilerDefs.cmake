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
#//      cmake macro that adds compiler definitions for some environment
#//      variables
#//
#//
#// Author List:
#//      Boris Grube          TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


# make some environment variables accessible via predefined macro variables
macro(add_env_compiler_defs)

	execute_process(COMMAND hostname
		OUTPUT_VARIABLE HOSTNAME
		RESULT_VARIABLE _HOSTNAME_RETURN
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	if(_HOSTNAME_RETURN)
		set(HOSTNAME "")
	endif()
	unset(_HOSTNAME_RETURN)
	if(DEBUG_OUTPUT)
		message(STATUS "Adding definitions for internal variables:")
	endif()
	foreach(_CMAKEVAR
			"CMAKE_HOST_SYSTEM_NAME"
			"CMAKE_HOST_SYSTEM_PROCESSOR"
			"CMAKE_HOST_SYSTEM_VERSION"
			"NMB_CPU_CORES"
			"HOSTNAME"
			"CMAKE_SOURCE_DIR"
			"CMAKE_BUILD_TYPE"
			"Boost_LIBRARY_VERSION"
			"Boost_INCLUDE_DIRS"
			"Libconfig_VERSION"
			"Libconfig_DIR"
			"ROOTSYS"
			"CUDA_VERSION"
			"CUDA_LIB_DIRS"
			"PYTHONLIBS_VERSION_STRING"
			"PYTHON_INCLUDE_DIRS"
			)
		if(_CMAKEVAR)
			if(DEBUG_OUTPUT)
				message(STATUS "        ${_CMAKEVAR} = ${${_CMAKEVAR}}")
			endif()
			add_definitions(-D'${_CMAKEVAR}=\"${${_CMAKEVAR}}\"')
		else()
			add_definitions(-D'${_CMAKEVAR}=\"\"')
		endif()
	endforeach()
	unset(_CMAKEVAR)
	if(DEBUG_OUTPUT)
		message(STATUS "Adding definitions for environment variables:")
	endif()
	foreach(_ENVVAR
			"USER"
			)
		set(_ENVVARVAL $ENV{${_ENVVAR}})
		if(_ENVVARVAL)
			if(DEBUG_OUTPUT)
				message(STATUS "        ${_ENVVAR} = ${_ENVVARVAL}")
			endif()
			add_definitions(-D'${_ENVVAR}=\"${_ENVVARVAL}\"')
		else()
			add_definitions(-D'${_ENVVAR}=\"\"')
		endif()
		unset(_ENVVARVAL)
	endforeach()
	unset(_ENVVAR)

endmacro(add_env_compiler_defs)
