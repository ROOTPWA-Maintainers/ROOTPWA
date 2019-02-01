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
#//      cmake module for finding Python 2.x installation
#//      replaces the official FindPYTHON that comes with Cmake which for several reasons is unusable
#//      requires executable of Python interpreter 'python' to be in PATH
#//
#//      following variables are defined:
#//      PYTHONINTERP_FOUND        - Was the Python executable found
#//      PYTHON_EXECUTABLE         - path to the Python interpreter
#//      PYTHON_VERSION_STRING     - Python version found e.g. 2.5.2
#//      PYTHON_VERSION_MAJOR      - Python major version found e.g. 2
#//      PYTHON_VERSION_MINOR      - Python minor version found e.g. 5
#//      PYTHON_VERSION_PATCH      - Python patch version found e.g. 2
#//      PYTHONLIBS_FOUND          - have the Python libs been found
#//      PYTHON_LIBRARIES          - path to the Python library
#//      PYTHON_INCLUDE_DIRS       - path to includes
#//      PYTHONLIBS_VERSION_STRING - version of the Python libs found
#//
#//      Example usage:
#//          find_package(PYTHON 2.7 REQUIRED)
#//
#//
#// Author List:
#//      Boris Grube          TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


set(PYTHON_FOUND        FALSE)
set(PYTHONINTERP_FOUND  FALSE)
set(PYTHONLIBS_FOUND    FALSE)
set(PYTHON_ERROR_REASON "")


# find Python interpreter
set(_PYTHON_EXEC_NAMES "python")  # allows to define search order for multiple python executable names
foreach(_PYTHON_EXEC_NAME ${_PYTHON_EXEC_NAMES})
	find_program(PYTHON_EXECUTABLE ${_PYTHON_EXEC_NAME})
	if(PYTHON_EXECUTABLE)
		break()
	endif()
endforeach()
if(NOT PYTHON_EXECUTABLE)
	set(PYTHON_ERROR_REASON "${PYTHON_ERROR_REASON} Cannot find executable of Python interpreter with name(s) '${_PYTHON_EXEC_NAMES}' in path. Make sure Python is setup correctly.")
else()
	set(PYTHONINTERP_FOUND TRUE)
endif()
unset(_PYTHON_EXEC_NAMES)


if(PYTHONINTERP_FOUND)

	# get version number
	execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print(sys.version_info[0])"
		OUTPUT_VARIABLE PYTHON_VERSION_MAJOR
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print(sys.version_info[1])"
		OUTPUT_VARIABLE PYTHON_VERSION_MINOR
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print(sys.version_info[2])"
		OUTPUT_VARIABLE PYTHON_VERSION_PATCH
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	set(PYTHON_VERSION_STRING "${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}.${PYTHON_VERSION_PATCH}")

	# protect against Python 3
	if(NOT PYTHON_VERSION_MAJOR EQUAL "2")
		set(PYTHONINTERP_FOUND FALSE)
		set(PYTHONLIBS_FOUND   FALSE)
		set(PYTHON_ERROR_REASON "${PYTHON_ERROR_REASON} Python 2 interpreter required. Found only version ${PYTHON_VERSION_STRING}.")
	endif()

	# compare version
	if(PYTHON_FIND_VERSION_EXACT)
		if(NOT PYTHON_VERSION_STRING VERSION_EQUAL PYTHON_FIND_VERSION)
			set(PYTHONINTERP_FOUND FALSE)
			set(PYTHONLIBS_FOUND   FALSE)
			set(PYTHON_ERROR_REASON "${PYTHON_ERROR_REASON} Python interpreter version ${PYTHON_VERSION_STRING} does not match requested version ${PYTHON_FIND_VERSION}.")
		endif()
	else()
		if(PYTHON_VERSION_STRING VERSION_LESS PYTHON_FIND_VERSION)
			set(PYTHONINTERP_FOUND FALSE)
			set(PYTHONLIBS_FOUND   FALSE)
			set(PYTHON_ERROR_REASON "${PYTHON_ERROR_REASON} Python interpreter version ${PYTHON_VERSION_STRING} is lower than requested version ${PYTHON_FIND_VERSION}.")
		endif()
	endif()

	if(PYTHON_VERSION_STRING VERSION_LESS 3.2)
		# get name of shared library
		execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print('python' + sysconfig.get_config_var('VERSION'))"
			OUTPUT_VARIABLE _PYTHON_LIBRARY_NAME
			OUTPUT_STRIP_TRAILING_WHITESPACE)
		# get path of shared library
		execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_config_var('LIBPL'))"
			OUTPUT_VARIABLE _PYTHON_LIBRARY_DIR
			OUTPUT_STRIP_TRAILING_WHITESPACE)
	else()
		# get name of shared library
		execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; import sysconfig; print('python' + sysconfig.get_config_var('VERSION') + sys.abiflags)"
			OUTPUT_VARIABLE _PYTHON_LIBRARY_NAME
			OUTPUT_STRIP_TRAILING_WHITESPACE)
		# get path of shared library
		execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_config_var('LIBPL'))"
			OUTPUT_VARIABLE _PYTHON_LIBRARY_DIR
			OUTPUT_STRIP_TRAILING_WHITESPACE)
	endif()
	# find shared library
	find_library(PYTHON_LIBRARIES
		NAMES ${_PYTHON_LIBRARY_NAME}
		PATHS ${_PYTHON_LIBRARY_DIR}
		NO_DEFAULT_PATH)
	if(NOT PYTHON_LIBRARIES)
		set(PYTHON_ERROR_REASON "${PYTHON_ERROR_REASON} Cannot find Python shared library '${_PYTHON_LIBRARY_NAME}' in '${_PYTHON_LIBRARY_DIR}'. Make sure Python is setup correctly.")
	else()
		set(PYTHONLIBS_FOUND TRUE)
	endif()
	unset(_PYTHON_LIBRARY_NAME)
	unset(_PYTHON_LIBRARY_DIR)


	if(PYTHONLIBS_FOUND)

		# get include directories
		if(PYTHON_VERSION_STRING VERSION_LESS 3.2)
			execute_process(
				COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print('{0};{1}'.format(sysconfig.get_python_inc(), sysconfig.get_python_inc(plat_specific=True)))"
				OUTPUT_VARIABLE PYTHON_INCLUDE_DIRS
				OUTPUT_STRIP_TRAILING_WHITESPACE)
		else()
			execute_process(
				COMMAND ${PYTHON_EXECUTABLE} -c "import sysconfig; print('{};{}'.format(sysconfig.get_path('include', 'posix_prefix'), sysconfig.get_path('platinclude', 'posix_prefix')))"
				OUTPUT_VARIABLE PYTHON_INCLUDE_DIRS
				OUTPUT_STRIP_TRAILING_WHITESPACE)
		endif()
		list(REMOVE_DUPLICATES PYTHON_INCLUDE_DIRS)
		foreach(_PYTHON_INCLUDE_DIR ${PYTHON_INCLUDE_DIRS})
			if(NOT EXISTS "${_PYTHON_INCLUDE_DIR}")
				set(PYTHONLIBS_FOUND FALSE)
				set(PYTHON_ERROR_REASON "${PYTHON_ERROR_REASON} Python include directory '${_PYTHON_INCLUDE_DIR}' does not exist.")
			else()
				# get library version from patchlevel.h
				if(EXISTS "${_PYTHON_INCLUDE_DIR}/patchlevel.h")
					# get line with version string
					file(STRINGS
						"${_PYTHON_INCLUDE_DIR}/patchlevel.h"
						_PYTHONLIBS_VERSION_STRING
						REGEX "^#define[ \t]+PY_VERSION[ \t]+\"[^\"]+\""
						)
					# parse version line
					string(REGEX REPLACE "^#define[ \t]+PY_VERSION[ \t]+\"([^\"]+)\".*" "\\1"
						PYTHONLIBS_VERSION_STRING "${_PYTHONLIBS_VERSION_STRING}")
					unset(_PYTHONLIBS_VERSION_STRING)
				endif()
			endif()
		endforeach()
		unset(_PYTHON_INCLUDE_DIR)
		if(NOT PYTHONLIBS_VERSION_STRING VERSION_EQUAL PYTHON_VERSION_STRING)
			set(PYTHONLIBS_FOUND FALSE)
			set(PYTHON_ERROR_REASON "${PYTHON_ERROR_REASON} Python library version ${PYTHONLIBS_VERSION_STRING} does not match Python interpreter version ${PYTHON_VERSION_STRING}.")
		endif()

	endif()

endif()  # PYTHONINTERP_FOUND


# make variables changeable
mark_as_advanced(
	PYTHON_LIBRARIES
	PYTHON_INCLUDE_DIRS
	)


# report result
if(PYTHONINTERP_FOUND)
	message(STATUS "Found Python interpreter version ${PYTHON_VERSION_STRING} at '${PYTHON_EXECUTABLE}'.")
else()
	if(PYTHON_FIND_REQUIRED)
		message(FATAL_ERROR "Unable to find requested Python interpreter:${PYTHON_ERROR_REASON}")
	else()
		if(NOT PYTHON_FIND_QUIETLY)
			message(STATUS "Python interpreter version ${PYTHON_FIND_VERSION}+ was not found:${PYTHON_ERROR_REASON}")
		endif()
	endif()
endif()
if(PYTHONLIBS_FOUND)
	message(STATUS "Found Python libraries version ${PYTHONLIBS_VERSION_STRING}.")
	message(STATUS "Using Python libraries '${PYTHON_LIBRARIES}'.")
	message(STATUS "Using Python include directory '${PYTHON_INCLUDE_DIRS}'.")
else()
	if(PYTHON_FIND_REQUIRED)
		message(FATAL_ERROR "Unable to find requested Python libraries:${PYTHON_ERROR_REASON}")
	else()
		if(NOT PYTHON_FIND_QUIETLY)
			message(STATUS "Python library version ${PYTHON_FIND_VERSION}+ was not found:${PYTHON_ERROR_REASON}")
		endif()
	endif()
endif()
if(PYTHONINTERP_FOUND AND PYTHONLIBS_FOUND)
	set(PYTHON_FOUND TRUE)
endif()
