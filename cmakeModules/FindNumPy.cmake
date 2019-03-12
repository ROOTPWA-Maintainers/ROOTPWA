#///////////////////////////////////////////////////////////////////////////
#//
#//    Copyright 2016
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
#//      cmake module for finding NumPy libraries and include files
#//
#//      following variables are defined:
#//      NumPy_FOUND       - indicates whether NumPy was found
#//      NumPy_VERSION     - version of NumPy
#//      NumPy_INCLUDE_DIR - include directory for NumPy
#//
#//      Example usage:
#//          find_package(NumPy 1.8 Optional)
#//
#//
#// Author List:
#//      Sebastian Uhl        TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


set(NumPy_FOUND        TRUE)
set(NumPy_ERROR_REASON "")

set(NumPy_VERSION     NOTFOUND)
set(NumPy_INCLUDE_DIR NOTFOUND)


# check for Python
if(NOT PYTHON_FOUND)
	set(NumPy_FOUND FALSE)
	set(NumPy_ERROR_REASON "${NumPy_ERROR_REASON} Did not find Python. Python needs to be set up prior to NumPy.")
endif()
if(NOT PYTHON_EXECUTABLE)
	set(NumPy_FOUND FALSE)
	set(NumPy_ERROR_REASON "${NumPy_ERROR_REASON} Did not find executable of Python interpreter. Python needs to be set up prior to NumPy.")
endif()


if(PYTHON_EXECUTABLE)
	# get version
	execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import numpy; print numpy.__version__"
		 RESULT_VARIABLE _NumPy_IMPORT_SUCCESS
		 OUTPUT_VARIABLE _NumPy_VERSION_RAW
		 OUTPUT_STRIP_TRAILING_WHITESPACE
		 ERROR_QUIET)
	if(_NumPy_IMPORT_SUCCESS EQUAL 0)
		# version extracted successfully
		set(NumPy_VERSION "${_NumPy_VERSION_RAW}")

		# get include path
		execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import numpy; print numpy.get_include()"
			 OUTPUT_VARIABLE _NumPy_INCLUDE_DIR_RAW
			 OUTPUT_STRIP_TRAILING_WHITESPACE
			 ERROR_QUIET)

		# check that this directory really contains the include files
		set(_NumPy_INCLUDE_FILE_NAME "numpy/numpyconfig.h")
		find_file(_NumPy_INCLUDE_FILE
			NAMES ${_NumPy_INCLUDE_FILE_NAME}
			PATHS "${_NumPy_INCLUDE_DIR_RAW}"
			NO_DEFAULT_PATH)
		if(_NumPy_INCLUDE_FILE)
			set(NumPy_INCLUDE_DIR "${_NumPy_INCLUDE_DIR_RAW}")
		else()
			set(NumPy_FOUND FALSE)
			set(NumPy_ERROR_REASON "${NumPy_ERROR_REASON} The NumPy module was found, but the include directory '${_NumPy_INCLUDE_DIR_RAW}' does not contain the required file '${_NumPy_INCLUDE_FILE_NAME}'.")
		endif()

		unset(_NumPy_INCLUDE_DIR_RAW)
		unset(_NumPy_INCLUDE_FILE_NAME)
		unset(_NumPy_INCLUDE_FILE)
	endif()
	unset(_NumPy_IMPORT_SUCCESS)
	unset(_NumPy_VERSION_RAW)
endif()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NumPy
	FOUND_VAR NumPy_FOUND
	REQUIRED_VARS NumPy_INCLUDE_DIR NumPy_VERSION
	VERSION_VAR NumPy_VERSION)
# additional reporting
if(NOT NumPy_FOUND)
	message(STATUS "Unable to find requested NumPy installation:${NumPy_ERROR_REASON}")
endif()


# hide variables from normal GUI
mark_as_advanced(
	NumPy_VERSION
	NumPy_INCLUDE_DIR
	)


if(NOT NumPy_FOUND)
	unset(NumPy_VERSION)
	unset(NumPy_INCLUDE_DIR)
endif()
