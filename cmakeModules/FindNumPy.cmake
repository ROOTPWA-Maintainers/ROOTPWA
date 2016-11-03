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
#//      cmake module for finding numpy libraries and include files
#//
#//      following variables are defined:
#//      NUMPY_FOUND              - numpy found
#//      NUMPY_INCLUDE_DIR        - include directory for numpy
#//      NUMPY_VERSION            - version of numpy
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


unset(NUMPY_FOUND)
unset(NUMPY_INCLUDE_DIR)
unset(NUMPY_VERSION)

# do not attempt to find numpy if Python was not found
if(PYTHON_EXECUTABLE)
	# get version
	execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import numpy; print numpy.__version__"
	                RESULT_VARIABLE _NUMPY_IMPORT_SUCCESS
	                OUTPUT_VARIABLE _NUMPY_VERSION_RAW
	                OUTPUT_STRIP_TRAILING_WHITESPACE
	                ERROR_QUIET
	               )
	if(_NUMPY_IMPORT_SUCCESS EQUAL 0)
		# version extracted successfully
		set(NUMPY_VERSION "${_NUMPY_VERSION_RAW}")

		# get include path
		execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import numpy; print numpy.get_include()"
		                OUTPUT_VARIABLE _NUMPY_INCLUDE_DIR_RAW
		                OUTPUT_STRIP_TRAILING_WHITESPACE
		                ERROR_QUIET
		               )

		# check that this directory really contains the include files
		set(_NUMPY_INCLUDE_FILE_NAME "numpy/numpyconfig.h")
		find_file(_NUMPY_INCLUDE_FILE
		          NAMES ${_NUMPY_INCLUDE_FILE_NAME}
		          PATHS ${_NUMPY_INCLUDE_DIR_RAW}
		          NO_DEFAULT_PATH)
		if(_NUMPY_INCLUDE_FILE)
			set(NUMPY_INCLUDE_DIR "${_NUMPY_INCLUDE_DIR_RAW}")
		else()
			message(WARNING "The numpy module was found, but the include directory '${_NUMPY_INCLUDE_DIR_RAW}' does not contain the required file '${_NUMPY_INCLUDE_FILE_NAME}'.")
		endif()

		unset(_NUMPY_INCLUDE_DIR_RAW)
		unset(_NUMPY_INCLUDE_FILE)
		unset(_NUMPY_INCLUDE_FILE_NAME)
	endif()
	unset(_NUMPY_IMPORT_SUCCESS)
	unset(_NUMPY_VERSION_RAW)
else()
	message(WARNING "Python needs to be set up prior to numpy.")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NUMPY
                                  REQUIRED_VARS NUMPY_INCLUDE_DIR NUMPY_VERSION
                                  VERSION_VAR NUMPY_VERSION
                                 )
