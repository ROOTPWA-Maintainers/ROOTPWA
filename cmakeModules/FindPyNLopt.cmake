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
#//      cmake module for finding nlopt libraries and include files
#//
#//      following variables are defined:
#//      PYNLOPT_FOUND              - nlopt found
#//      PYNLOPT_VERSION            - version of nlopt
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


set(PYNLOPT_FOUND FALSE)
unset(PYNLOPT_FOUND)
unset(PYNLOPT_VERSION)
unset(PYNLOPT_INCLUDE_DIR)

# do not attempt to find nlopt if Python was not found
if(PYTHON_EXECUTABLE)
	if(NLopt_INCLUDE_DIR)
		# get version
		execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import nlopt; print nlopt.__version__"
						RESULT_VARIABLE _PYNLOPT_IMPORT_SUCCESS
						OUTPUT_VARIABLE _PYNLOPT_VERSION_RAW
						OUTPUT_STRIP_TRAILING_WHITESPACE
						ERROR_QUIET
					   )
		if(_PYNLOPT_IMPORT_SUCCESS EQUAL 0)
			# version extracted successfully
			set(PYNLOPT_FOUND TRUE)
			set(PYNLOPT_VERSION "${_PYNLOPT_VERSION_RAW}")
			set(PYNLOPT_INCLUDE_DIR "${NLopt_INCLUDE_DIR}")

		endif()
		unset(_PYNLOPT_IMPORT_SUCCESS)
		unset(_PYNLOPT_VERSION_RAW)
	else()
		message(WARNING "NLopt needs to be set up prior to its Python interface")
	endif()
else()
	message(WARNING "Python needs to be set up prior to numpy.")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PYNLOPT
                                  REQUIRED_VARS PYNLOPT_INCLUDE_DIR PYNLOPT_VERSION
                                  VERSION_VAR PYNLOPT_VERSION
                                 )
