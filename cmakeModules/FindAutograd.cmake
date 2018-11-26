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
#//      cmake module for finding autograd library
#//
#//      following variables are defined:
#//      AUTOGRAD_FOUND              - autograd found
#//      AUTOGRAD_VERSION            - version of autograd
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


set(AUTOGRAD_FOUND FALSE)
unset(AUTOGRAD_FOUND)
unset(AUTOGRAD_VERSION)

# do not attempt to find autograd if Python was not found
if(PYTHON_EXECUTABLE)
	# get version
	execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import autograd; print '1.2'"
					RESULT_VARIABLE _AUTOGRAD_IMPORT_SUCCESS
					OUTPUT_VARIABLE _AUTOGRAD_VERSION_RAW
					OUTPUT_STRIP_TRAILING_WHITESPACE
					ERROR_QUIET
				   )
	if(_AUTOGRAD_IMPORT_SUCCESS EQUAL 0)
		# version extracted successfully
		set(AUTOGRAD_FOUND TRUE)
		set(AUTOGRAD_VERSION "${_AUTOGRAD_VERSION_RAW}")

	endif()
	unset(_AUTOGRAD_IMPORT_SUCCESS)
	unset(_AUTOGRAD_VERSION_RAW)
else()
	message(WARNING "Python needs to be set up prior to autograd.")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AUTOGRAD
                                  REQUIRED_VARS AUTOGRAD_VERSION
                                  VERSION_VAR AUTOGRAD_VERSION
                                 )
