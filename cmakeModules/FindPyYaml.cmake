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
#//      cmake module for finding PyYaml libraries and include files
#//
#//      following variables are defined:
#//      PyYaml_FOUND       - indicates whether PyYaml was found
#//      PyYaml_VERSION     - version of PyYaml
#//
#//      Example usage:
#//          find_package(PyYaml 1.8 Optional)
#//
#//
#// Author List:
#//      Sebastian Uhl        TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


set(PyYaml_FOUND        TRUE)
set(PyYaml_ERROR_REASON "")

set(PyYaml_VERSION     NOTFOUND)


# check for Python
if(NOT PYTHON_FOUND)
	set(PyYaml_FOUND FALSE)
	set(PyYaml_ERROR_REASON "${PyYaml_ERROR_REASON} Did not find Python. Python needs to be set up prior to PyYaml.")
endif()
if(NOT PYTHON_EXECUTABLE)
	set(PyYaml_FOUND FALSE)
	set(PyYaml_ERROR_REASON "${PyYaml_ERROR_REASON} Did not find executable of Python interpreter. Python needs to be set up prior to PyYaml.")
endif()


if(PYTHON_EXECUTABLE)
	# get version
	execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import yaml; print yaml.__version__"
		 RESULT_VARIABLE _PyYaml_IMPORT_SUCCESS
		 OUTPUT_VARIABLE _PyYaml_VERSION_RAW
		 OUTPUT_STRIP_TRAILING_WHITESPACE
		 ERROR_QUIET)
	if(_PyYaml_IMPORT_SUCCESS EQUAL 0)
		# version extracted successfully
		set(PyYaml_VERSION "${_PyYaml_VERSION_RAW}")
	endif()
	unset(_PyYaml_IMPORT_SUCCESS)
	unset(_PyYaml_VERSION_RAW)
endif()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PyYaml
	FOUND_VAR PyYaml_FOUND
	REQUIRED_VARS PyYaml_VERSION
	VERSION_VAR PyYaml_VERSION)
# additional reporting
if(NOT PyYaml_FOUND)
	message(STATUS "Unable to find requested PyYaml installation:${PyYaml_ERROR_REASON}")
endif()


# hide variables from normal GUI
mark_as_advanced(
	PyYaml_VERSION
	)


if(NOT PyYaml_FOUND)
	unset(PyYaml_VERSION)
endif()
