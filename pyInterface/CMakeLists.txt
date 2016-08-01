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
#//      build file for amplitude library
#//
#//
#// Author List:
#//      Boris Grube          TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


message_setup_this_dir()


add_subdirectory(bindings)

# 'pyRootPwa' module
set(RPWA_PYTHON_PACKAGE_FILES
	__init__.py
	_amplitude.py
	_config.py
	_fileManager.py
	_fit.py
	_integrals.py
	_integralsOnTheFly.py
	_likelihood.py
	core.py
	utils/__init__.py
	utils/_binningMapFromArgList.py
	utils/_binning.py
	utils/_fileUtils.py
	utils/_fitTreeUtils.py
	utils/_printingUtils.py
	utils/_progressBar.py
	utils/_root.py
	utils/_silencer.py
	utils/_waveDescThresUtils.py
)

set(RPWA_PYTHON_PACKAGE_BYTECODES)
set(RPWA_PYTHON_PACKAGE_DIRS)
set(RPWA_PYTHON_PACKAGE_LINKS)
foreach(_FILE ${RPWA_PYTHON_PACKAGE_FILES})
	get_filename_component(_DIR "${_FILE}" DIRECTORY)
	# list of directories to create
	list(APPEND RPWA_PYTHON_PACKAGE_DIRS "${PYTHON_LIBRARY_OUTPUT_PATH}/pyRootPwa/${_DIR}")
	add_custom_command(
		OUTPUT ${PYTHON_LIBRARY_OUTPUT_PATH}/pyRootPwa/${_DIR}
		COMMAND ${CMAKE_COMMAND} -E make_directory ${PYTHON_LIBRARY_OUTPUT_PATH}/pyRootPwa/${_DIR}
	)
	# list of files to link
	list(APPEND RPWA_PYTHON_PACKAGE_LINKS "${PYTHON_LIBRARY_OUTPUT_PATH}/pyRootPwa/${_FILE}")
	add_custom_command(
		OUTPUT ${PYTHON_LIBRARY_OUTPUT_PATH}/pyRootPwa/${_FILE}
		COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/module/${_FILE} ${PYTHON_LIBRARY_OUTPUT_PATH}/pyRootPwa/${_FILE}
		DEPENDS pyRootPwaPackageDirs
	)
	# list of files containing bytecode
	list(APPEND RPWA_PYTHON_PACKAGE_BYTECODES "${PYTHON_LIBRARY_OUTPUT_PATH}/pyRootPwa/${_FILE}c")
	add_custom_command(
		OUTPUT ${PYTHON_LIBRARY_OUTPUT_PATH}/pyRootPwa/${_FILE}c
		COMMAND ${PYTHON_EXECUTABLE} -m py_compile ${PYTHON_LIBRARY_OUTPUT_PATH}/pyRootPwa/${_FILE}
		DEPENDS pyRootPwaPackageLinks module/${_FILE}
	)
endforeach()
unset(_DIR)
unset(_FILE)

# target to create the module directory structure
add_custom_target(
	pyRootPwaPackageDirs ALL
	DEPENDS ${RPWA_PYTHON_PACKAGE_DIRS}
)
# target to link the python files
add_custom_target(
	pyRootPwaPackageLinks ALL
	DEPENDS pyRootPwaPackageDirs ${RPWA_PYTHON_PACKAGE_LINKS}
)
# target to create python bytecode
add_custom_target(
	pyRootPwaPackageBytecodes ALL
	DEPENDS pyRootPwaPackageLinks ${RPWA_PYTHON_PACKAGE_BYTECODES}
)


# link the shared library with the python bindings
add_custom_target(
	libRootPwaPyLink ALL
	COMMAND ${CMAKE_COMMAND} -E create_symlink ${LIBRARY_OUTPUT_PATH}/${CMAKE_SHARED_LIBRARY_PREFIX}RootPwaPy${CMAKE_SHARED_LIBRARY_SUFFIX} ${PYTHON_LIBRARY_OUTPUT_PATH}/pyRootPwa/libRootPwaPy.so
	DEPENDS RootPwaPy pyRootPwaPackageDirs
)


# scripts
set(RPWA_PYTHON_SCRIPTS_FILES
	calcAmplitudes.py
	calcCovMatrixForFitResult.py
	calcIntegrals.py
	convertEventFile.py
	convertEvtToTree.py
	convertTreeToEvt.py
	createFileManager.py
	deWeight.py
	eigenvectorLikelihoodSlices.py
	genPseudoData.py
	generatePhaseSpace.py
	likelihoodPointCalculator.py
	plotAngles.py
	printMetadata.py
	pwaFit.py
	pwaNloptFit.py
)

set(RPWA_PYTHON_SCRIPTS_LINKS)
foreach(_FILE ${RPWA_PYTHON_SCRIPTS_FILES})
	get_filename_component(_FILENAME "${_FILE}" NAME_WE)
	# list of files to link
	list(APPEND RPWA_PYTHON_SCRIPTS_LINKS "${EXECUTABLE_OUTPUT_PATH}/${_FILENAME}")
	add_custom_command(
		OUTPUT ${EXECUTABLE_OUTPUT_PATH}/${_FILENAME}
		COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${_FILE} ${EXECUTABLE_OUTPUT_PATH}/${_FILENAME}
	)
endforeach()
unset(_FILE)
unset(_FILENAME)

# target to link the python files
add_custom_target(
	pyRootPwaScriptsLinks ALL
	DEPENDS ${RPWA_PYTHON_SCRIPTS_LINKS}
)

# pylint tests
if(PYLINT_EXECUTABLE)
	# check 'pyRootPwa' module
	add_test(pyLintPyRootPwa ${PYLINT_EXECUTABLE} --reports=n --rcfile=${CMAKE_CURRENT_SOURCE_DIR}/pylintrc pyRootPwa)

	# check scripts
	foreach(_FILE ${RPWA_PYTHON_SCRIPTS_FILES})
		get_filename_component(_FILENAME ${_FILE} NAME_WE)
		add_test(NAME pyLint${_FILENAME} COMMAND ${PYLINT_EXECUTABLE} --reports=n --rcfile=${CMAKE_CURRENT_SOURCE_DIR}/pylintrc ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${_FILE})
	endforeach()
	unset(_FILE)
	unset(_FILENAME)
endif()
