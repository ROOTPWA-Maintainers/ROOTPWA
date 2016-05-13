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


set(PYUTILS_SUBDIR pyUtils)

set(DECAYAMPLITUDE_SUBDIR decayAmplitude)
set(GENERATORS_SUBDIR generators)
set(HIGHLEVELINTERFACE_SUBDIR highLevelInterface)
set(NBODYPHASESPACE_SUBDIR nBodyPhaseSpace)
set(PARTICLEDATA_SUBDIR particleData)
set(PARTIALWAVEFIT_SUBDIR partialWaveFit)
set(STORAGEFORMATS_SUBDIR storageFormats)
set(UTILITIES_SUBDIR utilities)


# set include directories
include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
	${PYUTILS_SUBDIR}
	${RPWA_DECAYAMPLITUDE_INCLUDE_DIR}
	${DECAYAMPLITUDE_SUBDIR}
	${RPWA_GENERATORS_INCLUDE_DIR}
	${GENERATORS_SUBDIR}
	${RPWA_HIGHLEVELINTERFACE_INCLUDE_DIR}
	${HIGHLEVELINTERFACE_SUBDIR}
	${RPWA_NBODYPHASESPACE_INCLUDE_DIR}
	${NBODYPHASESPACE_SUBDIR}
	${RPWA_PARTIALWAVEFIT_INCLUDE_DIR}
	${PARTIALWAVEFIT_SUBDIR}
	${RPWA_PARTICLEDATA_INCLUDE_DIR}
	${PARTICLEDATA_SUBDIR}
	${RPWA_STORAGEFORMATS_INCLUDE_DIR}
	${STORAGEFORMATS_SUBDIR}
	${RPWA_UTILITIES_INCLUDE_DIR}
	${UTILITIES_SUBDIR}
	SYSTEM
	${Boost_INCLUDE_DIRS}
	${Libconfig_INCLUDE_DIR}
	${PYTHON_INCLUDE_DIRS}
	${ROOT_INCLUDE_DIR}
	)


# source files that are compiled into library
set(SOURCES
	rootPwaPy.cc
	${PYUTILS_SUBDIR}/rootConverters_py.cc
	${PYUTILS_SUBDIR}/stlContainers_py.cc
	${DECAYAMPLITUDE_SUBDIR}/ampIntegralMatrix_py.cc
	${DECAYAMPLITUDE_SUBDIR}/ampIntegralMatrixMetadata_py.cc
	${DECAYAMPLITUDE_SUBDIR}/decayTopology_py.cc
	${DECAYAMPLITUDE_SUBDIR}/diffractiveDissVertex_py.cc
	${DECAYAMPLITUDE_SUBDIR}/fsVertex_py.cc
	${DECAYAMPLITUDE_SUBDIR}/interactionVertex_py.cc
	${DECAYAMPLITUDE_SUBDIR}/isobarAmplitude_py.cc
	${DECAYAMPLITUDE_SUBDIR}/isobarCanonicalAmplitude_py.cc
	${DECAYAMPLITUDE_SUBDIR}/isobarHelicityAmplitude_py.cc
	${DECAYAMPLITUDE_SUBDIR}/isobarDecayTopology_py.cc
	${DECAYAMPLITUDE_SUBDIR}/isobarDecayVertex_py.cc
	${DECAYAMPLITUDE_SUBDIR}/massDependence_py.cc
	${DECAYAMPLITUDE_SUBDIR}/phaseSpaceIntegral_py.cc
	${DECAYAMPLITUDE_SUBDIR}/productionVertex_py.cc
	${DECAYAMPLITUDE_SUBDIR}/waveDescription_py.cc
	${GENERATORS_SUBDIR}/beamAndVertexGenerator_py.cc
	${GENERATORS_SUBDIR}/generator_py.cc
	${GENERATORS_SUBDIR}/generatorManager_py.cc
	${GENERATORS_SUBDIR}/generatorParameters_py.cc
	${GENERATORS_SUBDIR}/generatorPickerFunctions_py.cc
	${HIGHLEVELINTERFACE_SUBDIR}/calcAmplitude_py.cc
	${HIGHLEVELINTERFACE_SUBDIR}/pwaFit_py.cc
	${NBODYPHASESPACE_SUBDIR}/nBodyPhaseSpaceGenerator_py.cc
	${NBODYPHASESPACE_SUBDIR}/nBodyPhaseSpaceKinematics_py.cc
	${NBODYPHASESPACE_SUBDIR}/randomNumberGenerator_py.cc
	${PARTIALWAVEFIT_SUBDIR}/complexMatrix_py.cc
	${PARTIALWAVEFIT_SUBDIR}/fitResult_py.cc
	${PARTIALWAVEFIT_SUBDIR}/partialWaveFitHelper_py.cc
	${PARTIALWAVEFIT_SUBDIR}/pwaLikelihood_py.cc
	${PARTICLEDATA_SUBDIR}/particle_py.cc
	${PARTICLEDATA_SUBDIR}/particleDataTable_py.cc
	${PARTICLEDATA_SUBDIR}/particleProperties_py.cc
	${STORAGEFORMATS_SUBDIR}/amplitudeFileWriter_py.cc
	${STORAGEFORMATS_SUBDIR}/amplitudeMetadata_py.cc
	${STORAGEFORMATS_SUBDIR}/amplitudeTreeLeaf_py.cc
	${STORAGEFORMATS_SUBDIR}/eventFileWriter_py.cc
	${STORAGEFORMATS_SUBDIR}/eventMetadata_py.cc
	${STORAGEFORMATS_SUBDIR}/hashCalculator_py.cc
	${UTILITIES_SUBDIR}/physUtils_py.cc
	${UTILITIES_SUBDIR}/reportingUtilsEnvironment_py.cc
	${HIGHLEVELINTERFACE_SUBDIR}/getMassShapes_py.cc
	)
if(USE_NLOPT)
	LIST(APPEND SOURCES ${HIGHLEVELINTERFACE_SUBDIR}/pwaNloptFit_py.cc)
endif()


# library
set(THIS_LIB "RootPwaPy")
make_shared_library(
	"${THIS_LIB}"
	"${SOURCES}"
	"${PYTHON_LIBRARIES}"
	"${Boost_PYTHON_LIBRARY}"
	"${ROOT_LIBS}"
	"${Libconfig_LIBS}"
	"${RPWA_DECAYAMPLITUDE_LIB}"
	"${RPWA_GENERATORS_LIB}"
	"${RPWA_HIGHLEVELINTERFACE_LIB}"
	"${RPWA_NBODYPHASESPACE_LIB}"
	"${RPWA_PARTICLEDATA_LIB}"
	"${RPWA_STORAGEFORMATS_LIB}"
	"${RPWA_UTILITIES_LIB}"
	)
if(USE_NLOPT)
	target_link_libraries(${THIS_LIB} "${NLopt_LIBS}")
endif()


# 'pyRootPwa' package
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
		COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/package/${_FILE} ${PYTHON_LIBRARY_OUTPUT_PATH}/pyRootPwa/${_FILE}
		DEPENDS pyRootPwaPackageDirs
	)
	# list of files containing bytecode
	list(APPEND RPWA_PYTHON_PACKAGE_BYTECODES "${PYTHON_LIBRARY_OUTPUT_PATH}/pyRootPwa/${_FILE}c")
	add_custom_command(
		OUTPUT ${PYTHON_LIBRARY_OUTPUT_PATH}/pyRootPwa/${_FILE}c
		COMMAND ${PYTHON_EXECUTABLE} -m py_compile ${PYTHON_LIBRARY_OUTPUT_PATH}/pyRootPwa/${_FILE}
		DEPENDS pyRootPwaPackageLinks package/${_FILE}
	)
endforeach()
unset(_DIR)
unset(_FILE)

# target to create the package directory structure
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
#	amplitudeBenchmark.py    # do not link tests
#	compareAmplitudeFiles.py # do not link tests
	convertEvtToTree.py
	convertTreeToEvt.py
	deWeight.py
	genPseudoData.py
	genpw.py
#	integralBenchmark.py     # do not link tests
	plotAngles.py
#	testFileManager.py       # do not link tests
#	testPyRootPwa.py         # do not link tests
	userInterface/calcAmplitudes.py
	userInterface/calcCovMatrixForFitResult.py
	userInterface/calcIntegrals.py
	userInterface/createFileManager.py
	userInterface/eigenvectorLikelihoodSlices.py
	userInterface/likelihoodPointCalculator.py
	userInterface/pwaFit.py
	userInterface/pwaNloptFit.py
)

set(RPWA_PYTHON_SCRIPTS_LINKS)
foreach(_FILE ${RPWA_PYTHON_SCRIPTS_FILES})
	get_filename_component(_FILENAME "${_FILE}" NAME)
	# list of files to link
	list(APPEND RPWA_PYTHON_SCRIPTS_LINKS "${EXECUTABLE_OUTPUT_PATH}/${_FILENAME}")
	add_custom_command(
		OUTPUT ${EXECUTABLE_OUTPUT_PATH}/${_FILENAME}
		COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/${_FILE} ${EXECUTABLE_OUTPUT_PATH}/${_FILENAME}
	)
endforeach()
unset(_FILE)
unset(_FILENAME)

# target to link the python files
add_custom_target(
	pyRootPwaScriptsLinks ALL
	DEPENDS ${RPWA_PYTHON_SCRIPTS_LINKS}
)


# test functionality of 'pyRootPwa' package
add_test(testPyRootPwa ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/pyInterface/testPyRootPwa.py)


# pylint tests
if(PYLINT_EXECUTABLE)
	# check 'pyRootPwa' package
	add_test(pyLintPyRootPwa ${PYLINT_EXECUTABLE} --reports=n --rcfile=${CMAKE_CURRENT_SOURCE_DIR}/pylintrc pyRootPwa)

	# check scripts
	foreach(_FILE ${RPWA_PYTHON_SCRIPTS_FILES})
		get_filename_component(_FILENAME ${_FILE} NAME_WE)
		add_test(NAME pyLint${_FILENAME} COMMAND ${PYLINT_EXECUTABLE} --reports=n --rcfile=${CMAKE_CURRENT_SOURCE_DIR}/pylintrc ${CMAKE_CURRENT_SOURCE_DIR}/${_FILE})
	endforeach()
	unset(_FILE)
	unset(_FILENAME)
endif()
