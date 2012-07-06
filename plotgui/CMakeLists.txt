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
#// File and Version Information:
#// $Rev::                             $: revision of last commit
#// $Author::                          $: author of last commit
#// $Date::                            $: date of last commit
#//
#// Description:
#//      build file for plotgui 
#//
#//
#// Author List:
#//      Stephan Schmeing          TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------

# set include directories
set(INCLUDE_DIRECTORIES
	${uprootedBoost_INCLUDE_DIRS}
	${Boost_INCLUDE_DIRS}
	${ROOT_INCLUDE_DIR}
	$ENV{QTROOTSYSDIR}/include
	${ROOTPWA_INCLUDE_DIR}
	${UTILITIES_INCLUDE_DIR}
	${AMPLITUDE_INCLUDE_DIR}
	${CMAKE_CURRENT_BINARY_DIR} #needed for Qt uis
	)
include_directories(${INCLUDE_DIRECTORIES})

# setup Qt
FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})

# set form files from Qt (.ui)
set(FORMS
	Qt-DesignerFiles/GuiPwaMain.ui
	)

# prozessing Qt Forms (.ui) with uic
QT4_WRAP_UI(FORMS_HEADERS ${FORMS})

# set header files that contain Q_OBJECT macro
set(HEADERS
	GuiPwaMain.h
	GuiWaveTreeModel.h
	)

# prozessing Q_OBJECT headers with moc
QT4_WRAP_CPP(HEADERS_MOC ${HEADERS})

# source files that are compiled into executable
set(SOURCES
	CompassPwaFileBase.cc
	CompassPwaFileFitResults.cc
	CompassPwaFileLoader.cc
	CompassPwaFileNormIntegrals.cc
	CompassPwaFileObject.cc
	CompassPwaFilePhaseSpaceIntegrals.cc
	GuiPwaMain.cc
	GuiWaveTreeModel.cc
	GuiWaveTreeModelItem.cc
	RootPwaDataObject.cc
	${HEADERS_MOC}
	${FORMS_HEADERS}
	)

# library
set(QTROOT_LIBDIR $ENV{QTROOTSYSDIR}/lib)
set(QTROOT_LIBS
	${QTROOT_LIBDIR}/libGQt.so
	${QTROOT_LIBDIR}/libQGLViewer.so
	${QTROOT_LIBDIR}/libQtGed.so
	${QTROOT_LIBDIR}/libQtImage.so
	${QTROOT_LIBDIR}/libQtRootGui.so
	${QTROOT_LIBDIR}/libQtRoot.so
	${QTROOT_LIBDIR}/libQtThread.so
	${QTROOT_LIBDIR}/libRQTGL.so
	)

set(THIS_LIB "RootPwaPlotGui")
make_shared_library("${THIS_LIB}" "${SOURCES}"
	"${BOOST_LIBS}"
	"${ROOT_LIBS}"
	"${QTROOT_LIBS}"
	"${ROOTPWA_LIBS}"
	"${QT_LIBRARIES}"
	)

# executables
make_executable(TestCompassPwaFileparsers	main_TestCompassPwaFileparsers.cc   ${THIS_LIB})
make_executable(CompassPwaFilesToRoot		main_CompassPwaFilesToRootTree.cc  	${THIS_LIB})
make_executable(GuiPwa						main_GuiPwa.cc						${THIS_LIB})
