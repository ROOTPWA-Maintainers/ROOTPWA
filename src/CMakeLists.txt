#///////////////////////////////////////////////////////////////////////////
#//
#//    Copyright 2009 Sebastian Neubert
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


PROJECT(rootpwa)

ENABLE_TESTING()

INCLUDE_DIRECTORIES($ENV{ROOTSYS}/include
	            $ENV{PWA2000}/include
		    .)

SET(SOURCES TFitBin.cc
	    TCMatrix.cc
	    TPWALikelihood.h
	    TPWALikelihoodC.cc
	    TPWAAmp.cc
	    TPWAPhaseConstraint.cc
	    TLogMultiGaus.cc
	    TCovEllipse.cc
	    mcmc/TMCMCMeta.cc
	    utilities.cc
	    rootpwaDict.cc)
	    
SET(MYFLAGS "-g -O3 -Wall -Werror")

SET(CMAKE_CXX_FLAGS ${MYFLAGS})
#SET(CMAKE_LD_FLAGS "${CMAKE_LD_FLAGS} -g")

ADD_LIBRARY(rootpwa SHARED ${SOURCES})

#add ROOT lib
SET(ROOT_LIB $ENV{ROOTSYS}/lib)
SET(ROOT_BIN $ENV{ROOTSYS}/bin)
SET(ROOT_LIBS ${ROOT_LIB}/libCore.so 
	      ${ROOT_LIB}/libCint.so
	      ${ROOT_LIB}/libRIO.so		
	      ${ROOT_LIB}/libNet.so
	      ${ROOT_LIB}/libHist.so
	      ${ROOT_LIB}/libGraf.so
	      ${ROOT_LIB}/libGraf3d.so 
	      ${ROOT_LIB}/libGpad.so
	      ${ROOT_LIB}/libTree.so
	      ${ROOT_LIB}/libRint.so
	      ${ROOT_LIB}/libPostscript.so
	      ${ROOT_LIB}/libMatrix.so 
	      ${ROOT_LIB}/libPhysics.so
	      ${ROOT_LIB}/libMathCore.so
	      ${ROOT_LIB}/libMinuit2.so
	      libdl.so
	      libm.so)

TARGET_LINK_LIBRARIES(rootpwa ${ROOT_LIBS})



#add PWA lib
SET(PWA_LIB $ENV{PWA2000}/lib)
SET(PWA_LIBS ${PWA_LIB}/libpp.a
 	     ${PWA_LIB}/libpwa.a)


#-------------------------------------------------
#------ Root Dictionary --------------------------
#-------------------------------------------------

FIND_PROGRAM(ROOT_CINT_EXECUTABLE
  NAMES rootcint
  PATHS ${ROOT_BIN}
  NO_DEFAULT_PATH
 )

ADD_CUSTOM_COMMAND(OUTPUT rootpwaDict.h rootpwaDict.cc
      COMMAND ${ROOT_CINT_EXECUTABLE}
      ARGS -f rootpwaDict.cc -c -DHAVE_CONFIG_H ${INCLUDE_DIRECTORIES} TFitBin.h TCMatrix.h TCovEllipse.h mcmc/TMCMCMeta.h rootpwaLinkDef.h)


#-------------------------------------------------
#------- Main executables ------------------------
#-------------------------------------------------

ADD_EXECUTABLE(listwaves listwaves.cc)
TARGET_LINK_LIBRARIES(listwaves rootpwa)
TARGET_LINK_LIBRARIES(listwaves ${PWA_LIBS})
TARGET_LINK_LIBRARIES(listwaves ${ROOT_LIBS})

ADD_EXECUTABLE(pwafit pwafit.cc)
TARGET_LINK_LIBRARIES(pwafit rootpwa)
TARGET_LINK_LIBRARIES(pwafit ${PWA_LIBS})
TARGET_LINK_LIBRARIES(pwafit ${ROOT_LIBS})

ADD_EXECUTABLE(pwaConstraintFit pwaConstraintFit.cc)
TARGET_LINK_LIBRARIES(pwaConstraintFit rootpwa)
TARGET_LINK_LIBRARIES(pwaConstraintFit ${PWA_LIBS})
TARGET_LINK_LIBRARIES(pwaConstraintFit ${ROOT_LIBS})

ADD_EXECUTABLE(testPWALikeli testPWALikeli.cc)
TARGET_LINK_LIBRARIES(testPWALikeli rootpwa)
TARGET_LINK_LIBRARIES(testPWALikeli ${PWA_LIBS})
TARGET_LINK_LIBRARIES(testPWALikeli ${ROOT_LIBS})

ADD_EXECUTABLE(testPWAConstraint testPWAConstraint.cc)
TARGET_LINK_LIBRARIES(testPWAConstraint rootpwa)
TARGET_LINK_LIBRARIES(testPWAConstraint ${PWA_LIBS})
TARGET_LINK_LIBRARIES(testPWAConstraint ${ROOT_LIBS})

ADD_EXECUTABLE(testPWAAmp testPWAAmp.cc)
TARGET_LINK_LIBRARIES(testPWAAmp rootpwa)
TARGET_LINK_LIBRARIES(testPWAAmp ${PWA_LIBS})
TARGET_LINK_LIBRARIES(testPWAAmp ${ROOT_LIBS})

#ADD_EXECUTABLE(hamMCpwa mcmc/hamMCpwa.cc)
#TARGET_LINK_LIBRARIES(hamMCpwa rootpwa)
#TARGET_LINK_LIBRARIES(hamMCpwa ${PWA_LIBS})
#TARGET_LINK_LIBRARIES(hamMCpwa ${ROOT_LIBS})

#ADD_EXECUTABLE(hamMCtest mcmc/hamMCtest.cc)
#TARGET_LINK_LIBRARIES(hamMCtest rootpwa)
#TARGET_LINK_LIBRARIES(hamMCtest ${PWA_LIBS})
#TARGET_LINK_LIBRARIES(hamMCtest ${ROOT_LIBS})

#ADD_EXECUTABLE(mcmc2fitbin mcmc/mcmc2fitbin.cc)
#TARGET_LINK_LIBRARIES(mcmc2fitbin rootpwa)
#TARGET_LINK_LIBRARIES(mcmc2fitbin ${PWA_LIBS})
#TARGET_LINK_LIBRARIES(mcmc2fitbin ${ROOT_LIBS})

ADD_EXECUTABLE(selector selector.cc)
TARGET_LINK_LIBRARIES(selector rootpwa)
TARGET_LINK_LIBRARIES(selector ${PWA_LIBS})
TARGET_LINK_LIBRARIES(selector ${ROOT_LIBS})


#-------------------------------------------------
#------- Install Commands ------------------------
#-------------------------------------------------

INSTALL(FILES librootpwa.so DESTINATION lib)
INSTALL(FILES pwafit pwaConstraintFit selector DESTINATION bin PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE)


#-------------------------------------------------
#------- Test Executables ------------------------
#-------------------------------------------------

ADD_TEST(testBuild     ${CMAKE_BUILD_TOOL})
ADD_TEST(testPWAAmp    ${EXECUTABLE_OUTPUT_PATH}/testPWAAmp)
ADD_TEST(testPWALikeli ${EXECUTABLE_OUTPUT_PATH}/testPWALikeli.sh)
ADD_TEST(testPWAConstrain ${EXECUTABLE_OUTPUT_PATH}/testPWAConstraint.sh)
ADD_TEST(testPWAFit    ${EXECUTABLE_OUTPUT_PATH}/testPWAFit.sh)
