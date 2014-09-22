///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      tells rootcint for which classes to generate method interface stubs
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifdef __CINT__


#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;


#pragma link C++ class rpwa::waveDescription+;
// data model evolution rule that triggers parsing of key file string
// whenever waveDescription is read from file
// see http://root.cern.ch/root/html/io/DataModelEvolution.html
// and http://indico.cern.ch/contributionDisplay.py?contribId=210&sessionId=59&confId=35523
#pragma read sourceClass="rpwa::waveDescription" version="[1-]"	  \
	targetClass="rpwa::waveDescription" \
	source="" target="" \
	code="{ newObj->parseKeyFileContent(newObj->keyFileContent()); }"


#pragma link C++ class std::vector<std::complex<double> >+;
#pragma link C++ class std::vector<std::string>+;
#pragma link C++ class std::vector<rpwa::waveDescription>+;
#pragma link C++ class rpwa::ampIntegralMatrix-;


#endif
