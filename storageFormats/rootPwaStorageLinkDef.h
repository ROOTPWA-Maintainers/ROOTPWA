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


#pragma link C++ class rpwa::eventMetadata+;

#ifdef USE_STD_COMPLEX_TREE_LEAFS
#pragma link C++ class std::vector<std::complex<double> >+;
#pragma link C++ class std::vector<std::string>+;
#pragma link C++ class rpwa::amplitudeTreeLeaf+;
#pragma read sourceClass="rpwa::amplitudeTreeLeaf" version="[1-]"	  \
	targetClass="rpwa::amplitudeTreeLeaf" \
	source="" target="" \
	code="{ newObj->rebuildSubAmpLabelMap(); }"
#endif


#endif
