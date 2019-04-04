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
#pragma link C++ class rpwa::amplitudeMetadata+;


#pragma read                                                                                        \
        sourceClass="rpwa::eventMetadata"                                                           \
        source="std::map<std::string, std::pair<double, double>> _binningMap"                       \
        version="[-2]"                                                                              \
        targetClass="rpwa::eventMetadata"                                                           \
        target="_multibinBoundaries"                                                                \
        code="{                                                                                     \
{                                                                                                   \
    _multibinBoundaries = onfile._binningMap;                                                       \
}                                                                                                   \
              }";


#pragma read                                                                                        \
        sourceClass="rpwa::eventMetadata"                                                           \
        source="std::string _userString; std::vector<std::string>  _additionalSavedVariableLabels"  \
        version="[-3]"                                                                              \
        targetClass="rpwa::eventMetadata"                                                           \
        target="_auxString, _additionalTreeVariableNames"                                           \
        code="{                                                                                     \
{                                                                                                   \
    _auxString = onfile._userString;                                                                \
    _additionalTreeVariableNames = onfile._additionalSavedVariableLabels;                           \
}                                                                                                   \
              }";


#pragma read                                                                                        \
        sourceClass="rpwa::eventMetadata"                                                           \
        source=""                                                                                   \
        version="[-4]"                                                                              \
        targetClass="rpwa::eventMetadata"                                                           \
        target="_datasetLabel, _datasetDescription"                                                 \
        code="{                                                                                     \
{                                                                                                   \
    _datasetLabel = \"not-defined\";                                                                \
    _datasetDescription = \"Data-set not defined. Default data-set label for backwards compatibility\";\
}                                                                                                   \
              }";


#pragma link C++ class std::vector<std::complex<double> >+;
#pragma link C++ class std::vector<std::string>+;
#pragma link C++ class rpwa::amplitudeTreeLeaf+;
#pragma read sourceClass="rpwa::amplitudeTreeLeaf" version="[1-]" \
	targetClass="rpwa::amplitudeTreeLeaf" \
	source="" target="" \
	code="{ newObj->rebuildSubAmpLabelMap(); }"


#endif
