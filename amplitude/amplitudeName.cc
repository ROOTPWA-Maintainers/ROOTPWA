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
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      base class that encapsulates naming scheme for production and
//      decay amplitudes
//
//      in the general case the intensity is parameterized by
//
//      I = sum_j sum_k sum_l |sum_i T_i^{j l} A_i^{j k}|^2
//
//      where
//      T_i^{j l} - transition amplitude
//      A_i^{j k} - decay amplitude
//      i         - coherent sum index common for T and A (e.g. I^G J^PC M [isobar decay chain])
//      j         - incoherent sum index common for T and A (e.g. reflectivity)
//      k         - incoherent sum index for A only (e.g. FS-particle helicity)
//      l         - incoherent sum index for T only (e.g. rank)
//
//      so both T and A are defined by 3 sets of quantum numbers
//      the total amplitude name thus consists of 3 substrings
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "amplitudeName.h"

  
using namespace std;
using namespace boost;
using namespace rpwa;


ClassImp(amplitudeName);


bool amplitudeName::_debug = false;


amplitudeName::amplitudeName()
	: rpwa::waveName  (),
		_incohAmpQnLabel("")
{
	//amplitudeName::Class()->IgnoreTObjectStreamer();  // don't store TObject's fBits and fUniqueID
}


amplitudeName::amplitudeName(const rpwa::waveName& name,
                             const string&         incohAmpQnLabel)
	: rpwa::waveName  (name),
	  _incohAmpQnLabel(incohAmpQnLabel)
{ }


amplitudeName::amplitudeName(const isobarAmplitudePtr& amp,
                             const string&             incohAmpQnLabel)
	: rpwa::waveName  (amp),
	  _incohAmpQnLabel(incohAmpQnLabel)
{ }


amplitudeName::~amplitudeName()
{ }


void
amplitudeName::clear()
{
	waveName::clear();
	_incohAmpQnLabel = "";
}


amplitudeName&
amplitudeName::operator =(const amplitudeName& ampName)
{
	if (this != &ampName) {
		waveName::operator =(ampName);
		_incohAmpQnLabel   = ampName._incohAmpQnLabel;
	}
	return *this;
}
