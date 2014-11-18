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
//      TTree leaf persistency storage class for amplitude information
//      needed by fit program
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <algorithm>

#include "TClass.h"

#include "reportingUtils.hpp"
#include "amplitudeTreeLeaf.h"


using namespace std;
using namespace rpwa;


ClassImp(amplitudeTreeLeaf);


bool amplitudeTreeLeaf::_debug = false;


amplitudeTreeLeaf::amplitudeTreeLeaf()
	: TObject           (),
	  _incohSubAmps     (1, 0),
	  _incohSubAmpLabels(),
	  _labelToIndexMap  ()
{
	amplitudeTreeLeaf::Class()->IgnoreTObjectStreamer();  // don't store TObject's fBits and fUniqueID
}


amplitudeTreeLeaf::~amplitudeTreeLeaf()
{ }


void amplitudeTreeLeaf::clear()
{
	_incohSubAmps.clear     ();
	_incohSubAmps.resize    (1, 0);
	_incohSubAmpLabels.clear();
	_labelToIndexMap.clear  ();
}


amplitudeTreeLeaf&
amplitudeTreeLeaf::operator =(const amplitudeTreeLeaf& amp)
{
	if (this != &amp) {
		TObject::operator =(amp);
		_incohSubAmps      = amp._incohSubAmps;
		_incohSubAmpLabels = amp._incohSubAmpLabels;
		_labelToIndexMap   = amp._labelToIndexMap;
	}
	return *this;
}


amplitudeTreeLeaf&
amplitudeTreeLeaf::operator +=(const amplitudeTreeLeaf& amp)
{
	if (_incohSubAmpLabels != amp._incohSubAmpLabels) {
		printErr << "cannot add " << *this << endl
		         << "and " << amp << endl
		         << "because the two amplitudes have different incoherent sub-amplitudes. "
		         << "Aborting..." << endl;
		throw;
	}
	for (unsigned int i = 0; i < nmbIncohSubAmps(); ++i)
		_incohSubAmps[i] += amp.incohSubAmp(i);
	return *this;
}


amplitudeTreeLeaf&
amplitudeTreeLeaf::operator -=(const amplitudeTreeLeaf& amp)
{
	if (_incohSubAmpLabels != amp._incohSubAmpLabels) {
		printErr << "cannot subtract " << amp << endl
		         << "from " << *this << endl
		         << "because the two amplitudes have different incoherent sub-amplitudes. "
		         << "Aborting..." << endl;
		throw;
	}
	for (unsigned int i = 0; i < nmbIncohSubAmps(); ++i)
		_incohSubAmps[i] -= amp.incohSubAmp(i);
	return *this;
}


bool
amplitudeTreeLeaf::containsIncohSubAmp(const string& subAmpLabel) const
{
	labelToIndexMapIterator entry = _labelToIndexMap.find(subAmpLabel);
	if (entry == _labelToIndexMap.end())
		return false;
	return true;
}


unsigned int
amplitudeTreeLeaf::incohSubAmpIndex(const string& subAmpLabel) const
{
	labelToIndexMapIterator entry = _labelToIndexMap.find(subAmpLabel);
	if (entry == _labelToIndexMap.end()) {
		printErr << "cannot find subamp '" << subAmpLabel << "'. Aborting..." << endl;
		throw;
	}
	return entry->second;
}


const string&
amplitudeTreeLeaf::incohSubAmpName(const unsigned int subAmpIndex) const
{
	if (subAmpIndex < _incohSubAmpLabels.size())
		return _incohSubAmpLabels[subAmpIndex];
	else {
		printErr << "subamp index " << subAmpIndex << " is out of range. Aborting..." << endl;
		throw;
	}
}


void
amplitudeTreeLeaf::defineIncohSubAmps(const vector<string>& subAmpLabels)
{
	const unsigned int nmbSubAmps = subAmpLabels.size();
	if (nmbSubAmps < 2)
		printWarn << "vector with subamp labels '" << subAmpLabels << "' has less than 2 entries. "
		          << "using single amplitude value." << endl;
	// make sure labels are unique
	vector<string> uniqueLabels = subAmpLabels;
	sort(uniqueLabels.begin(), uniqueLabels.end());
	uniqueLabels.erase(unique(uniqueLabels.begin(), uniqueLabels.end()), uniqueLabels.end());
	if (uniqueLabels.size() < nmbSubAmps) {
		printErr << "vector with subamp labels '" << subAmpLabels << "' contains dublicate entries. "
		         << "Aborting..." << endl;
		throw;
	}
	_incohSubAmps.resize(nmbSubAmps, 0);
	_incohSubAmpLabels = subAmpLabels;
	rebuildSubAmpLabelMap();
}


ostream&
amplitudeTreeLeaf::print(ostream& out) const
{
	out << "amplitude tree leaf:"  << endl;
	const unsigned int nmbSubAmps = nmbIncohSubAmps();
	if (nmbSubAmps > 1) {
		out << "    number of incoherent sub-amps ... " << nmbSubAmps << endl;
		out << "    amplitude value(s):" << endl;
		for (unsigned int i = 0; i < nmbSubAmps; ++i)
			out << "        [" << setw(4) << i << "] = '" << _incohSubAmpLabels[i]
			    << "' = " << maxPrecisionDouble(_incohSubAmps[i]) << endl;
	} else if (nmbSubAmps == 1)
		out << "    amplitude ... " << _incohSubAmps[0] << endl;
	else
		out << "    no amplitude value" << endl;
	return out;
}


void
amplitudeTreeLeaf::rebuildSubAmpLabelMap()
{
	_labelToIndexMap.clear();
	for (unsigned int i = 0; i < _incohSubAmpLabels.size(); ++i)
		_labelToIndexMap[_incohSubAmpLabels[i]] = i;
}
