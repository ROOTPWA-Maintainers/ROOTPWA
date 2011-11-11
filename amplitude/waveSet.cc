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
//      handles wave sets for fitting
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "libconfig.h++"

#include "libConfigUtils.hpp"
#include "waveSet.h"

  
using namespace std;
using namespace libconfig;
using namespace rpwa;


#ifdef USE_STD_COMPLEX_TREE_LEAFS
ClassImp(waveSet);
#endif

bool waveSet::_debug = false;


waveSet::waveSet()
	: TObject           (),
	  decayAmpTreeNames (),
	  decayAmpMassRanges()
{
	//waveSet::Class()->IgnoreTObjectStreamer();  // don't store TObject's fBits and fUniqueID
}


waveSet::~waveSet()
{ }


void
waveSet::clear()
{
	decayAmpTreeNames.clear ();
	decayAmpMassRanges.clear();
}


waveSet&
waveSet::operator =(const waveSet& set)
{
	if (this != &set) {
		TObject::operator  =(set);
		decayAmpTreeNames  = set.decayAmpTreeNames;
		decayAmpMassRanges = set.decayAmpMassRanges;
	}
	return *this;
}


bool
waveSet::parseWaveSetFile(const string& waveSetFileName)
{
	libconfig::Config config;
	if (not parseLibConfigFile(waveSetFileName, config, _debug)) {
		printWarn << "problems reading wave set file '" << waveSetFileName << "'. "
		          << "cannot construct wave set." << endl;
		return false;
	}

	const Setting& configRoot = config.getRoot();
	printInfo << "constructing wave set from file '" << waveSetFileName << "'" << endl;

	// find waveSet group
	const Setting* configWaveSet = findLibConfigGroup(configRoot, "waveSet", true);
	if (not configWaveSet) {
		printWarn << "cannot construct wave set." << endl;
		return false;
	}

	// find list of decay amplitude specifications
	const Setting* configDecayAmps = findLibConfigList(*configWaveSet, "decayAmplitudes", true);
	if (not configDecayAmps) {
		printWarn << "cannot construct wave set." << endl;
		return false;
	}

	// read decay amplitude list
	const int nmbEntries = configDecayAmps->getLength();
	decayAmpTreeNames.clear ();
	decayAmpMassRanges.clear();
	decayAmpTreeNames.resize (nmbEntries);
	decayAmpMassRanges.resize(nmbEntries);
	bool success = true;
	for (int i = 0; i < nmbEntries; ++i) {
		// get tree name
		string treeName = "";
		if (not (*configDecayAmps)[i].lookupValue("treeName", treeName)) {
			printWarn << "entry [" << i << "] in decay amplitude list does not specify a tree name.";
			success = false;
		}
		decayAmpTreeNames[i] = treeName;
		// get mass range
		pair<double, double> massRange(0, numeric_limits<double>::infinity());
		const Setting*       configMassRange = findLibConfigArray((*configDecayAmps)[i], "massRange", false);
		if (configMassRange) {
			if (configMassRange->getLength() == 1)
				massRange.first  = (*configMassRange)[0];
			else if (configMassRange->getLength() == 2) {
				massRange.first  = (*configMassRange)[0];
				massRange.second = (*configMassRange)[1];
			} else {
				printWarn << "cannot read mass range from entry [" << i << "] in decay amplitude list. "
				          << "array length is neither 1 nor 2." << endl;
				success = false;
			}
		}
		decayAmpMassRanges[i] = massRange;
	}	

	if (success)
		printSucc << "constructed wave set from file '" << waveSetFileName << "'" << endl;
	else
		printSucc << "problems constructing wave set from file '" << waveSetFileName << "'" << endl;
	return success;
}


ostream&
waveSet::print(ostream& out) const
{
	out << "wave set" << endl;
	for (unsigned int i = 0; i < decayAmpTreeNames.size(); ++i)
		out << "    [" << i << "]: tree name = '" << decayAmpTreeNames[i] << "', mass range = "
		    << "[" << decayAmpMassRanges[i].first << ", " << decayAmpMassRanges[i].second << "] MeV/c^2"
		    << endl;
	return out;
}
