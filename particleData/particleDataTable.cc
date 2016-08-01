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
//      singleton class that manages all particle data
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <fstream>
#include <iomanip>
#include <iterator>

#include <boost/assign/list_inserter.hpp>

#include "libConfigUtils.hpp"
#include "reportingUtils.hpp"
#include "particleDataTable.h"
#include "libconfig.h++"


using namespace std;
using namespace rpwa;
using namespace libconfig;
using namespace boost::bimaps;


namespace {

	typedef bimap<string, unsigned int> nameGeantIdBimap;
	nameGeantIdBimap initNameGeantIdTranslator()
	{
		nameGeantIdBimap translator;
		boost::assign::insert(translator)
			(nameGeantIdBimap::value_type("gamma",       1))
			(nameGeantIdBimap::value_type("e+",          2))
			(nameGeantIdBimap::value_type("e-",          3))
			(nameGeantIdBimap::value_type("mu+",         5))
			(nameGeantIdBimap::value_type("mu-",         6))
			(nameGeantIdBimap::value_type("pi0",         7))
			(nameGeantIdBimap::value_type("pi+",         8))
			(nameGeantIdBimap::value_type("pi-",         9))
			(nameGeantIdBimap::value_type("K_L",        10))
			(nameGeantIdBimap::value_type("K+",         11))
			(nameGeantIdBimap::value_type("K-",         12))
			(nameGeantIdBimap::value_type("n",          13))
			(nameGeantIdBimap::value_type("p+",         14))
			(nameGeantIdBimap::value_type("pbar-",      15))
			(nameGeantIdBimap::value_type("K_S",        16))
			(nameGeantIdBimap::value_type("eta",        17))
			(nameGeantIdBimap::value_type("Lambda",     18))
			(nameGeantIdBimap::value_type("nbar",       25))
			(nameGeantIdBimap::value_type("Lambdabar",  26))
			(nameGeantIdBimap::value_type("rho(770)0",  57))
			(nameGeantIdBimap::value_type("rho(770)+",  58))
			(nameGeantIdBimap::value_type("rho(770)-",  59))
			(nameGeantIdBimap::value_type("omega(782)", 60))
			(nameGeantIdBimap::value_type("eta'(982)",  61))
			(nameGeantIdBimap::value_type("phi(1020)",  62));
		return translator;
	}

}


particleDataTable               particleDataTable::_instance;
map<string, particleProperties> particleDataTable::_dataTable;
bool                            particleDataTable::_debug = false;
nameGeantIdBimap                particleDataTable::_nameGeantIdMap = initNameGeantIdTranslator();


string particleDataTable::particleNameFromGeantId(const int id) {
	nameGeantIdBimap::right_const_iterator it = _nameGeantIdMap.right.find(id);
	if (it == _nameGeantIdMap.right.end()) {
		printErr << id << " is unknown GEANT particle ID. returning particle name 'unknown'." << endl;
		return "unknown";
	}
	// make sure charge is correctly put into particle name
	int charge;
	const string bareName = particleProperties::chargeFromName(it->second, charge);
	return particleProperties::nameWithCharge(bareName, charge);
}


void particleDataTable::geantIdAndChargeFromParticleName(const string& name,
                                                         int&          id,
                                                         int&          charge)
{
	id = 0;
	const string bareName = particleProperties::chargeFromName(name, charge);
	nameGeantIdBimap::left_const_iterator it = _nameGeantIdMap.left.find(name);
	if(it == _nameGeantIdMap.left.end()) {
		// try again with charge stripped from name
		it = _nameGeantIdMap.left.find(bareName);
	}
	if(it == _nameGeantIdMap.left.end()) {
		printErr << "particle '" << name << "' is unknown. returning GEANT particle ID 0." << endl;
		return;
	}
	id = it->second;
}


unsigned int particleDataTable::geantIdFromParticleName(const std::string& name) {
	int retval;
	int tmp;
	geantIdAndChargeFromParticleName(name, retval, tmp);
	return retval;
}


bool
particleDataTable::isInTable(const string& partName)
{
	if (_dataTable.find(partName) == _dataTable.end()) {
		return false;
	} else
		return true;
}


const particleProperties*
particleDataTable::entry(const string& partName,
                         const bool    warnIfNotExistent)
{
	iterator i = _dataTable.find(partName);
	if (i == _dataTable.end()) {
		if (warnIfNotExistent)
			printWarn << "could not find entry for particle '" << partName << "'" << endl;
		return 0;
	} else
		return &(i->second);
}


vector<const particleProperties*>
particleDataTable::entriesMatching(const particleProperties&            prototype,
                                   const string&                        sel,
                                   const double                         minMass,
                                   const double                         minMassWidthFactor,
                                   const vector<string>&                whiteList,
                                   const vector<string>&                blackList,
                                   const particleProperties::decayMode& decay,
                                   const bool&                          forceDecayCheck)
{
	const pair<particleProperties, string> selector(prototype, sel);
	vector<const particleProperties*>      matchingEntries;
	for (iterator i = _dataTable.begin(); i != _dataTable.end(); ++i) {
		if (i->second != selector)
			continue;
		const particleProperties& partProp = i->second;
		// limit isobar mass, if minMass > 0
		if ((minMass > 0) and (partProp.mass() + minMassWidthFactor * partProp.width() < minMass))
		{
			if(_debug)printDebug << partProp.name() << " not in mass window " << flush;
			continue;
		}
		// apply white list
		bool whiteListMatch = (whiteList.size() == 0) ? true : false;
		for (size_t j = 0; j < whiteList.size(); ++j)
			if ((partProp.name() == whiteList[j]) or (partProp.bareName() == whiteList[j])) {
				whiteListMatch = true;
				break;
			}
		if (not whiteListMatch){
			if(_debug)printDebug << partProp.name() << " not in whitelist " << endl;
			continue;
		}
		// apply black list
		bool blackListMatch = false;
		for (size_t j = 0; j < blackList.size(); ++j)
			if ((partProp.name() == blackList[j]) or (partProp.bareName() == blackList[j])) {
				blackListMatch = true;
				break;
			}
		if (blackListMatch) {
			if(_debug)printDebug << partProp.name() << " on blacklist " << endl;
			continue;
		}
		// apply list of decays
		bool decayMatch = true;
		if ((decay._daughters.size() > 0) and partProp.nmbDecays() > 0) {
			if (not partProp.hasDecay(decay))
				decayMatch = false;
		} else if (forceDecayCheck)
			decayMatch = false;
		if (not decayMatch) {
			if (_debug)
				printDebug << partProp.name() << " does not have a decay into " << decay << endl;
			continue;
		}

		if (_debug) {
			printDebug << "found entry " << partProp.name() << " matching " << prototype
				<< " and '" << sel << "'" << flush;
			if (minMass > 0)
				cout << " with mass > " << minMass - minMassWidthFactor * partProp.width() << " GeV";
			if (whiteList.size() > 0)
				cout << " ; in white list";
			if (blackList.size() > 0)
				cout << " ; not in black list";
			if (decayMatch)
				cout << " ; with allowed decay into " << decay << endl;
		}
		matchingEntries.push_back(&partProp);
	}
	return matchingEntries;
}


bool
particleDataTable::addEntry(const particleProperties& partProp)
{
	const string name = partProp.name();
	iterator     i    = _dataTable.find(name);
	if (i != _dataTable.end()) {
		printWarn << "trying to add entry for particle '" << name << "' "
		          << "which already exists in table"     << endl
		          << "    existing entry: " << i->second << endl
		          << "    conflicts with: " << partProp  << endl
		          << "    entry was not added to table." << endl;
		return false;
	} else {
		_dataTable[name] = partProp;
		if (_debug)
			printDebug << "added entry for '" << name << "' into particle data table" << endl;
		return true;
	}
}


ostream&
particleDataTable::print(ostream& out)
{
	unsigned int countEntries = 0;
	for (iterator i = begin(); i != end(); ++i) {
		++countEntries;
		out << "entry " << setw(3) << countEntries << ": " << i->second << endl;
	}
	return out;
}


ostream&
particleDataTable::dump(ostream& out)
{
	for (iterator i = begin(); i != end(); ++i) {
		i->second.dump(out);
		out << endl;
	}
	return out;
}


bool
particleDataTable::readFile(const string& fileName)
{
	printInfo << "reading particle data from file '" << fileName << "'" << endl;
	ifstream file(fileName.c_str());
	if (not file or not file.good()) {
		printWarn << "cannot open file '" << fileName << "'" << endl;
		return false;
	}
	return read(file);
}


bool
particleDataTable::read(istream& in)
{
	if (not in or not in.good()) {
		printWarn << "cannot read from input stream" << endl;
		return false;
	}
	if (_debug)
		printDebug << "data table has " << nmbEntries() << " entries (before reading)" << endl;
	unsigned int countEntries = 0;
	while (in.good()) {
		particleProperties partProp;
		if (in >> partProp) {
			if (addEntry(partProp))
				++countEntries;
			if (not partProp.isItsOwnAntiPart() and addEntry(partProp.antiPartProperties()))
				++countEntries;
		}
	}
	printSucc << "read " << countEntries << " new entries into particle data table" << endl;
	if (_debug)
		cout << "    data table has " << nmbEntries() << " entries (after reading)" << endl;
	return true;
}


bool
particleDataTable::readDecayModeFile(const string& fileName)
{
	printInfo << "reading particle decay modes from file '" << fileName << "'" << endl;
	// parse decay mode file
  Config file;
	if (not parseLibConfigFile(fileName, file, _debug)) {
		printWarn << "problems reading file with decay modes '" << fileName << "'. "
		          << "cannot set decay modes." << endl;
		return false;
	}
	const Setting& root = file.getRoot();
	// find list of particle decays
	const Setting* partDecayList = findLibConfigList(root, "particleDecayList");
	if (not partDecayList) {
		printWarn << "cannot find 'particleDecayList' list. cannot set decay modes." << endl;
		return false;
	}
  const unsigned int nmbEntries = partDecayList->getLength();
  if (_debug)
	  printDebug << "found particle decay list with " << nmbEntries << " entries in file "
	             << "'" << fileName << "'" << endl;
  // loop over decay list entries and add decays to particleDataTable
  unsigned int countEntries = 0;
  for (unsigned int entryIndex = 0; entryIndex < nmbEntries; ++entryIndex) {
	  const Setting& partDecayEntry = (*partDecayList)[entryIndex];
    string         partName;
    if (not partDecayEntry.lookupValue("name", partName)) {
	    printWarn << "particle decay list entry [" << entryIndex << "] does not hava a 'name' field. "
	              << "ignoring entry." << endl;
	    continue;
    }
    if (_debug)
	    printDebug << "setting decays for '" << partName << "':" << endl;
    // get list of decay modes for this particle
    const Setting* partDecays = findLibConfigList(partDecayEntry, "decays");
    if (not partDecays) {
	    printWarn << "cannot find 'decays' list for " << partName
	              << " particle decay (list entry [" << entryIndex << "]). "
	              << "ignoring entry." << endl;
	    continue;
    }
    // lookup particle in database
    map<string, particleProperties>::iterator dataTableEntry = _dataTable.find(partName);
    if (dataTableEntry == _dataTable.end()) {
	    printWarn << "could not find particle " << partName << " in data table. "
	              << "ignoring entry." << endl;
	    continue;
    }
    particleProperties& particleProp = dataTableEntry->second;
    // loop over decay modes for this particle
    const unsigned int nmbDecays = partDecays->getLength();
    for (unsigned int decayIndex = 0; decayIndex < nmbDecays; ++decayIndex) {
      // get string array of decay daughters
	    const Setting* decayDaughters = findLibConfigArray((*partDecays)[decayIndex], "products");
	    if (not decayDaughters) {
		    printWarn << "cannot find 'products' array for " << partName
		              << " particle decay list entry [" << decayIndex << "]. "
		              << "ignoring entry." << endl;
		    continue;
	    }
	    if (_debug)
		    cout << "    -> ";
	    // loop over decay products
      const unsigned int nmbDaughters = decayDaughters->getLength();
      multiset<string>   daughters;
      for (unsigned int daughterIndex = 0; daughterIndex < nmbDaughters; ++daughterIndex) {
	      const string daughterName = (*decayDaughters)[daughterIndex];
	      daughters.insert(daughterName);
	      if (_debug)
		      cout << daughterName << ((daughterIndex < nmbDaughters - 1) ? "  " : "");
      }
      if(_debug)
	      cout << endl;
      particleProperties::decayMode decay(daughters);
      (*partDecays)[decayIndex].lookupValue("L", decay._L);
      (*partDecays)[decayIndex].lookupValue("S", decay._S);

			const Setting* configMassDependencies = findLibConfigArray((*partDecays)[decayIndex], "massDeps", false);
			if (configMassDependencies) {
				vector<string> massDependencies;
				if (configMassDependencies->getLength() > 0 && (*configMassDependencies)[0].getType() == Setting::TypeString) {
					for (int j = 0; j < configMassDependencies->getLength(); ++j)
						massDependencies.push_back((*configMassDependencies)[j]);
				} else if (configMassDependencies->getLength() > 0) {
					printWarn << "mass dependencies setting for isobar '" << partName << "' is not an array of strings. "
					          << "not using special mass dependencies for this isobar." << endl;
				}
				decay._massDependencies = massDependencies;
			}

      particleProp.addDecayMode(decay);
      ++countEntries;
    }
  }  // end loop over decay list entries

  printSucc << "read " << countEntries << " particle decays into particle data table" << endl;
  return true;
}
