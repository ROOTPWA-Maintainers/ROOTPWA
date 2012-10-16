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

#include "libConfigUtils.hpp"
#include "reportingUtils.hpp"
#include "particleDataTable.h"
#include "libconfig.h++"

	
using namespace std;
using namespace rpwa;
using namespace libconfig;


particleDataTable               particleDataTable::_instance;
map<string, particleProperties> particleDataTable::_dataTable;
bool                            particleDataTable::_debug = false;


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
particleDataTable::entriesMatching(const particleProperties& prototype,
                                   const string&             sel,
                                   const double              minMass,
                                   const double              minMassWidthFactor,
                                   const vector<string>&     whiteList,
                                   const vector<string>&     blackList,
                                   const multiset<string>&   decayProducts,
                                   const bool&               forceDecayCheck)
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
		if ((decayProducts.size() > 0) and partProp.nmbDecays() > 0) {
			if (not partProp.hasDecay(decayProducts))
				decayMatch = false;
		} else if (forceDecayCheck)
			decayMatch = false;
		if (not decayMatch) {
			if (_debug) {
				printDebug << partProp.name() << " does not have a decay into ";
				copy(decayProducts.begin(), decayProducts.end(), ostream_iterator<string>(cout, "  "));
				cout << endl;
			}
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
				cout << " ; with allowed decay into " ;
			copy(decayProducts.begin(), decayProducts.end(), ostream_iterator<string>(cout, "  "));
			cout << endl;
		}      
		matchingEntries.push_back(&partProp);
	}
	return matchingEntries;
}


bool
particleDataTable::addEntry(const particleProperties& partProp)
{
	const string name = partProp.name();
	iterator i    = _dataTable.find(name);
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
  for (unsigned int i = 0; i < nmbEntries; ++i) {
	  const Setting& partDecayEntry = (*partDecayList)[i];
    string         partName;
    if (not partDecayEntry.lookupValue("name", partName)) {
	    printWarn << "particle decay list entry [" << i << "] does not hava a 'name' field. "
	              << "ignoring entry." << endl;
	    continue;
    }
    if (_debug)
	    printDebug << "setting decays for '" << partName << "':" << endl;
    // get list of decay modes for this particle
    const Setting* partDecays = findLibConfigList(partDecayEntry, "decays");
    if (not partDecays) {
	    printWarn << "cannot find 'decays' list for " << partName
	              << " particle decay (list entry [" << i << "]). "
	              << "ignoring entry." << endl;
	    continue;
    }
    // lookup particle in database
    map<string, particleProperties>::iterator i = _dataTable.find(partName);
    if (i == _dataTable.end()) {
	    printWarn << "could not find particle " << partName << " in data table. "
	              << "ignoring entry." << endl;
	    continue;
    }
    particleProperties& particleProp = i->second;
    // loop over decay modes for this particle
    const unsigned int nmbDecays = partDecays->getLength();
    for (unsigned int j = 0; j < nmbDecays; ++j) {
      // get string array of decay daughters
	    const Setting* decayDaughters = findLibConfigArray((*partDecays)[j], "products");
	    if (not decayDaughters) {
		    printWarn << "cannot find 'products' array for " << partName
		              << " particle decay list entry [" << j << "]. "
		              << "ignoring entry." << endl;
		    continue;
	    }
	    if (_debug)
		    cout << "    -> ";
	    // loop over decay products
      const unsigned int nmbDaughters = decayDaughters->getLength();
      multiset<string>   daughters;
      for (unsigned int k = 0; k < nmbDaughters; ++k) {
	      const string daughterName = (*decayDaughters)[k];
	      daughters.insert(daughterName);
	      if (_debug)
		      cout << daughterName << ((k < nmbDaughters - 1) ? "  " : "");
      }
      if(_debug)
	      cout << endl;
      particleProp.addDecayMode(daughters);
      ++countEntries;
    }
  }  // end loop over decay list entries

  printSucc << "read " << countEntries << " particle decays into particle data table" << endl;
  return true;
}
