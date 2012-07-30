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
	dataIterator i = _dataTable.find(partName);
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
				   const set<string>&        decayproducts,
				   const bool& forceDecayCheck)
{
	const pair<particleProperties, string> selector(prototype, sel);
	vector<const particleProperties*>      matchingEntries;
	for (dataIterator i = _dataTable.begin(); i != _dataTable.end(); ++i) {
		if (i->second != selector)
			continue;
		// limit isobar mass, if minMass > 0
		if ((minMass > 0) and (i->second.mass() + minMassWidthFactor * i->second.width() < minMass))
		  { printDebug << i->second.name() << " not in mass window " << flush;
			continue;
		  }
		// apply white list
		bool whiteListMatch = (whiteList.size() == 0) ? true : false;
		for (size_t j = 0; j < whiteList.size(); ++j)
			if ((i->second.name() == whiteList[j]) or (i->second.bareName() == whiteList[j])) {
				whiteListMatch = true;
				break;
			}
		if (not whiteListMatch)
		  { printDebug << i->second.name() << " not in whitelist " << endl;
			continue;
		  }
		// apply black list
		bool blackListMatch = false;
		for (size_t j = 0; j < blackList.size(); ++j)
			if ((i->second.name() == blackList[j]) or (i->second.bareName() == blackList[j])) {
				blackListMatch = true;
				break;
			}
		if (blackListMatch)
		  { printDebug << i->second.name() << " on blacklist " << endl;
			continue;
		  }
		// apply list of decays
		bool decaymatch = true;
		// check if there are decays defined at all
		if(decayproducts.size()>0 && i->second.nDecays()>0){
		  if(!i->second.hasDecay(decayproducts))decaymatch=false;
		}
		else if(forceDecayCheck)decaymatch = false;

		if (!decaymatch)
		  { printDebug << i->second.name() << " does not have a decay into ";
		     std::copy(decayproducts.begin(), decayproducts.end(), std::ostream_iterator<string>(std::cout, " "));
		     cout << endl;
		     continue;
		  }


		if (_debug) {
			printDebug << "found entry " << i->second.name() << " matching " << prototype
			           << " and '" << sel << "'" << flush;
			if (minMass > 0)
				cout << " with mass > " << minMass - minMassWidthFactor * i->second.width() << " GeV";
			if (whiteList.size() > 0)
				cout << " ; in white list";
			if (blackList.size() > 0)
				cout << " ; not in black list";
			if(decaymatch)
			  cout << " ; with allowed decay " ;
			 std::copy(decayproducts.begin(), decayproducts.end(), std::ostream_iterator<string>(std::cout, " "));
		         cout << endl;
		}      
		matchingEntries.push_back(&(i->second));
	}
	return matchingEntries;
}


bool
particleDataTable::addEntry(const particleProperties& partProp)
{
	const string name = partProp.name();
	dataIterator i    = _dataTable.find(name);
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
	for (dataIterator i = begin(); i != end(); ++i) {
		++countEntries;
		out << "entry " << setw(3) << countEntries << ": " << i->second << endl;
	}
	return out;
}


ostream&
particleDataTable::dump(ostream& out)
{
	for (dataIterator i = begin(); i != end(); ++i) {
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
particleDataTable::readDecayFile(const string& fileName)
{
  libconfig::Config config;
  bool result=false;
  parseLibConfigFile(fileName,config,result);
    
  const libconfig::Setting& particles=config.lookup("particles");
  unsigned int npart=particles.getLength();
  if(_debug)printInfo << npart << " Particles found in decay file." << endl;

  // loop through decay entries and add decays to particleDataTable
  for(unsigned int ipart=0; ipart<npart;++ipart){
    const libconfig::Setting& part=particles[ipart];
    std::string name;part.lookupValue("name",name);
    // lookup particle in database
    // convert to modifiable entry
    particleProperties* particleProp=const_cast<particleProperties*>(entry(name));
    if(particleProp==NULL)continue;
    if(_debug){
      printInfo << name << endl;
      printInfo << "decays into: "<< endl;
    }

    // get and loop over decay modesfor this particle
    const libconfig::Setting& decays=part["decays"];
    unsigned int ndec=decays.getLength();
    for(unsigned int idec=0;idec<ndec;++idec){
      // get list of decay products
      const libconfig::Setting& prod=decays[idec]["products"];
      unsigned int n=prod.getLength();
      set<string> daughters;
      // loop over decay products
      for(unsigned int i=0;i<n;++i){
	string prodname=prod[i];
        if(_debug)cout << prodname << " ";
	daughters.insert(prodname);
      }// end loop over decay products
      if(_debug)cout << endl;
      particleProp->addDecayMode(daughters);
    }// end loop over decay modes of this particle
   }// end loop over entries in decaytable

  return result;
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
