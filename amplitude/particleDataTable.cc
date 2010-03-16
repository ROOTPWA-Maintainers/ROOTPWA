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

#include "utilities.h"
#include "particleDataTable.h"

	
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
particleDataTable::entry(const string& partName)
{
  dataIterator i = _dataTable.find(partName);
  if (i == _dataTable.end()) {
    printWarn << "could not find entry for particle '" << partName << "'" << endl;
    return 0;
  } else
    return &(i->second);
}



    
std::vector<const particleProperties*> 
particleDataTable::entrylist(const particleProperties& prototype, const std::string& opt){
  pair<particleProperties,string> selector(prototype,opt);
  std::vector<const particleProperties*> matchingentries;
  dataIterator i = _dataTable.begin();
  while(i != _dataTable.end()) {
    if( i->second == selector)matchingentries.push_back(&(i->second));
    ++i;
  } 
  return matchingentries;
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
	      << "    conflicts with: " << partProp  << endl;
    return false;
  } else {
    _dataTable[name] = partProp;
    if (_debug)
      printInfo << "added entry for '" << name << "' into particle data table" << endl;
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
  if (!file || !file.good()) {
    printWarn << "cannot open file '" << fileName << "'" << endl;
    return false;
  }
  return read(file);
}


bool
particleDataTable::read(istream& in)
{
  if (!in || !in.good()) {
    printWarn << "cannot read from input stream" << endl;
    return false;
  }
  if (_debug)
    printInfo << "data table has " << nmbEntries() << " entries (before reading)" << endl;
  unsigned int countEntries = 0;
  while (in.good()) {
    particleProperties partProp;
    if (in >> partProp)
      if(addEntry(partProp))
	++countEntries;
  }
  printInfo << "successfully read " << countEntries << " new entries into particle data table" << endl;
  if (_debug)
    cout << "    data table has " << nmbEntries() << " entries (after reading)" << endl;
  return true;
}
