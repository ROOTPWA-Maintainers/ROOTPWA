///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009
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


#ifndef PARTICLEDATATABLE_H
#define PARTICLEDATATABLE_H


#include <string>

#include "particleProperties.h"


namespace rpwa {


  class particleDataTable {

  public:

    static particleDataTable& instance() { return _instance; }

    static const particleProperties* operator [] (const string& partName);   ///< access properties by particle name
    typedef iterator = std::map<std::string, particleData>::const_iterator;
    static iterator begin() { return _dataTable.begin(); } 
    static iterator end()   { return _dataTable.end();   } 
    static unsigend int nmbEntries() const { return _dataTable.size(); }

    static void print(std::ostream& out) const;  ///< prints particle data in human-readable form
    static void dump (std::ostream& out) const;  ///< dumps particle properties in format of data file
    friend std::ostream& operator << (std::ostream&             out,
				      const particleProperties& partProp);

    static bool readFile(const string& fileName = "./particleDataTable.txt");  ///< reads in particle data from file
    static bool read(std::istream& in);  ///< reads whitespace separated properties from stream
    friend std::istream& operator << (std::istream&       in,
				      particleProperties& partProp);

    static void reset() { _dataTable.clear(); }


  private:

    particleDataTable () { }
    ~particleDataTable() { }
    particleDataTable (const particleDataTable&);
    particleDataTable& operator = (const particleDataTable&);
  

    static particleDataTable                   _instance;   ///< singleton instance
    static std::map<std::string, particleData> _dataTable;  ///< map with particle data
	
  };


} // namespace rpwa
	

#endif  // PARTICLEDATATABLE_H
