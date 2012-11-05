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


#ifndef PARTICLEDATATABLE_H
#define PARTICLEDATATABLE_H


#include <string>
#include <vector>
#include <map>
#include <set>


#include "particleProperties.h"


namespace rpwa {


	class particleDataTable {

	public:

		static particleDataTable& instance() { return _instance; }  ///< get singleton instance

		static bool isInTable(const std::string& partName);  ///< returns, whether particle has a table entry

		static const particleProperties* entry(const std::string& partName,
		                                       const bool         warnIfNotExistent = true);  ///< access properties by particle name
	
		static bool addEntry(const particleProperties& partProp);  ///< adds entry to particle data table

		static std::vector<const particleProperties*>
		entriesMatching(const particleProperties&         prototype,
		                const std::string&                sel,
		                const double                      minMass            = 0,
		                const double                      minMassWidthFactor = 0,
		                const std::vector<std::string>&   whiteList          = std::vector<std::string>(),
		                const std::vector<std::string>&   blackList          = std::vector<std::string>(),
		                const std::multiset<std::string>& decayProducts      = std::multiset<std::string>(),
		                const bool&                       forceDecayCheck    = true);  ///< returns entries that have the same quantum numbers as prototype property; quantum numbers to be compared are selected by sel string; if minMass > 0 the isobar mass is limited; checks for allowed decays if they are defined; decay checks can be forced, then particles which have no specified decays will be discarded

		static unsigned int nmbEntries() { return _dataTable.size(); }  ///< returns number of entries in particle data table

		typedef std::map<std::string, particleProperties>::const_iterator iterator;
		static iterator begin() { return _dataTable.begin(); }  ///< returns iterator pointing at first entry of particle data table
		static iterator end()   { return _dataTable.end();   }  ///< returns iterator pointing after last entry of particle data table

		static std::ostream& print(std::ostream& out);  ///< prints particle data in human-readable form
		static std::ostream& dump (std::ostream& out);  ///< dumps particle properties in format of data file

		static bool readFile(const std::string& fileName = "./particleDataTable.txt");  ///< reads in particle data from file
		static bool read(std::istream& in);  ///< reads whitespace separated properties from stream

		static bool readDecayModeFile(const std::string& fileName);  ///< reads in decay modes for list of particles from file; requires particle properties


		static void clear() { _dataTable.clear(); }  ///< deletes all entries in particle data table

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag
     

	private:

		particleDataTable () { }
		~particleDataTable() { }
		particleDataTable (const particleDataTable&);
		particleDataTable& operator =(const particleDataTable&);

		static particleDataTable                         _instance;   ///< singleton instance
		static std::map<std::string, particleProperties> _dataTable;  ///< map with particle data

		static bool _debug;  ///< if set to true, debug messages are printed

	};

  
	inline
	std::ostream&
	operator <<(std::ostream&            out,
	            const particleDataTable& dataTable)
	{
		return dataTable.print(out);
	}


	inline
	std::istream&
	operator >>(std::istream&      in,
	            particleDataTable& dataTable)
	{
		dataTable.read(in);
		return in;
	}


} // namespace rpwa
	

#endif  // PARTICLEDATATABLE_H
