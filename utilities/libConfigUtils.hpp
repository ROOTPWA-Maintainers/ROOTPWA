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
//      collection of useful libConfig routines
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef LIBCONFIGUTILS_H
#define LIBCONFIGUTILS_H


#include "libconfig.h++"

#include "reportingUtils.hpp"


namespace rpwa {

	
	inline
	bool
	parseLibConfigFile(const std::string& libConfigFileName,
	                   libconfig::Config& config,
	                   const bool         debug)  ///< lets config object parse a libConfig file
	{
		if (debug)
			printDebug << "parsing libConfig file '" << libConfigFileName << "'" << std::endl;
		try {
			config.readFile(libConfigFileName.c_str());
		} catch (const libconfig::FileIOException& ioEx) {
			printWarn << "I/O error while reading libConfig file "
			          << "'" << libConfigFileName << "'" << std::endl;
			return false;
		} catch (const libconfig::ParseException&  parseEx) {
			printWarn << "parse error in '" << parseEx.getFile() << "' line " << parseEx.getLine()
			          << ": " << parseEx.getError() << std::endl;
			return false;
		}
		return true;
	}


	inline
	bool
	parseLibConfigString(const std::string& libConfigString,
	                     libconfig::Config& config,
	                     const bool         debug)  ///< lets config object parse a libConfig string
	{
		if (debug)
			printDebug << "parsing libConfig string" << std::endl;
		try {
			config.readString(libConfigString);
		} catch (const libconfig::ParseException&  parseEx) {
			printWarn << "parse error in line " << parseEx.getLine() << " of libConfig string: "
			          << parseEx.getError() << std::endl;
			return false;
		}
		return true;
	}


	inline
	const libconfig::Setting*
	findLibConfigGroup(const libconfig::Setting& parent,
	                   const std::string&        groupName,
	                   const bool                mustExist = true)  ///< finds field in setting and makes sure it is a group
	{
		// find field
		if (not parent.exists(groupName)) {
			if (mustExist)
				printWarn << "cannot find field '" << groupName << "' in setting "
				          << "'" << parent.getPath() << "' in file "
				          << "'" << parent.getSourceFile() << "'" << std::endl;
			return 0;
		}
		const libconfig::Setting& group = parent[groupName];
		// check that it is a group
		if (not group.isGroup()) {
			printWarn << "field '" << groupName << "' in setting '" << parent.getPath() << "' "
			          << "in file '" << parent.getSourceFile() << "' is not a group" << std::endl;
			return 0;
		}
		return &group;
	}


	inline
	const libconfig::Setting*
	findLibConfigGroupAtIndex(const libconfig::Setting& parent,
	                          const int                 index,
	                          const bool                mustExist = true)  ///< finds field in group, array, or list and makes sure it is a group
	{
		// find field
		const libconfig::Setting* group = 0;
		try {
			group = &(parent[index]);
		} catch (const libconfig::SettingNotFoundException&) {
			if (mustExist)
				printWarn << "no entry at index [" << index << "] in setting "
				          << "'" << parent.getPath() << "' in file "
				          << "'" << parent.getSourceFile() << "'" << std::endl;
			return 0;
		}
		// check that it is a group
		if (not (*group).isGroup()) {
			printWarn << "entry [" << index << "] in setting '" << parent.getPath() << "' "
			          << "in file '" << parent.getSourceFile() << "' is not a group" << std::endl;
			return 0;
		}
		return group;
	}


	inline
	const libconfig::Setting*
	findLibConfigList(const libconfig::Setting& parent,
	                  const std::string&        listName,
	                  const bool                mustExist = true)  ///< finds field in setting and makes sure it is a non-empty list
	{
		// find field
		if (not parent.exists(listName)) {
			if (mustExist)
				printWarn << "cannot find field '" << listName << "' in setting '" << parent.getPath() << "' "
				          << "in file '" << parent.getSourceFile() << "'" << std::endl;
			return 0;
		}
		const libconfig::Setting& list = parent[listName];
		// check that it is a list
		if (not list.isList()) {
			printWarn << "field '" << listName << "' in setting '" << parent.getPath() << "' "
			          << "in file '" << parent.getSourceFile() << "' is not a list. "
			          << "check that braces are correct." << std::endl;
			return 0;
		}
		// check that it is not empty
		if (list.getLength() < 1) {
			printWarn << "list '" << listName << "' in setting '" << parent.getPath() << "' "
			          << "in file '" << parent.getSourceFile() << "' is empty" << std::endl;
			return 0;
		}
		return &list;
	}


	inline
	const libconfig::Setting*
	findLibConfigListAtIndex(const libconfig::Setting& parent,
	                         const int                 index,
	                         const bool                mustExist = true)  ///< finds field in group, array, or list and makes sure it is a non-empty list
	{
		// find field
		const libconfig::Setting* list = 0;
		try {
			list = &(parent[index]);
		} catch (const libconfig::SettingNotFoundException& notFound) {
			printWarn << "no entry at index [" << index << "] in setting '" << parent.getPath() << "' "
			          << "in file '" << parent.getSourceFile() << "'" << std::endl;
			return 0;
		}
		// check that it is a list
		if (not (*list).isList()) {
			printWarn << "entry [" << index << "] in setting '" << parent.getPath() << "' "
			          << "in file '" << parent.getSourceFile() << "' is not a list. "
			          << "check that braces are correct." << std::endl;
			return 0;
		}
		
		// check that it is not empty
		if ((*list).getLength() < 1) {
			printWarn << "list at [" << index << "] in setting '" << parent.getPath() << "' "
			          << "in file '" << parent.getSourceFile() << "' is empty" << std::endl;
			return 0;
		}
		return list;
	}

	
	inline
	const libconfig::Setting*
	findLibConfigArray(const libconfig::Setting& parent,
	                   const std::string&        arrayName,
	                   const bool                mustExist = true)  ///< finds field in setting and makes sure it is a non-empty array
	{
		// find field
		if (not parent.exists(arrayName)) {
			if (mustExist)
				printWarn << "cannot find field '" << arrayName << "' in setting "
				          << "'" << parent.getPath() << "' in file "
				          << "'" << parent.getSourceFile() << "'" << std::endl;
			return 0;
		}
		const libconfig::Setting& array = parent[arrayName];
		// check that it is a list
		if (not array.isArray()) {
			printWarn << "field '" << arrayName << "' in setting '" << parent.getPath() << "' "
			          << "in file '" << parent.getSourceFile() << "' is not an array. "
			          << "check that braces are correct." << std::endl;
			return 0;
		}
		// check that it is not empty
		if (array.getLength() < 1) {
			printWarn << "array '" << arrayName << "' in setting '" << parent.getPath() << "' "
			          << "in file '" << parent.getSourceFile() << "' is empty" << std::endl;
			return 0;
		}
		return &array;
	}

	
	inline
	const libconfig::Setting*
	findLibConfigArrayAtIndex(const libconfig::Setting& parent,
	                          const int                 index,
	                          const bool                mustExist = true)  ///< finds field in group, array, or list and makes sure it is a non-empty array
	{
		// find field
		const libconfig::Setting* array = 0;
		try {
			array = &(parent[index]);
		} catch (const libconfig::SettingNotFoundException& NotFound) {
			printWarn << "no entry at index [" << index << "] in setting '" << parent.getPath() << "' "
			          << "in file '" << parent.getSourceFile() << "'" << std::endl;
			return 0;
		}
		// check that it is a list
		if (not (*array).isArray()) {
			printWarn << "entry [" << index << "] in setting '" << parent.getPath() << "' "
			          << "in file '" << parent.getSourceFile() << "' is not an array. "
			          << "check that braces are correct." << std::endl;
			return 0;
		}
		// check that it is not empty
		if ((*array).getLength() < 1) {
			printWarn << "array at [" << index << "] in setting '" << parent.getPath() << "' "
			          << "in file '" << parent.getSourceFile() << "' is empty" << std::endl;
			return 0;
		}
		return array;
	}
	
	
} // namespace rpwa


#endif  // LIBCONFIGUTILS_H
