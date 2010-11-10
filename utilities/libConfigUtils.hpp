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
// $Rev:: 462                         $: revision of last commit
// $Author:: bgrube                   $: author of last commit
// $Date:: 2010-11-03 19:15:35 +0100 #$: date of last commit
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

#include "utilities.h"


namespace rpwa {


	inline
	const libconfig::Setting*
	findGroup(const libconfig::Setting& parent,
	          const std::string&        groupName,
	          const bool                mustExist = true)  ///< finds field in keyfile and makes sure it is a group
	{
		// find field
		if (not parent.exists(groupName)) {
			if (mustExist)
				printWarn << "cannot find '" << groupName << "' field in '" << parent.getPath() << "' "
				          << "of key file '" << parent.getSourceFile() << "'" << std::endl;
			return 0;
		}
		const libconfig::Setting& groupKey = parent[groupName];
		// check that it is a group
		if (not groupKey.isGroup()) {
			printWarn << "'" << groupName << "' field in '" << parent.getPath() << "' "
			          << "of key file '" << parent.getSourceFile() << "' is not a group" << std::endl;
			return 0;
		}
		return &groupKey;
	}


	inline
	const libconfig::Setting*
	findList(const libconfig::Setting& parent,
	         const std::string&        listName,
	         const bool                mustExist = true)  ///< finds field in keyfile and makes sure it is a non-empty list
	{
		// find field
		if (not parent.exists(listName)) {
			if (mustExist)
				printWarn << "cannot find '" << listName << "' field in '" << parent.getPath() << "' "
				          << "of key file '" << parent.getSourceFile() << "'" << std::endl;
			return 0;
		}
		const libconfig::Setting& listKey = parent[listName];
		// check that it is a list
		if (not listKey.isList()) {
			printWarn << "'" << listName << "' field in '" << parent.getPath() << "' "
			          << "of key file '" << parent.getSourceFile() << "' is not a list. "
			          << "check that braces are correct." << std::endl;
			return 0;
		}
		// check that it is not empty
		if (listKey.getLength() < 1) {
			printWarn << "list '" << listName << "' in '" << parent.getPath() << "' "
			          << "of key file '" << parent.getSourceFile() << "' is empty" << std::endl;
			return 0;
		}
		return &listKey;
	}

	
	inline
	const libconfig::Setting*
	findArray(const libconfig::Setting& parent,
	          const std::string&        arrayName,
	          const bool                mustExist = true)  ///< finds field in keyfile and makes sure it is a non-empty array
	{
		// find field
		if (not parent.exists(arrayName)) {
			if (mustExist)
				printWarn << "cannot find '" << arrayName << "' field in '" << parent.getPath() << "' "
				          << "of key file '" << parent.getSourceFile() << "'" << std::endl;
			return 0;
		}
		const libconfig::Setting& arrayKey = parent[arrayName];
		// check that it is a list
		if (not arrayKey.isArray()) {
			printWarn << "'" << arrayName << "' field in '" << parent.getPath() << "' "
			          << "of key file '" << parent.getSourceFile() << "' is not an array. "
			          << "check that braces are correct." << std::endl;
			return 0;
		}
		// check that it is not empty
		if (arrayKey.getLength() < 1) {
			printWarn << "array '" << arrayName << "' in '" << parent.getPath() << "' "
			          << "of key file '" << parent.getSourceFile() << "' is empty" << std::endl;
			return 0;
		}
		return &arrayKey;
	}
  
  
} // namespace rpwa


#endif  // LIBCONFIGUTILS_H
