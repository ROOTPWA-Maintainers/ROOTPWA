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


#include<map>
#include "libconfig.h++"

#include "reportingUtils.hpp"


namespace rpwa {
	typedef boost::shared_ptr<libconfig::Config> configPtr;

	inline
	std::string
	getLibConfigSourceFilePath(const libconfig::Setting& setting)
	{
		std::string sourceFile = "";
		const char* sourceFilePtr = setting.getSourceFile();
		if(sourceFilePtr) {
			sourceFile = sourceFilePtr;
		}
		return sourceFile;
	}


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
	libconfig::Setting*
	findLibConfigGroup(const libconfig::Setting& parent,
	                   const std::string&        groupName,
	                   const bool                mustExist = true)  ///< finds field in keyfile and makes sure it is a group
	{
		// find field
		std::string sourceFile = getLibConfigSourceFilePath(parent);
		if (not parent.exists(groupName)) {
			if (mustExist)
				printWarn << "cannot find '" << groupName << "' field in '" << parent.getPath() << "' "
				          << "of key file '" << sourceFile << "'" << std::endl;
			return 0;
		}
		libconfig::Setting& groupKey = parent[groupName];
		// check that it is a group
		if (not groupKey.isGroup()) {
			printWarn << "'" << groupName << "' field in '" << parent.getPath() << "' "
			          << "of key file '" << sourceFile << "' is not a group" << std::endl;
			return 0;
		}
		return &groupKey;
	}


	inline
	const libconfig::Setting*
	findLibConfigList(const libconfig::Setting& parent,
	                  const std::string&        listName,
	                  const bool                mustExist = true)  ///< finds field in keyfile and makes sure it is a non-empty list
	{
		// find field
		std::string sourceFile = getLibConfigSourceFilePath(parent);
		if (not parent.exists(listName)) {
			if (mustExist)
				printWarn << "cannot find '" << listName << "' field in '" << parent.getPath() << "' "
				          << "of key file '" << sourceFile << "'" << std::endl;
			return 0;
		}
		const libconfig::Setting& listKey = parent[listName.c_str()];
		// check that it is a list
		if (not listKey.isList()) {
			printWarn << "'" << listName << "' field in '" << parent.getPath() << "' "
			          << "of key file '" << sourceFile << "' is not a list. "
			          << "check that braces are correct." << std::endl;
			return 0;
		}
		// check that it is not empty
		if (listKey.getLength() < 1) {
			printWarn << "list '" << listName << "' in '" << parent.getPath() << "' "
			          << "of key file '" << sourceFile << "' is empty" << std::endl;
			return 0;
		}
		return &listKey;
	}


	inline
	const libconfig::Setting*
	findLibConfigArray(const libconfig::Setting& parent,
	                   const std::string&        arrayName,
	                   const bool                mustExist = true)  ///< finds field in keyfile and makes sure it is a non-empty array
	{
		// find field
		std::string sourceFile = getLibConfigSourceFilePath(parent);
		if (not parent.exists(arrayName)) {
			if (mustExist)
				printWarn << "cannot find '" << arrayName << "' field in '" << parent.getPath() << "' "
				          << "of key file '" << sourceFile << "'" << std::endl;
			return 0;
		}
		const libconfig::Setting& arrayKey = parent[arrayName.c_str()];
		// check that it is a list
		if (not arrayKey.isArray()) {
			printWarn << "'" << arrayName << "' field in '" << parent.getPath() << "' "
			          << "of key file '" << sourceFile << "' is not an array. "
			          << "check that braces are correct." << std::endl;
			return 0;
		}
		// check that it is not empty
		if (arrayKey.getLength() < 1) {
			printWarn << "array '" << arrayName << "' in '" << parent.getPath() << "' "
			          << "of key file '" << sourceFile << "' is empty" << std::endl;
			return 0;
		}
		return &arrayKey;
	}


	inline
	const libconfig::Setting*
	findLibConfigGroup(const libconfig::Setting& parent,
	                   const int                 groupIndex,
	                   const bool                mustExist = true)  ///< finds field in keyfile and makes sure it is a group
	{
		// find field

		const libconfig::Setting* groupKey = 0;
		try {
			groupKey = &(parent[groupIndex]);
		} catch (const libconfig::SettingNotFoundException& NotFound) {
			if (mustExist)
				printWarn << "cannot find group[" << groupIndex << "] in '" << parent.getPath() << "' "
				          << "of key file '" << parent.getSourceFile() << "'" << std::endl;
			return 0;
		}

		// check that it is a group
		if (not (*groupKey).isGroup()) {
			printWarn << "setting[" << groupIndex << "] in '" << parent.getPath() << "' "
			          << "of key file '" << parent.getSourceFile() << "' is not a group" << std::endl;
			return 0;
		}
		return groupKey;
	}


	inline
	const libconfig::Setting*
	findLibConfigList(const libconfig::Setting& parent,
	                  const int                 listIndex,
	                  const bool                mustExist = true)  ///< finds field in keyfile and makes sure it is a non-empty list
	{
		// find field

		const libconfig::Setting* listKey = 0;
		try {
			listKey = &(parent[listIndex]);
		} catch (const libconfig::SettingNotFoundException& NotFound) {
			if (mustExist)
				printWarn << "cannot find list[" << listIndex << "] in '" << parent.getPath() << "' "
				          << "of key file '" << parent.getSourceFile() << "'" << std::endl;
			return 0;
		}
		// check that it is a list
		if (not (*listKey).isList()) {
			printWarn << "setting[" << listIndex << "] in '" << parent.getPath() << "' "
			          << "of key file '" << parent.getSourceFile() << "' is not a list. "
			          << "check that braces are correct." << std::endl;
			return 0;
		}

		// check that it is not empty
		if ((*listKey).getLength() < 1) {
			printWarn << "list[" << listIndex << "] in '" << parent.getPath() << "' "
			          << "of key file '" << parent.getSourceFile() << "' is empty" << std::endl;
			return 0;
		}
		return listKey;
	}


	inline
	const libconfig::Setting*
	findLibConfigArray(const libconfig::Setting& parent,
	                   const int                 arrayIndex,
	                   const bool                mustExist = true)  ///< finds field in keyfile and makes sure it is a non-empty array
	{
		// find field

		const libconfig::Setting* arrayKey = 0;
		try {
			arrayKey = &(parent[arrayIndex]);
		} catch (const libconfig::SettingNotFoundException& NotFound) {
			if (mustExist)
				printWarn << "cannot find array[" << arrayIndex << "] in '" << parent.getPath() << "' "
				          << "of key file '" << parent.getSourceFile() << "'" << std::endl;
			return 0;
		}

		// check that it is a list
		if (not (*arrayKey).isArray()) {
			printWarn << "setting[" << arrayIndex << "] in '" << parent.getPath() << "' "
			          << "of key file '" << parent.getSourceFile() << "' is not an array. "
			          << "check that braces are correct." << std::endl;
			return 0;
		}
		// check that it is not empty
		if ((*arrayKey).getLength() < 1) {
			printWarn << "array[" << arrayIndex << "] in '" << parent.getPath() << "' "
			          << "of key file '" << parent.getSourceFile() << "' is empty" << std::endl;
			return 0;
		}
		return arrayKey;
	}


	inline
	bool
	checkIfAllVariablesAreThere(const libconfig::Setting* group,
	                            const std::map<std::string, libconfig::Setting::Type>& nameMap)
	{
		for(std::map<std::string, libconfig::Setting::Type>::const_iterator it = nameMap.begin(); it != nameMap.end(); ++it) {
			std::string name = it->first;
			libconfig::Setting::Type type = it->second;
			if(not group->exists(name)) {
				printWarn << "'" << name << "' not found "
				          << "in '" << (group->getName() ? group->getName() : "???") << "' section." << std::endl;
				return false;
			}
			if((*group)[name.c_str()].getType() != type) {
				printWarn << "'" << name << "' in section '" << (group->getName() ? group->getName() : "???") << "' has wrong type." << std::endl;
				return false;
			}
		}
		return true;
	}

	inline
	std::string
	getConfigString(const configPtr cfgIn)
	{
		int pipeFileDescriptors[2];
		if (pipe(pipeFileDescriptors) == -1) {
			printErr << "failed to create pipe. cannot write keys." << std::endl;
			return "";
		};
		FILE* pipeWriteEnd = fdopen(pipeFileDescriptors[1], "wt");
		if (!pipeWriteEnd) {
			printErr << "could not open write end of pipe. cannot write keys." << std::endl;
			return "";
		};
		cfgIn->write(pipeWriteEnd);
		fclose(pipeWriteEnd);
		// read keys from pipe
		char         buf;
		unsigned int countChar = 0;
		std::stringstream out;
		while (read(pipeFileDescriptors[0], &buf, 1) > 0) {
			out << buf;
			++countChar;
		}
		close(pipeFileDescriptors[0]);
		if (countChar == 0) {
			printWarn << "nothing was written" << std::endl;
		};
		return out.str();
	};

	inline
	configPtr
	copyConfig(const configPtr cfgIn)
	{
		configPtr cfg(new libconfig::Config);
		cfg->readString(getConfigString(cfgIn));
		return cfg;
	};




} // namespace rpwa


#endif  // LIBCONFIGUTILS_H
