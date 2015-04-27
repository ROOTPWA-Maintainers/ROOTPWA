///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2015 Sebastian Uhl (TUM)
//
//    This file is part of ROOTPWA
//
//    ROOTPWA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ROOTPWA is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ROOTPWA. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      collection of useful routines around yaml-cpp library
//
//-------------------------------------------------------------------------


#ifndef YAMLCPPUTILS_H
#define YAMLCPPUTILS_H


#include <yaml-cpp/yaml.h>

#include "reportingUtils.hpp"


namespace rpwa {


	namespace YamlCppUtils {


		enum Type {
			TypeBoolean,
			TypeFloat,
			TypeInt,
			TypeSequence,
			TypeString
		};


		inline
		std::string
		TypeToStr(const Type type) {
			switch (type) {
				case TypeBoolean:  return "boolean";
				case TypeFloat:    return "floating point number";
				case TypeInt:      return "integer";
				case TypeSequence: return "YAML sequence";
				case TypeString:   return "string";
			}
			return "unknown";
		}


		inline
		bool
		parseYamlFile(const std::string& yamlFileName,
		              YAML::Node&        config,
		              const bool         debug)  ///< reads a yaml file into the config object
		{
			if (debug)
				printDebug << "parsing yaml file '" << yamlFileName << "'" << std::endl;
			try {
				config = YAML::LoadFile(yamlFileName);
			} catch (const YAML::BadFile& badFileEx) {
				printWarn << "I/O error while reading yaml file "
				          << "'" << yamlFileName << "': " << badFileEx.msg << std::endl;
				return false;
			} catch (const YAML::ParserException& parserEx) {
				printWarn << "parse error in '" << yamlFileName << "' line " << parserEx.mark.line+1
				          << ", column " << parserEx.mark.column+1 << ": " << parserEx.msg << std::endl;
				return false;
			}
			return true;
		}


		inline
		bool
		checkVariableType(const YAML::Node&        node,
		                  const YamlCppUtils::Type type)
		{
			switch (type) {
				case YamlCppUtils::TypeBoolean:
					try {
						node.as<bool>();
						return true;
					} catch (const YAML::TypedBadConversion<bool>&) {
						return false;
					}
				case YamlCppUtils::TypeFloat:
					try {
						node.as<double>();
						return true;
					} catch (const YAML::TypedBadConversion<double>&) {
						return false;
					}
				case YamlCppUtils::TypeInt:
					try {
						node.as<int>();
						return true;
					} catch (const YAML::TypedBadConversion<int>&) {
						return false;
					}
				case YamlCppUtils::TypeSequence:
					return node.IsSequence();
				case YamlCppUtils::TypeString:
					try {
						node.as<std::string>();
						return true;
					} catch (const YAML::TypedBadConversion<std::string>&) {
						return false;
					}
			}
			return false;
		}


		inline
		bool
		checkIfAllVariablesAreThere(const YAML::Node&                                parent,
		                            const std::map<std::string, YamlCppUtils::Type>& nameMap)
		{
			for(std::map<std::string, YamlCppUtils::Type>::const_iterator it = nameMap.begin(); it != nameMap.end(); ++it) {
				const std::string& name = it->first;

				const YAML::Node& node = parent[name];
				if(not node) {
					printWarn << "required variable '" << name << "' not found." << std::endl;
					return false;
				}

				const YamlCppUtils::Type type = it->second;
				if(not checkVariableType(node, type)) {
					printWarn << "type of variable '" << name << "' is not correct (expected " << YamlCppUtils::TypeToStr(type) << ")." << std::endl;
					return false;
				}
			}
			return true;
		}


	}  // namespace YamlCppUtils


} // namespace rpwa


#endif  // YAMLCPPUTILS_H
