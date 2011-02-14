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
// $Rev:: 503                         $: revision of last commit
// $Author:: bgrube                   $: author of last commit
// $Date:: 2010-11-22 12:28:07 +0100 #$: date of last commit
//
// Description:
//      some common routines for file (name) handling
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef FILEUTILS_HPP
#define FILEUTILS_HPP


#include <string>
#include <iostream>
#include <vector>
#ifndef __CINT__
#include <glob.h>
#else
struct glob_t;
int glob(const char *,
         int,
         int(*)(const char*, int),
         glob_t*);
void globfree(glob_t *);
#endif 


namespace rpwa {


	inline
	std::string
	directoryFromPath(const std::string& path)
	{
	  return path.substr(0, path.find_last_of('/') + 1);
  }


	inline
	std::string
	fileNameFromPath(const std::string& path)
	{
		return path.substr(path.find_last_of('/') + 1);
  }


	inline
	std::string
	fileNameNoExtFromPath(const std::string& path)
	{
		const std::string fileName = fileNameFromPath(path);
		return fileName.substr(0, fileName.find_last_of('.'));
  }


	inline
	std::string
	extensionFromPath(const std::string& path)
	{
		const std::string fileName = fileNameFromPath(path);
		return fileName.substr(fileName.find_last_of('.') + 1);
  }


	inline
	std::string
	changeFileExtension(const std::string& path,
	                    const std::string& newExt = "")
  {
	  const std::string directory     = directoryFromPath    (path);
	  const std::string fileNameNoExt = fileNameNoExtFromPath(path);
	  if (newExt != "")
		  return directory + fileNameNoExt + ((newExt[0] == '.') ? "" : ".") + newExt;
	  else
		  return directory + fileNameNoExt;
  }


	inline
	std::streampos
	fileSize(std::istream& file)
	{
		if (not file)
			return 0;
		const std::streampos currentPos = file.tellg();
		file.seekg(0, std::ios::end);
		const std::streampos size = file.tellg();
		file.seekg(currentPos);
		return size;
	}


	// expands glob pattern into list of file names
	inline
	std::vector<std::string>
	filesMatchingGlobPattern(const std::string& globPattern,
	                         const bool         sortList = false)
	{
		std::vector<std::string> fileList;
		glob_t globBuffer;
		if (sortList)
			glob(globPattern.c_str(), 0,           NULL, &globBuffer);
		else
			glob(globPattern.c_str(), GLOB_NOSORT, NULL, &globBuffer);
		for (unsigned int i = 0; i < globBuffer.gl_pathc; ++i)
			fileList.push_back(globBuffer.gl_pathv[i]);
		globfree(&globBuffer);
		return fileList;
	} 


}  // namespace rpwa


#endif  // FILEUTILS_HPP
