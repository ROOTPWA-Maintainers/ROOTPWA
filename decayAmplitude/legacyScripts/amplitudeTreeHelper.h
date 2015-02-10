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
//      helper functions that convert between standard binary PWA2000
//      .amp files and the new ROOT tree format
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <string>


class TTree;
class TChain;


namespace rpwa {


	bool fillTreeFromAmp(const std::string& inFileName,
	                     TTree&             outTree,
	                     const long int     maxNmbEvents  = -1,
	                     const std::string& ampLeafName   = "amplitude",
	                     const long int     treeCacheSize = 1000000,  // 1 MByte ROOT tree read cache
	                     const bool         debug         = false);


	bool writeAmpFromTree(TChain&            inTree,
	                      const std::string& outFileName,
	                      const long int     maxNmbEvents = -1,
	                      const std::string& ampLeafName  = "amplitude",
	                      const bool         debug        = false);


}  // namespace rpwa
