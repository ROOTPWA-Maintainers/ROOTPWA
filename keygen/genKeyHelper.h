///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
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
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#ifndef KEYGENHELPER_HH
#define KEYGENHELPER_HH


#include <string>


namespace rpwa {


  class waveKey;


  bool testKeyFile(const std::string& keyFileName,                            // file name of key file under test
		   const int          refl,                                   // reflectivity of wave
		   const std::string& dataFileName     = "./testEvents.evt",  // file with test data in .evt format
		   const std::string& pdgTableFileName = "./pdgTable.txt",    // path to PDG table file
		   const double       precision        = 1e-12,               // warn threshold for comparison |1 - amp. / refl. amp.|
		   const bool         printAmp         = false);              // if true amplitude values are printed in addition to ratios amp. / refl. amp.
  
  
  void generateKeyFile(const rpwa::waveKey& wave,                                   // complete isobar decay spec
		       const std::string&   srcMacroFileName,                       // path to macro that will be copied to <wave name>.C
		       const bool           testKey          = true,                // if true amplitude behavior under reflectivity is tested
		       const std::string&   dataFileName     = "./testEvents.evt",  // file with test data in .evt format
		       const std::string&   pdgTableFileName = "./pdgTable.txt");   // path to PDG table file
  
  
}


#endif  // KEYGENHELPER_HH


//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
