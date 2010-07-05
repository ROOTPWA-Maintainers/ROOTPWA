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
//      example macro to create 5pi key file for gamp     
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#include <string>

#include "particleKey.h"
#include "waveKey.h"

#include "genKeyHelper.h"


using namespace std;
using namespace rpwa;


void
genKey(const bool    testKey          = true,
       const string& dataFileName     = "./testEvents.evt",  // file with test data in .evt format
       const string& pdgTableFileName = "./pdgTable.txt")
{
  // define final state particles
  particleKey p1("pi-");
  particleKey p2("pi+");
  particleKey p3("pi-");
  particleKey p4("pi+");
  particleKey p5("pi-");
  
  // define isobars: (name, daughter1, daughter2, L, S, mass dependence)
  particleKey i11("sigma",    &p1, &p2,  0, 0, "amp_ves");
  particleKey i1 ("a1(1269)", &p4, &i11, 1);
  particleKey i2 ("f1(1285)", &p3, &i1,  1, 1);
  const int J    =  2;
  const int P    = +1;
  const int M    =  1;
  const int refl = +1;
  const int L    =  1; 
  const int S    =  1;

  // define X system
  particleKey X("X", &p5, &i2, L, S);
  waveKey     wave(&X, J, P, M, refl);

  // test and generate key file
  const string thisFilePath = __FILE__;
  generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
}
