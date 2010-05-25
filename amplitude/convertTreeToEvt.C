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
//      program that converts standrad ASCII PWA2000 .evt files into
//      new ROOT tree format
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <fstream>
#include <string>
#include <cassert>
#include <iomanip>
#include <limits>
#include <algorithm>

#include "TChain.h"
#include "TClonesArray.h"
#include "TObjString.h"
#include "TVector3.h"

#include "pputil.h"
#include "utilities.h"
#include "particleDataTable.h"
#include "evtTreeHelper.h"


using namespace std;
using namespace rpwa;


bool
convertTreeToEvt(const string&  inFileNamePattern     = "testEvents.root",
		 const string&  outFileName           = "testTree.evt",
		 const long int maxNmbEvents          = -1,
		 const string&  pdgTableFileName      = "./particleDataTable.txt",
		 const string&  inTreeName            = "rootPwaEvtTree",
		 const string&  leafNameIsPartNames   = "initialStateNames",
		 const string&  leafNameIsPartMomenta = "initialStateMomenta",
		 const string&  leafNameFsPartNames   = "finalStateNames",
		 const string&  leafNameFsPartMomenta = "finalStateMomenta",
		 const bool     debug                 = false)
{
  // open input file
  printInfo << "opening input file(s) '" << inFileNamePattern << "'" << endl;
  TChain chain(inTreeName.c_str());
  if (chain.Add(inFileNamePattern.c_str()) < 1) {
    printWarn << "no events in input file(s) '" << inFileNamePattern << "'" << endl;
    return false;
  }
  chain.GetListOfFiles()->ls();

  // create output file
  printInfo << "creating output file '" << outFileName << "'" << endl;
  ofstream outFile(outFileName.c_str());
  if (!outFile) {
    printWarn << "cannot open output file '" << outFileName << "'" << endl;
    return false;
  }

  rpwa::particleDataTable& pdt = rpwa::particleDataTable::instance();
  pdt.readFile(pdgTableFileName);

  const bool success = writeEvtFromTree(chain, outFile, maxNmbEvents, inTreeName,
					leafNameIsPartNames, leafNameIsPartMomenta,
					leafNameFsPartNames, leafNameFsPartMomenta, debug);

  if (success)
    printInfo << "wrote events to file '" << outFileName << "'" << endl;
  else
    printWarn << "problems processing events" << endl;
  return success;
}
