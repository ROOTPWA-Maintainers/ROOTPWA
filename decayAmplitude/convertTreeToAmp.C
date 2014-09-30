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
#include <iomanip>
#include <limits>
#include <algorithm>

#include "TChain.h"
#include "TClonesArray.h"
#include "TObjString.h"
#include "TVector3.h"

#include "reportingUtils.hpp"
#include "amplitudeTreeHelper.h"


using namespace std;
using namespace rpwa;


bool
convertTreeToAmp(const string&  inFileNamePattern = "testAmpTree.root",
                 const string&  inTreeName        = "1-1++1+rho770_21_pi-.amp",
                 const string&  outFileName       = "1-1++1+rho770_21_pi-.test.amp",
                 const long int maxNmbEvents      = -1,
                 const string&  ampLeafName       = "amplitude",
                 const bool     debug             = false)
{
	// open input file
	printInfo << "opening input file(s) '" << inFileNamePattern << "'" << endl;
	TChain chain(inTreeName.c_str());
	if (chain.Add(inFileNamePattern.c_str()) < 1) {
		printWarn << "no events in input file(s) '" << inFileNamePattern << "'" << endl;
		return false;
	}
	chain.GetListOfFiles()->ls();
	printInfo << "reading from tree '" << inTreeName << "'" << endl;

	// doit
	const bool success = writeAmpFromTree(chain, outFileName, maxNmbEvents, ampLeafName, debug);

	if (success)
		printInfo << "wrote events to file '" << outFileName << "'" << endl;
	else
		printWarn << "problems processing events" << endl;
	return success;
}
