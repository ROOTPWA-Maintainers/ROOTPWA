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

#include "TTree.h"
#include "TClonesArray.h"
#include "TObjString.h"
#include "TVector3.h"

#include "reportingUtils.hpp"
#ifndef __CINT__
#include "particleDataTable.h"
#endif
#include "evtTreeHelper.h"


using namespace std;
using namespace rpwa;


bool
convertTreeToEvt(const string&  inFileNamePattern        = "testEvents.root",
                 const string&  outFileName              = "testEvents.evt",
                 const string&  pdgTableFileName         = "./particleDataTable.txt",
                 const long int maxNmbEvents             = -1,
                 const string&  inTreeName               = "rootPwaEvtTree",
                 const string&  prodKinPartNamesObjName  = "prodKinParticles",
                 const string&  prodKinMomentaLeafName   = "prodKinMomenta",
                 const string&  decayKinPartNamesObjName = "decayKinParticles",
                 const string&  decayKinMomentaLeafName  = "decayKinMomenta",
                 const bool     debug                    = false)
{
	// open input files
	TTree*        inTree            = 0;
	TClonesArray* prodKinPartNames  = 0;
	TClonesArray* decayKinPartNames = 0;
	{
		vector<string> rootFileNames;
		rootFileNames.push_back(inFileNamePattern);
		vector<string> evtFileNames;
		vector<TTree*> inTrees;
		if (not openRootEvtFiles(inTrees, prodKinPartNames, decayKinPartNames,
		                         rootFileNames, evtFileNames,
		                         inTreeName, prodKinPartNamesObjName, prodKinMomentaLeafName,
		                         decayKinPartNamesObjName, decayKinMomentaLeafName, debug)) {
			printErr << "problems opening input file(s). exiting." << endl;
			return false;
		}
		inTree = inTrees[0];
	}

	// create output file
	printInfo << "creating output file '" << outFileName << "'" << endl;
	ofstream outFile(outFileName.c_str());
	if (!outFile) {
		printWarn << "cannot open output file '" << outFileName << "'. exiting." << endl;
		return false;
	}

	rpwa::particleDataTable& pdt = rpwa::particleDataTable::instance();
	pdt.readFile(pdgTableFileName);

	// doit
	const bool success = writeEvtFromTree(*inTree, outFile, *prodKinPartNames, *decayKinPartNames,
	                                      maxNmbEvents, inTreeName,
	                                      prodKinMomentaLeafName, decayKinMomentaLeafName, debug);

	if (success)
		printSucc << "wrote events to file '" << outFileName << "'" << endl;
	else
		printWarn << "problems processing events" << endl;
	return success;
}
