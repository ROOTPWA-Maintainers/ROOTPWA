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
#include <sstream>
#include <string>
#include <cassert>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TObjString.h"
#include "TVector3.h"

#include "pputil.h"
#include "utilities.h"
#include "evtTreeHelper.h"


using namespace std;
using namespace rpwa;


bool
convertEvtToTree(//const string&  evtFileName               = "1500.1540.3pi.evt",
                 const string&  evtFileName               = "1720.1840.5pi.evt",
                 const string&  outFileName               = "testEvents.root",
                 //const string&  evtFileName               = "testTree.evt",
                 //const string&  outFileName               = "testEvents2.root",
                 const long int maxNmbEvents              = -1,
                 const string&  outTreeName               = "rootPwaEvtTree",
                 const string&  prodKinParticlesLeafName  = "prodKinParticles",
                 const string&  prodKinMomentaLeafName    = "prodKinMomenta",
                 const string&  decayKinParticlesLeafName = "decayKinParticles",
                 const string&  decayKinMomentaLeafName   = "decayKinMomenta",
                 const string&  targetParticleName        = "p+",
                 const bool     debug                     = false)
{
	// open input file
	printInfo << "opening input file '" << evtFileName << "'" << endl;
	ifstream evtFile(evtFileName.c_str());
	if (not evtFile or not evtFile.good()) {
		printWarn << "cannot open input file '" << evtFileName << "'" << endl;
		return false;
	}

	// create output file
	printInfo << "creating output file '" << outFileName << "'" << endl;
	TFile* outFile = TFile::Open(outFileName.c_str(), "RECREATE");
	if (not outFile) {
		printErr << "cannot open output file '" << outFileName << "'" << endl;
		return false;
	}

	// create tree
	TTree* tree = new TTree(outTreeName.c_str(), outTreeName.c_str());
	if (not tree) {
		printErr << "problems creating tree '" << outTreeName << "' "
		         << "in file '" << outFileName << "'" << endl;
		return false;
	}

	// doit
	const bool success = fillTreeFromEvt(evtFile, *tree, maxNmbEvents,
	                                     prodKinParticlesLeafName,  prodKinMomentaLeafName,
	                                     decayKinParticlesLeafName, decayKinMomentaLeafName,
	                                     targetParticleName, debug);
	tree->Write();
	tree->OptimizeBaskets(10000000, 1.1, "d");

	outFile->Close();
	if (success)
		printInfo << "wrote events to file '" << outFileName << "'" << endl;
	else
		printWarn << "problems processing events" << endl;
	return success;
}
