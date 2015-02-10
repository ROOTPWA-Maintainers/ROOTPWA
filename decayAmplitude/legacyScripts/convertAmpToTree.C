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
//      program that converts standard binary PWA2000 .amp files into
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

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TObjString.h"
#include "TVector3.h"

#include "reportingUtils.hpp"
#include "amplitudeTreeHelper.h"


using namespace std;
using namespace rpwa;


bool
convertAmpToTree(const string&  ampFileName  = "1-1++1+rho770_21_pi-.amp",
                 const string&  outFileName  = "testAmpTree.root",
                 const long int maxNmbEvents = -1,
                 const string&  ampLeafName  = "amplitude",
                 const bool     debug        = false)
{
	// create output file
	printInfo << "creating output file '" << outFileName << "'" << endl;
	TFile* outFile = TFile::Open(outFileName.c_str(), "RECREATE");
	if (not outFile) {
		printErr << "cannot open output file '" << outFileName << "'" << endl;
		return false;
	}
	outFile->SetCompressionLevel(9);

	// create tree
	TTree* outTree = new TTree();
	if (not outTree) {
		printErr << "problems creating tree in file '" << outFileName << "'" << endl;
		return false;
	}

	// doit
	const bool success = fillTreeFromAmp(ampFileName, *outTree, maxNmbEvents, ampLeafName, debug);
	outTree->Write();

	outFile->Close();
	if (success)
		printInfo << "wrote events to file '" << outFileName << "'" << endl;
	else
		printWarn << "problems processing events" << endl;
	return success;
}
