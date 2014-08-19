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


#include <fstream>
#include <cassert>

#include <boost/progress.hpp>

#include "TTree.h"
#include "TChain.h"

#include "reportingUtils.hpp"
#include "fileUtils.hpp"
#include "amplitudeTreeLeaf.h"
#include "amplitudeTreeHelper.h"


using namespace std;
using namespace boost;


namespace rpwa {


	bool
	fillTreeFromAmp(const string&  inFileName,
	                TTree&         outTree,
	                const long int maxNmbEvents,
	                const string&  ampLeafName,
	                const long int treeCacheSize,
	                const bool)
	{
		// open input file
		printInfo << "opening input file '" << inFileName << "'" << endl;
		ifstream inFile(inFileName.c_str());
		if (not inFile or not inFile.good()) {
			printErr << "cannot open amplitude file '" << inFileName << "'." << endl;
			return false;
		}

		// create leaf variable and connect it to tree branch
		amplitudeTreeLeaf* ampTreeLeaf = new amplitudeTreeLeaf();
		const string       ampTreeName = fileNameFromPath(inFileName);
		printInfo << "writing to tree '" << ampTreeName << "'" << endl;
		outTree.SetName (ampTreeName.c_str());
		outTree.SetTitle(ampTreeName.c_str());
		const int splitLevel = 99;
		const int bufSize    = 256000;
		outTree.Branch(ampLeafName.c_str(), &ampTreeLeaf, bufSize, splitLevel);

		// loop over events and fill tree
		printInfo << "writing amplitudes..." << endl;
		long int         countEvents = 0;
		streampos        fileLength  = fileSize(inFile);
		streampos        lastPos     = inFile.tellg();
		progress_display progressIndicator(fileLength, cout, "");
		complex<double>  amp;
		while (inFile.read((char*)&amp, sizeof(complex<double>))) {
			ampTreeLeaf->setAmp(amp);
			outTree.Fill();
			++countEvents;
			progressIndicator += inFile.tellg() - lastPos;
			lastPos = inFile.tellg();
			if ((maxNmbEvents > 0) and (countEvents >= maxNmbEvents))
				break;
		}

		printInfo << "optimizing tree" << endl;
		//outTree.Print();
		outTree.OptimizeBaskets(treeCacheSize, 1, "d");
		//outTree.Print();

		printInfo << "wrote amplitudes for " << countEvents << " events to tree "
		          << "'" << outTree.GetName() << "'" << endl;
		return true;
	}


	bool
	writeAmpFromTree(TChain&            inTree,
	                 const std::string& outFileName,
	                 const long int     maxNmbEvents,
	                 const std::string& ampLeafName,
	                 const bool)
	{
		// create output file
		printInfo << "creating output file '" << outFileName << "'" << endl;
		ofstream outFile(outFileName.c_str());
		if (not outFile) {
			printWarn << "cannot create amplitude file '" << outFileName << "'." << endl;
			return false;
		}

		// connect leaf variable to tree branch
		amplitudeTreeLeaf* ampTreeLeaf = 0;
		inTree.SetBranchAddress(ampLeafName.c_str(), &ampTreeLeaf);

		// loop over events
		printInfo << "writing amplitudes..." << endl;
		bool             success       = true;
		const long int   nmbEventsTree = inTree.GetEntries();
		const long int   nmbEvents     = ((maxNmbEvents > 0) ? min(maxNmbEvents, nmbEventsTree)
		                                  : nmbEventsTree);
		progress_display progressIndicator(nmbEvents, cout, "");
		for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
			++progressIndicator;

			if (inTree.LoadTree(eventIndex) < 0)
				break;
			inTree.GetEntry(eventIndex);
			assert(ampTreeLeaf);

			if (ampTreeLeaf->nmbIncohSubAmps() == 0) {
				printWarn << "amplitude leaf for event " << eventIndex << " is empty. skipping." << endl;
				success = false;
				continue;
			}
			if (ampTreeLeaf->nmbIncohSubAmps() > 1)
				printWarn << "amplitude leaf for event " << eventIndex
				          << " has a size of " << ampTreeLeaf->nmbIncohSubAmps()
				          << " writing only first amplitude."<< endl;

			complex<double> amp = ampTreeLeaf->amp();
			outFile.write((char*)(&amp), sizeof(complex<double>));
		}

		outFile.close();
		printInfo << "wrote amplitudes for " << nmbEvents << " events to output file "
		          << "'" << outFileName << "'" << endl;
		return success;
	}


}  // namespace rpwa
