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
//      basic test program for wave set class
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "TROOT.h"
#include "TTree.h"

#include "reportingUtils.hpp"
#include "amplitudeTreeLeaf.h"
#include "particleDataTable.h"
#include "waveSet.h"


using namespace std;
using namespace rpwa;


int
main(int argc, char** argv)
{
	waveSet::setDebug(true);

	printCompilerInfo();
	printSvnVersion  ();

	if (1) {
		const string waveSetFileName = "testWaveSet.waveset";
		waveSet      set;
		set.buildWaveSet(waveSetFileName);
		printDebug << set;

		vector<string> fileNames(3);
		fileNames[0] = "testAmp1.root";
		fileNames[1] = "testAmp2.root";
		fileNames[2] = "testAmp3.root";
		set.setDecayAmpFileNames(fileNames);
		set.getDecayAmplitudeTrees();

		gROOT->ProcessLine("#include <complex>");
		const string   ampLeafName   = "decayAmp";
		const long int treeCacheSize = 1000000;
		const vector<TTree*>& ampTrees = set.decayAmpTrees();
		for (unsigned int i = 0; i < ampTrees.size(); ++i) {
			amplitudeTreeLeaf* ampLeaf = 0;
			ampTrees[i]->SetBranchAddress(ampLeafName.c_str(), &ampLeaf);
			ampTrees[i]->SetCacheSize(treeCacheSize);
			ampTrees[i]->AddBranchToCache(ampLeafName.c_str(), true);
			ampTrees[i]->StopCacheLearningPhase();
			for (unsigned int j = 0; j < ampTrees[i]->GetEntriesFast(); ++j) {
				ampTrees[i]->GetEntry(j);
				if (j < 5)
					printDebug << "read event " << j << " from tree '" << ampTrees[i]->GetName()
					           << "': " << *ampLeaf;
			}
			ampTrees[i]->PrintCacheStats();
			cout << endl;			
		}

		particleDataTable::readFile("./particleDataTable.txt");
		set.constructDecayAmps();
		const vector<isobarAmplitudePtr>& amps = set.decayAmps();
		for (unsigned int i = 0; i < amps.size(); ++i)
			printDebug << "[" << i << "] = " << *(amps[i]) << endl;
	}

}
