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
#include "TFile.h"

#include "reportingUtils.hpp"
#include "amplitudeTreeLeaf.h"
#include "particleDataTable.h"
#include "waveSet.h"


using namespace std;
using namespace rpwa;


int
main(int argc, char** argv)
{
	using rpwa::cout;
	printCompilerInfo();
	printLibraryInfo ();
	printSvnVersion  ();
	cout << endl;

	waveSet::setDebug(true);

	if (1) {
		const string waveSetFileName = "testWaveSet.waveset";
		waveSet      set;
		set.parseWaveSetFile(waveSetFileName);

		printDebug << set;

		vector<string> fileNames(3);
		fileNames[0] = "testAmp1.root";
		fileNames[1] = "testAmp2.root";
		fileNames[2] = "testAmp3.root";
		set.setDecayAmpFileNames(fileNames);
		set.getDecayAmpTrees();
		set.getWaveDescs();

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
		// const vector<isobarAmplitudePtr>& amps = set.decayAmps();
		// for (unsigned int i = 0; i < amps.size(); ++i)
		// 	printDebug << "amplitude[" << i << "] = " << *(amps[i]) << endl;

		set.constructWaveNames();
		const vector<waveName>& waveNames = set.waveNames();
		for (unsigned int i = 0; i < waveNames.size(); ++i)
			printDebug << "wave name[" << i << "] = " << waveNames[i] << endl;
		
		printDebug << set;

		const string waveSetName = "testWaveSet";
		{
			printDebug << "write test" << endl;
			TFile* f = TFile::Open("testWaveSet.root", "RECREATE");
			set.Write(waveSetName.c_str());
			f->Close();
			delete f;
		}
		set.clear();
		printDebug << "after clear(): " << set;
		{
			TFile*   f    = TFile::Open("testWaveSet.root", "READ");
			waveSet* set2 = 0;
			f->GetObject(waveSetName.c_str(), set2);
			if (not set2)
				printErr << "cannot find wave set '" << waveSetName << "'" << endl;
			else {
				printInfo << "wave set read from file: " << *set2;
			}
			f->Close();
			delete f;
		}
		
	}

}
