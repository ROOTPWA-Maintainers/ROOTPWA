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
//      test program for amplitude persistency class
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <boost/assign/list_of.hpp>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

#include "reportingUtils.hpp"
#include "amplitudeTreeLeaf.h"


using namespace std;
using namespace boost::assign;
using namespace rpwa;


int
main(int argc, char** argv)
{
	printCompilerInfo();
	printSvnVersion();

#ifdef USE_STD_COMPLEX_TREE_LEAFS
	
	const unsigned int nmbEvents       = 1000000;
	const unsigned int nmbIncohSubAmps = 3;
	
	// force loading predefined std::complex dictionary
	// see http://root.cern.ch/phpBB3/viewtopic.php?f=5&t=9618&p=50164
	gROOT->ProcessLine("#include <complex>");
	
	if (1) {
		TFile*               outFile      = TFile::Open("testAmplitudeTree.root", "RECREATE");
		amplitudeTreeLeaf*   ampLeaf      = new amplitudeTreeLeaf();
		const vector<string> subAmpLabels = list_of("lambda=-1")("lambda=0")("lambda=+1");
		TTree*               tree         = new TTree("test", "test");

		waveDescription waveDesc;
		waveDesc.parseKeyFile("testWaveDescription.key");
		tree->GetUserInfo()->AddFirst(&waveDesc);

		gRandom->SetSeed(123456789);
		tree->Branch("amp", &ampLeaf, 256000, 99);
		for (unsigned int i = 0; i < nmbEvents; ++i) {
			ampLeaf->clear();
			ampLeaf->defineIncohSubAmps(subAmpLabels);
			for (unsigned int j = 0; j < nmbIncohSubAmps; ++j) {
				// ampLeaf->setIncohSubAmp(complex<double>(i, -(double)j), j);
				ampLeaf->setIncohSubAmp(complex<double>(gRandom->Rndm(), gRandom->Rndm()), j);
				//ampLeaf->_waveDesc = &waveDesc;
			}
			tree->Fill();
			if (i < 5)
				cout << "written event " << i << ": " << *ampLeaf;
		}
		cout << endl;
		tree->Print();
 		tree->OptimizeBaskets(1000000, 1, "d");
		tree->Write();
		tree->Print();
		outFile->Close();
		for (unsigned int i = 0; i < subAmpLabels.size(); ++i)
			cout << subAmpLabels[i] << ": ["  << ampLeaf->incohSubAmpIndex(subAmpLabels[i]) << "]" << endl;
		cout << endl;
	}

	if (1) {
		TFile* inFile = TFile::Open("testAmplitudeTree.root", "READ");

		gRandom->SetSeed(123456789);
		TTree* tree;
		inFile->GetObject("test", tree);
		tree->GetUserInfo()->Print();
		tree->GetUserInfo()->ls();

		TList*           userInfo = tree->GetUserInfo();
		waveDescription* waveDesc = static_cast<waveDescription*>(userInfo->FindObject("rpwa::waveDescription"));
		waveDesc->printKeyFileContents();
		
		amplitudeTreeLeaf* ampLeaf = 0;
		tree->SetBranchAddress("amp", &ampLeaf);
		tree->SetCacheSize(1000000);
		tree->AddBranchToCache("amp",  true);
		tree->StopCacheLearningPhase();
		for (unsigned int i = 0; i < tree->GetEntriesFast(); ++i) {
			tree->GetEntry(i);
			if (i < 5)
				cout << "read event " << i << ": " << *ampLeaf;
			for (unsigned int j = 0; j < ampLeaf->nmbIncohSubAmps(); ++j) {
				const complex<double> amp      = ampLeaf->incohSubAmp(j);
				const complex<double> expected = complex<double>(gRandom->Rndm(), gRandom->Rndm());
				if (amp != expected)
					printWarn << "mismatch: read " << amp << " expected " << expected << endl;
			}
		}
		tree->PrintCacheStats();
		cout << endl;
		
		// test arthmetic functions
		printInfo << "original: " << *ampLeaf << endl;
		amplitudeTreeLeaf ampLeaf2(*ampLeaf);
		printInfo << "copy: "<< ampLeaf2 << endl;
		if (ampLeaf2 != *ampLeaf)
			printErr << "problem with assignment" << endl;
		amplitudeTreeLeaf ampLeaf3 = 0.1 * (*ampLeaf) + 0.9 * ampLeaf2;
		printInfo << "arithmetic test: " << ampLeaf3 << endl;
		if (ampLeaf3 != *ampLeaf)
			printErr << "problem with arithmetic" << endl;
	}

#endif  // USE_STD_COMPLEX_TREE_LEAFS

}
