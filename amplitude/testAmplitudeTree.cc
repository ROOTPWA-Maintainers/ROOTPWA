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


#include "TFile.h"
#include "TTree.h"

#include "svnVersion.h"
#include "amplitudeTreeLeaf.h"


using namespace std;
using namespace rpwa;


int
main(int argc, char** argv)
{
	printCompilerInfo();
	printSvnVersion();

#if AMPLITUDETREELEAF_ENABLED
	
	const unsigned int nmbEvents                 = 10;
	const unsigned int nmbFsSpinProjCombinations = 3;
	
	if (1) {
		TFile*             outFile = TFile::Open("testAmplitudeTree.root", "RECREATE");
		TTree*             tree    = new TTree("test", "test");
		amplitudeTreeLeaf* ampLeaf = new amplitudeTreeLeaf();
		tree->Branch("amp", &ampLeaf);
		// ampLeaf->setNmbFsSpinProjCombinations(nmbFsSpinProjCombinations);
		ampLeaf->_foo.resize(nmbFsSpinProjCombinations, 0);
		for (unsigned int i = 0; i < nmbEvents; ++i) {
			for (unsigned int j = 0; j < nmbFsSpinProjCombinations; ++j) {
				// ampLeaf->setSubAmp(complex<double>(i, -(double)j), j);
				ampLeaf->_foo[j] = (double)i - j;
			}
			tree->Fill();
			cout << "event " << i << ": ";
			// for (unsigned int j = 0; j < ampLeaf->nmbFsSpinProjCombinations(); ++j)
			// 	cout << ampLeaf->subAmp(j) << "   ";
			// cout << "|   ";
			for (unsigned int j = 0; j < nmbFsSpinProjCombinations; ++j)
				//for (unsigned int j = 0; j < ampLeaf->nmbFsSpinProjCombinations(); ++j)
				cout << ampLeaf->_foo[j] << "   ";
			cout << endl;
		}
		outFile->Write();
		outFile->Close();
	}

	if (1) {
		TFile* inFile = TFile::Open("testAmplitudeTree.root", "READ");
		TTree* tree;
		inFile->GetObject("test", tree);
		amplitudeTreeLeaf* ampLeaf = 0;
		tree->SetBranchAddress("amp", &ampLeaf);
		for (unsigned int i = 0; i < tree->GetEntriesFast(); ++i) {
			tree->GetEntry(i);
			cout << "event " << i << ": ";
			// for (unsigned int j = 0; j < ampLeaf->nmbFsSpinProjCombinations(); ++j)
			// 	cout << ampLeaf->subAmp(j) << "   ";
			// cout << "|   ";
			for (unsigned int j = 0; j < nmbFsSpinProjCombinations; ++j)
				//for (unsigned int j = 0; j < ampLeaf->nmbFsSpinProjCombinations(); ++j)
				cout << ampLeaf->_foo[j] << "   ";
			cout << endl;
		}
	}

#endif  // AMPLITUDETREELEAF_ENABLED

}
