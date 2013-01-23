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
	printGitHash();

#ifdef USE_STD_COMPLEX_TREE_LEAFS

	const unsigned int nmbEvents       = 1000000;
	const unsigned int nmbIncohSubAmps = 3;
	gRandom->SetSeed(123456789);

	// force loading predefined std::complex dictionary
	// see http://root.cern.ch/phpBB3/viewtopic.php?f=5&t=9618&p=50164
	gROOT->ProcessLine("#include <complex>");

	if (1) {
		TFile*               outFile      = TFile::Open("testAmplitudeTree.root", "RECREATE");
		amplitudeTreeLeaf*   ampLeaf      = new amplitudeTreeLeaf();
		const vector<string> subAmpLabels = list_of("lambda=-1")("lambda=0")("lambda=+1");
		TTree*               tree         = new TTree("test", "test");

		tree->Branch("amp", &ampLeaf);
		for (unsigned int i = 0; i < nmbEvents; ++i) {
			ampLeaf->clear();
			ampLeaf->defineIncohSubAmps(subAmpLabels);
			for (unsigned int j = 0; j < nmbIncohSubAmps; ++j)
				// ampLeaf->setIncohSubAmp(complex<double>(i, -(double)j), j);
				ampLeaf->setIncohSubAmp(complex<double>(gRandom->Rndm(), gRandom->Rndm()), j);
			tree->Fill();
			if (i < 5)
				cout << "written event " << i << ": " << *ampLeaf;
		}
		cout << endl;
		tree->Write();
		outFile->Close();
		for (unsigned int i = 0; i < subAmpLabels.size(); ++i)
			cout << subAmpLabels[i] << ": ["  << ampLeaf->incohSubAmpIndex(subAmpLabels[i]) << "]" << endl;
		cout << endl;
	}

	if (1) {
		TFile* inFile = TFile::Open("testAmplitudeTree.root", "READ");
		TTree* tree;
		inFile->GetObject("test", tree);
		amplitudeTreeLeaf* ampLeaf = 0;
		tree->SetBranchAddress("amp", &ampLeaf);
		for (unsigned int i = 0; i < tree->GetEntriesFast(); ++i) {
			tree->GetEntry(i);
			if (i < 5)
				cout << "read event " << i << ": " << *ampLeaf;
		}
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
