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
//      example macro that reads uDST tree and generates mass bin
//      directory structure with .root files in ROOTPWA format
//      parts that depend on the analyzed final state and the uDST
//      tree structure are marked accordingly
//
//      to run this macro type
//      root rootlogon.C fillUdstDataIntoMassBins_example.C+
//
//
// Author List:
//      Stefan Pflueger          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>

#include <boost/progress.hpp>

#include "TStopwatch.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "reportingUtils.hpp"


using namespace std;
using namespace boost;
using namespace rpwa;


// generic function that creates mass bin directory structure and
// returns array of event files and trees, one for each mass bin
bool
createMassBinFiles(vector<TFile*>&    pwaFiles,
                   vector<TTree*>&    pwaTrees,
                   const string&      dirBaseName  = "/tmp/",
                   const unsigned int nmbMassBins  = 50,
                   const double       massBinWidth = 50,   // [MeV/c^2]
                   const double       massRangeMin = 500,  // [MeV/c^2]
                   const string&      treeName     = "rootPwaEvtTree")
{
	printInfo << "creating mass bin directories/files in '" << dirBaseName << "':" << endl;
	// cleanup
	for (unsigned int i = 0; i < pwaTrees.size(); ++i)
		if (pwaTrees[i])
			delete pwaTrees[i];
	pwaTrees.clear();
	for (unsigned int i = 0; i < pwaFiles.size(); ++i)
		if (pwaFiles[i]) {
			pwaFiles[i]->Close();
			delete pwaFiles[i];
		}
	pwaFiles.clear();
	pwaFiles.resize(nmbMassBins, 0);
	pwaTrees.resize(nmbMassBins, 0);
	bool success = true;
	for (unsigned int i = 0; i < nmbMassBins; ++i) {
		const double binMin = massRangeMin + i * massBinWidth;
		const double binMax = binMin + massBinWidth;
		// create mass bin directory
		stringstream n;
		n << binMin << "." << binMax;
		const string dirName = dirBaseName + "/" + n.str();
		mkdir(dirName.c_str(), S_IRWXU | S_IRWXG); // create directory read/writable by owner and group
		// create output file
		const string pwaFileName = dirName + "/" + n.str() + ".root";
		TFile*       pwaFile     = TFile::Open(pwaFileName.c_str(), "RECREATE");
		if (not pwaFile or pwaFile->IsZombie()) {
			printWarn << "problems creating file '" << pwaFileName << "'" << endl;
			success = false;
		} else {
			pwaFiles[i] = pwaFile;
			pwaTrees[i] = new TTree(treeName.c_str(), treeName.c_str());
			if (not pwaTrees[i]) {
				printWarn << "problems creating tree '" << treeName << "' " << "in file "
				          << "'" << pwaFileName << "'" << endl;
				success = false;
			} else {
				pwaTrees[i]->SetDirectory(pwaFile);
				printSucc << "created tree '" << treeName << "' in file " << "'" << pwaFileName
				          << "'" << endl;
			}
		}
	}
	return success;
}


// fills event data into correct mass bin
bool
writeEvent(vector<TTree*>&       pwaTrees,
           const TLorentzVector& beamLv,
           // !!! <channel-dependent part> !!!
           const TLorentzVector& piZero0,
           const TLorentzVector& piZero1,
           const TLorentzVector& piMinus,
           // !!! </channel-dependent part> !!!
           const double          XMass,                          // [GeV/c^2]
           const unsigned int    nmbMassBins             = 50,
           const double          massBinWidth            = 50,   // [MeV/c^2]
           const double          massRangeMin            = 500,  // [MeV/c^2]
           const string&         prodKinMomentaLeafName  = "prodKinMomenta",
           const string&         decayKinMomentaLeafName = "decayKinMomenta",
           const bool            debug                   = false)
{

	const double mass = 1000 * XMass; // convert from GeV/c^2 to MeV/c^2
	// make sure that mass is in range
	if ((mass < massRangeMin) or (mass > (massRangeMin + nmbMassBins * massBinWidth)))
		return false;
	const unsigned int bin = (unsigned int) ((mass - massRangeMin) / massBinWidth);
	if (not pwaTrees[bin]) {
		printWarn << "null pointer for tree for mass bin [" << massRangeMin + bin * massBinWidth << ", "
		          << massRangeMin + (bin + 1) * massBinWidth << "]" << endl;
		return false;
	}
	// fill tree
	if (debug)
		printDebug << "filling tree for bin " << bin << " = ["
		           << massRangeMin + bin * massBinWidth << ", "
		           << massRangeMin + (bin + 1) * massBinWidth << "] MeV/c^2" << endl;

	// create tree leafs
	static TClonesArray* prodKinMomenta  = new TClonesArray("TVector3");
	static TClonesArray* decayKinMomenta = new TClonesArray("TVector3");

	// connect leaf variables to tree branches or create branches, if they don't exist yet
	TTree* outTree = pwaTrees[bin];
	if (outTree->SetBranchAddress(prodKinMomentaLeafName.c_str(), &prodKinMomenta) < 0) {
		printWarn << "could not connect variable to branch '" << prodKinMomentaLeafName << "'. "
		          << "skipping." << endl;
		return false;
	}
	if (outTree->SetBranchAddress(decayKinMomentaLeafName.c_str(), &decayKinMomenta) < 0) {
		printWarn << "could not connect variable to branch '" << prodKinMomentaLeafName << "'. "
		          << "skipping." << endl;
		return false;
	}

	// clear arrays
	prodKinMomenta->Clear ();
	decayKinMomenta->Clear();

	// set leaf variables
	// beam particle
	new ((*prodKinMomenta)[0]) TVector3(beamLv.Vect());

	// !!! <channel-dependent part> !!!
	// for target particle elastic scattering is assumed
	// outgoing hadrons
	new ((*decayKinMomenta)[0]) TVector3(piZero0.Vect());
	new ((*decayKinMomenta)[1]) TVector3(piZero1.Vect());
	new ((*decayKinMomenta)[2]) TVector3(piMinus.Vect());
	// !!! </channel-dependent part> !!!

	// fill tree
	outTree->Fill();
	return true;
}


void
fillUdstDataIntoMassBins_example(const string&      inFileNamePattern = "fillUdstDataIntoMassBins_example.root",
                                 const string&      dirName           = "./test",
                                 const long int     maxNmbEvents      = -1,
                                 const unsigned int nmbMassBins       = 50,
                                 const double       massBinWidth      = 40,   // [MeV/c^2]
                                 const double       massRangeMin      = 500,  // [MeV/c^2]
                                 const string&      uDstTreeName      = "pwaDataTree",
                                 const string&      pwaTreeName       = "rootPwaEvtTree",
                                 const long int     treeCacheSize     = 25000000,  // 25 MByte ROOT tree read cache
                                 const bool         debug             = false)
{
	const string prodKinPartNamesObjName  = "prodKinParticles";
	const string prodKinMomentaLeafName   = "prodKinMomenta";
	const string decayKinPartNamesObjName = "decayKinParticles";
	const string decayKinMomentaLeafName  = "decayKinMomenta";

	TStopwatch timer;
	timer.Start();

	printInfo << "reading uDST file(s) '" << inFileNamePattern << "'" << endl
	          << "    writing " << nmbMassBins << " mass bins in mass interval "
	          << "[" << massRangeMin << ", " << massRangeMin + nmbMassBins * massBinWidth << "] "
	          << "MeV/c^2 to '" << dirName << "'" << endl
	          << "    reading uDST data from tree '" << uDstTreeName << "'" << endl
	          << "    writing PWA data to tree '" << pwaTreeName << "'" << endl;

	// create chain and connect tree leaf variables to branches
	TChain uDstChain(uDstTreeName.c_str());
	if (uDstChain.Add(inFileNamePattern.c_str()) < 1)
		printWarn << "no events in uDST input file(s) '" << inFileNamePattern << "'" << endl;
	const long int nmbEventsUdstChain = uDstChain.GetEntries();
	uDstChain.GetListOfFiles()->ls();

	// !!! <channel-dependent part> !!!
	// connect tree leafs
	TLorentzVector* photons[4] = {0, 0, 0, 0};
	TLorentzVector* piMinus    = 0;
	TLorentzVector* beam       = 0;
	TLorentzVector* recoil     = 0;
	uDstChain.SetBranchAddress("gamma1", &(photons[0]));
	uDstChain.SetBranchAddress("gamma2", &(photons[1]));
	uDstChain.SetBranchAddress("gamma3", &(photons[2]));
	uDstChain.SetBranchAddress("gamma4", &(photons[3]));
	uDstChain.SetBranchAddress("pi_out", &piMinus);
	uDstChain.SetBranchAddress("pi_in",  &beam);
	uDstChain.SetBranchAddress("proton", &recoil);
	uDstChain.SetCacheSize(treeCacheSize);
	uDstChain.AddBranchToCache("gamma1", true);
	uDstChain.AddBranchToCache("gamma2", true);
	uDstChain.AddBranchToCache("gamma3", true);
	uDstChain.AddBranchToCache("gamma4", true);
	uDstChain.AddBranchToCache("pi_out", true);
	uDstChain.AddBranchToCache("pi_in",  true);
	uDstChain.AddBranchToCache("proton", true);
	uDstChain.StopCacheLearningPhase();
	// !!! </channel-dependent part> !!!

	// create directories and .root files
	vector<TFile*> pwaDataFiles;
	vector<TTree*> pwaDataTrees;
	if (not createMassBinFiles(pwaDataFiles, pwaDataTrees, dirName, nmbMassBins, massBinWidth,
	                           massRangeMin, pwaTreeName)) {
		printErr << "there were problems creating the mass bin directories/files. Aborting..." << endl;
		return;
	}
	printSucc << "created " << pwaDataFiles.size() << " directories/files" << endl;

	// write arrays with production and decay particle names to root files
	{
		TClonesArray prodKinPartNames ("TObjString", 1);
		TClonesArray decayKinPartNames("TObjString", 3);

		// !!! <channel-dependent part> !!!
		new (prodKinPartNames [0]) TObjString("pi-"); // beam particle
		new (decayKinPartNames[0]) TObjString("pi0");
		new (decayKinPartNames[1]) TObjString("pi0");
		new (decayKinPartNames[2]) TObjString("pi-");
		// !!! </channel-dependent part> !!!

		for (unsigned int i = 0; i < pwaDataFiles.size(); ++i) {
			pwaDataFiles[i]->cd();
			prodKinPartNames.Write (prodKinPartNamesObjName.c_str (), TObject::kSingleKey);
			decayKinPartNames.Write(decayKinPartNamesObjName.c_str(), TObject::kSingleKey);
		}
		printSucc << "wrote particle name arrays to all files. "
		          << "beam = 'pi-', decay = {'pi0', 'pi0', 'pi-'}." << endl;
	}

	// create tree leafs
	{
		TClonesArray* prodKinMomenta  = new TClonesArray("TVector3");
		TClonesArray* decayKinMomenta = new TClonesArray("TVector3");
		const int     splitLevel      = 99;
		const int     bufSize         = 256000;
		for (unsigned int i = 0; i < pwaDataTrees.size(); ++i) {
			pwaDataTrees[i]->Branch(prodKinMomentaLeafName.c_str(),  "TClonesArray", &prodKinMomenta,
			                        bufSize, splitLevel);
			pwaDataTrees[i]->Branch(decayKinMomentaLeafName.c_str(), "TClonesArray", &decayKinMomenta,
			                        bufSize, splitLevel);
		}
		printSucc << "created branches for all trees" << endl;
	}

	// loop over events
	const long int nmbEvents = ((maxNmbEvents > 0) ? maxNmbEvents : nmbEventsUdstChain);
	printInfo << "writing events into mass bin files" << endl
	          << "    looping over " << nmbEvents << " tree entries" << endl;
	unsigned long int countEvWritten = 0;
	progress_display  progressIndicator(nmbEvents, cout, "");
	for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
		++progressIndicator;
		if ((uDstChain.LoadTree(eventIndex) < 0) or (uDstChain.GetEntry(eventIndex) == 0)) {
			printWarn << "error reading event " << eventIndex << " from tree. skipping." << endl;
			continue;
		}

		// !!! <channel-dependent part> !!!
		// now just make some minor calculations before
		// construct two pi0s
		const TLorentzVector piZeros[2] = {*(photons[0]) + *(photons[1]), *(photons[2]) + *(photons[3])};
		// construct intermediate state X
		const TLorentzVector X = piZeros[0] + piZeros[1] + *piMinus;
		// calculate t'
		const double t      = (*beam - X) * (*beam - X);
		const double tPrime = fabs(t) - fabs((X.M() * X.M() - beam->M() * beam->M())*(X.M() * X.M() - beam->M() * beam->M())
		                                     / (4 * (beam->Vect() * beam->Vect())));
		// cut on t'
		if ((tPrime < 0.1) or (tPrime > 1.0))
			continue;
		// write out PWA data
		const double         fsEnergy    = X.E() + recoil->E() - recoil->M();  // measured total energy of final state
		const double         scaleFactor = fsEnergy / beam->E();
		const TLorentzVector beamScaled(beam->Vect() * scaleFactor, fsEnergy);
		if (writeEvent(pwaDataTrees, beamScaled, piZeros[0], piZeros[1], *piMinus,
		               X.M(), nmbMassBins, massBinWidth, massRangeMin,
		               prodKinMomentaLeafName, decayKinMomentaLeafName, debug))
			++countEvWritten;
		// !!! </channel-dependent part> !!!
	}

	// write trees
	long unsigned int countTreeEvents = 0;
	for (unsigned int i = 0; i < nmbMassBins; ++i) {
		pwaDataTrees[i]->GetCurrentFile()->Write();
		long unsigned int nmbEvents = pwaDataTrees[i]->GetEntries();
		printSucc << "written " << setw(10) << nmbEvents << " events to file " << "'"
		          << pwaDataTrees[i]->GetCurrentFile()->GetName() << "'" << endl;
		countTreeEvents += nmbEvents;
		pwaDataTrees[i]->GetCurrentFile()->Close();
	}
	pwaDataFiles.clear();
	pwaDataTrees.clear();

	printInfo << "wrote " << min(countEvWritten, countTreeEvents)
	          << " out of " << nmbEvents << " events" <<endl;
	timer.Stop();
	printInfo << "this job consumed: ";
	timer.Print();
}
