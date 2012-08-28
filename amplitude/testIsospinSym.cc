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
// $Rev:: 932                         $: revision of last commit
// $Author:: bgrube                   $: author of last commit
// $Date:: 2012-08-22 21:38:53 +0200 #$: date of last commit
//
// Description:
//      test program for amplitude isospin symmetrization
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
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TStopwatch.h"

#include "keyfile.h"
#include "event.h"

#include "particleDataTable.h"
#include "waveDescription.h"
#include "isobarHelicityAmplitude.h"
#include "evtTreeHelper.h"


extern particleDataTable PDGtable;


using namespace std;
using namespace boost;
using namespace rpwa;


long int
calcNewAmps(const string&             rootInFileName,
            const string&             keyFileName,
            vector<complex<double> >& amps,
            const long int            maxNmbEvents)
{
	waveDescription    waveDesc;
	isobarAmplitudePtr amp;
	if (not waveDesc.parseKeyFile(keyFileName) or not waveDesc.constructAmplitude(amp)) {
		printErr << "problems constructing amplitude. returning 0." << endl;
		return 0;
	}
	isobarDecayTopologyPtr topo = amp->decayTopology();
	printInfo << *amp;
	amp->init();
			
	// open input file
	vector<TTree*>       inTrees;
	TClonesArray*        prodKinPartNames  = 0;
	TClonesArray*        decayKinPartNames = 0;
	const vector<string> rootFileNames(1, rootInFileName);
	const vector<string> evtFileNames;
	const string&        inTreeName               = "rootPwaEvtTree";
	const string&        prodKinPartNamesObjName  = "prodKinParticles";
	const string&        prodKinMomentaLeafName   = "prodKinMomenta";
	const string&        decayKinPartNamesObjName = "decayKinParticles";
	const string&        decayKinMomentaLeafName  = "decayKinMomenta";
	if (not openRootEvtFiles(inTrees, prodKinPartNames, decayKinPartNames,
	                         rootFileNames, evtFileNames,
	                         inTreeName, prodKinPartNamesObjName, prodKinMomentaLeafName,
	                         decayKinPartNamesObjName, decayKinMomentaLeafName, true)) {
		printErr << "problems opening input files. returning 0." << endl;
		return 0;
	}
	const long int nmbEventsChain = inTrees[0]->GetEntries();

	// create branch pointers and leaf variables
	TBranch*      prodKinMomentaBr  = 0;
	TBranch*      decayKinMomentaBr = 0;
	TClonesArray* prodKinMomenta    = 0;
	TClonesArray* decayKinMomenta   = 0;
  
	// connect leaf variables to tree branches
	inTrees[0]->SetBranchAddress(prodKinMomentaLeafName.c_str(),  &prodKinMomenta,  &prodKinMomentaBr );
	inTrees[0]->SetBranchAddress(decayKinMomentaLeafName.c_str(), &decayKinMomenta, &decayKinMomentaBr);

	// loop over events
	if (not topo->initKinematicsData(*prodKinPartNames, *decayKinPartNames)) {
		printErr << "problems initializing input data. returning 0." << endl;
		return 0;
	}
	const long int   nmbEvents = ((maxNmbEvents > 0) ? min(maxNmbEvents, nmbEventsChain)
	                              : nmbEventsChain);
	progress_display progressIndicator(nmbEvents);
	for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
		++progressIndicator;
    
		if (inTrees[0]->LoadTree(eventIndex) < 0)
			break;
		// read only required branches
		prodKinMomentaBr->GetEntry (eventIndex);
		decayKinMomentaBr->GetEntry(eventIndex);
	      
		if (not prodKinMomenta or not decayKinMomenta) {
			printWarn << "at least one data array is null pointer: "
			          << "prodKinMomenta = "    << prodKinMomenta    << ", "
			          << "decayKinMomenta = "   << decayKinMomenta   << ". "
			          << "skipping event." << endl;
			continue;
		}
	      
		if (topo->readKinematicsData(*prodKinMomenta, *decayKinMomenta)) {
			amps.push_back((*amp)());
			if ((amps.back().real() == 0) or (amps.back().imag() == 0))
				printWarn << "event " << eventIndex << ": " << amps.back() << endl;
			// topo->productionVertex()->productionAmp();
		}
	}
	return nmbEvents;
}


long int
calcPwa2kAmps(const string&             evtInFileName,
              const string&             keyFileName,
              vector<complex<double> >& amps,
              const long int            maxNmbEvents)
{
	ifstream eventData(evtInFileName.c_str());
	keyfile  key;
	event    ev;
	key.open(keyFileName);
	ev.setIOVersion(1);
	long int       countEvent  = 0;
	const long int logInterval = 1000;
	while ((countEvent < maxNmbEvents) and (not (eventData >> ev).eof())) {
		complex<double> amp;
		key.run(ev, amp, true);
		amps.push_back(amp);
		key.rewind();
		++countEvent;
		if (countEvent % logInterval == 0)
			cout << "    " << countEvent << endl;
	}
	return countEvent;
}


int
main(int argc, char** argv)
{
	printCompilerInfo();
	printSvnVersion();
	
	//const long int maxNmbEvents = 10000;
	const long int maxNmbEvents = 1;

	rpwa::particleDataTable& pdt = rpwa::particleDataTable::instance();
	pdt.readFile();
	PDGtable.initialize("../keyfiles/key5pi/pdgTable.txt");
	TStopwatch timer;

	// waves with isospin symmetrization

	// rel. delta = (1.4319815006796074e-09, 1.4120124973618052e-09)
	// rms 1.51e-9, 1.91e-9
	// phi = 0, R = 1 ---> 1 / sqrt(2) * (a1 + a2)
	// const string newKeyFileName = "test5pi/charly/sym/1-1+00+rho1700=a11260-=rho770_01_pi-_01_pi+_01_pi-.key";
	const string newKeyFileName = "test5pi/charly/sym/1-1+00+rho1700=a11260+=rho770_01_pi+_01_pi-_01_pi-.key";
	const string pwa2kKeyFileName[2] = {
		"test5pi/sebastian/sym/1-1++0+pi-_01_rho1700=a11269=pi+_0_rho770_0_pi-.key",
		"test5pi/sebastian/sym/1-1++0+pi-_01_rho1700=a11269=pi-_0_rho770_0_pi+.key"
		// "test5pi/sebastian/sym/1-1++0+pi-_01_rho1700=a11269=pi+_0_rho770_0_pi-_noBose.key",
		// "test5pi/sebastian/sym/1-1++0+pi-_01_rho1700=a11269=pi-_0_rho770_0_pi+_noBose.key"
	};
	const double phi   = 0;
	const double ratio = 1;

	const string evtInFileName  = "test5pi/1900.1960.genbod.regen.evt";
	const string rootInFileName = "test5pi/1900.1960.genbod.root";
	// const string evtInFileName  = "test5pi/oneEvent.evt";
	// const string rootInFileName = "test5pi/oneEvent.root";

	//decayTopology::setDebug          (true);
	//isobarDecayTopology::setDebug    (true);
	//massDependence::setDebug         (true);
	//isobarAmplitude::setDebug        (true);
	//isobarHelicityAmplitude::setDebug(true);

	timer.Reset();
	timer.Start();
	vector<complex<double> > newAmps;
	calcNewAmps(rootInFileName, newKeyFileName, newAmps, maxNmbEvents);
	timer.Stop();
	printSucc << "read " << newAmps.size() << " events from file(s) "
	          << "'" << rootInFileName << "' and calculated amplitudes" << endl;
	cout << "needed ";
	printInfo << "newAmps[0] = " << maxPrecisionDouble(newAmps[0]) << endl;
	timer.Print();

	timer.Reset();
	timer.Start();
	vector<complex<double> > pwa2kAmps;
	{
		vector<complex<double> > amps[2];
		calcPwa2kAmps(evtInFileName, pwa2kKeyFileName[0], amps[0], maxNmbEvents);
		calcPwa2kAmps(evtInFileName, pwa2kKeyFileName[1], amps[1], maxNmbEvents);
		const unsigned int nmbAmps = amps[0].size();
		assert(amps[1].size() == nmbAmps);
		pwa2kAmps.resize(nmbAmps);
		complex<double> phase(ratio * cos(phi), ratio * sin(phi));
		printInfo << "PWA2000 term[0] = " << maxPrecisionDouble(amps[0][0])
		          << ", term[1] = " << maxPrecisionDouble(amps[1][0])
		          << ", phase = " << phase << endl;
		for (unsigned int i = 0; i < nmbAmps; ++i)
			pwa2kAmps[i] = (1 / sqrt(1 + ratio)) * (amps[0][i] + phase * amps[1][i]);
	}
	timer.Stop();
	printSucc << "read " << pwa2kAmps.size() << " events from file(s) "
	          << "'" << evtInFileName << "' and calculated amplitudes" << endl;
	cout << "needed ";
	timer.Print();
	printInfo << "newAmps[0] = " << maxPrecisionDouble(newAmps[0]) << " vs. pwa2kAmps[0] = "
	          << maxPrecisionDouble(pwa2kAmps[0]) << ", abs. delta = "
	          << maxPrecisionDouble(newAmps[0] - pwa2kAmps[0]) << ", rel. delta = "
	          << "(" << maxPrecision((newAmps[0].real() - pwa2kAmps[0].real()) / newAmps[0].real())
	          << ", " << maxPrecision((newAmps[0].imag() - pwa2kAmps[0].imag()) / newAmps[0].imag())
	          << ")" << endl;
	
	if (1) {
		const string outFileName = "testAmplitudeDiff.root";
		printInfo << "writing comparison plots to " << outFileName << endl;
		TFile* f              = TFile::Open(outFileName.c_str(), "RECREATE");
		TH1D*  hMyAmpsReal    = new TH1D("hMyAmpsReal",    "hMyAmpsReal;Event Number;#Rgothic[Amplitude]",    newAmps.size(),    -0.5, newAmps.size()    - 0.5);
		TH1D*  hMyAmpsImag    = new TH1D("hMyAmpsImag",    "hMyAmpsImag;Event Number;#Jgothic[Amplitude]",    newAmps.size(),    -0.5, newAmps.size()    - 0.5);
		TH1D*  hPwa2kAmpsReal = new TH1D("hPwa2kAmpsReal", "hPwa2kAmpsReal;Event Number;#Rgothic[Amplitude]", pwa2kAmps.size(), -0.5, pwa2kAmps.size() - 0.5);
		TH1D*  hPwa2kAmpsImag = new TH1D("hPwa2kAmpsImag", "hPwa2kAmpsImag;Event Number;#Jgothic[Amplitude]", pwa2kAmps.size(), -0.5, pwa2kAmps.size() - 0.5);
		TH1D*  hDiffReal      = new TH1D("hDiffReal", "hDiffReal;#Rgothic[Amplitude] Difference;Count", 100000, -1e-7, 1e-7);
		TH1D*  hDiffImag      = new TH1D("hDiffImag", "hDiffImag;#Jgothic[Amplitude] Difference;Count", 100000, -1e-7, 1e-7);
		TH2D*  hCorrReal      = new TH2D("hCorrReal", "hCorrReal;#Rgothic[My Amp];#Rgothic[PWA2000 Amp]", 1000, -2, 2, 1000, -2, 2);
		TH2D*  hCorrImag      = new TH2D("hCorrImag", "hCorrImag;#Jgothic[My Amp];#Jgothic[PWA2000 Amp]", 1000, -2, 2, 1000, -2, 2);
		for (unsigned int i = 0; i < newAmps.size(); ++i) {
			hMyAmpsReal->SetBinContent   (i + 1, newAmps[i].real());
			hMyAmpsImag->SetBinContent   (i + 1, newAmps[i].imag());
			hPwa2kAmpsReal->SetBinContent(i + 1, pwa2kAmps[i].real());
			hPwa2kAmpsImag->SetBinContent(i + 1, pwa2kAmps[i].imag());
			hDiffReal->Fill((pwa2kAmps[i].real() - newAmps[i].real()) / pwa2kAmps[i].real());
			hDiffImag->Fill((pwa2kAmps[i].imag() - newAmps[i].imag()) / pwa2kAmps[i].imag());
			// hDiffReal->Fill(pwa2kAmps[i].real() - newAmps[i].real());
			// hDiffImag->Fill(pwa2kAmps[i].imag() - newAmps[i].imag());
			hCorrReal->Fill(newAmps[i].real(), pwa2kAmps[i].real());
			hCorrImag->Fill(newAmps[i].imag(), pwa2kAmps[i].imag());
		}
		f->Write();
		f->Close();
	}
}
