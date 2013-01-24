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

#ifdef USE_PWA2000
#include "keyfile.h"
#include "event.h"
extern particleDataTable PDGtable;
#endif

#include "mathUtils.hpp"
#include "fileUtils.hpp"
#include "particleDataTable.h"
#include "waveDescription.h"
#include "isobarHelicityAmplitude.h"
#include "evtTreeHelper.h"


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
#ifdef USE_PWA2000
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
#else
	return 0;
	printWarn << "code disabled, because compilation of PWA2000 is disabled" << endl;
#endif
}


struct symInfo {
	string waveNames[2];
	string symWaveName;
	double phi;
	double ratio;
};


vector<symInfo>
readDatFile(const string& fileName,
            const bool    debug = false)
{
	ifstream        datFile(fileName.c_str());
	vector<symInfo> symInfos;
	while (datFile.good()) {
		symInfo s;
		getline(datFile, s.waveNames[0]);
		getline(datFile, s.waveNames[1]);
		getline(datFile, s.symWaveName);
		string phi;
		getline(datFile, phi);
		if ((phi == "pi") or (phi == "3.1416") or (phi == "1"))
			s.phi = rpwa::pi;
		else
			s.phi = atof(phi.c_str());
		string ratio;
		getline(datFile, ratio);
		s.ratio = atof(ratio.c_str());
		symInfos.push_back(s);

		if (debug)
			printDebug << "read " << symInfos.size() << ". symInfo from file '" << fileName << "'" << endl
			           << "    name1   = '" << s.waveNames[0] << "'" << endl
			           << "    name1   = '" << s.waveNames[1] << "'" << endl
			           << "    symName = '" << s.symWaveName  << "'" << endl
			           << "    phi     = "  << maxPrecision(s.phi)   << endl
			           << "    ratio   = "  << maxPrecision(s.ratio) << endl;

	}
	return symInfos;
}


int
main(int argc, char** argv)
{
	printCompilerInfo();
	printGitHash();

	const long int maxNmbEvents = 1000;
	//const long int maxNmbEvents = 1;

	rpwa::particleDataTable& pdt = rpwa::particleDataTable::instance();
	pdt.readFile();
	TStopwatch timer;

	// waves with isospin symmetrization
	// rel. delta = (1.4396529281641782e-09, 1.4018652315552070e-09)
	// rms 4.88e-10, 1.56e-9
	// const string newKeyFileName = "test5pi/charly/sym/1-1+00+rho1700=a11260-=rho770_01_pi-_01_pi+_01_pi-.key";
	// const string newKeyFileName = "test5pi/charly/sym/1-1+00+rho1700=a11260+=rho770_01_pi+_01_pi-_01_pi-.key";
	// rel. delta = (1.3916463537100299e-09, 1.3553984601344382e-09)
	// rms 5.65e-10, 1.89e-9
	const string newKeyFileName = "test5pi/charly/sym/1-1+00+f11285=a11260-=rho770_01_pi-_11_pi+_11_pi-.key";
	// const string newKeyFileName = "test5pi/charly/sym/1-1+00+f11285=a11260+=rho770_01_pi+_11_pi-_11_pi-.key";
	// rel. delta = (1.4267325099126538e-09, 1.4274825765880905e-09)
	// rms 1.50e-9, 3.18e-9
	// const string newKeyFileName = "test5pi/charly/sym/1-2-00+f21270=a11260-=rho770_01_pi-_11_pi+_02_pi-.key";
	// rel. delta = (1.3479408700334261e-09, 1.3596883619015898e-09)
	// rms 2.48e-9, 6.67e-10
	// const string newKeyFileName = "test5pi/charly/sym/1-2-00+f21270=a11260-=rho770_01_pi-_11_pi+_22_pi-.key";
	const string pwa2kKeyFileName[2] = {
		// "test5pi/sebastian/sym/1-1++0+pi-_01_rho1700=a11269=pi+_0_rho770_0_pi-.key",
		// "test5pi/sebastian/sym/1-1++0+pi-_01_rho1700=a11269=pi-_0_rho770_0_pi+.key"
		// "test5pi/sebastian/sym/1-1++0+pi-_01_rho1700=a11269=pi+_0_rho770_0_pi-_noBose.key",
		// "test5pi/sebastian/sym/1-1++0+pi-_01_rho1700=a11269=pi-_0_rho770_0_pi+_noBose.key"
		"test5pi/sebastian/sym/1-1++0+pi-_11_f11285=pi-_11_a11269=pi+_0_rho770.key",
		"test5pi/sebastian/sym/1-1++0+pi-_11_f11285=pi+_11_a11269=pi-_0_rho770.key"
		// "test5pi/sebastian/sym/1-2-+0+pi-_02_f21270=pi-_1_a11269=pi+_0_rho770.key",
		// "test5pi/sebastian/sym/1-2-+0+pi-_02_f21270=pi+_1_a11269=pi-_0_rho770.key"
		// "test5pi/sebastian/sym/1-2-+0+pi-_22_f21270=pi-_1_a11269=pi+_0_rho770.key",
		// "test5pi/sebastian/sym/1-2-+0+pi-_22_f21270=pi+_1_a11269=pi-_0_rho770.key"
	};

	// read symmetrization info from .dat files
	vector<string>  datFiles = filesMatchingGlobPattern("test5pi/sebastian/sym/*.dat");
	vector<symInfo> symInfos;
	for (unsigned int i = 0; i < datFiles.size(); ++i) {
		vector<symInfo> s = readDatFile(datFiles[i]);
		symInfos.insert(symInfos.end(), s.begin(), s.end());
	}
	// get symmetrization parameters from .dat files
	double phi   = 0;
	double ratio = 0;
	for (unsigned int i = 0; i < symInfos.size(); ++i) {
		const string waveName = fileNameNoExtFromPath(pwa2kKeyFileName[0]) + ".amp";
		if (   (symInfos[i].waveNames[0] == waveName)
		    or (symInfos[i].waveNames[1] == waveName)) {
			phi   = symInfos[i].phi;
			ratio = symInfos[i].ratio;
			break;
		}
	}
	printInfo << "using phi = " << maxPrecision(phi) << " and ratio = " << maxPrecision(ratio)
	          << " for symmetrization" << endl;

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

#ifdef USE_PWA2000
	PDGtable.initialize("../../keyfiles/key5pi/pdgTable.txt");

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
	{
		double relDiff[2] = {(pwa2kAmps[0].real() - newAmps[0].real()) / pwa2kAmps[0].real(),
		                     (pwa2kAmps[0].imag() - newAmps[0].imag()) / pwa2kAmps[0].imag()};
		for (unsigned int i = 0; i < 2; ++i) {
			if (relDiff[i] < 1)
				relDiff[i] += 2;
			else if (relDiff[i] > 1)
				relDiff[i] -= 2;
		}
		printInfo << "newAmps[0] = " << maxPrecisionDouble(newAmps[0]) << " vs. pwa2kAmps[0] = "
		          << maxPrecisionDouble(pwa2kAmps[0]) << ", abs. delta = "
		          << maxPrecisionDouble(newAmps[0] - pwa2kAmps[0]) << ", rel. delta = "
		          << "(" << maxPrecision(relDiff[0]) << ", " << maxPrecision(relDiff[1]) << ")" << endl;
	}

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
			double relDiff[2] = {(pwa2kAmps[i].real() - newAmps[i].real()) / pwa2kAmps[i].real(),
			                     (pwa2kAmps[i].imag() - newAmps[i].imag()) / pwa2kAmps[i].imag()};
			for (unsigned int j = 0; j < 2; ++j) {
				if (relDiff[j] < 1)
					relDiff[j] += 2;
				else if (relDiff[j] > 1)
					relDiff[j] -= 2;
			}
			hMyAmpsReal->SetBinContent   (i + 1, newAmps[i].real());
			hMyAmpsImag->SetBinContent   (i + 1, newAmps[i].imag());
			hPwa2kAmpsReal->SetBinContent(i + 1, pwa2kAmps[i].real());
			hPwa2kAmpsImag->SetBinContent(i + 1, pwa2kAmps[i].imag());
			hDiffReal->Fill(relDiff[0]);
			hDiffImag->Fill(relDiff[1]);
			// hDiffReal->Fill(pwa2kAmps[i].real() - newAmps[i].real());
			// hDiffImag->Fill(pwa2kAmps[i].imag() - newAmps[i].imag());
			hCorrReal->Fill(newAmps[i].real(), pwa2kAmps[i].real());
			hCorrImag->Fill(newAmps[i].imag(), pwa2kAmps[i].imag());
		}
		f->Write();
		f->Close();
	}
#else
		printWarn << "code disabled, because compilation of PWA2000 is disabled" << endl;
#endif
}
