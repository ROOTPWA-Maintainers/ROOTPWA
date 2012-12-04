///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
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
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////

/** @brief Simple partial wave event generator (homogeneous in m)
 */


#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cstdlib>
#include <string>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include <event.h>

#include <boost/progress.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TClonesArray.h"

#include "libconfig.h++"

#include "breitWignerProductionAmp.h"
#include "diffractivePhaseSpace.h"
#include "fileUtils.hpp"
#include "particleDataTable.h"
#include "partialWaveWeight.h"
#include "productionAmp.h"
#include "fileUtils.hpp"
#include "TFitBin.h"

#include "primaryVertexGen.h"

#include "randomNumberGenerator.h"
#include "generatorManager.h"

using namespace boost;
using namespace libconfig;
using namespace rpwa;
using namespace std;

extern ::particleDataTable PDGtable;


void printUsage(char* prog, int errCode = 0)
{
	cerr << "usage:" << endl
	     << prog
	     << " -n # [-a # -m # -M # -B # -s #] -o <file> -p <file> -w <file> -k <path> -i <file> -r <file>" << endl
	     << "    where:" << endl
	     << "        -n #       (max) number of events to generate (default: 100)" << endl
	     << "        -a #       (max) number of attempts to do (default: infinity)" << endl
	     << "        -m #       maxWeight FEATURE DISABLED" << endl
	     << "        -o <file>  ROOT output file (if not specified, generated automatically)" << endl
	     << "        -p         path to particle data table file (default: ./particleDataTable.txt)" << endl
	     << "        -w <file>  wavelist file (contains production amplitudes) FEATURE DISABLED" << endl
	     << "        -w <file.root>  to use TFitBin tree as input FEATURE PERMANENTLY REMOVED" << endl
	     << "        -c <0/1>   if 1 a comgeant eventfile (.fort.26) is written with same naming as the root file (default 0)" << endl
	     << "        -k <path>  path to keyfile directory (all keyfiles have to be there) FEATURE DISABLED" << endl
	     << "        -i <file>  integral file FEATURE DISABLED" << endl
	     << "        -r <file>  reaction config file" << endl
	     << "        -s #   set seed " << endl
	     << "        -M #   lower boundary of mass range in MeV (overwrites values from config file) " << endl
	     << "        -B #   width of mass bin in MeV" << endl
	     << endl
	     << "A comment regarding the disabled features: these options have been taken out "
	     << "for the time being. If you want to get them back, check SVN revision 1072." << endl
	     << endl;
	exit(errCode);
}


int main(int argc, char** argv)
{

	unsigned int nEvents = 100;
	unsigned int maxAttempts = 0;
	string outputFileName = ""; // either given by option or generated automatically by mass range
	string outputEvtFileName = "";
	string outputWhtFileName = "";
	string outputComgeantFileName = "";
//	string integralsFileName;
//	bool hasIntegralsFile = false;
	string pdgFileName = "./particleDataTable.txt";
//	string wavelistFileName; // format: name Re Im
//	string pathToKeyfiles = "./";
	string reactionFile;
//	double maxWeight = 0;
	int seed = 123456;
	bool seedSet = false;
	int massLower = 0;
	int massBinWidth = 0;
	bool overwriteMass = false;
	bool writeComgeantOut = false;

	int c;
	while ((c = getopt(argc, argv, "n:a:o:p:w:k:i:r:m:s:M:B:hc")) != -1) {
		switch (c) {
			case 'n':
				nEvents = atoi(optarg);
				break;
			case 'a':
				maxAttempts = atoi(optarg);
				break;
			case 's':
				seed = atoi(optarg);
				seedSet = true;
				break;
			case 'o':
				outputEvtFileName = optarg;
				break;
			case 'p':
				pdgFileName = optarg;
				break;
			case 'w':
				printErr << "this feature has been removed for the time being. "
				         << "If you want it back, check SVN revision 1072." << endl;
				exit(10);
//				wavelistFileName = optarg;
				break;
			case 'i':
				printErr << "this feature has been removed for the time being. "
				         << "If you want it back, check SVN revision 1072." << endl;
				exit(10);
//				integralsFileName = optarg;
//				hasIntegralsFile = true;
				break;
			case 'k':
				printErr << "this feature has been removed for the time being. "
				         << "If you want it back, check SVN revision 1072." << endl;
				exit(10);
//				pathToKeyfiles = optarg;
				break;
			case 'r':
				reactionFile = optarg;
				break;
			case 'm':
				printErr << "this feature has been removed for the time being. "
				         << "If you want it back, check SVN revision 1072." << endl;
				exit(10);
//				maxWeight = atof(optarg);
				break;
			case 'c':
				writeComgeantOut = true;
				break;
			case 'M':
				massLower = atoi(optarg);
				overwriteMass = true;
				break;
			case 'B':
				massBinWidth = atoi(optarg);
				overwriteMass = true;
				break;

			case 'h':
				printUsage(argv[0]);
				break;
			default:
				printUsage(argv[0], 5);
				break;
		}
	}

	if(not seedSet) {
		printInfo << "Setting random seed to " << seed << endl;
	}
	randomNumberGenerator::instance()->setSeed(seed);

	rpwa::particleDataTable::readFile(pdgFileName);
	generatorManager generatorMgr;
	if(not generatorMgr.readReactionFile(reactionFile)) {
		printErr << "could not read reaction file. Aborting..." << endl;
		exit(1);
	}
	if(not generatorMgr.initializeGenerator()) {
		printErr << "could not initialize generator. Aborting..." << endl;
		exit(1);
	}

	if(overwriteMass) {
		generatorMgr.overrideMassRange(massLower / 1000., (massLower + massBinWidth) / 1000.);
	}

	if(outputEvtFileName == "") {
		stringstream fileName;
		fileName << massLower << "." << massLower + massBinWidth << ".genbod.evt";
		outputEvtFileName = fileName.str();
	}
	ofstream outputEvtFile(outputEvtFileName.c_str());
	printInfo << "output event file: " << outputEvtFileName << endl;

	ofstream outputComgeantFile;
	if(writeComgeantOut) {
		outputComgeantFileName = changeFileExtension(outputEvtFileName, ".fort.26");
		printInfo << "output comgeant file: " << outputComgeantFileName << endl;
		outputComgeantFile.open(outputComgeantFileName.c_str());
	}

	progress_display* progressIndicator = new progress_display(nEvents, cout, "");

	unsigned int attempts = 0;
	for(unsigned int i = 0; i < nEvents; ++i) {

		attempts += generatorMgr.event();
		const generator& gen = generatorMgr.getGenerator();
		rpwa::particle beam = gen.getGeneratedBeam();
		std::vector<rpwa::particle> finalState = gen.getGeneratedFinalState();
		generator::convertEventToAscii(outputEvtFile, beam, finalState);
		if(writeComgeantOut) {
			rpwa::particle recoil = gen.getGeneratedRecoil();
			TVector3 vertex = gen.getGeneratedVertex();
			generator::convertEventToComgeant(outputComgeantFile, beam, recoil, vertex, finalState, false);
		}
		++(*progressIndicator);

	}

	outputEvtFile.close();
	if(writeComgeantOut) {
		outputComgeantFile.close();
	}

	printSucc << "generated " << nEvents << " events." << endl;
	printInfo << "attempts: " << attempts << endl;
	printInfo << "efficiency: " << 100 * ((double)nEvents / attempts) << "%" << endl;

	exit(255);

// -------------------------------------------------------------------------------

	PDGtable.initialize();

//	partialWaveWeight weighter;
	Config reactConf;
	reactConf.readFile(reactionFile.c_str());

	// variable that need to get initialized either by input options
	// or the config file
	// will be stored in the tree later
//	double weight, impweight;
	TClonesArray* momenta = new TClonesArray("TLorentzVector");
	TLorentzVector beam;
	TVector3 vertex;
//	double tPrime = 0.;
//	int qBeam;
	vector<int> charges; // array of charges
/*
	// generate the filename automatically if not specified
	if (outputFileName == "") {
		stringstream _filename;
		_filename << massLower << "." << massLower+massBinWidth << ".genbod.root";
		printInfo << "Guessed filename '" << _filename.str() << "'." << endl;
		outputFileName = _filename.str();
	}
	outputEvtFileName = changeFileExtension(outputFileName, ".evt");
	outputWhtFileName = changeFileExtension(outputFileName, ".wht");
	outputComgeantFileName = changeFileExtension(outputFileName, ".fort.26");
*/
	// now create the root file to store the events
//	TFile* outputFile = TFile::Open(outputFileName.c_str(), "RECREATE");
//	TH1D* weightsHistogram = new TH1D("hWeights", "PW Weights", 100, 0, 100);
//	TTree* outputTree = new TTree("pwevents", "pwevents");
/*
	outputTree->Branch("weight", &weight, "weight/d");
	outputTree->Branch("impweight", &impweight, "impweight/d");
	outputTree->Branch("p", &momenta);
	outputTree->Branch("beam", &beam);
	outputTree->Branch("vertex", &vertex);
	outputTree->Branch("q", &charges);
	outputTree->Branch("qbeam", &qBeam,"qbeam/I");
	outputTree->Branch("tprime", &tPrime);
*/
	//string theta_file= reactConf.lookup("finalstate.theta_file");
/*
	if(not seedSet) {
		printWarn << "Seed not set on command line, using default seed '"
		          << seed << "'." << endl;
	} else {
		printInfo << "Random seed: " << seed << endl;
	}
*/
	diffractivePhaseSpace diffPS;
/*
	double importanceMass;
	double importanceWidth;

	if(reactConf.lookupValue("importance.mass",importanceMass) &&
	   reactConf.lookupValue("importance.width",importanceWidth))
	{
		int act = reactConf.lookup("importance.active");
		if(act == 1) {
//			diffPS.setImportanceBW(importanceMass,importanceWidth);
		}
	}
*/
	const Setting& root = reactConf.getRoot();
	const Setting& fspart = root["finalstate"]["particles"];
	int nParticles = fspart.getLength();
/*
	// see if we have a resonance in this wave
	map<string, breitWignerProductionAmp*> bwAmps;
	if(reactConf.exists("resonances")) {
		const Setting &bws = root["resonances"]["breitwigners"];
		// loop through breitwigners
		int nbw = bws.getLength();
		printInfo << "Found " << nbw << " the following Breit-Wigners in config:" << endl;
		for(int ibw = 0; ibw < nbw; ++ibw) {
			const Setting &bw = bws[ibw];
			string jpcme;
			double mass, width;
			double cRe, cIm;
			bw.lookupValue("jpcme",       jpcme);
			bw.lookupValue("mass",        mass);
			bw.lookupValue("width",       width);
			bw.lookupValue("coupling_Re", cRe);
			bw.lookupValue("coupling_Im", cIm);
			complex<double> coupl(cRe, cIm);
			cout << "    JPCME = " << jpcme << ", mass = " << mass << " GeV/c^2, "
			     << "width = " << width << " GeV/c^2, coupling = " << coupl << endl;
			bwAmps[jpcme] = new breitWignerProductionAmp(mass, width, coupl);
		}
	}

	// check if TFitBin is used as input
	if(extensionFromPath(wavelistFileName) == "root") {
		printInfo << "Using TFitBin as input." << endl;
		TFile* fitResultsFile = TFile::Open(wavelistFileName.c_str(),"READ");
		TFitBin* fitBin = NULL;
		if(!fitResultsFile || fitResultsFile->IsZombie()) {
			printErr << "Cannot open start fit results file '" << wavelistFileName << "'." << endl;
			exit(1);
		}
		// get tree with start values
		TTree* tree;
		fitResultsFile->GetObject("pwa", tree);
		if(!tree) {
			printErr << "Cannot find TFitBin tree 'pwa' in file '" << wavelistFileName << "'." << endl;
			exit(1);
		} else {
			fitBin = new TFitBin();
			tree->SetBranchAddress("fitbin", &fitBin);
			// find entry which is closest to mass bin center
			unsigned int iBest = 0;
			double mBest = 0;
			for(unsigned int i = 0; i < tree->GetEntriesFast(); ++i) {
				tree->GetEntry(i);
				if(fabs(binCenter - fitBin->mass()) <= fabs(binCenter - mBest)) {
					iBest = i;
					mBest = fitBin->mass();
				}
			}  // end loop over TFitBins
			printInfo << "Using data from Mass bin with m = " << mBest << endl;
			tree->GetEntry(iBest);
			// write wavelist file for generator
			string tmpname = tmpnam(NULL);
			ofstream tmpfile(tmpname.c_str());
			fitBin->printAmpsGenPW(tmpfile);
			tmpfile.close();
			wavelistFileName=tmpname;
		}
	} // end root file given

	// read input wavelist and amplitudes
	printInfo << "Reading wavelist file '" << wavelistFileName << "':" << endl;
	ifstream wavefile(wavelistFileName.c_str());
	while(wavefile.good()) {
		TString wavename;
		double RE, IM;
		int rank = 0;
		int refl = 0;
		wavefile >> wavename >> RE >> IM;

		if(wavename.Contains("flat") || wavename.Length() < 2) {
			continue;
		}
		if(RE == 0 && IM == 0) {
			continue;
		}
		// check if there is rank information
		if(wavename(0) == 'V') {
			// we multiply rank by two to make space for refl+- production vectors
			rank = 2 * atoi(wavename(1,1).Data());
			// check reflecitivity to sort into correct production vector
			refl = wavename(9)=='+' ? 0 : 1;
			wavename = wavename(3, wavename.Length());
		}

		TString jpcme = wavename(2,5);

		wavename.ReplaceAll(".amp", ".key");
		wavename.Prepend(pathToKeyfiles.c_str());

		complex<double> amp(RE, IM);
		printInfo << "Got wave: " << wavename << " " << amp << ": r=" << rank/2
		          << " eps=" << refl << " qn=" << jpcme << endl;
		wavefile.ignore(256, '\n');

		productionAmp* prodAmp;
		if(bwAmps[jpcme.Data()] != NULL) {
			prodAmp = bwAmps[jpcme.Data()];
			printInfo << "Using BW for " << jpcme << endl;
			// production vector index: rank+refl
			weighter.addWave(wavename.Data(), prodAmp, amp, rank+refl);
		} else {
			prodAmp = new productionAmp(amp);
			weighter.addWave(wavename.Data(), prodAmp, complex<double>(1, 0), rank+refl);
		}
	}

	if(hasIntegralsFile){
		// read integral files
		ifstream integralsFile(integralsFileName.c_str());
		while(integralsFile.good()) {
			string filename;
			double mass;
			integralsFile >> filename >> mass;
			weighter.loadIntegrals(filename, mass);
			integralsFile.ignore(256, '\n');
		} // loop over integralfile
	}// endif hasIntegraFile

	double maxweight = -1;
*/
//	difPS.setVerbose(true);
//	ofstream outputEvtFile(outputEvtFileName.c_str());
//	ofstream outputComgeantFile;
//	if(writeComgeantOut) {
//		outputComgeantFile.open(outputComgeantFileName.c_str());
//	}
//	ofstream outputWhtFile(outputWhtFileName.c_str());
//	outputWhtFile << setprecision(10);

	{ // Block to avoid "i" and "attempts" in global namespace.

		unsigned int attempts = 0;
		unsigned int i = 0;
		progress_display* progressIndicator = new progress_display(nEvents, cout, "");

		for(i = 0; i < nEvents; ++i)
		{

			++attempts;

			if(maxAttempts > 0 && attempts > maxAttempts) {
				break;
			}

			momenta->Delete(); // clear output array

			stringstream str;
			if (writeComgeantOut) {
//				diffPS.event(str, outputComgeantFile);
			} else {
//				diffPS.event(str);
			}
//			impweight = diffPS.impWeight();

			for(int ip = 0; ip < nParticles; ++ip){
				new ((*momenta)[ip]) TLorentzVector(diffPS.getGeneratedFinalState()[ip].lzVec());
			}

			beam = diffPS.getGeneratedBeam().lzVec();
			vertex = diffPS.getGeneratedVertex();
//			tPrime = diffPS.tPrime();

			// calculate weight
//			event e;
//			e.setIOVersion(1);

//			str >> e;
//			outputEvtFile << e;
//			outputWhtFile << impweight << endl;

			//cerr << e << endl;
/*
			weight = weighter.weight(e);
			if(weight > maxweight) {
				maxweight = weight;
			}
//			weightsHistogram->Fill(weight);
*/
/*			if(maxWeight > 0) { // do weighting
				cout << weight << endl;
				//if(weight>maxWeight)maxWeight=weight;
				if(randomNumberGenerator::instance()->getGenerator()->Uniform() > weight / maxWeight) {
					continue;
				}
			}
*/
//			outputTree->Fill();
			++(*progressIndicator);

		} // end event loop

		printInfo << "Generated " << i << " Events." << endl;
		printInfo << "Attempts: " << attempts << endl;
		printInfo << "Efficiency: " << (double)i / (double)attempts << endl;
//		printInfo << "Maximum weight found: " << maxweight << endl;

	} // Block to avoid "i" and "attempts" in global namespace.
/*
	outputFile->cd();
	weightsHistogram->Write();
	outputTree->Write();
	outputFile->Close();
	outputEvtFile.close();
	outputWhtFile.close();
*/
	if(writeComgeantOut) {
		outputComgeantFile.close();
	}

	return 0;

}

