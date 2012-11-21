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
	     << "        -m #       maxWeight" << endl
	     << "        -o <file>  ROOT output file (if not specified, generated automatically)" << endl
	     << "        -p         path to particle data table file (default: ./particleDataTable.txt)" << endl
	     << "        -w <file>  wavelist file (contains production amplitudes)" << endl
	     << "        -w <file.root>  to use TFitBin tree as input" << endl
	     << "        -c <0/1>   if 1 a comgeant eventfile (.fort.26) is written with same naming as the root file (default 0)" << endl
	     << "        -k <path>  path to keyfile directory (all keyfiles have to be there)" << endl
	     << "        -i <file>  integral file" << endl
	     << "        -r <file>  reaction config file" << endl
	     << "        -s #   set seed " << endl
	     << "        -M #   lower boundary of mass range in MeV (overwrites values from config file) " << endl
	     << "        -B #   width of mass bin in MeV" << endl
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
	string integralsFileName;
	bool hasIntegralsFile = false;
	string pdgFileName = "./particleDataTable.txt";
	string wavelistFileName; // format: name Re Im
	string pathToKeyfiles = "./";
	string reactionFile;
	double maxWeight = 0;
	int seed = 123456;
	bool seedSet = false;
	int massLower = 0;
	int massBinWidth = 0;
	bool overwriteMass = false;
	bool writeComgeantOut = false;

	int c;
	while ((c = getopt(argc, argv, "n:a:o:p:w:k:i:r:m:s:M:B:h:c")) != -1) {
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
				outputFileName = optarg;
				break;
			case 'p':
				pdgFileName = optarg;
				break;
			case 'w':
				wavelistFileName = optarg;
				break;
			case 'i':
				integralsFileName = optarg;
				hasIntegralsFile = true;
				break;
			case 'k':
				pathToKeyfiles = optarg;
				break;
			case 'r':
				reactionFile = optarg;
				break;
			case 'm':
				maxWeight = atof(optarg);
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
				printUsage(argv[0]);
				break;
		}
	}

	PDGtable.initialize();

	partialWaveWeight weighter;
	Config reactConf;
//	reactConf.readFile(reactionFile.c_str());

	rpwa::particleDataTable::readFile(pdgFileName);
	generatorManager generatorMgr;
	generatorMgr.readReactionFile(reactionFile);




	// variable that need to get initialized either by input options
	// or the config file
	// will be stored in the tree later
	double weight, impweight;
	TClonesArray* momenta = new TClonesArray("TLorentzVector");
	TLorentzVector beam;
	TVector3 vertex;
	double tPrime = 0.;
	int qBeam;
	vector<int> charges; // array of charges

	double beamMom = reactConf.lookup("beam.momentum");
	double beamMomSigma = reactConf.lookup("beam.sigma_momentum");
	double beamPartMass = 0.13957018;
	if(reactConf.exists("beam.mass")){
		beamPartMass = reactConf.lookup("beam.mass");
	} else {
		printWarn << "Beam particle mass not found in config file. Using '"
				  << setprecision(8) << beamPartMass << "'." << endl;
	}
	double beamDxDz = reactConf.lookup("beam.DxDz");
	double beamDxDzSigma = reactConf.lookup("beam.sigma_DxDz");
	double beamDyDz = reactConf.lookup("beam.DyDz");
	double beamDyDzSigma = reactConf.lookup("beam.sigma_DyDz");

	double targetZPosition = reactConf.lookup("target.pos.z");
	double targetLength = reactConf.lookup("target.length");
	double targetRadius = reactConf.lookup("target.radius");
	double massRecoil = reactConf.lookup("target.mrecoil");

	double tPrimeMin = 0.;
	double tPrimeMax = numeric_limits<double>::max();
	reactConf.lookupValue("finalstate.t_min", tPrimeMin);
	reactConf.lookupValue("finalstate.t_max", tPrimeMax);

	double minimumFinalStateMass = reactConf.lookup("finalstate.mass_min");
	double maximumFinalStateMass = reactConf.lookup("finalstate.mass_max");
	if(overwriteMass) {
		minimumFinalStateMass = massLower / 1000.0;
		maximumFinalStateMass = (massLower + massBinWidth) / 1000.0;
	}
	// array of tslopes even when only one is existing
	double* tSlopes = NULL;
	double* finalStateInvariantMasses  = NULL;
	int numberOftSlopes = 1;
	if(reactConf.lookup("finalstate.t_slope").isArray()) {
		numberOftSlopes = reactConf.lookup("finalstate.t_slope").getLength();
		if(reactConf.lookup("finalstate.inv_m").getLength() != numberOftSlopes) {
			printErr << " Error: please check number of t' values and the corresponding invariant masses in the Configuration File! " << endl;
			exit(1);
		}
		tSlopes = new double[numberOftSlopes];
		finalStateInvariantMasses  = new double[numberOftSlopes];
		printInfo << "Found array of t' slopes. Reading " << numberOftSlopes << "values...";
		for(int i = 0; i < numberOftSlopes; ++i) {
			tSlopes[i] = reactConf.lookup("finalstate.t_slope")[i];
			finalStateInvariantMasses[i] = reactConf.lookup("finalstate.inv_m")[i];
		}
		cout << "done." << endl;
	} else {
		tSlopes = new double[1];
		tSlopes[0] = reactConf.lookup("finalstate.t_slope");
		printInfo << "Found one t' slope: " << tSlopes[0] << "." << endl;
	}
	double binCenter = 500 * (minimumFinalStateMass + maximumFinalStateMass);

  // check whether to use a primary vertex generator as requested by the config file
	primaryVertexGen* primaryVtxGen = NULL;
	string histFilenamePrimVertex = "";
	if(reactConf.lookupValue("primvertex.histfilename", histFilenamePrimVertex)) {
		primaryVtxGen = new primaryVertexGen(
		    histFilenamePrimVertex,
		    beamPartMass,
		    beamMom,
		    beamMomSigma
		    );
		if(!primaryVtxGen->check()) {
			printErr << "Could not load histogram file '" << histFilenamePrimVertex << "with beam properties" << endl;
			delete primaryVtxGen;
			primaryVtxGen = NULL;
		}
	}

	if(!reactConf.lookupValue("beam.charge", qBeam)) {
		printWarn << "Beam charge not found in config file. Setting it to '-1'." << endl;
		qBeam = -1;
	}

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

	// now create the root file to store the events
	TFile* outputFile = TFile::Open(outputFileName.c_str(), "RECREATE");
	TH1D* weightsHistogram = new TH1D("hWeights", "PW Weights", 100, 0, 100);
	TTree* outputTree = new TTree("pwevents", "pwevents");

	outputTree->Branch("weight", &weight, "weight/d");
	outputTree->Branch("impweight", &impweight, "impweight/d");
	outputTree->Branch("p", &momenta);
	outputTree->Branch("beam", &beam);
	outputTree->Branch("vertex", &vertex);
	outputTree->Branch("q", &charges);
	outputTree->Branch("qbeam", &qBeam,"qbeam/I");
	outputTree->Branch("tprime", &tPrime);

	//string theta_file= reactConf.lookup("finalstate.theta_file");

	if(not seedSet) {
		printWarn << "Seed not set on command line, using default seed '"
		          << seed << "'." << endl;
	} else {
		printInfo << "Random seed: " << seed << endl;
	}

	diffractivePhaseSpace diffPS;
	diffPS.setSeed(seed);
	diffPS.setBeam(beamMom, beamMomSigma, beamDxDz, beamDxDzSigma, beamDyDz, beamDyDzSigma);
	diffPS.setTarget(targetZPosition, targetLength, targetRadius, massRecoil);
	diffPS.setTPrimeSlope(tSlopes, finalStateInvariantMasses, numberOftSlopes);
	diffPS.setMassRange(minimumFinalStateMass, maximumFinalStateMass);
	diffPS.setPrimaryVertexGen(primaryVtxGen);
	if(tPrimeMin >= 0.) {
		diffPS.setTPrimeMin(tPrimeMin);
	} else {
		printErr << "t_min (minimum t') must be positive, found '" << tPrimeMin << "'." << endl;
		exit(1);
	}
	if(tPrimeMax >= 0.) {
		diffPS.setTPrimeMax(tPrimeMax);
	} else {
		printErr << "t_max (maximum t') must be positive, found '" << tPrimeMax << "'." << endl;
		exit(1);
	}

	double importanceMass;
	double importanceWidth;

	if(reactConf.lookupValue("importance.mass",importanceMass) &&
	   reactConf.lookupValue("importance.width",importanceWidth))
	{
		int act = reactConf.lookup("importance.active");
		if(act == 1) {
			diffPS.setImportanceBW(importanceMass,importanceWidth);
		}
	}

	const Setting& root = reactConf.getRoot();
	const Setting& fspart = root["finalstate"]["particles"];
	int nParticles = fspart.getLength();

	for(int ifs = 0; ifs < nParticles; ++ifs) {
		const Setting &part = fspart[ifs];
		int id;
		part.lookupValue("g3id", id);
		int myq;
		part.lookupValue("charge", myq);
		double m;
		part.lookupValue("mass", m);
		charges.push_back(myq);
		diffPS.addDecayProduct(particleInfo(id, myq, m));
	}

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

	//difPS.setVerbose(true);
	ofstream outputEvtFile(outputEvtFileName.c_str());
	ofstream outputComgeantFile;
	if(writeComgeantOut) {
		outputComgeantFile.open(outputComgeantFileName.c_str());
	}
	ofstream outputWhtFile(outputWhtFileName.c_str());
	outputWhtFile << setprecision(10);

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
				diffPS.event(str, outputComgeantFile);
			} else {
				diffPS.event(str);
			}
			impweight = diffPS.impWeight();

			for(int ip = 0; ip < nParticles; ++ip){
				new ((*momenta)[ip]) TLorentzVector(*diffPS.decay(ip));
			}

			beam = *diffPS.beam();
			vertex = *diffPS.vertex();
			tPrime = diffPS.tPrime();

			// calculate weight
			event e;
			e.setIOVersion(1);

			str >> e;
			outputEvtFile << e;
			outputWhtFile << impweight << endl;

			//cerr << e << endl;

			weight = weighter.weight(e);
			if(weight > maxweight) {
				maxweight = weight;
			}
			weightsHistogram->Fill(weight);

			if(maxWeight > 0) { // do weighting
				cout << weight << endl;
				//if(weight>maxWeight)maxWeight=weight;
				if(gRandom->Uniform() > weight / maxWeight) {
					continue;
				}
			}

			outputTree->Fill();
			++(*progressIndicator);

		} // end event loop

		printInfo << "Generated " << i << " Events." << endl;
		printInfo << "Attempts: " << attempts << endl;
		printInfo << "Efficiency: " << (double)i / (double)attempts << endl;
		printInfo << "Maximum weight found: " << maxweight << endl;

	} // Block to avoid "i" and "attempts" in global namespace.

	outputFile->cd();
	weightsHistogram->Write();
	outputTree->Write();
	outputFile->Close();
	outputEvtFile.close();
	outputWhtFile.close();
	if(writeComgeantOut) {
		outputComgeantFile.close();
	}

	delete [] tSlopes;
	delete [] finalStateInvariantMasses;

	return 0;

}

