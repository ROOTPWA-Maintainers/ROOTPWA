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
//-----------------------------------------------------------
//
// Description:
//      fitting program for rootpwa
//      minimizes pwaLikelihood function
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#include <iostream>

#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>

#include "ampIntegralMatrix.h"
#include "fileUtils.hpp"
#include "fitResult.h"
#include "reportingUtilsEnvironment.h"

using namespace std;
using namespace rpwa;


ampIntegralMatrix readIntegralMatrix(const string& intFileName, const string& integralTKeyName = "integral") {
	printInfo << "loading normalization integral from '" << intFileName << "'" << endl;
	ampIntegralMatrix integral;
	const string normIntFileExt  = extensionFromPath(intFileName);
#ifdef USE_STD_COMPLEX_TREE_LEAFS
	if (normIntFileExt == "root") {
		TFile* intFile  = TFile::Open(intFileName.c_str(), "READ");
		if (not intFile or intFile->IsZombie()) {
			printErr << "could open normalization integral file '" << intFileName << "'. "
			         << "aborting." << endl;
			throw;
		}
		ampIntegralMatrix* integralPtr = 0;
		intFile->GetObject(integralTKeyName.c_str(), integralPtr);
		if (not integralPtr) {
			printErr << "cannot find integral object in TKey '" << integralTKeyName << "' in file "
			         << "'" << intFileName << "'. aborting." << endl;
			throw;
		}
		integral = *integralPtr;
		intFile->Close();
	} else
#endif  // USE_STD_COMPLEX_TREE_LEAFS
		if(normIntFileExt == "int") {
			integral.readAscii(intFileName);
	} else {
		printErr << "unknown file type '" << intFileName << "'. "
		         << "only .int and .root files are supported. aborting." << endl;
		throw;
	}
	return integral;
}



void
usage(const string& progName,
      const int     errCode = 0)
{

	cerr << "adds normalization integrals from fit result file" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " -i inuptFile -o outputFile -n PSIntFile [-a ACCIntFile -f]" << endl
	     << "    where:" << endl
	     << "        -i         input fit result file" << endl
	     << "        -o         output fit result file" << endl
	     << "        -n         phase space integral file to be added" << endl
	     << "        -a         accepted integral file to be added" << endl
	     << "        -f         force integral replacement (default: false)" << endl
	     << "        -h         print help" << endl
	     << endl;
	exit(errCode);
}


int
main(int    argc,
     char** argv)
{
	printCompilerInfo();
	printLibraryInfo ();
	printGitHash     ();
	cout << endl;

	// force loading predefined std::complex dictionary
	// see http://root.cern.ch/phpBB3/viewtopic.php?f=5&t=9618&p=50164
	gROOT->ProcessLine("#include <complex>");

	const string treeName = "pwa";
	const string branchName = "fitResult_v2";

	const string progName            = argv[0];
	string inputFileName = "";
	string outputFileName = "";
	string phaseSpaceIntegralFileName = "";
	string acceptedIntegralFileName = "";
	bool forceIntegralReplacement = false;
	extern char* optarg;
	// extern int optind;
	int c;
	while ((c = getopt(argc, argv, "i:o:n:a:fh")) != -1)
		switch (c) {
		case 'i':
			inputFileName = optarg;
			break;
		case 'o':
			outputFileName = optarg;
			break;
		case 'n':
			phaseSpaceIntegralFileName = optarg;
			break;
		case 'a':
			acceptedIntegralFileName = optarg;
			break;
		case 'f':
			forceIntegralReplacement = true;
			break;
		case 'h':
			usage(progName);
			break;
		}

	TFile* inputFile = TFile::Open(inputFileName.c_str(), "READ");
	if(not inputFile) {
		printErr << "could not open input file '" << inputFileName << "'. Aborting..." << endl;
		return 1;
	}

	const ampIntegralMatrix providedNormIntegral = readIntegralMatrix(phaseSpaceIntegralFileName);
	ampIntegralMatrix providedAccIntegral = readIntegralMatrix(acceptedIntegralFileName);

	TFile* outputFile = TFile::Open(outputFileName.c_str(), "NEW");
	if(not outputFile) {
		printErr << "could not open output file '" << outputFileName << "'. Aborting..." << endl;
		return 1;
	}

	TTree* inResultTree = 0;
	inputFile->GetObject(treeName.c_str(), inResultTree);
	if(not inResultTree) {
		printErr << "could not find input tree with name '" << treeName << "' in input file '" << inputFileName << "'. Aborting..." << endl;
		return 1;
	}

	fitResult* inResult = 0;
	inResultTree->SetBranchAddress(branchName.c_str(), &inResult);

	outputFile->cd();
	TTree* outResultTree = new TTree(treeName.c_str(), treeName.c_str());
	fitResult* outResult = 0;
	outResultTree->Branch(branchName.c_str(), &outResult);

	for(long i = 0; i < inResultTree->GetEntries(); ++i) {
		inResultTree->GetEntry(i);
		const unsigned int& nmbEvents = inResult->nmbEvents();
		const unsigned int& normNmbEvents = inResult->normNmbEvents();
		const double& massBinCenter = inResult->massBinCenter();
		const double& logLikelihood = inResult->logLikelihood();
		const int& rank = inResult->rank();
		const vector<TComplex> prodAmpsTComplex = inResult->prodAmps();
		const unsigned int nmbProdAmps = prodAmpsTComplex.size();
		vector<complex<double> > prodAmps(nmbProdAmps);
		for(unsigned int i = 0; i < nmbProdAmps; ++i) {
			prodAmps[i] = complex<double>(prodAmpsTComplex[i].Re(), prodAmpsTComplex[i].Im());
		}
		const vector<string>& prodAmpNames = inResult->prodAmpNames();
		const TMatrixT<double> fitParCovMatrix = inResult->fitParCovMatrix();
		const vector<pair<int, int> >& fitParCovMatrixIndices = inResult->fitParCovIndices();
		complexMatrix normIntegral = inResult->normIntegralMatrix();
		complexMatrix accIntegral = inResult->acceptedNormIntegralMatrix();
		vector<double> phaseSpaceIntegral = inResult->phaseSpaceIntegralVector();
		const bool& converged = inResult->converged();
		const bool& hasHesse = inResult->hasHessian();

		if(normIntegral.nCols() != 0 or normIntegral.nRows() != 0 or
		   accIntegral.nCols() != 0 or accIntegral.nRows() != 0 or
		   not phaseSpaceIntegral.empty())
		{
			printWarn << "input file already has integral matrices.";
			if(forceIntegralReplacement) {
				cout << "Replacing anyway..." << endl;
			} else {
				cout << "Aborting..." << endl;
				return 1;
			}
		}
		{
			const unsigned int nmbWaves = inResult->nmbWaves();
			normIntegral = complexMatrix(nmbWaves, nmbWaves);
			accIntegral = complexMatrix(nmbWaves, nmbWaves);
			phaseSpaceIntegral = vector<double>(nmbWaves, 0.);
			const double totalAcceptance = (double)providedAccIntegral.nmbEvents() / (double)providedNormIntegral.nmbEvents();
			providedAccIntegral.setNmbEvents(providedNormIntegral.nmbEvents());
			for(unsigned int waveIndex_i = 0; waveIndex_i < nmbWaves - 1; ++waveIndex_i) {
				const string& waveName_i(inResult->waveName(waveIndex_i).Data());
				if(waveName_i == "flat") {
					printErr << "encountered flat wave prematurely. Aborting..." << endl;
					return 1;
				}
				for(unsigned int waveIndex_j = 0; waveIndex_j < nmbWaves - 1; ++waveIndex_j) {
					const string& waveName_j(inResult->waveName(waveIndex_j).Data());
					if(waveName_j == "flat") {
						printErr << "encountered flat wave prematurely. Aborting..." << endl;
						return 1;
					}
					normIntegral.set(waveIndex_i, waveIndex_j, providedNormIntegral.element(waveName_i, waveName_j));
					accIntegral.set(waveIndex_i, waveIndex_j, providedAccIntegral.element(waveName_i, waveName_j));
				}
			}
			for(unsigned int i = 0; i < phaseSpaceIntegral.size() - 1; ++i) {
				phaseSpaceIntegral[i] = sqrt(normIntegral.get(i, i).real());
			}
			phaseSpaceIntegral[nmbWaves-1] = 1.;
			for(unsigned int waveIndex_i = 0; waveIndex_i < nmbWaves - 1; ++waveIndex_i) {
				for(unsigned int waveIndex_j = 0; waveIndex_j < nmbWaves - 1; ++waveIndex_j) {
					if(waveIndex_i == waveIndex_j) {
						normIntegral.set(waveIndex_i, waveIndex_j, 1.);
					} else {
						normIntegral.set(waveIndex_i, waveIndex_j, normIntegral.get(waveIndex_i, waveIndex_j) / (phaseSpaceIntegral[waveIndex_i] * phaseSpaceIntegral[waveIndex_j]));
					}
					accIntegral.set(waveIndex_i, waveIndex_j, accIntegral.get(waveIndex_i, waveIndex_j) / (phaseSpaceIntegral[waveIndex_i] * phaseSpaceIntegral[waveIndex_j]));
				}
			}
			for (unsigned int i = 0; i < normIntegral.nCols(); ++i) {
				normIntegral.set(nmbWaves-1, i, complex<double>(0., 0.));
				normIntegral.set(i, nmbWaves-1, complex<double>(0., 0.));
				accIntegral.set(nmbWaves-1, i, complex<double>(0., 0.));
				accIntegral.set(i, nmbWaves-1, complex<double>(0., 0.));
			}
			normIntegral.set(nmbWaves-1, nmbWaves-1, complex<double>(1., 0.));
			accIntegral.set (nmbWaves-1, nmbWaves-1, complex<double>(totalAcceptance, 0.));
		}

		outResult->reset();
		outResult->fill(nmbEvents,
		                normNmbEvents,
		                massBinCenter,
		                logLikelihood,
		                rank,
		                prodAmps,
		                prodAmpNames,
		                fitParCovMatrix,
		                fitParCovMatrixIndices,
		                normIntegral,  // contains the sqrt of the integral matrix diagonal elements!!!
		                accIntegral,
		                phaseSpaceIntegral,
		                converged,
		                hasHesse);
		outResultTree->Fill();

	}

	outputFile->cd();
	outResultTree->Write();
	outputFile->Close();

	if(outResult) {
		delete outResult;
	}

	return 0;
}
