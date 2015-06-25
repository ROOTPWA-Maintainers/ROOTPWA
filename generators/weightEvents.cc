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
//      calculates weights (= intensities) from amplitude and fit result
//      files and writes them to output file
//
//      limitations:
//
//      * the parsing of wave names to extract quantum numbers is not
//        implemented in an error safe way
//
//      * a better handling of wave quantum numbers would allow to
//        have a more flexible way of calculating weights for model
//        components; at the moment the weights for individual waves
//        and the weights for positive and negative reflectivity
//        totals are hard-coded
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#include <vector>

#include <boost/progress.hpp>

#include <TFile.h>
#include <TROOT.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TTree.h>

#include "ampIntegralMatrix.h"
#include "amplitudeTreeLeaf.h"
#include "fitResult.h"
#include "fileUtils.hpp"
#include "partialWaveFitHelper.h"
#include "reportingUtilsEnvironment.h"


using namespace std;
using namespace boost;
using namespace rpwa;


unsigned long
openRootAmpFiles(const string&               ampDirName,
                 const vector<string>&       waveNames,
                 vector<TTree*>&             ampRootTrees,
                 vector<amplitudeTreeLeaf*>& ampRootLeafs,
                 const string&               ampLeafName = "amplitude")
{
	ampRootTrees.clear();
	ampRootLeafs.assign(waveNames.size(), NULL);
	unsigned long nmbAmpValues = 0;
	for (size_t iWave = 0; iWave < waveNames.size(); ++iWave) {
		// no amplitude file for the flat wave
		if (waveNames[iWave] == "flat") {
			ampRootTrees.push_back(NULL);
			continue;
		}

		// try to open amplitude file
		const string ampFilePath = changeFileExtension(ampDirName + "/" + waveNames[iWave], ".root");
		printInfo << "opening amplitude file '" << ampFilePath << "'" << endl;
		TFile* ampFile = TFile::Open(ampFilePath.c_str(), "READ");
		if (not ampFile or ampFile->IsZombie()) {
			printErr << "cannot open amplitude file '" << ampFilePath << "'. skipping wave." << endl;
			continue;
		}

		// find amplitude tree
		TTree* ampTree = NULL;
		const string ampTreeName = changeFileExtension(waveNames[iWave], ".amp");
		ampFile->GetObject(ampTreeName.c_str(), ampTree);
		if (not ampTree) {
			printErr << "cannot find tree '" << ampTreeName << "' in file "
			         << "'" << ampFilePath << "'. skipping wave." << endl;
			continue;
		}

		// check that all tree have the same number of entries
		const unsigned long nmbEntries = ampTree->GetEntriesFast();
		if (nmbEntries == 0) {
			printErr << "amplitude tree '" << ampTree->GetName() << "' has zero entries. "
			         << "skipping wave." << endl;
			continue;
		}
		if (nmbAmpValues == 0)
			nmbAmpValues = nmbEntries;
		else if (nmbEntries != nmbAmpValues) {
			printErr << "amplitude tree '" << ampTree->GetName() << "' has different number "
			         << "of entries than previous tree. skipping wave." << endl;
			continue;
		}

		// connect tree leaf
		ampRootTrees.push_back(ampTree);
		ampTree->SetBranchAddress(ampLeafName.c_str(), &ampRootLeafs[iWave]);
	}

	return nmbAmpValues;
}


void
usage(const string& progName,
      const int     errCode = 0)
{
	cerr << "weights (phase space) events according to production amplitudes given by fit result" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " [-o output file -s -w fit-result file -n # of samples "
	     << "-i integral file -d amplitude directory -R] "
	     << "-m mass [-b mass bin width -t tree name -v -h]" << endl
	     << "    where:" << endl
	     << "        -o file    ROOT output file (default: './genpw.root')"<< endl
	     << "        -s         write out weights for each single wave (caution: this vastly increase the size of the output file)" << endl
	     << "        -w file    fit-result file containing the fitResult tree to be used as input (default: './fitresult.root')"<< endl
	     << "        -n #       if > 1, additional production amplitudes are generated according to covariances (default: 1)"<< endl
	     << "        -i file    integral file (default: './norm.root')"<< endl
	     << "        -d dir     path to directory with decay amplitude files (default: '.')" << endl
	     << "        -m #       central mass of mass bin [MeV/c^2]"<< endl
	     << "        -b #       width of mass bin [MeV/c^2] (default: 60 MeV/c^2)"<< endl
	     << "        -t name    name of tree in output file (default: rootPwaWeightTree)" << endl
	     << "        -v         verbose; print debug output (default: false)" << endl
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

#if ROOT_VERSION_CODE < ROOT_VERSION(6, 0, 0)
	// force loading predefined std::complex dictionary
	// see http://root.cern.ch/phpBB3/viewtopic.php?f=5&t=9618&p=50164
	gROOT->ProcessLine("#include <complex>");
#endif

	// parse command line options
	const string   progName                 = argv[0];
	string         outFileName              = "./genpw.root";
	string         fitResultFileName        = "./fitresult.root";
	const string   fitResultTreeName        = "pwa";
	const string   fitResultLeafName        = "fitResult_v2";
	string         intFileName              = "./norm.root";
	string         ampDirName               = ".";
	unsigned int   nmbProdAmpSamples        = 1;
	bool           writeSingleWaveWeights   = false;
	double         massBinCenter            = 0;                       // [MeV/c^2]
	double         massBinWidth             = 60;                      // [MeV/c^2]
	string         outTreeName              = "rootPwaWeightTree";
	bool           debug                    = false;

	int c;
	while ((c = getopt(argc, argv, "o:sw:n:i:d:m:b:t:vh")) != -1) {
		switch (c) {
		case 'o':
			outFileName = optarg;
			break;
		case 's':
			writeSingleWaveWeights = true;
			break;
		case 'w':
			fitResultFileName = optarg;
			break;
		case 'n':
			nmbProdAmpSamples = atoi(optarg);
			break;
		case 'i':
			intFileName = optarg;
			break;
		case 'd':
			ampDirName = optarg;
			break;
		case 'm':
			massBinCenter = atof(optarg);
			break;
		case 'b':
			massBinWidth = atof(optarg);
			break;
		case 't':
			outTreeName = optarg;
			break;
		case 'v':
			debug = true;
			break;

		case 'h':
		default:
			usage(progName);
			break;
		}
	}

	// the mass bin has to be specified
	if (massBinCenter == 0) {
		printErr << "central mass of mass bin to be processed has to be different from 0. Aborting..." << endl;
		usage(progName, 1);
	}

	// get input file names
	if (optind != argc) {
		printErr << "additional unhandled options encountered. Aborting..." << endl;
		usage(progName, 1);
	}

	// load integrals
	TFile* intFile  = TFile::Open(intFileName.c_str(), "READ");
	if (not intFile or intFile->IsZombie()) {
		printErr << "cannot open normalization integral file '" << intFileName << "'. "
		         << "Aborting..." << endl;
		exit(1);
	}
	ampIntegralMatrix* integral = 0;
	intFile->GetObject(ampIntegralMatrix::integralObjectName.c_str(), integral);
	if (not integral) {
		printErr << "cannot find integral object in TKey '" << ampIntegralMatrix::integralObjectName << "' in file "
		         << "'" << intFileName << "'. Aborting..." << endl;
		exit(1);
	}

	// load production amplitudes and wave information
	vector<string>                    waveNames;                    // wave names [wave index]
	vector<unsigned int>              waveIndex;                    // wave index [production amplitude index]
	                                                                // mapping from production amplitude to wave
	vector<string>                    prodAmpNames;                 // names of production amplitudes [production amplitude index]
	vector<vector<complex<double> > > prodAmps(nmbProdAmpSamples);  // production amplitudes [sample index][production amplitude index]
	vector<int>                       reflectivities;               // [production amplitude index]
	vector<int>                       ranks;                        // [production amplitude index]
	int                               maxRank = 0;
	{
		// open fit-result file
		TFile* fitResultFile = TFile::Open(fitResultFileName.c_str(), "READ");
		if (not fitResultFile or fitResultFile->IsZombie()) {
			printErr << "cannot open fit-result file '" << fitResultFileName << "'. Aborting..." << endl;
			exit(1);
		}
		// find fit-result tree
		TTree* fitResultTree;
		fitResultFile->GetObject(fitResultTreeName.c_str(), fitResultTree);
		if (not fitResultTree) {
			printErr << "cannot find fit-result tree '" << fitResultTreeName << "' in file '"
			         << fitResultFileName << "'. Aborting..." << endl;
			exit(1);
		}
		// loop over fit results and the ones which are closest to the given mass-bin center
		// if there are more than one such fit results, take the one which has the best likelihood
		printInfo << "searching for the best fit result with mass-bin center closest to m = "
		          << massBinCenter << " MeV/c^2." << endl;
		fitResult* resultBest = NULL;
		{
			unsigned int indexBest         = 0;
			double       massBinCenterBest = 0;
			double       logLikeBest       = 0;
			fitResult*   result            = NULL;
			fitResultTree->SetBranchAddress(fitResultLeafName.c_str(), &result);
			for (unsigned int iTree = 0; iTree < fitResultTree->GetEntries(); ++iTree) {
				fitResultTree->GetEntry(iTree);
				const double mass    = result->massBinCenter();
				const double logLike = result->logLikelihood();  // its actually -ln L
				if ((iTree == 0) or (fabs(massBinCenter - mass) < fabs(massBinCenter - massBinCenterBest))) {
					// this fit result is closer to the given mass-bin center
					indexBest         = iTree;
					massBinCenterBest = mass;
					logLikeBest       = logLike;
				} else if (fabs(massBinCenter - mass) == fabs(massBinCenter - massBinCenterBest)) {
					// this fit result is as close to the given mass-bin center as before
					// so pick the one with better likelihood
					if (logLike < logLikeBest) {  // smaller values are better
						indexBest         = iTree;
						massBinCenterBest = mass;
						logLikeBest       = logLike;
					}
				}
				if (debug)
					printDebug << "tree index [" << iTree << "] with m = " << mass
					           << " MeV/c^2 and log-likelihood " << logLike
					           << ((indexBest == iTree) ? " (best candidate so far)." : ".") << endl;
			}  // end loop over fit results

			// check that mass-bin center of found fit result lies in given bin width
			if (   (massBinCenterBest < (massBinCenter - massBinWidth / 2))
			       or (massBinCenterBest > (massBinCenter + massBinWidth / 2))) {
				printErr << "no fit found for mass region "  << massBinCenter << " +- "
				         << massBinWidth / 2. << " MeV/c^2" << endl;
				exit(1);
			}
			if (debug)
				printDebug << "found best matching fit result centered at m = " << massBinCenterBest
				           << " MeV/c^2 at tree index [" << indexBest << "]." << endl;
			fitResultTree->GetEntry(indexBest);
			resultBest = result;
		}

		// resize wave names to number of waves in fit result
		waveNames.resize(resultBest->nmbWaves(), "");

		// read production amplitudes, wave names, reflectivities, and rank
		// if nmbProdAmpSamples > 1 vary production amplitudes according to covariances
		for (unsigned int iSample = 0; iSample < nmbProdAmpSamples; ++iSample) {

			fitResult* result = NULL;
			if (iSample == 0) {
				result = resultBest;
			} else {
				result = resultBest->variedProdAmps();
			}

			for (unsigned int iProdAmp = 0; iProdAmp < result->nmbProdAmps(); ++iProdAmp) {
				const string prodAmpName = result->prodAmpName(iProdAmp).Data();
				const string waveName    = result->waveNameForProdAmp(iProdAmp).Data();
				const int iWave          = result->waveIndex(waveName);

				// extract rank and reflectivity from wave name
				const int rank = result->rankOfProdAmp(iProdAmp);
				const int refl = partialWaveFitHelper::getReflectivity(waveName);
				// read production amplitude
				const complex<double> prodAmp = result->prodAmp(iProdAmp);
				prodAmps[iSample].push_back(prodAmp);
				// for first fit result store wave info
				if (iSample == 0) {
					waveNames[iWave] = waveName;
					waveIndex.push_back     (iWave);
					prodAmpNames.push_back  (prodAmpName);
					reflectivities.push_back(refl);
					ranks.push_back         (rank);
					if (maxRank < rank)
						maxRank = rank;
				}

				if (debug)
					printDebug << "read production amplitude '" << prodAmpName << "'"
					           << " [" << iProdAmp << "] = " << prodAmp
					           << " for wave '" << waveName << "'; rank = " << rank
					           << ", reflectivity = " << refl << endl;
			}

			if (iSample > 0)
				delete result;

			++maxRank;
			printInfo << "rank of fit is " << maxRank << endl;
		}  // end loop over variations of production amplitudes
	}

	const unsigned int nmbWaves = waveNames.size();
	const unsigned int nmbProdAmps = prodAmpNames.size();

	// open decay amplitude files
	vector<TTree*>             ampRootTrees;
	vector<amplitudeTreeLeaf*> ampRootLeafs;
	const unsigned long        nmbEvents   = openRootAmpFiles(ampDirName, waveNames, ampRootTrees, ampRootLeafs);
	// test that an amplitude file was opened for each wave
	// note that ampRootLeafs cannot be used for this check
	if (waveNames.size() != ampRootTrees.size()) {
		printErr << "error opening ROOT amplitude files." << endl;
		exit(1);
	}
	if (nmbEvents == 0) {
		printErr << "could not read any amplitude from files." << endl;
		exit(1);
	}

	// create output file and tree
	TFile* outFile = TFile::Open(outFileName.c_str(), "RECREATE");
	TTree* outTree = new TTree(outTreeName.c_str(), outTreeName.c_str());
	// leaf variables
	double         weight;
	double         weightPosRef;
	double         weightNegRef;
	double         weightFlat;
	vector<double> weightWaves(nmbWaves);                    // branches will take pointer to elements
	vector<double> weightProdAmps(nmbProdAmps);              // branches will take pointer to elements
	vector<double> weightProdAmpSamples(nmbProdAmpSamples);  // branches will take pointer to elements
	// book branches
	outTree->Branch("weight",       &weight,       "weight/D");
	outTree->Branch("weightPosRef", &weightPosRef, "weightPosRef/D");
	outTree->Branch("weightNegRef", &weightNegRef, "weightNegRef/D");
	outTree->Branch("weightFlat",   &weightFlat,   "weightFlat/D");
	if (writeSingleWaveWeights) {
		// create weight branches for each individual wave
		for (unsigned int iWave = 0; iWave < nmbWaves; ++iWave) {
			TString weightName("weightWave_");
			weightName += waveNames[iWave];
			outTree->Branch(weightName.Data(), &weightWaves[iWave], (weightName + "/D").Data());
		}
		// if not a rank-1 fit, also create weights for each rank
		if (maxRank > 1) {
			for (unsigned int iProdAmp = 0; iProdAmp < nmbProdAmps; ++iProdAmp) {
				TString weightName("weightProdAmp_");
				weightName += prodAmpNames[iProdAmp];
				outTree->Branch(weightName.Data(), &weightProdAmps[iProdAmp], (weightName + "/D").Data());
			}
		}
	}
	// create branches for the weights calculated from the varied production amplitudes
	// weightProdAmpSamples[0] == weight
	for (unsigned int iSample = 0; iSample < nmbProdAmpSamples; ++iSample) {
		TString weightName("W");
		weightName += iSample;
		outTree->Branch(weightName.Data(), &weightProdAmpSamples[iSample], (weightName + "/D").Data());
	}

	// read data from tree(s) and calculate weight for each event
	TStopwatch timer;
	timer.Reset();
	timer.Start();

	// loop over amplitudes and calculate weight
	progress_display progressIndicator(nmbEvents, cout, "");
	for (unsigned long iEvent = 0; iEvent < nmbEvents; ++iEvent) {
		++progressIndicator;

		// read decay amplitudes for this event
		vector<complex<double> > decayAmps(nmbWaves);
			for (unsigned int iWave = 0; iWave < nmbWaves; ++iWave) {
			if (not ampRootTrees[iWave])  // e.g. flat wave
				decayAmps[iWave] = complex<double>(0);
			else {
				ampRootTrees[iWave]->GetEntry(iEvent);
				assert(ampRootLeafs[iWave]->nmbIncohSubAmps() == 1);
				decayAmps[iWave] = ampRootLeafs[iWave]->incohSubAmp(0);
			}
		}

		// calculate weight for each sample of production amplitudes
		for (unsigned int iSample = 0; iSample < nmbProdAmpSamples; ++iSample) {

			vector<complex<double> > posReflAmpSums(maxRank, 0);  // positive-reflectivity amplitude sums [rank]
			vector<complex<double> > negReflAmpSums(maxRank, 0);  // negative-reflectivity amplitude sums [rank]
			complex<double>          flatAmp = 0;

			for (unsigned int iProdAmp = 0; iProdAmp < nmbProdAmps; ++iProdAmp) {
				const unsigned int    iWave     = waveIndex        [iProdAmp];
				const complex<double> decayAmp  = decayAmps        [iWave];
				const complex<double> prodAmp   = prodAmps[iSample][iProdAmp];
				const string          waveName  = waveNames        [iWave];
				const double          nmbEvents = integral->nmbEvents();
				if (waveName == "flat") {
					flatAmp = prodAmp / sqrt(nmbEvents);
					if (iSample == 0)  // set weights of individual production amplitudes
						weightProdAmps[iProdAmp] = norm(flatAmp);
				} else {
					const complex<double> amp = decayAmp * prodAmp / sqrt(integral->element(waveName, waveName).real() * nmbEvents);
					if (iSample == 0)  // set weights of individual production amplitudes
						weightProdAmps[iProdAmp] = norm(amp);
					if (reflectivities[iProdAmp] == +1)
						posReflAmpSums[ranks[iProdAmp]] += amp;
					else if (reflectivities[iProdAmp] == -1)
						negReflAmpSums[ranks[iProdAmp]] += amp;
				}
			}  // end loop over waves

			// incoherent sum over rank
			double sampleWeightPosRef = 0;
			double sampleWeightNegRef = 0;
			for (int iRank = 0; iRank < maxRank; ++iRank) {
				sampleWeightPosRef += norm(posReflAmpSums[iRank]);
				sampleWeightNegRef += norm(negReflAmpSums[iRank]);
			}
			// total weight is incoherent sum of the two reflectivities and the flat wave
			double sampleWeightFlat = norm(flatAmp);
			double sampleWeight     = sampleWeightPosRef + sampleWeightNegRef + sampleWeightFlat;

			if (iSample == 0) {
				for (unsigned int iWave = 0; iWave < nmbWaves; ++iWave)
					weightWaves[iWave] = 0;
				// in the end the following corresponds to a sum over ranks,
				// which is an incoherent sum
				for (unsigned int iProdAmp = 0; iProdAmp < nmbProdAmps; ++iProdAmp)
					weightWaves[waveIndex[iProdAmp]] += weightProdAmps[iProdAmp];

				weightPosRef = sampleWeightPosRef;
				weightNegRef = sampleWeightNegRef;
				weightFlat   = sampleWeightFlat;
				weight       = sampleWeight;
			}

			weightProdAmpSamples[iSample] = sampleWeight;

		}  // end loop over production-amplitude samples
		outTree->Fill();
	}

	printSucc << "calculated weight for " << nmbEvents << " events" << endl;
	timer.Stop();
	printInfo << "this job consumed: ";
	timer.Print();

	outFile->cd();
	outTree->Write();
	outFile->Close();

	// cleanup
	intFile->Close();
	delete intFile;

	ampRootTrees.clear();
	ampRootLeafs.clear();

	prodAmps.clear();

	return 0;
}
