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
//      calculates weights (= intensities) from .amp files and writes
//      them to output file
//
//      limitations:
//
//      * at the moment the event data files are used to define the
//        looping over the events. in the case where the data are not
//        copied to the output file this is unnecessary. instead the
//        looping should be based on the amplitude files which are the
//        sole input needed to calculate the weight(s). this change
//        should be implemented along with the support for .root
//        amplitude files. also the use of friend trees should be
//        investigated in order to avoid copying of event data.
//
//      * the code reading the production amplitudes and calculating
//        the weigt is correct only for rank-1
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


#include <unistd.h>
#include <vector>

#include <boost/progress.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TStopwatch.h"

#include "ampIntegralMatrix.h"
#include "fitResult.h"
#include "fileUtils.hpp"
#include "evtTreeHelper.h"


using namespace std;
using namespace boost;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{
	cerr << "weights (phase space) events according to production amplitudes given by fit result" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " [-o output file -e -s -w fit-result file -n # of samples -i integral file -d amplitude directory] -m mass [-b mass bin width] "
	     << "-t tree name -l leaf names -v -h] "
	     << "input data file(s) (.evt or .root format)" << endl
	     << "    where:" << endl
	     << "        -o file    ROOT output file (default: './genpw.root')"<< endl
	     << "        -e         do _not_ copy event data to output file"<< endl
	     << "        -s         write out weights for each single wave (caution: this vastly increase the size of the output file)" << endl
	     << "        -w file    fit-result file containing the fitResult tree to be used as input (default: './fitresult.root')"<< endl
	     << "        -n #       if > 1, additional production amplitudes are generated according to covariances (default: 1)"<< endl
	     << "        -i file    integral file (default: './norm.int')"<< endl
	     << "        -d dir     path to directory with decay amplitude files (default: '.')" << endl
	     << "        -m #       central mass of mass bin [MeV/c^2]"<< endl
	     << "        -b #       width of mass bin [MeV/c^2] (default: 60 MeV/c^2)"<< endl
	     << "        -t name    name of tree in ROOT data files (default: rootPwaEvtTree)" << endl
	     << "        -l names   semicolon separated object/leaf names in input data" << endl
	     << "                   (default: 'prodKinParticles;prodKinMomenta;decayKinParticles;decayKinMomenta')" << endl
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

	// parse command line options
	const string   progName          = argv[0];
	string         outFileName       = "./genpw.root";
	bool           doCopyEventData   = true;
	string         fitResultFileName = "./fitresult.root";
	const string   fitResultTreeName = "pwa";
	const string   fitResultLeafName = "fitResult_v2";
	string         intFileName       = "./norm.int";
	string         ampDirName        = ".";
	unsigned int   nmbProdAmpSamples = 1;
	bool           writeSingleWaveWeights        = false;
	double         massBinCenter     = 0;   // [MeV/c^2]
	double         massBinWidth      = 60;  // [MeV/c^2]
	string         inTreeName        = "rootPwaEvtTree";
	string         leafNames         = "prodKinParticles;prodKinMomenta;decayKinParticles;decayKinMomenta";
	bool           debug             = false;
	const long int treeCacheSize     = 1000000;  // 1 MByte ROOT tree read cache size

	int c;
	while ((c = getopt(argc, argv, "o:esw:n:i:d:m:b:t:l:vh")) != -1) {
		switch (c) {
		case 'o':
			outFileName = optarg;
			break;
		case 'e':
			doCopyEventData = false;
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
			inTreeName = optarg;
			break;
		case 'l':
			leafNames = optarg;
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
		printErr << "central mass of mass bin to be processed has to be different from 0. aborting." << endl;
		usage(progName, 1);
	}

	// get input file names
	if (optind >= argc) {
		printErr << "you need to specify at least one data file to process. aborting." << endl;
		usage(progName, 1);
	}
	vector<string> rootInFileNames;
	vector<string> evtInFileNames;
	while (optind < argc) {
		const string fileName = argv[optind++];
		const string fileExt  = extensionFromPath(fileName);
		if (fileExt == "root")
			rootInFileNames.push_back(fileName);
		else if (fileExt == "evt")
			evtInFileNames.push_back(fileName);
		else
			printWarn << "input file '" << fileName << "' is neither a .root nor a .evt file. "
			          << "skipping." << endl;
	}
	if ((rootInFileNames.size() == 0) and (evtInFileNames.size() == 0)) {
		printErr << "none of the specified input files is a .root or .evt file. aborting.";
		usage(progName, 1);
	}

	// get object and leaf names for event data
	string prodKinPartNamesObjName,  prodKinMomentaLeafName;
	string decayKinPartNamesObjName, decayKinMomentaLeafName;
	parseLeafAndObjNames(leafNames, prodKinPartNamesObjName, prodKinMomentaLeafName,
	                     decayKinPartNamesObjName, decayKinMomentaLeafName);

	// open .root and .evt files for event data
	vector<TTree*> inTrees;
	TClonesArray*  prodKinPartNames  = 0;
	TClonesArray*  decayKinPartNames = 0;
	if (not openRootEvtFiles(inTrees, prodKinPartNames, decayKinPartNames,
	                         rootInFileNames, evtInFileNames,
	                         inTreeName, prodKinPartNamesObjName, prodKinMomentaLeafName,
	                         decayKinPartNamesObjName, decayKinMomentaLeafName, debug)) {
		printErr << "problems opening input file(s). aborting." << endl;
		exit(1);
	}

	// load integrals
	ampIntegralMatrix normInt;
	if (not normInt.readAscii(intFileName)) {
		printErr << "cannot read normalization integral from file '"
		         << intFileName << "'. aborting." << endl;
		exit(1);
	}

	// load production amplitudes and wave information
	vector<vector<complex<double> > > prodAmps(nmbProdAmpSamples);  // production amplitudes [sample index][wave index]
	vector<string>                    waveNames;                    // wave names [wave index]
	vector<int>                       reflectivities;               // [wave index]
	vector<int>                       ranks;                        // [wave index]
	int                               maxRank = 0;
	{
		// open fit-result file
		TFile* fitResultFile = TFile::Open(fitResultFileName.c_str(), "READ");
		if (not fitResultFile or fitResultFile->IsZombie()) {
			printErr << "cannot open fit-result file '" << fitResultFileName << "'. aborting." << endl;
			exit(1);
		}
		// find fit-result tree
		TTree* fitResultTree;
		fitResultFile->GetObject(fitResultTreeName.c_str(), fitResultTree);
		if (not fitResultTree) {
			printErr << "cannot find fit-result tree '" << fitResultTreeName << "' in file '"
			         << fitResultFileName << "'. aborting." << endl;
			exit(1);
		}
		// loop over fit results and the ones which are closest to the given mass-bin center
		// if there are more than one such fit results, take the one which has the best likelihood
		printInfo << "searching for the best fit result with mass-bin center closest to m = "
		          << massBinCenter << " MeV/c^2." << endl;
		fitResult* resultBest = 0;
		{
			unsigned int indexBest         = 0;
			double       massBinCenterBest = 0;
			double       logLikeBest       = 0;
			fitResult*   result            = 0;
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
			}  // end loop over fit results

			// check that mass-bin center of found fit result lies in given bin width
			if (   (massBinCenterBest < (massBinCenter - massBinWidth / 2))
			       or (massBinCenterBest > (massBinCenter + massBinWidth / 2))) {
				printErr << "no fit found for mass region "  << massBinCenter << " +- "
				         << massBinWidth / 2. << " MeV/c^2" << endl;
				exit(1);
			}
			printInfo << "found best matching fit result centered at m = " << massBinCenterBest
			          << " MeV/c^2 at tree index [" << indexBest << "]." << endl;
			fitResultTree->GetEntry(indexBest);
			resultBest = result;
		}

		// read production amplitudes, wave names, reflectivities, and rank
		// if nmbProdAmpSamples > 1 vary production amplitudes according to covariances
		for (unsigned int iSample = 0; iSample < nmbProdAmpSamples; ++iSample) {

			fitResult* result = 0;
			if (iSample == 0) {
				result = resultBest;
			} else {
				result = resultBest->variedProdAmps();
			}

			for (unsigned int iWave = 0; iWave < result->nmbWaves(); ++iWave) {
				TString  waveName = result->waveName(iWave);
				if (waveName.Length() < 2)
					continue;
				// extract rank and reflectivity from wave name
				int rank = 0;
				int refl = 0;
				if (waveName(0) == 'V') {
					rank     = atoi(waveName(1, 1).Data());
					// check reflectivity to sort into correct production vector
					refl     = ((waveName(9) == '+') ? +1 : -1);
					waveName = waveName(3, waveName.Length());
				} else if (waveName != "flat")
					refl = ((waveName(6) == '+') ? +1 : -1);
				// read production amplitude
				const complex<double> prodAmp = result->prodAmp(iWave);
				if (prodAmp == 0.)
					continue;
				prodAmps[iSample].push_back(prodAmp);
				// for first fit result store wave info
				if (iSample == 0) {
					waveNames.push_back     (waveName.Data());
					reflectivities.push_back(refl);
					ranks.push_back         (rank);
					if (maxRank < rank)
						maxRank = rank;
				}

				printDebug << "read production amplitude [" << iWave << "] = " << prodAmp
				           << " for wave '" << waveName << "'; rank = " << rank
				           << ", reflectivity = " << refl << endl;
			}

			if (iSample > 0)
				delete result;

			++maxRank;
			printDebug << "rank of fit is " << maxRank << endl;
		}  // end loop over variations of production amplitudes
	}

	const unsigned int nmbWaves = waveNames.size();

	// open decay amplitude files
	vector<ifstream*> ampFiles(nmbWaves, 0);
	for (unsigned int iWave = 0; iWave < nmbWaves; ++iWave) {
		if (waveNames[iWave] == "flat") {
			ampFiles.push_back(0);
			continue;
		}
		const string ampFilePath = ampDirName + "/" + waveNames[iWave];
		printInfo << "opening amplitude file '" << ampFilePath << "'" << endl;
		ifstream* ampFile = new ifstream(ampFilePath.c_str());
		if (not ampFile->good()) {
			printErr << "cannot open amplitude file '" << ampFilePath << "'. aborting." << endl;
			exit(1);
		}
		ampFiles[iWave] = ampFile;
	}

	// create output file and tree
	TFile* outFile = TFile::Open(outFileName.c_str(), "RECREATE");
	TTree* outTree = new TTree(inTreeName.c_str(), inTreeName.c_str());
	// leaf variables
	double         weight;
	double         weightPosRef;
	double         weightNegRef;
	double         weightFlat;
	vector<double> weights(nmbWaves);                        // branches will take pointer to elements
	vector<double> weightProdAmpSamples(nmbProdAmpSamples);  // branches will take pointer to elements
	// book branches
	outTree->Branch("weight",       &weight,       "weight/D");
	outTree->Branch("weightPosRef", &weightPosRef, "weightPosRef/D");
	outTree->Branch("weightNegRef", &weightNegRef, "weightNegRef/D");
	outTree->Branch("weightFlat",   &weightFlat,   "weightFlat/D");
	if (writeSingleWaveWeights) {
		// create weight branches for each individual wave
		for (unsigned int iWave = 0; iWave < nmbWaves; ++iWave) {
			TString weightName("Wintens_");
			weightName += waveNames[iWave];
			outTree->Branch(weightName.Data(), &weights[iWave], (weightName + "/D").Data());
		}
	}
	// create branches for the weights calculated from the varied production amplitudes
	// weightProdAmpSamples[0] == weight
	for (unsigned int iSample = 0; iSample < nmbProdAmpSamples; ++iSample) {
		TString weightName("W");
		weightName += iSample;
		outTree->Branch(weightName.Data(), &weightProdAmpSamples[iSample], (weightName + "/D").Data());
	}
	TBranch* outProdKinMomentaBr  = 0;
	TBranch* outDecayKinMomentaBr = 0;

	// read data from tree(s) and calculate weight for each event
	TStopwatch timer;
	timer.Reset();
	timer.Start();
	unsigned int eventCounter = 0;
	for (unsigned int iTree = 0; iTree < inTrees.size(); ++iTree) {
		printInfo << "processing ";
		if ((rootInFileNames.size() > 0) and (iTree == 0))
			cout << "chain of .root files";
		else
			cout << ".evt tree[" << ((rootInFileNames.size() > 0) ? iTree : iTree + 1) << "]";
		cout << endl;

		// create leaf variables and branch pointers
		TClonesArray* prodKinMomenta    = 0;
		TClonesArray* decayKinMomenta   = 0;
		TBranch*      prodKinMomentaBr  = 0;
		TBranch*      decayKinMomentaBr = 0;

		if (doCopyEventData) {
			// connect leaf variables to tree branches
			inTrees[iTree]->SetBranchAddress(prodKinMomentaLeafName.c_str(),  &prodKinMomenta,  &prodKinMomentaBr );
			inTrees[iTree]->SetBranchAddress(decayKinMomentaLeafName.c_str(), &decayKinMomenta, &decayKinMomentaBr);
			inTrees[iTree]->SetCacheSize(treeCacheSize);
			inTrees[iTree]->AddBranchToCache(prodKinMomentaLeafName.c_str(),  true);
			inTrees[iTree]->AddBranchToCache(decayKinMomentaLeafName.c_str(), true);
			inTrees[iTree]->StopCacheLearningPhase();
			// connect TClonesArrays of input tree to those of output tree
			const int splitLevel = 99;
			const int bufSize    = 256000;
			if (not outProdKinMomentaBr)
				outProdKinMomentaBr  = outTree->Branch(prodKinMomentaLeafName.c_str(),  "TClonesArray", &prodKinMomenta,  bufSize, splitLevel);
			else
				outProdKinMomentaBr->SetAddress(&prodKinMomenta);
			if (not outDecayKinMomentaBr)
				outDecayKinMomentaBr = outTree->Branch(decayKinMomentaLeafName.c_str(), "TClonesArray", &decayKinMomenta, bufSize, splitLevel);
			else
				outDecayKinMomentaBr->SetAddress(&decayKinMomenta);
		}

		// loop over data events and calculate weight
		const long int   nmbEventsTree = inTrees[iTree]->GetEntries();
		const long int   maxNmbEvents  = 0;
		const long int   nmbEvents     = ((maxNmbEvents > 0) ? min(maxNmbEvents, nmbEventsTree) : nmbEventsTree);
		progress_display progressIndicator(nmbEvents, cout, "");
		for (long int iEvent = 0; iEvent < nmbEvents; ++iEvent) {
			++eventCounter;
			++progressIndicator;

			if (inTrees[iTree]->LoadTree(iEvent) < 0)
				break;
			if (doCopyEventData) {
				// read only required branches
				prodKinMomentaBr->GetEntry (iEvent);
				decayKinMomentaBr->GetEntry(iEvent);
			}

			// read decay amplitudes for this event
			vector<complex<double> > decayAmps(nmbWaves);
			for (unsigned int iWave = 0; iWave < nmbWaves; ++iWave) {
				if (not ampFiles[iWave])  // e.g. flat wave
					decayAmps[iWave] = complex<double>(0);
				else {
					complex<double> decayAmp;
					ampFiles[iWave]->read((char*) &decayAmp, sizeof(complex<double>));
					decayAmps[iWave] = decayAmp;
				}
			}

			// calculate weight for each sample of production amplitudes
			for (unsigned int iSample = 0; iSample < nmbProdAmpSamples; ++iSample) {

				vector<complex<double> > posReflAmpSums(maxRank, 0);  // positive-reflectivity amplitude sums [rank]
				vector<complex<double> > negReflAmpSums(maxRank, 0);  // negative-reflectivity amplitude sums [rank]
				complex<double>          flatAmp = 0;

				for (unsigned int iWave = 0; iWave < nmbWaves; ++iWave) {
					const complex<double> decayAmp  = decayAmps        [iWave];
					const complex<double> prodAmp   = prodAmps[iSample][iWave];
					const string          waveName  = waveNames        [iWave];
					const double          nmbEvents = normInt.nmbEvents();
					if (waveName == "flat") {
						flatAmp = prodAmp / sqrt(nmbEvents);
						if (iSample == 0)  // set weights of individual waves
							weights[iWave] = norm(flatAmp);
					} else {
						const complex<double> amp = decayAmp * prodAmp / sqrt(normInt.element(waveName, waveName).real() * nmbEvents);
						if (iSample == 0)  // set weights of individual waves
							weights[iWave] = norm(amp);
						if (reflectivities[iWave] == +1)
							posReflAmpSums[ranks[iWave]] += amp;
						else if (reflectivities[iWave] == -1)
							negReflAmpSums[ranks[iWave]] += amp;
					}
				}  // end loop over waves

				// incoherent sum over rank and reflectivities
				weightPosRef = 0;
				weightNegRef = 0;
				for (int iRank = 0; iRank < maxRank; ++iRank) {
					weightPosRef += norm(posReflAmpSums[iRank]);
					weightNegRef += norm(negReflAmpSums[iRank]);
				}
				weightFlat = norm(flatAmp);
				// weight is incoherent sum of the two reflectivities and the flat wave
				weight = weightPosRef + weightNegRef + weightFlat;

				weightProdAmpSamples[iSample] = weight;

			}  // end loop over production-amplitude samples
			outTree->Fill();
		}

		inTrees[iTree]->PrintCacheStats();
	}
	printSucc << "calculated weight for " << eventCounter << " events" << endl;
	timer.Stop();
	printInfo << "this job consumed: ";
	timer.Print();

	outFile->cd();
	prodKinPartNames ->Write(prodKinPartNamesObjName.c_str(),  TObject::kSingleKey);
	decayKinPartNames->Write(decayKinPartNamesObjName.c_str(), TObject::kSingleKey);
	outTree->Write();
	outFile->Close();

	// cleanup
	for (unsigned int iWave = 0; iWave < nmbWaves; ++iWave) {
		if (ampFiles[iWave]) {
			ampFiles[iWave]->close();
			delete ampFiles[iWave];
		}
	}
	ampFiles.clear();
	prodAmps.clear();

	return 0;
}
