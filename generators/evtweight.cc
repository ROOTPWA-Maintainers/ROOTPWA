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
//      Calculate weight of event list from amp files
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#include <cstdlib>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <unistd.h>

#include <boost/progress.hpp>

#include "TClonesArray.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TTree.h"

#include "integral.h"
#include "fitResult.h"
#include "event.h"
#include "fileUtils.hpp"
#include "../amplitude/particle.h"
#include "evtTreeHelper.h"
#include "particleDataTable.h"


using namespace std;
using namespace boost;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{
	cerr << "weight (phase space) events with fit result" << endl
         << endl
         << "usage:" << endl
	     << progName
	     << " [-o output file -w fit result file -i integral file -d amplitude directory -n samples] -m mass [-b mass bin width] "
	     << "-t tree name -l leaf names -v -h] "
         << "input data file(s) (.evt or .root format)" << endl
	     << "    where:" << endl
	     << "        -o file    ROOT output file (default: './genpw.root')"<< endl
	     << "        -w file    file containing the TFitBin tree to be used as input (default: './fitresult.root')"<< endl
	     << "        -i file    integral file (default: './norm.int')"<< endl
	     << "        -d dir     path to directory with decay amplitude files (default: '.')" << endl
	     << "        -n #       number of samples of model parameters (default: 1)"<< endl
	     << "        -m #       central mass of mass bin [MeV/c^2]"<< endl
	     << "        -b #       width of mass bin [MeV/c^2] (default: 60 MeV/c^2)"<< endl
	     << "        -p file    path to particle data table file (default: ./particleDataTable.txt)" << endl
	     << "        -t name    name of tree in ROOT data files (default: rootPwaEvtTree)" << endl
	     << "        -l names   semicolon separated object/leaf names in input data" << endl
             << "                   (default: 'prodKinParticles;prodKinMomenta;decayKinParticles;decayKinMomenta')" << endl
	     << "        -v         verbose; print debug output (default: false)" << endl
	     << "        -h         print help" << endl
	     << endl;
	exit(errCode);
}


int
getReflectivity(const TString& waveName)
{
	int refl = 0;
	unsigned int reflIndex = 6;  // position of reflectivity in wave
	// check whether it is parameter or wave name
	if(waveName[0] == 'V') {
		reflIndex = 9;
	}
	if(waveName[reflIndex] == '-') {
		refl = -1;
	} else if(waveName[reflIndex] == '+') {
		refl = +1;
	} else {
		printErr << "Cannot parse parameter/wave name '" << waveName << "'. Cannot not determine reflectivity. Aborting." << endl;
		throw;
	}
	return refl;
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
	const string   progName       = argv[0];
	string         outFileName    = "./genpw.root";
	string         fitFileName    = "./fitresult.root";
	string         intFileName    = "./norm.int";
	string         ampDirName     = ".";
	unsigned int   nmbSamples     = 1;
	double         massBinCenter  = 0.;
	double         massBinWidth   = 60.;
	string         pdgFileName    = "./particleDataTable.txt";
	string         inTreeName     = "rootPwaEvtTree";
	string         leafNames      = "prodKinParticles;prodKinMomenta;decayKinParticles;decayKinMomenta";
	bool           debug          = false;
	const long int treeCacheSize  = 1000000;  // 1 MByte ROOT tree read cache size

	int c;
	while ((c = getopt(argc, argv, "o:w:i:d:m:n:b:p:t:l:vh")) != -1) {
		switch (c) {
		case 'o':
			outFileName = optarg;
			break;
		case 'w':
			fitFileName = optarg;
			break;
		case 'i':
			intFileName = optarg;
			break;
		case 'd':
			ampDirName = optarg;
			break;
		case 'n':
			nmbSamples = atoi(optarg);
			break;
		case 'm':
			massBinCenter = atof(optarg);
			break;
		case 'b':
			massBinWidth = atof(optarg);
			break;
		case 'p':
			pdgFileName = optarg;
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
	if (massBinCenter == 0.) {
		printErr << "central mass of mass bin to be processed has to be specified. aborting." << endl;
		usage(progName, 1);
	}

        // get input file names
        if (optind >= argc) {
                printErr << "you need to specify at least one data file to process. aborting." << endl;
                usage(progName, 1);
        }
        vector<string> rootFileNames;
        vector<string> evtFileNames;
        while (optind < argc) {
                const string fileName = argv[optind++];
                const string fileExt  = extensionFromPath(fileName);
                if (fileExt == "root")
                        rootFileNames.push_back(fileName);
                else if (fileExt == "evt")
                        evtFileNames.push_back(fileName);
                else
                        printWarn << "input file '" << fileName << "' is neither a .root nor a .evt file. "
                                  << "skipping." << endl;
        }
        if ((rootFileNames.size() == 0) and (evtFileNames.size() == 0)) {
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
                                 rootFileNames, evtFileNames,
                                 inTreeName, prodKinPartNamesObjName, prodKinMomentaLeafName,
                                 decayKinPartNamesObjName, decayKinMomentaLeafName, debug)) {
                printErr << "problems opening input file(s). aborting." << endl;
                exit(1);
        }

	// initialize particle data table
	rpwa::particleDataTable::readFile(pdgFileName);
  
	TFile* outfile = TFile::Open(outFileName.c_str(), "RECREATE");
	TH1D* hWeights = new TH1D("hWeights", "PW Weights", 100, 0, 100);
	TTree* outtree = new TTree("pwevents", "pwevents");
	double weight;
	TClonesArray* p = new TClonesArray("TLorentzVector");
	TLorentzVector beam;
	double qbeam;
	std::vector<int> q; // array of charges

	outtree->Branch("weight", &weight, "weight/D");
	outtree->Branch("p", &p);
	outtree->Branch("beam", &beam);
	outtree->Branch("q", &q);
	outtree->Branch("qbeam", &qbeam, "qbeam/i");

	TBranch* outProdKinMomentaBr  = 0;
	TBranch* outDecayKinMomentaBr = 0;

	// load integrals ---------------------------------------------------
	integral normInt;
	ifstream intFile(intFileName.c_str());
	if(!intFile) {
		printErr << "Cannot open file '"
		         << intFileName << "'. Exiting." << endl;
		throw;
	}
	// !!! integral.scan() performs no error checks!
	normInt.scan(intFile);
	intFile.close();

	// load production amplitudes ------------------------------------------
	// read TFitResult is used as input
	TFile* fitresults = TFile::Open(fitFileName.c_str(), "READ");
	fitResult* Bin = NULL;
	if(!fitresults || fitresults->IsZombie()) {
		cerr << "Cannot open start fit results file " << fitFileName << endl;
		return 1;
	}
	// get tree with start values
	vector<vector<complex<double> >*> prodAmps(nmbSamples); //production amplitudes
	vector<string> waveNames;
	vector<int> reflectivities;
	vector<int> ms;
	vector<int> ranks;
	int maxrank = 0;

	bool hasfit = true;
	TTree* tree;
	fitresults->GetObject("pwa", tree);
	if(!tree) {
		cerr << "Cannot find fitbin tree '"<< "pwa" << "' "<< endl;
		exit(1);
	} else {
		Bin = new fitResult();
		tree->SetBranchAddress("fitResult_v2", &Bin);
		// find entry which is closest to mass bin center
		// and has the lowest likelihood
		printInfo << "Searching for best fit in mass bin centered at " << massBinCenter << " MeV/c^2." << std::endl;
		unsigned int iBest = 0;
		double mBest = 0;
		double loglike = std::numeric_limits<double>::max();
		for(unsigned int i = 0; i < tree->GetEntries(); ++i) {
			tree->GetEntry(i);
			if (fabs(massBinCenter - Bin->massBinCenter()) < fabs(massBinCenter - mBest)) {
				// this mass bin is closer to mass bin searched for
				iBest = i;
				mBest = Bin->massBinCenter();
				loglike = Bin->logLikelihood();
			} else if (fabs(massBinCenter - Bin->massBinCenter()) <= fabs(massBinCenter - mBest)) {
				// this mass bin is close as before, also compare likelihood
				if (Bin->logLikelihood() < loglike) {
					iBest = i;
					mBest = Bin->massBinCenter();
					loglike = Bin->logLikelihood();
				}
			}
		}  // end loop over TFitBins


		if((mBest < (massBinCenter-massBinWidth / 2.)) || (mBest > (massBinCenter + massBinWidth / 2.))) {
			cerr << "No fit found for Mass bin m=" << massBinCenter << " (mBest=" << mBest << ", massBinCenter=" << massBinCenter << ", massBinWidth=" << massBinWidth << ")" << endl;
			Bin->reset();
			hasfit = false;
			return 1;
		} else {
			printInfo << "Found best matching mass bin to be centered at " << mBest << " MeV/c^2, index " << iBest << "." << endl;
			tree->GetEntry(iBest);
		}
		// write wavelist files for generator
		// string tmpname("genamps.txt");
		// create slightly modified versions of the model according to covariances
		for(unsigned int isamples = 0; isamples < nmbSamples; ++isamples) {
			prodAmps[isamples] = new vector<complex<double> >();

			fitResult* fit = NULL;
			if(isamples == 0) {
				fit = Bin;
			} else {
				fit = Bin->cloneVar();
			}

			for (unsigned int iwaves = 0; iwaves < fit->nmbWaves(); iwaves++) {
				TString wavename = fit->waveName(iwaves);
				const double RE = fit->prodAmp(iwaves).real();
				const double IM = fit->prodAmp(iwaves).imag();

				int rank = 0;
				int refl = 0;
				int m = 0;

				if(wavename.Length() < 2) {
					continue;
				}
				if((RE == 0) && (IM == 0)) {
					continue;
				}
				// check if there is rank information
				if(wavename(0) == 'V') {
					// we multiply rank by to to make space for refl+- production vectors
					rank = 2 * atoi(wavename(1, 1).Data());
					// check reflecitivity to sort into correct production vector
					refl = wavename(9)=='+' ? 0 : 1;
					m = wavename(8)=='0' ? 0 : 1;
					//cerr << wavename(9) << endl;
					wavename = wavename(3, wavename.Length());
				} else if(wavename != "flat") {
					refl = wavename(6)=='+' ? 0 : 1;
					m = wavename(5)=='0' ? 0 : 1;
					//cerr << wavename(6) << endl;
				}

				std::complex<double> amp(RE, IM);
				prodAmps[isamples]->push_back(amp);
				cerr << wavename << " " << amp << " r=" << rank/2
				     << " eps=" << refl
				     << " m="   << m << endl;
				// for first file store info on waves
				if(isamples == 0) {
					waveNames.push_back(wavename.Data());
					reflectivities.push_back(refl);
					ms.push_back(m);
					ranks.push_back(rank);
					if(maxrank < rank) {
						maxrank = rank;
					}
				}
			}

			if (isamples != 0) {
				delete fit;
			}

			cerr << "Rank of fit was:" << maxrank + 1 << endl;
		} // end loop over model versions
	}

	unsigned int nmbWaves = waveNames.size();
	// TODO reserve list of wheight branches!

	vector<ifstream*> ampfiles;
	// reserve vector beforehand because Branch
	// will take a pointer onto the elements
	//vector<double> weights((nmbWaves+1)*nmbWaves/2);
	vector<double> weights(nmbWaves);
	//unsigned int wcount=0;
	//create wheight vectors for individual intensities and interference terms
	// for(unsigned int iw=0;iw<nmbWaves;++iw){
	//   for(unsigned int jw=iw;jw<nmbWaves;++jw){
	//     TString weightname("W_");
	//     if(iw==jw)weightname+=waveNames[iw];
	//     else weightname+=waveNames[iw] +"_"+ waveNames[jw];

	//     outtree->Branch(weightname.Data(),&weights[wcount++],(weightname+"/d").Data());
	//   }
	// }
	for(unsigned int iw = 0; iw < nmbWaves; ++iw) {
		TString weightname("Wintens_");
		weightname += waveNames[iw];
		outtree->Branch(weightname.Data(), &weights[iw], (weightname + "/D").Data());
	}

	// create branches for the weights of the different model variants
	vector<double> modelweights(nmbSamples);
	for(unsigned int isamples = 0; isamples < nmbSamples; ++isamples) {
		TString weightname("W");
		weightname += isamples;
		outtree->Branch(weightname.Data(), &modelweights[isamples], (weightname + "/D").Data());
	}

	// open decay amplitude files --------------------------------------------
	for(unsigned int iw = 0; iw < nmbWaves; ++iw) {
		if (waveNames[iw] == "flat") {
			ampfiles.push_back(NULL);
			continue;
		}
		string ampFilePath = ampDirName + "/" + waveNames[iw];
		printInfo << "loading amplitude data from '" << ampFilePath << "'" << endl;
		ifstream* ampfile = new ifstream(ampFilePath.c_str());
		if (!ampfile->good()) {
			printErr << "cannot open amplitude file '" << ampFilePath << "'. aborting." << endl;
			exit(1);
		}
		ampfiles.push_back(ampfile);
	}


	// read data from tree(s) and calculate weight
	TStopwatch timer;
	timer.Reset();
	timer.Start();
	unsigned int counter = 0;
	for (unsigned int i = 0; i < inTrees.size(); ++i) {
		printInfo << "processing ";
		if ((rootFileNames.size() > 0) and (i == 0)) 
			cout << "chain of .root files";
		else
			cout << ".evt tree[" << ((rootFileNames.size() > 0) ? i : i + 1) << "]";
		cout << endl;

		// create branch pointers and leaf variables
		TBranch*      prodKinMomentaBr  = 0;
		TBranch*      decayKinMomentaBr = 0;
		TClonesArray* prodKinMomenta    = 0;
		TClonesArray* decayKinMomenta   = 0;

		// connect leaf variables to tree branches
		inTrees[i]->SetBranchAddress(prodKinMomentaLeafName.c_str(),  &prodKinMomenta,  &prodKinMomentaBr );
		inTrees[i]->SetBranchAddress(decayKinMomentaLeafName.c_str(), &decayKinMomenta, &decayKinMomentaBr);
		inTrees[i]->SetCacheSize(treeCacheSize);
		inTrees[i]->AddBranchToCache(prodKinMomentaLeafName.c_str(),  true);
		inTrees[i]->AddBranchToCache(decayKinMomentaLeafName.c_str(), true);
		inTrees[i]->StopCacheLearningPhase();

		// connect TClonesArrays from input to output tree
		if (outProdKinMomentaBr == NULL) {
			const int splitLevel = 99;
			const int bufSize    = 256000;
			outProdKinMomentaBr  = outtree->Branch(prodKinMomentaLeafName.c_str(),  "TClonesArray", &prodKinMomenta,  bufSize, splitLevel);
		} else {
			outProdKinMomentaBr->SetAddress(&prodKinMomenta);
		}
		if (outDecayKinMomentaBr == NULL) {
			const int splitLevel = 99;
			const int bufSize    = 256000;
			outDecayKinMomentaBr = outtree->Branch(decayKinMomentaLeafName.c_str(), "TClonesArray", &decayKinMomenta, bufSize, splitLevel);
		} else {
			outDecayKinMomentaBr->SetAddress(&decayKinMomenta);
		}

		const long int   nmbEventsTree = inTrees[i]->GetEntries();
		const long int   maxNmbEvents  = 0;
		const long int   nmbEvents     = ((maxNmbEvents > 0) ? min(maxNmbEvents, nmbEventsTree)
		                                  : nmbEventsTree);
		progress_display progressIndicator(nmbEvents, cout, "");
		for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
			++counter;
			++progressIndicator;

			if (inTrees[i]->LoadTree(eventIndex) < 0)
				break;
			// read only required branches
			prodKinMomentaBr->GetEntry (eventIndex);
			decayKinMomentaBr->GetEntry(eventIndex);

			// read decay amps for this event
			vector<complex<double> > decayamps(nmbWaves);
			for(unsigned int iw = 0; iw < nmbWaves; ++iw) {
				if (ampfiles[iw] == NULL) { // e.g. flat wave
					decayamps[iw] = complex<double>(0., 0.);
				} else {
					complex<double> decayamp;
					ampfiles[iw]->read((char*) &decayamp, sizeof(complex<double> ));
					decayamps[iw] = decayamp;
				}
			}

			// weighting - do this for each model-sample
			for(unsigned int isample=0; isample<nmbSamples; ++isample){

				vector<complex<double> > posm0amps(maxrank + 1); // positive refl vector m=0
				vector<complex<double> > posm1amps(maxrank + 1); // positive refl vector m=1

				vector<complex<double> > negm0amps(maxrank + 1); // negative refl vector m=0
				vector<complex<double> > negm1amps(maxrank + 1); // negative refl vector m=1

				complex<double> flatamp;

				for(unsigned int iw = 0; iw < nmbWaves; ++iw) {
					complex<double> decayamp = decayamps[iw];
					string w1 = waveNames[iw];
					if(w1 == "flat"){
						flatamp = prodAmps[isample]->at(iw);
						continue;
					}
					//cerr << w1 << "  " << decayamp << endl;
					double nrm = sqrt(normInt.val(w1, w1).real());
					complex<double> amp = decayamp / nrm * prodAmps[isample]->at(iw);
					if(isample == 0) {// fill wheights of individual waves
						weights[iw] = norm(amp);
					}

					if (reflectivities[iw] == 0) {
						if(ms[iw] == 0) {
							posm0amps[ranks[iw]] += amp;
						} else if(ms[iw] == 1) {
							posm1amps[ranks[iw]] += amp;
						}
					} else if(reflectivities[iw] == 1) {
						if(ms[iw] == 0) {
							negm0amps[ranks[iw]] += amp;
						} else if(ms[iw] == 1) {
							negm1amps[ranks[iw]] += amp;
						}
					}
				}
				// end loop over waves

				// incoherent sum over rank and diffrerent reflectivities
				weight = 0;
				if(hasfit) {
					for(int ir = 0; ir < maxrank + 1; ++ir) {
						weight += norm(posm0amps[ir] + posm1amps[ir]);
						//weight += norm(posm1amps[ir]);
						weight += norm(negm0amps[ir] + negm1amps[ir]);
						// weight += norm(negm1amps[ir]);
					}
				}
				// add flat
				weight += norm(flatamp);


				if(isample == 0) {
					hWeights->Fill(weight);
				}
				modelweights[isample] = weight;

			}// end loop over model samples
			outtree->Fill();
		}

		inTrees[i]->PrintCacheStats();
	}
	printSucc << "calculated weight for " << counter << " events" << endl;
	timer.Stop();
	printInfo << "this job consumed: ";
	timer.Print();

	outfile->cd();
	hWeights->Write();
	outtree->Write();
	outfile->Close();

	for(unsigned int iw = 0; iw < nmbWaves; ++iw) {
		if (ampfiles[iw] != NULL) {
			ampfiles[iw]->close();
			delete ampfiles[iw];
		}
	}
	ampfiles.clear();
	for(unsigned int isamples=0; isamples < nmbSamples; ++isamples) {
		delete prodAmps[isamples];
	}
	prodAmps.clear();

	return 0;
}
