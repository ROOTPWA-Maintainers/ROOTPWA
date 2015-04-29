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
#include <string>
#include <unistd.h>

#include <RVersion.h>
#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>

#include <amplitudeMetadata.h>
#include <fitResult.h>
#include <pwaFit.h>
#include <reportingUtils.hpp>
#include <reportingUtilsEnvironment.h>


using namespace std;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{

	cerr << "performs PWA fit for given mass bin and list of waves" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " -l # -u # -w wavelist [-d amplitude directory -o outfile -S start value file -s seed -n normfile"
	     << " -a normfile -A # normalisation events -r rank -H -q -h]" << endl
	     << "    where:" << endl
	     << "        -l #       lower edge of mass bin [MeV/c^2]" << endl
	     << "        -u #       upper edge of mass bin [MeV/c^2]" << endl
	     << "        -w file    path to wavelist file" << endl
	     << "        -d dir     path to directory with decay amplitude files (default: '.')" << endl
	     << "        -o file    path to output file (default: 'fitresult.root')" << endl
	     << "        -S file    path to file with start values (default: none; highest priority)" << endl
	     << "        -s #       seed for random start values (default: 0)" << endl
	     << "        -n file    path to normalization integral file (default: 'norm.int')" << endl
	     << "        -a file    path to acceptance integral file (default: 'norm.int')" << endl
	     << "        -A #       number of input events to normalize acceptance to (default: use number of events from normalization integral file)" << endl
	     << "        -r #       rank of spin density matrix (default: 1)" << endl
	     << "        -H         check analytical Hessian eigenvalues (default: false)" << endl
	     << "        -q         run quietly (default: false)" << endl
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

	// ---------------------------------------------------------------------------
	// internal parameters
	const string progName      = argv[0];
	const string valTreeName   = "pwa";
	const string valBranchName = "fitResult_v2";

	// ---------------------------------------------------------------------------
	// parse command line options
	double       massBinMin       = 0.;                     // [MeV/c^2]
	double       massBinMax       = 0.;                     // [MeV/c^2]
	string       waveListFileName = "";                     // wavelist filename
	string       ampDirName       = ".";                    // decay amplitude directory name
	string       outFileName      = "fitresult.root";       // output filename
	string       startValFileName = "";                     // file with start values
	int          startValSeed     = 0;
	string       normIntFileName  = "";                     // file with normalization integrals
	string       accIntFileName   = "";                     // file with acceptance integrals
	unsigned int numbAccEvents    = 0;                      // number of events used for acceptance integrals
	unsigned int rank             = 1;                      // rank of fit
	bool         checkHessian     = false;                  // if true checks analytical Hessian eigenvalues
	bool         verbose          = true;

	extern char* optarg;
	// extern int optind;
	int c;
	while ((c = getopt(argc, argv, "l:u:w:d:o:S:s:n:a:A:r:Hqh")) != -1)
		switch (c) {
		case 'l':
			massBinMin = atof(optarg);
			break;
		case 'u':
			massBinMax = atof(optarg);
			break;
		case 'w':
			waveListFileName = optarg;
			break;
		case 'd':
			ampDirName = optarg;
			break;
		case 'o':
			outFileName = optarg;
			break;
		case 'S':
			startValFileName = optarg;
			break;
		case 's':
			startValSeed = atoi(optarg);
			break;
		case 'n':
			normIntFileName = optarg;
			break;
		case 'a':
			accIntFileName = optarg;
			break;
		case 'A':
			numbAccEvents = atoi(optarg);
			break;
		case 'r':
			rank = atoi(optarg);
			break;
		case 'H':
			checkHessian = true;
			break;
		case 'q':
			verbose = false;
			break;
		case 'h':
			usage(progName);
			break;
		}
	if (normIntFileName.length() <= 1) {
		normIntFileName = "norm.int";
		printWarn << "using default normalization integral file '" << normIntFileName << "'" << endl;
	}
	if (accIntFileName.length() <= 1) {
		accIntFileName = "norm.int";
		printWarn << "using default acceptance normalization integral file "
		          << "'" << accIntFileName << "'" << endl;
	}
	if (waveListFileName.length() <= 1) {
		printErr << "no wavelist file specified. Aborting..." << endl;
		usage(progName, 1);
	}

	// read wave list: wave names and thresholds
	printInfo << "reading amplitude names and thresholds from wave list file '" << waveListFileName << "'." << endl;
	ifstream waveListFile(waveListFileName.c_str());
	if (not waveListFile) {
		printErr << "cannot open file '" << waveListFileName << "'. Aborting..." << endl;
		exit(1);
	}

	//extract wave names and thresholds
	vector<string> waveNames;
	vector<double> waveThresholds;

	string       line;
	unsigned int lineNmb = 0;
	while (getline(waveListFile, line)) {
		if (line[0] == '#')  // comments start with #
			continue;
		stringstream lineStream;
		lineStream.str(line);
		string waveName;
		if (lineStream >> waveName) {
			double threshold;
			// !!! it would be safer to make the threshold value in the wave list file mandatory
			if (not (lineStream >> threshold))
				threshold = 0;
			if (verbose)
				printDebug << "reading line " << setw(3) << lineNmb + 1 << ": " << waveName<< ", "
				           << "threshold = " << setw(4) << threshold << " MeV/c^2" << endl;
			waveNames.push_back(waveName);
			waveThresholds.push_back(threshold);
		} else
			printWarn << "cannot parse line '" << line << "' in wave list file '" << waveListFileName << "'" << endl;
		++lineNmb;
	}
	waveListFile.close();

	TFile* normIntFile = TFile::Open(normIntFileName.c_str(), "READ");
	if (not normIntFile or normIntFile->IsZombie()) {
		printErr << "could not open normalization integral file '" << normIntFileName << "'. "
		         << "Aborting..." << endl;
		exit(1);
	}
	ampIntegralMatrix* normMatrix = NULL;
	normIntFile->GetObject(ampIntegralMatrix::integralObjectName.c_str(), normMatrix);
	if (not normMatrix) {
		printErr << "cannot find integral object in TKey '" << ampIntegralMatrix::integralObjectName << "' in file "
		         << "'" << normIntFileName << "'. Aborting..." << endl;
		exit(1);
	}

	TFile* accIntFile = TFile::Open(accIntFileName.c_str(), "READ");
	if (not accIntFile or accIntFile->IsZombie()) {
		printErr << "could not open acceptance integral file '" << accIntFileName << "'. "
		         << "Aborting..." << endl;
		exit(1);
	}
	ampIntegralMatrix* accMatrix = NULL;
	accIntFile->GetObject(ampIntegralMatrix::integralObjectName.c_str(), accMatrix);
	if (not accMatrix) {
		printErr << "cannot find integral object in TKey '" << ampIntegralMatrix::integralObjectName << "' in file "
		         << "'" << accIntFileName << "'. Aborting..." << endl;
		exit(1);
	}

	// open amplitude files
	map<string, TTree*> ampTrees;
	for (size_t i=0; i<waveNames.size(); ++i) {
		const string waveName    = waveNames[i];
		const string ampFileName = waveName + ".root";
		TFile* file = TFile::Open(ampFileName.c_str());
		if (file == NULL || file->IsZombie()) {
			printErr << "amplitude file '" << ampFileName << "' cannot be opened. Aborting..." << endl;
			exit(1);
		}

		const amplitudeMetadata* ampMeta = amplitudeMetadata::readAmplitudeFile(file, waveName);
		if (ampMeta == NULL) {
			printErr << "cannot read event data from event file '" << ampFileName << "'. Aborting..." << endl;
			exit(1);
		}

		ampTrees[waveName] = ampMeta->amplitudeTree();
	}

	fitResultPtr result = hli::pwaFit(ampTrees,
	                                  *normMatrix,
	                                  *accMatrix,
	                                  waveNames,
	                                  waveThresholds,
	                                  massBinMin,
	                                  massBinMax,
	                                  startValSeed,
	                                  startValFileName,
	                                  numbAccEvents,
	                                  rank,
	                                  checkHessian,
	                                  verbose);

	// ---------------------------------------------------------------------------
	// write out result
	printInfo << "writing result to '" << outFileName << "'" << endl;
	{
		// open output file and create tree for writing
		TFile* outFile = new TFile(outFileName.c_str(), "UPDATE");
		if ((not outFile) or outFile->IsZombie())
			printWarn << "cannot open output file '" << outFileName << "'. "
			          << "no results will be written." << endl;
		else {
			// check whether output tree already exists
			TTree*     tree;
			fitResult* resultPtr = result.get();
			outFile->GetObject(valTreeName.c_str(), tree);
			if (not tree) {
				printInfo << "file '" << outFileName << "' is empty. "
				          << "creating new tree '" << valTreeName << "' for PWA result." << endl;
				tree = new TTree(valTreeName.c_str(), valTreeName.c_str());
				tree->Branch(valBranchName.c_str(), &resultPtr);
			} else {
				tree->SetBranchAddress(valBranchName.c_str(), &resultPtr);
			}

			// write result to file
			tree->Fill();
			tree->Write("", TObject::kOverwrite);
			outFile->Close();
		}
	}

	return 0;
}
