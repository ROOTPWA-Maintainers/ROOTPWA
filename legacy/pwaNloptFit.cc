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
#include <particleDataTable.h>
#include <pwaNloptFit.h>
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
	     << " -a normfile -A # normalisation events -r rank -C -P width -H -c -z -p PDG file -q -h]" << endl
	     << "    where:" << endl
	     << "        -l #       lower edge of mass bin [MeV/c^2]" << endl
	     << "        -u #       upper edge of mass bin [MeV/c^2]" << endl
	     << "        -w file    path to wavelist file" << endl
	     << "        -d dir     path to directory with decay amplitude files (default: '.')" << endl
	     << "        -o file    path to output file (default: 'fitresult.root')" << endl
	     << "        -S file    path to file with start values (default: none; highest priority)" << endl
	     << "        -s #       seed for random start values (default: 0)" << endl
	     << "        -n file    path to normalization integral file (default: 'norm.root')" << endl
	     << "        -a file    path to acceptance integral file (default: 'norm.root')" << endl
	     << "        -A #       number of input events to normalize acceptance to (default: use number of events from normalization integral file)" << endl
	     << "        -r #       rank of spin density matrix (default: 1)" << endl
	     << "        -C         use half-Cauchy priors (default: false)" << endl
	     << "        -P #       width of half-Cauchy priors (default: 0.5)" << endl
	     << "        -H         check analytical Hessian eigenvalues (default: false)" << endl
#ifdef USE_CUDA
	     << "        -c         enable CUDA acceleration (default: off)" << endl
#else
	     << "        -c         enable CUDA acceleration [not supported by your platform]" << endl
#endif
	     << "        -z         save space by not saving integral and covariance matrices (default: false)" << endl
	     << "        -p file    path to particle data table file (default: ./particleDataTable.txt)" << endl
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
	bool         cauchy           = false;
	double       cauchyWidth      = 0.5;
	bool         checkHessian     = false;                  // if true checks analytical Hessian eigenvalues
	bool         cudaEnabled      = false;                  // if true CUDA kernels are activated
	bool         saveSpace        = false;
	string       pdgFileName      = "./particleDataTable.txt";
	bool         verbose          = true;

	extern char* optarg;
	// extern int optind;
	int c;
	while ((c = getopt(argc, argv, "l:u:w:d:o:S:s:n:a:A:r:CP:Hczp:qh")) != -1)
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
		case 'C':
			cauchy = true;
			break;
		case 'P':
			cauchyWidth = atof(optarg);
			break;
		case 'H':
			checkHessian = true;
			break;
		case 'c':
#ifdef USE_CUDA
			cudaEnabled = true;
#endif
			break;
		case 'z':
			saveSpace = true;
			break;
		case 'p':
			pdgFileName = optarg;
			break;
		case 'q':
			verbose = false;
			break;
		case 'h':
			usage(progName);
			break;
		}
	if (normIntFileName.length() == 0) {
		normIntFileName = "norm.root";
		printWarn << "using default normalization integral file '" << normIntFileName << "'" << endl;
	}
	if (accIntFileName.length() == 0) {
		accIntFileName = "norm.root";
		printWarn << "using default acceptance normalization integral file "
		          << "'" << accIntFileName << "'" << endl;
	}
	if (waveListFileName.length() == 0) {
		printErr << "no wavelist file specified. Aborting..." << endl;
		usage(progName, 1);
	}
	// report parameters
	printInfo << "running " << progName << " with the following parameters:" << endl;
	cout << "    mass bin [" <<  massBinMin << ", " <<  massBinMax << "] MeV/c^2" << endl
	     << "    path to wave list file ......................... '" << waveListFileName << "'" << endl
	     << "    path to amplitude directory .................... '" << ampDirName       << "'" << endl
	     << "    path to output file ............................ '" << outFileName      << "'" << endl
	     << "    path to file with normalization integral ....... '" << normIntFileName  << "'" << endl
	     << "    path to file with acceptance integral .......... '" << accIntFileName   << "'" << endl
	     << "    number of acceptance norm. events .............. "  << numbAccEvents           << endl
	     << "    rank of spin density matrix .................... "  << rank                    << endl
	     << "    using half-Cauchy priors........................ "  << yesNo(cauchy)           << endl;
	if (cauchy) {
		cout << "    width of cauchy priors.......................... "  << cauchyWidth << endl;
	}
	cout << "    CUDA acceleration .............................. "  << enDisabled(cudaEnabled) << endl
	     << "    quiet .......................................... "  << yesNo(not verbose)      << endl;

	// ---------------------------------------------------------------------------
	// initialize particle data table
	particleDataTable::readFile(pdgFileName);

	// ---------------------------------------------------------------------------
	// read wave list to get wave names and thresholds
	printInfo << "reading amplitude names and thresholds from wave list file '" << waveListFileName << "'." << endl;
	ifstream waveListFile(waveListFileName.c_str());
	if (not waveListFile) {
		printErr << "cannot open file '" << waveListFileName << "'. Aborting..." << endl;
		exit(1);
	}

	//extract wave names and thresholds
	vector<pwaLikelihood<complex<double> >::waveDescThresType> waveDescThres;

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

			const string waveFileName = ampDirName + "/" + waveName + ".root";
			TFile* file = TFile::Open(waveFileName.c_str());
			if (file == NULL || file->IsZombie()) {
				printErr << "amplitude file '" << waveFileName << "' cannot be opened. Aborting..." << endl;
				exit(1);
			}

			const amplitudeMetadata* ampMeta = amplitudeMetadata::readAmplitudeFile(file, waveName);
			if (ampMeta == NULL) {
				printErr << "cannot read event data from event file '" << waveFileName << "'. Aborting..." << endl;
				exit(1);
			}

			waveDescription waveDesc(ampMeta);

			waveDescThres.push_back(boost::tuples::make_tuple(waveName, waveDesc, threshold));

			file->Close();
			delete file;
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
	ampIntegralMatrix* normIntMatrix = NULL;
	normIntFile->GetObject(ampIntegralMatrix::integralObjectName.c_str(), normIntMatrix);
	if (not normIntMatrix) {
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
	ampIntegralMatrix* accIntMatrix = NULL;
	accIntFile->GetObject(ampIntegralMatrix::integralObjectName.c_str(), accIntMatrix);
	if (not accIntMatrix) {
		printErr << "cannot find integral object in TKey '" << ampIntegralMatrix::integralObjectName << "' in file "
		         << "'" << accIntFileName << "'. Aborting..." << endl;
		exit(1);
	}

	// ---------------------------------------------------------------------------
	// setup likelihood function
	const double massBinCenter = (massBinMin + massBinMax) / 2;
	printInfo << "creating and setting up likelihood function" << endl;
	pwaLikelihood<complex<double> > L;
	if (not verbose)
		L.setQuiet();
	L.useNormalizedAmps(true);
#ifdef USE_CUDA
	L.enableCuda(cudaEnabled);
#endif
	if (not L.init(waveDescThres, rank, massBinCenter)) {
		printErr << "error while initializing likelihood function. Aborting..." << endl;
		exit(1);
	}

	if (cauchy) {
		L.setPriorType(L.HALF_CAUCHY);
		L.setCauchyWidth(cauchyWidth);
	}

	if (not L.addNormIntegral(*normIntMatrix)) {
		printErr << "error while adding normalization integral to likelihood function. Aborting..." << endl;
		exit(1);
	}
	if (not L.addAccIntegral(*accIntMatrix, numbAccEvents)) {
		printErr << "error while adding acceptance integral to likelihood function. Aborting..." << endl;
		exit(1);
	}
	for (size_t i=0; i<waveDescThres.size(); ++i) {
		const string waveName     = boost::tuples::get<0>(waveDescThres[i]);

		const string waveFileName = ampDirName + "/" + waveName + ".root";
		TFile* file = TFile::Open(waveFileName.c_str());
		if (file == NULL || file->IsZombie()) {
			printErr << "amplitude file '" << waveFileName << "' cannot be opened. Aborting..." << endl;
			exit(1);
		}

		const amplitudeMetadata* ampMeta = amplitudeMetadata::readAmplitudeFile(file, waveName);
		if (ampMeta == NULL) {
			printErr << "cannot read event data from event file '" << waveFileName << "'. Aborting..." << endl;
			exit(1);
		}

		if (not L.addAmplitude(*ampMeta)) {
			printErr << "error while adding amplitude of wave '" << waveName << "' to likelihood function. Aborting..." << endl;
			exit(1);
		}

		file->Close();
		delete file;
	}
	if (not L.finishInit()) {
		printErr << "error while finished initialization of likelihood function. Aborting..." << endl;
		exit(1);
	}

	// ---------------------------------------------------------------------------
	// likelihood minimization
	fitResultPtr result = hli::pwaNloptFit(L,
	                                       massBinMin,
	                                       massBinMax,
	                                       startValSeed,
	                                       startValFileName,
	                                       checkHessian,
	                                       saveSpace,
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
