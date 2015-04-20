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
//      calculates integral matrix for set of amplitudes
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <stdlib.h>
#include <unistd.h>

#include "TROOT.h"
#include "TFile.h"

#include "reportingUtils.hpp"
#include "reportingUtilsEnvironment.h"
#include "fileUtils.hpp"
#include "ampIntegralMatrix.h"


using namespace std;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{
	cerr << "computes integral matrix for given amplitude files" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " [-o output file -i TKey name -n max. # of events -r max. # of events -w weight file -v -h] "
	     << "amplitude files" << endl
	     << "    where:" << endl
	     << "        -o path    path to output file (default: './norm.int')"                      << endl
	     << "        -i name    integral TKey name (only for .root format, default: 'integral')"  << endl
	     << "        -n #       maximum number of events to process (default: all)"               << endl
	     << "        -r #       number of events to renormalize to (default: no renormalization)" << endl
	     << "        -w path    path to MC weight file for de-weighting (default: none)"          << endl
	     << "        -v         verbose; print debug output (default: false)"                     << endl
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
	const string progName        = argv[0];
	string       outFileName     = "./norm.int";
	string       integralName    = "integral";
	unsigned int maxNmbEvents    = 0;
	unsigned int nmbEventsRenorm = 0;
	string       weightFileName  = "";
	bool         debug           = false;
	extern char* optarg;
	extern int   optind;
	int          c;
	while ((c = getopt(argc, argv, "o:i:n:r:w:vh")) != -1)
		switch (c) {
		case 'o':
			outFileName = optarg;
			break;
		case 'i':
			integralName = optarg;
			break;
		case 'n':
			maxNmbEvents = atoi(optarg);
			break;
		case 'r':
			nmbEventsRenorm = atoi(optarg);
			break;
		case 'w':
			weightFileName = optarg;
			break;
		case 'v':
			debug = true;
			break;
		case 'h':
		default:
			usage(progName);
		}

	// switch debug output
	ampIntegralMatrix::setDebug(debug);

	// get input file names
	if (optind >= argc) {
		printErr << "you need to specify at least one amplitude file to process. Aborting..." << endl;;
		usage(progName, 1);
	}
	vector<string> rootAmpFileNames;
	vector<string> binAmpFileNames;
	while (optind < argc) {
		const string fileName = argv[optind++];
		const string fileExt  = extensionFromPath(fileName);
		if (fileExt == "root")
			rootAmpFileNames.push_back(fileName);
		else if (fileExt == "amp")
			binAmpFileNames.push_back(fileName);
		else
			printWarn << "input file '" << fileName << "' is neither a .root nor a .amp file. "
			          << "skipping." << endl;
	}
	if ((rootAmpFileNames.size() == 0) and (binAmpFileNames.size() == 0)) {
		printErr << "none of the specified input files is a .root or .amp file. Aborting...";
		usage(progName, 1);
	}

	// calculate integral
	ampIntegralMatrix integral;
	integral.integrate(binAmpFileNames, rootAmpFileNames, maxNmbEvents, weightFileName);
	if (nmbEventsRenorm > 0)
		integral.renormalize(nmbEventsRenorm);

	// write out integral
	const string outFileExt = extensionFromPath(outFileName);
	if (outFileExt == "root") {
		TFile* outFile = TFile::Open(outFileName.c_str(), "RECREATE");
		if (not outFile) {
			printErr << "cannot open output file '" << outFileName << "'. Aborting..." << endl;
			exit(1);
		}
		const int nmbBytes = integral.Write(integralName.c_str());
		outFile->Close();
		if (nmbBytes == 0) {
			printErr << "problems writing integral to TKey '" << integralName << "' "
			         << "in file '" << outFileName << "'" << endl;
			exit(1);
		} else
			printSucc << "wrote integral to TKey '" << integralName << "' "
			          << "in file '" << outFileName << "'" << endl;
	} else if (outFileExt == "int")
		integral.writeAscii(outFileName);
	else {
		printErr << "output file '" << outFileName << "' should be either a .root or a .int file. "
		         << "Aborting...";
		exit(1);
	}

	return 0;
}
