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
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      Calculates integral matrix for set of amplitudes
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <stdlib.h>
#include <unistd.h>

#include "reportingUtils.hpp"
#include "fileUtils.hpp"
#include "normalizationIntegral.h"


using namespace std;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{
	cerr << "generates set of all allowed waves given the constraints" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " -o output file [-n max. # of events -r max. # of events -w weight file -v -h] "
	     << "amplitude files" << endl
	     << "    where:" << endl
	     << "        -o path    path to output file (default: './norm.int')"           << endl
	     << "        -n #       maximum number of events to process (default: all)"    << endl
	     << "        -r #       number of events to renormalize to (default: no done)" << endl
	     << "        -w path    path to MC weight file (default: none)"                << endl
	     << "        -v         verbose; print debug output (default: false)"          << endl
	     << "        -h         print help" << endl
	     << endl;
	exit(errCode);
}


int
main(int    argc,
     char** argv)
{
	printCompilerInfo();
	printSvnVersion();
	
	// parse command line options
	const string progName        = argv[0];
	string       outFileName     = "./norm.int";
	unsigned int maxNmbEvents    = 0;
	unsigned int nmbEventsRenorm = 0;
	string       weightFileName  = "";
	bool         debug           = false;
	extern char* optarg;
	extern int   optind;
	int          c;
	while ((c = getopt(argc, argv, "o:n:r:w:vh")) != -1)
		switch (c) {
		case 'o':
			outFileName = optarg;
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
	if (debug)
		normalizationIntegral::setDebug(true);
	
	// get input file names
	if (optind >= argc) {
		printErr << "you need to specify at least one amplitude file to process. aborting." << endl;;
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
		printErr << "none of the specified input files is a .root or .amp file. aborting.";
		usage(progName, 1);
	}

	// calculate integral
	normalizationIntegral integral;
	integral.integrate(binAmpFileNames);
	integral.writeAscii(outFileName);

	return 0;
}
