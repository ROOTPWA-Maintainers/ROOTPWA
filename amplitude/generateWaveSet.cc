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
//      generates a set of all waves allowed by convervation laws and
//      constrained by user specified parameters
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


// #include <iostream>
// #include <fstream>
#include <unistd.h>
// #include <vector>
// #include <complex>

#include "svnVersion.h"
// #include "particleDataTable.h"
// #include "keyFileParser.h"
#include "waveSetGenerator.h"


using namespace std;
using namespace boost;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{
	cerr << "generates set of all allowed waves given the constraints" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " -k template key file -o output directory [-v -h]" << endl
	     << "    where:" << endl
	     << "        -k file    path to template key file" << endl
	     << "        -o dir     path to directory where key files will be written (default: '.')" << endl
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
	printSvnVersion();
	
	// parse command line options
	const string progName    = argv[0];
	string       keyFileName = "";
	string       outDirName  = ".";
	bool         debug       = false;
	extern char* optarg;
	//extern int   optind;
	int          c;
	while ((c = getopt(argc, argv, "k:o:vh")) != -1)
		switch (c) {
		case 'k':
			keyFileName = optarg;
			break;
		case 'o':
			outDirName = optarg;
			break;
		case 'v':
			debug = true;
			break;
		case 'h':
		default:
			usage(progName);
		}

	return 0;
}
