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
//      reads in data files in .evt or tree format, calculates
//      amplitudes for each event based on given key file, and writes
//      out amplitudes in PWA2000 binary or ascii format
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <iostream>
#include <unistd.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cctype>

#include "TSystem.h"

#include "svnVersion.h"
#include "utilities.h"
#include "particleDataTable.h"
#include "keyFileParser.h"


using namespace std;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{
  cerr << "creates graphical representation of decay graph(s) specified in given key file(s)" << endl
       << "output files will have the same name as the corresponding key files" << endl
       << endl
       << "usage:" << endl
       << progName
       << " [-p PDG file -f output format -v -h] key file(s)" << endl
       << "    where:" << endl
       << "        -p file    path to particle data table file (default: ./particleDataTable.txt)" << endl
       << "        -f format  output format: dot, ps, eps (default), svg, fig, dia, png, gif" << endl
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
  const string progName     = argv[0];
  string       pdgFileName  = "./particleDataTable.txt";
  string       outFormat    = "eps";
  bool         debug        = false;
  extern char* optarg;
  extern int   optind;
  int          c;
  while ((c = getopt(argc, argv, "p:f:vh")) != -1)
    switch (c) {
    case 'p':
      pdgFileName = optarg;
      break;
    case 'f':
      outFormat = optarg;
      break;
    case 'v':
      debug = true;
      break;
    case 'h':
    default:
      usage(progName);
    }

  // get input file names
  if (optind >= argc) {
    printErr << "you need to specify at least one key file to process. aborting." << endl;;
    usage(progName, 1);
  }
  vector<string> keyFileNames;
  while (optind < argc) {
    const string fileName = argv[optind++];
    if (fileName.substr(fileName.length() - 4) == ".key")
      keyFileNames.push_back(fileName);
    else
	    printWarn << "input file '" << fileName << "' is not a .key file. "
	              << "skipping." << endl;
  }
  if (keyFileNames.size() == 0) {
    printErr << "none of the specified input files is a .key file. aborting.";
    usage(progName, 1);
  }

  // check output format
  transform(outFormat.begin(), outFormat.end(), outFormat.begin(), (int(*)(int))tolower);
  const string validFormats[] = {"dot", "ps", "eps", "svg", "fig", "dia", "png", "gif"};
  bool         isValidFormat  = false;
  for (unsigned int i = 0; i < sizeof(validFormats) / sizeof(validFormats[0]); ++i)
    if (outFormat == validFormats[i]) {
      isValidFormat = true;
      break;
    }
  if (not isValidFormat) {
    printErr << "requested format '" << outFormat << "' is not supported." << endl;
    usage(progName, 1);
  }

  // initialize particle data table
  particleDataTable& pdt = particleDataTable::instance();
  pdt.readFile(pdgFileName);

  // loop over key files
  unsigned int countSuccess = 0;
  for (unsigned int i = 0; i < keyFileNames.size(); ++i) {

    // parse key file and create decay topology
    keyFileParser&         parser = keyFileParser::instance();
    isobarDecayTopologyPtr decayTopo;
    if (not parser.parse(keyFileNames[i]) or not parser.constructDecayTopology(decayTopo)) {
      printErr << "problems constructing decay topology from key file '" << keyFileNames[i] << "'. "
               << "skipping." << endl;
      continue;
    }
    const string dotFileName = keyFileNames[i].substr(0, keyFileNames[i].length() - 4) + ".dot";
    if (debug)
	    printInfo << "writing graph to file '" << dotFileName << "'" << endl;
    if (not decayTopo->writeGraphViz(dotFileName)) {
	    printWarn << "there were problems writing graph to file '" << dotFileName << "'. "
	              << "skipping." << endl;
	    continue;
    }
    
    // convert file to output format
    const string outFileName = keyFileNames[i].substr(0, keyFileNames[i].length() - 4)
                               + "." + outFormat;
    if (debug)
	    printInfo << "converting graph to file '" << outFileName << "'" << endl;
    stringstream cmd;
    if (   (outFormat == "ps" ) or (outFormat == "eps")
        or (outFormat == "fig") or (outFormat == "dia")
        or (outFormat == "png") or (outFormat == "gif")
        or (outFormat == "svg")) {
	    cmd << "dot -T" << outFormat << " -o " << outFileName << " " << dotFileName;
	    if (debug)
		    printInfo << "executing command '" << cmd.str() << "'" << endl;
	    if (gSystem->Exec(cmd.str().c_str()) != 0)
		    printWarn << "command '" << cmd.str() << "' was not successful." << endl;
	    else
		    ++countSuccess;
	    // cleanup
	    cmd.str("");
	    cmd << "rm " << (debug ? "--verbose " : "") << dotFileName;
	    if (gSystem->Exec(cmd.str().c_str()) != 0)
		    printWarn << "command '" << cmd.str() << "' was not successful." << endl;
    }
  }

  printInfo << "successfully created " << countSuccess << " out of " << keyFileNames.size()
            << " diagram files" << endl;
  
  return 0;
}
