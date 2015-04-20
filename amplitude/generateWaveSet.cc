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
//      generates a set of all waves allowed by convervation laws and
//      constrained by user specified parameters
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <unistd.h>
#include <fstream>

#include "reportingUtilsEnvironment.h"
#include "particleDataTable.h"
#include "waveSetGenerator.h"
#include "waveDescription.h"


using namespace std;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{
	cerr << "generates set of all allowed waves given by template .key file" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " -k template key file -o output directory [-p PDG file -n -v -h]" << endl
	     << "    where:" << endl
	     << "        -k file    path to template key file" << endl
	     << "        -p file    path to particle data table file (default: ./particleDataTable.txt)" << endl
	     << "        -d file    path to decay config file (default: decay info will not be used)" << endl
	     << "                   isobars without defined decays will be allowed to decay to all possibilities" << endl
	     << "        -f         force decay check (only works with -d option); isobars without defined decay modes will be ignored" << endl
	     << "        -o dir     path to directory where key files will be written (default: '.')" << endl
	     << "        -t file    path to waveset LaTeX output file" << endl
	     << "        -n         use new key file name convention (default: false)" << endl
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
	const string progName                 = argv[0];
	unsigned int entriesPerPage           = 20;
	string       keyFileName              = "";
	string       pdgFileName              = "./particleDataTable.txt";
	string       decayFileName            = "";
	string       texFileName              = "waveset.tex";
	string       outDirName               = ".";
	bool         newKeyFileNameConvention = false;
	bool         debug                    = false;
	bool         forceDecayCheck          = false;
	extern char* optarg;
	int          c;
	while ((c = getopt(argc, argv, "k:p:d:o:t:fnvh")) != -1)
		switch (c) {
		case 'k':
			keyFileName = optarg;
			break;
		case 'p':
			pdgFileName = optarg;
			break;
		case 'd':
			decayFileName = optarg;
			break;
		case 't':
			texFileName = optarg;
			break;
		case 'o':
			outDirName = optarg;
			break;
		case 'n':
			newKeyFileNameConvention = true;
			break;
		case 'v':
			debug = true;
			break;
		case 'f':
			forceDecayCheck = true;
			break;
		case 'h':
		default:
			usage(progName);
		}

	// initialize particle data table
	particleDataTable::readFile(pdgFileName);
	bool useDecays = false;
	if (decayFileName != "") {
		if (particleDataTable::readDecayModeFile(decayFileName))
			useDecays = true;
		else {
			printErr << "could not read particle decay modes from file '" << decayFileName << "'. "
			         << "Aborting..." << endl;
			exit(1);
		}
	}
	if (debug)
		printDebug << "particle data table:" << endl << particleDataTable::instance();
	particleDataTable::setDebug(debug);

	// print latex header
	ofstream wavelistTeX;
	bool     doTeX = false;
	if (texFileName != "") {
		doTeX = true;
		wavelistTeX.open(texFileName.c_str());
		wavelistTeX << "\\documentclass[10pt,a4paper]{article}" << endl
		            << "\\usepackage{amsmath,amsthm,amssymb}"   << endl
		            << "\\def\\dst{\\displaystyle}"             << endl
		            << "\\def\\vsp{\\hbox{\\vrule height12.5pt depth3.5pt width0pt}}" << endl
		            << "\\def\\ells#1#2{\\Big[\\hskip-5pt\\vsp\\begin{array}{c}\\dst#1\\\\[-4pt]\\dst#2\\end{array}\\vsp\\hskip-5pt\\Big]}" << endl
		            << "\\begin{document}"                      << endl
		            << "\\begin{align*} \n \\begin{aligned}"    << endl;
	}


	printInfo << "generating wave set from template key file '" << keyFileName << "'" << endl;
	waveSetGenerator waveSetGen;
	waveSetGenerator::setDebug(debug);
	if (not waveSetGen.setWaveSetParameters(keyFileName)) {
		printErr << "could not initialize wave set generator. Aborting..." << endl;
		exit(1);
	}

	if (useDecays)
		waveSetGen.setForceDecayCheck(forceDecayCheck);
	printInfo << waveSetGen;
	waveSetGen.generateWaveSet();

	printInfo << "checking generated waves..." << endl;
	vector<isobarDecayTopology>& decayTopos            = waveSetGen.waveSet();
	unsigned int                 nmbInconsistentDecays = 0;
	for (unsigned int i = 0; i < decayTopos.size(); ++i) {
		bool isConsistent = decayTopos[i].checkTopology() and decayTopos[i].checkConsistency();
		cout << "    " << setw(4) << i << ": "
		     << waveDescription::waveNameFromTopology(decayTopos[i], newKeyFileNameConvention) << " ... ";
		if (isConsistent) {
		  cout << "okay" << endl;
		  if (doTeX) {
		     wavelistTeX << waveDescription::waveLaTeXFromTopology(decayTopos[i]) << "\\\\" << endl;
		     if ((i + 1) % entriesPerPage ==0 ) {
			     wavelistTeX << "\\end{aligned} \n \\end{align*}"
			                 << "\\pagebreak\n"
			                 << "\\begin{align*} \n \\begin{aligned}" << endl;
		     }
		  }
		}	else {
			cout << "problems!" << endl;
			++nmbInconsistentDecays;
		}
	}

	printInfo << "writing .key files for generated waves to directory '" << outDirName << "'" << endl;
	waveSetGen.writeKeyFiles(outDirName, newKeyFileNameConvention);

	printInfo << ((nmbInconsistentDecays == 0) ? "successfully " : "")
	          << "generated " << decayTopos.size() << " waves";
	if (nmbInconsistentDecays != 0)
		cout << ". " << nmbInconsistentDecays << " waves have problems.";
	cout<< endl;

	if (doTeX) {
	  wavelistTeX << "\\end{aligned} \n \\end{align*}" << endl;
	  wavelistTeX << "\\end{document}" << endl;
	  wavelistTeX.close();
	}
	return 0;
}
