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

/** @brief Simple partial wave event generator (homogeneous in m)
 */


#include<cstdlib>
#include<iostream>
#include<fstream>
#include<getopt.h>
#include<unistd.h>

#include <boost/progress.hpp>

#include<fileUtils.hpp>
#include<generator.h>
#include<generatorManager.h>
#include<particleDataTable.h>
#include<randomNumberGenerator.h>
#include<reportingUtils.hpp>


using namespace rpwa;
using namespace std;


void printUsage(char* prog, int errCode = 0)
{
	cerr << "usage:" << endl
	     << prog
	     << " -n # [-a # -m # -M # -B # -s #] -o <file> -p <file> -w <file> -k <path> -i <file> -r <file>" << endl
	     << "    where:" << endl
	     << "        -n #       (max) number of events to generate (default: 100)" << endl
	     << "        -a #       (max) number of attempts to do (default: infinity)" << endl
	     << "        -m #       maxWeight FEATURE DISABLED" << endl
	     << "        -o <file>  ASCII output file (if not specified, generated automatically)" << endl
	     << "        -p         path to particle data table file (default: ./particleDataTable.txt)" << endl
	     << "        -w <file>  wavelist file (contains production amplitudes) FEATURE DISABLED" << endl
	     << "        -w <file.root>  to use TFitBin tree as input FEATURE PERMANENTLY REMOVED" << endl
	     << "        -c         if 1 a comgeant eventfile (.fort.26) is written with same naming as the root file (default 0)" << endl
	     << "        -k <path>  path to keyfile directory (all keyfiles have to be there) FEATURE DISABLED" << endl
	     << "        -i <file>  integral file FEATURE DISABLED" << endl
	     << "        -r <file>  reaction config file" << endl
	     << "        -s #       set seed " << endl
	     << "        -M #       lower boundary of mass range in MeV (overwrites values from config file)" << endl
	     << "        -B #       width of mass bin in MeV" << endl
	     << endl
	     << "A comment regarding the disabled features: these options have been taken out\n"
	     << "for the time being. If you want to get them back, check GIT revision\n"
	     << "cb48b651809e1058ab740441c2b6bd8a1579d46d or use the _v1 branch (or one of its\n"
	     << "tags)." << endl
	     << endl;
	exit(errCode);
}


int main(int argc, char** argv)
{

	unsigned int nEvents = 100;
	unsigned int maxAttempts = 0;
	string outputEvtFileName = "";
	string outputWhtFileName = "";
	string outputComgeantFileName = "";
	string pdgFileName = "./particleDataTable.txt";
	string reactionFile;
	int seed = 123456;
	bool seedSet = false;
	int massLower = 0;
	int massBinWidth = 0;
	bool overrideMass = false;
	bool writeComgeantOut = false;

	int c;
	while ((c = getopt(argc, argv, "n:a:o:p:w:k:i:r:m:s:M:B:hc")) != -1) {
		switch (c) {
			case 'n':
				nEvents = atoi(optarg);
				break;
			case 'a':
				maxAttempts = atoi(optarg);
				break;
			case 's':
				seed = atoi(optarg);
				seedSet = true;
				break;
			case 'o':
				outputEvtFileName = optarg;
				break;
			case 'p':
				pdgFileName = optarg;
				break;
			case 'w':
				printErr << "this feature has been removed for the time being. "
				         << "If you want it back, check GIT revision "
				         << "cb48b651809e1058ab740441c2b6bd8a1579d46d." << endl;
				exit(10);
				break;
			case 'i':
				printErr << "this feature has been removed for the time being. "
				         << "If you want it back, check GIT revision "
				         << "cb48b651809e1058ab740441c2b6bd8a1579d46d." << endl;
				exit(10);
				break;
			case 'k':
				printErr << "this feature has been removed for the time being. "
				         << "If you want it back, check GIT revision "
				         << "cb48b651809e1058ab740441c2b6bd8a1579d46d." << endl;
				exit(10);
				break;
			case 'r':
				reactionFile = optarg;
				break;
			case 'm':
				printErr << "this feature has been removed for the time being. "
				         << "If you want it back, check GIT revision "
				         << "cb48b651809e1058ab740441c2b6bd8a1579d46d." << endl;
				exit(10);
				break;
			case 'c':
				writeComgeantOut = true;
				break;
			case 'M':
				massLower = atoi(optarg);
				overrideMass = true;
				break;
			case 'B':
				massBinWidth = atoi(optarg);
				overrideMass = true;
				break;

			case 'h':
				printUsage(argv[0]);
				break;
			default:
				printUsage(argv[0], 5);
				break;
		}
	}

	if(maxAttempts && (maxAttempts < nEvents)) {
		printWarn << "Maximum attempts is smaller than the number of events. Setting it to infinity." << endl;
		maxAttempts = 0;
	}

	if(overrideMass && not (massLower && massBinWidth)) {
		printErr << "'-M' and '-B' can only be set together. Aborting..." << endl;
		exit(2);
	}

	if(not seedSet) {
		printInfo << "Setting random seed to " << seed << endl;
	}
	randomNumberGenerator::instance()->setSeed(seed);

	particleDataTable::readFile(pdgFileName);
	generatorManager generatorMgr;
	if(not generatorMgr.readReactionFile(reactionFile)) {
		printErr << "could not read reaction file. Aborting..." << endl;
		exit(1);
	}
	if(not generatorMgr.initializeGenerator()) {
		printErr << "could not initialize generator. Aborting..." << endl;
		exit(1);
	}

	if(overrideMass) {
		generatorMgr.overrideMassRange(massLower / 1000., (massLower + massBinWidth) / 1000.);
	}

	if(outputEvtFileName == "") {
		stringstream fileName;
		fileName << massLower << "." << massLower + massBinWidth << ".genbod.evt";
		outputEvtFileName = fileName.str();
	}
	ofstream outputEvtFile(outputEvtFileName.c_str());
	printInfo << "output event file: " << outputEvtFileName << endl;

	ofstream outputComgeantFile;
	if(writeComgeantOut) {
		outputComgeantFileName = changeFileExtension(outputEvtFileName, ".fort.26");
		printInfo << "output comgeant file: " << outputComgeantFileName << endl;
		outputComgeantFile.open(outputComgeantFileName.c_str());
	}

	boost::progress_display* progressIndicator = new boost::progress_display(nEvents, cout, "");

	unsigned int attempts = 0;
	unsigned int eventsGenerated = 0;
	for(; eventsGenerated < nEvents; ++eventsGenerated) {

		attempts += generatorMgr.event();
		const generator& gen = generatorMgr.getGenerator();
		rpwa::particle beam = gen.getGeneratedBeam();
		std::vector<rpwa::particle> finalState = gen.getGeneratedFinalState();
		generator::convertEventToAscii(outputEvtFile, beam, finalState);
		if(writeComgeantOut) {
			rpwa::particle recoil = gen.getGeneratedRecoil();
			TVector3 vertex = gen.getGeneratedVertex();
			generator::convertEventToComgeant(outputComgeantFile, beam, recoil, vertex, finalState, false);
		}
		if(maxAttempts && (attempts > maxAttempts)) {
			printWarn << "reached maximum attempts. Aborting..." << endl;
			break;
		}
		++(*progressIndicator);

	}

	outputEvtFile.close();
	if(writeComgeantOut) {
		outputComgeantFile.close();
	}

	printSucc << "generated " << eventsGenerated << " events." << endl;
	printInfo << "attempts: " << attempts << endl;
	printInfo << "efficiency: " << setprecision(3) << 100 * ((double)eventsGenerated / attempts) << "%" << endl;

}

