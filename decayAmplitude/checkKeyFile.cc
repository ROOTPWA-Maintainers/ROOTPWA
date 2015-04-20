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
//      checks validity of amplitude specified by key file by
//      performing a number of consistency checks
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <iostream>
#include <fstream>
#include <unistd.h>
#include <vector>
#include <complex>

#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TStopwatch.h"

#include "reportingUtilsEnvironment.h"
#include "fileUtils.hpp"
#include "particleDataTable.h"
#include "evtTreeHelper.h"
#include "waveDescription.h"
#include "isobarHelicityAmplitude.h"
#include "isobarCanonicalAmplitude.h"


using namespace std;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{
	cerr << "verify that given .key files do not violate physics" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " [-d test data -n max. # of events -p PDG file -t tree name -e max. diff. -l leaf names -v -h] key file(s)" << endl
	     << "    where:" << endl
	     << "        -d file    path to file with test data (.evt or .root format)" << endl
	     << "        -n #       maximum number of events to read (default: all)" << endl
	     << "        -p file    path to particle data table file (default: ./particleDataTable.txt)" << endl
	     << "        -t name    name of tree in ROOT data files (default: rootPwaEvtTree)" << endl
	     << "        -e #       maximum deviation of amplitude ratios from 1 (default: 1E-6)" << endl
	     << "        -l names   semicolon separated object/leaf names in input data (default: 'prodKinParticles;prodKinMomenta;decayKinParticles;decayKinMomenta')" << endl
	     << "        -v         verbose; print debug output (default: false)" << endl
	     << "        -h         print help" << endl
	     << endl;
	exit(errCode);
}


bool testAmplitude(TTree*              inTree,
                   const string&       keyFileName,
                   vector<string>&     keyFileErrors,
                   const long int      maxNmbEvents            = -1,
                   const bool          debug                   = false,
                   const double        maxDelta                = 1e-6,
                   const TClonesArray* prodKinPartNames        = 0,
                   const TClonesArray* decayKinPartNames       = 0,
                   const string&       prodKinMomentaLeafName  = "prodKinMomenta",
                   const string&       decayKinMomentaLeafName = "decayKinMomenta")
{
	// parse key file and create decay topology and amplitude instances
	waveDescription        waveDesc;
	isobarDecayTopologyPtr decayTopo;
	if (   not waveDesc.parseKeyFile(keyFileName)
	    or not waveDesc.constructDecayTopology(decayTopo)) {
		printWarn << "problems constructing decay topology from key file '" << keyFileName << "'. "
		          << "skipping." << endl;
		keyFileErrors.push_back("parsing errors");
		return false;
	}

	printInfo << "checking decay topology" << endl;
	bool success = true;
	if (not decayTopo->checkTopology()) {
		keyFileErrors.push_back("problematic topology");
		success = false;
	}
	if (not decayTopo->checkConsistency()) {
		keyFileErrors.push_back("inconsistent decay");
		success = false;
	}
	if (not success)
		return false;

	if (not inTree or not prodKinPartNames or not decayKinPartNames)
		return true;  // no data; cannot check symmetry properties of amplitude

	printInfo << "calculating amplitudes for ";
	if (maxNmbEvents < 0)
		cout << "all";
	else
		cout << maxNmbEvents;
	cout << " events" << endl;

	// construct amplitude
	isobarAmplitudePtr amplitude;
	waveDesc.constructAmplitude(amplitude, decayTopo);


	// read data from tree and calculate amplitudes
	vector<complex<double> > ampValues;
	if (not processTree(*inTree, *prodKinPartNames, *decayKinPartNames,
	                    amplitude, ampValues, maxNmbEvents,
	                    prodKinMomentaLeafName, decayKinMomentaLeafName, false)) {
		printWarn << "problems reading tree" << endl;
		return false;
	}
	if (ampValues.size() < 1) {
		printWarn << "no amplitude were calculated" << endl;
		return false;
	}

	// calculate amplitudes for parity transformed decay daughters
	vector<complex<double> > ampSpaceInvValues;
	amplitude->enableSpaceInversion(true);
	if (not processTree(*inTree, *prodKinPartNames, *decayKinPartNames,
	                    amplitude, ampSpaceInvValues, maxNmbEvents,
	                    prodKinMomentaLeafName, decayKinMomentaLeafName, false)) {
		printWarn << "problems reading tree" << endl;
		return false;
	}

	// calculate amplitudes for decay daughters reflected through production plane
	vector<complex<double> > ampReflValues;
	amplitude->enableSpaceInversion(false);
	amplitude->enableReflection    (true);
	if (not processTree(*inTree, *prodKinPartNames, *decayKinPartNames,
	                    amplitude, ampReflValues, maxNmbEvents,
	                    prodKinMomentaLeafName, decayKinMomentaLeafName, false)) {
		printWarn << "problems reading tree" << endl;
		return false;
	}

	if (   (ampValues.size() != ampSpaceInvValues.size())
	    or (ampValues.size() != ampReflValues.size    ())) {
		printWarn << "different number of amplitudes for space inverted "
		          << "(" << ampSpaceInvValues.size() << "), reflected "
		          << "(" << ampReflValues.size() << "), and unmodified data "
		          << "(" << ampValues.size() << ")." << endl;
		return false;
	}

	printInfo << "checking symmetry properties of amplitudes" << endl;
	unsigned int countAmpZero               = 0;
	unsigned int countAmpRatioNotOk         = 0;
	unsigned int countSpaceInvEigenValNotOk = 0;
	unsigned int countReflEigenValNotOk     = 0;
	for (unsigned int i = 0; i < ampValues.size(); ++i) {
		// check that amplitude is non-zero
		bool ampZero = false;
		if (ampValues[i] == complex<double>(0, 0)) {
			ampZero = true;
			++countAmpZero;
		}
		const unsigned int nmbDigits = numeric_limits<double>::digits10 + 1;
		ostringstream s;
		s.precision(nmbDigits);
		s.setf(ios_base::scientific, ios_base::floatfield);
		if (debug) {
			printDebug << "amplitude [" << i << "]: " << endl;
			s << "        ampl.            = " << ampValues[i];
			if (ampZero)
				s << " <! zero amplitude";
			s << endl
			  << "        ampl. space inv. = " << ampSpaceInvValues[i] << endl
			  << "        ampl. refl.      = " << ampReflValues    [i] << endl;
			cout << s.str();
		}

		// check space inversion symmetry
		const complex<double> spaceInvRatio
			= complex<double>((ampSpaceInvValues[i].real() != 0) ?
			                  ampValues[i].real() / ampSpaceInvValues[i].real() : 0,
			                  (ampSpaceInvValues[i].imag() != 0) ?
			                  ampValues[i].imag() / ampSpaceInvValues[i].imag() : 0);
		bool spaceInvAmpRatioOk = true;
		if (   (fabs(spaceInvRatio.real()) - 1 > maxDelta)
		    or (fabs(spaceInvRatio.imag()) - 1 > maxDelta)) {
			spaceInvAmpRatioOk = false;
			++countAmpRatioNotOk;
		}
		const int spaceInvEigenValue = decayTopo->spaceInvEigenValue();
		bool      spaceInvEigenValOk = true;
		if (   (    (spaceInvRatio.real() != 0)
		        and (spaceInvRatio.real() / fabs(spaceInvRatio.real()) != spaceInvEigenValue))
		    or (    (spaceInvRatio.imag() != 0)
		        and (spaceInvRatio.imag() / fabs(spaceInvRatio.imag()) != spaceInvEigenValue))) {
			spaceInvEigenValOk = false;
			++countSpaceInvEigenValNotOk;
		}
		if (debug) {
			s.str("");
			s << "Re[ampl.] / Re[ampl. space inv.] = " << setw(23) << spaceInvRatio.real() << ", "
			  << "Im[ampl.] / Im[ampl. space inv.] = " << setw(23) << spaceInvRatio.imag();
			cout << "        " << s.str();
			if (not spaceInvAmpRatioOk)
				cout << " <! larger than " << maxDelta;
			if (not spaceInvEigenValOk)
				cout << " <! eigenvalue != " << spaceInvEigenValue;
			cout << endl;
		}

		// check reflection symmetry through production plane
		const complex<double> reflRatio
			= complex<double>((ampReflValues[i].real() != 0) ?
			                  ampValues[i].real() / ampReflValues[i].real() : 0,
			                  (ampReflValues[i].imag() != 0) ?
			                  ampValues[i].imag() / ampReflValues[i].imag() : 0);
		bool reflAmpRatioOk = true;
		if (   (fabs(reflRatio.real()) - 1 > maxDelta)
	      or (fabs(reflRatio.imag()) - 1 > maxDelta)) {
			reflAmpRatioOk = false;
			++countAmpRatioNotOk;
		}
		const int reflEigenValue = decayTopo->reflectionEigenValue();
		bool      reflEigenValOk = true;
		if (   (    (reflRatio.real() != 0)
		        and (reflRatio.real() / fabs(reflRatio.real()) != reflEigenValue))
		    or (    (reflRatio.imag() != 0)
		        and (reflRatio.imag() / fabs(reflRatio.imag()) != reflEigenValue))) {
			reflEigenValOk = false;
			++countReflEigenValNotOk;
		}
		if (debug) {
			s.str("");
			s << "Re[ampl.] / Re[ampl. refl.]      = " << setw(23) << reflRatio.real() << ", "
			  << "Im[ampl.] / Im[ampl. refl.]      = " << setw(23) << reflRatio.imag();
			cout << "        " << s.str();
			if (not reflAmpRatioOk)
				cout << " <! larger than " << maxDelta;
			if (not reflEigenValOk)
				cout << " <! eigenvalue != " << reflEigenValue;
			cout << endl;
		}
	}

	if (countAmpZero > 0) {
		keyFileErrors.push_back("zero amplitude");
		success = false;
	}
	if (countAmpRatioNotOk > 0) {
		stringstream s;
		s << "amplitude deviation in symmetry check larger than " << maxDelta;
		keyFileErrors.push_back(s.str());
		success = false;
	}
	if (countSpaceInvEigenValNotOk > 0) {
		keyFileErrors.push_back("wrong space inversion eigenvalue");
		success = false;
	}
	if (countReflEigenValNotOk > 0) {
		keyFileErrors.push_back("wrong eigenvalue for reflection through production plane");
		success = false;
	}
	return success;
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
	const string progName     = argv[0];
	string       dataFileName = "";
	long int     maxNmbEvents = -1;
	string       pdgFileName  = "./particleDataTable.txt";
	string       inTreeName   = "rootPwaEvtTree";
	double       maxDelta     = 1e-6;
	string       leafNames    = "prodKinParticles;prodKinMomenta;decayKinParticles;decayKinMomenta";
	bool         debug        = false;
	extern char* optarg;
	extern int   optind;
	int          c;
	while ((c = getopt(argc, argv, "d:n:p:t:e:l:r:vh")) != -1)
		switch (c) {
		case 'd':
			dataFileName = optarg;
			break;
		case 'n':
			maxNmbEvents = atol(optarg);
			break;
		case 'p':
			pdgFileName = optarg;
			break;
		case 't':
			inTreeName = optarg;
			break;
		case 'e':
			maxDelta = atof(optarg);
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
		}

	// set debug options
	if (debug) {
		waveDescription::setDebug(true);
		//particleProperties::setDebug(true);
		isobarHelicityAmplitude::setDebug(true);
		isobarCanonicalAmplitude::setDebug(true);
		//massDependence::setDebug(true);
	}

	// get key file names
	if (optind >= argc) {
		printErr << "you need to specify at least one key file to process. Aborting..." << endl;;
		usage(progName, 1);
	}
	vector<string> keyFileNames;
	while (optind < argc) {
		const string fileName = argv[optind++];
		keyFileNames.push_back(fileName);
	}

	// get object and leaf names for event data
	string prodKinPartNamesObjName,  prodKinMomentaLeafName;
	string decayKinPartNamesObjName, decayKinMomentaLeafName;
	parseLeafAndObjNames(leafNames, prodKinPartNamesObjName, prodKinMomentaLeafName,
	                     decayKinPartNamesObjName, decayKinMomentaLeafName);

	// open input file
	TTree*        inTree            = 0;
	TClonesArray* prodKinPartNames  = 0;
	TClonesArray* decayKinPartNames = 0;
	if (dataFileName == "")
		printInfo << "No test data file was specified (option -d). "
		          << "skipping symmetry checks of amplitude." << endl;
	else {
		vector<string> rootFileNames;
		vector<string> evtFileNames;
		const string   fileExt  = extensionFromPath(dataFileName);
		if (fileExt == "root")
			rootFileNames.push_back(dataFileName);
		else if (fileExt == "evt")
			evtFileNames.push_back(dataFileName);
		else {
			printErr << "specified input files is neither a .root nor a .evt file. Aborting...";
			usage(progName, 1);
		}
		vector<TTree*> inTrees;
		if (not openRootEvtFiles(inTrees, prodKinPartNames, decayKinPartNames,
		                         rootFileNames, evtFileNames,
		                         inTreeName, prodKinPartNamesObjName, prodKinMomentaLeafName,
		                         decayKinPartNamesObjName, decayKinMomentaLeafName, debug)) {
			printErr << "problems opening input file(s). Aborting..." << endl;
			exit(1);
		}
		inTree = inTrees[0];
	}

	// initialize particle data table
	particleDataTable::readFile(pdgFileName);

	// loop over key files
	map<string, vector<string> > keyFileErrors;  // maps error description to key files
	unsigned int                 countKeyFileErr = 0;
	for (unsigned int i = 0; i < keyFileNames.size(); ++i) {
		cout << endl;
		printInfo << "checking key file '" << keyFileNames[i] << "'" << endl;
		vector<string> errors;
		if (not testAmplitude(inTree, keyFileNames[i], errors, maxNmbEvents, debug, maxDelta,
		                      prodKinPartNames,  decayKinPartNames,
		                      prodKinMomentaLeafName, decayKinMomentaLeafName)) {
			++countKeyFileErr;
			printWarn << "key file '" << keyFileNames[i] << "' did not pass all tests. "
			          << "see summary below." << endl;
		} else
			printSucc << "key file '" << keyFileNames[i] << "' passed all tests" << endl;
		// collect errors
		for (unsigned int j = 0; j < errors.size(); ++j)
			keyFileErrors[errors[j]].push_back(keyFileNames[i]);
	}

	// clean up
	if (inTree)
		delete inTree;

	// report final result and set exit status accordingly
	cout << endl;
	if (countKeyFileErr == 0) {
		printSucc << "all " << keyFileNames.size() << " keyfile(s) passed all tests" << endl;
		return 0;
	}

	printInfo << keyFileNames.size() - countKeyFileErr << " of " << keyFileNames.size()
	          << " keyfile(s) passed all tests" << endl;
	if (keyFileErrors.size() > 0) {
		printInfo << countKeyFileErr << " problematic keyfile(s):" << endl;
		for (map<string, vector<string> >::const_iterator entry = keyFileErrors.begin();
		     entry != keyFileErrors.end(); ++entry) {
			cout << "        " << entry->first << ":" << endl;
			for (unsigned int i = 0; i < entry->second.size(); ++i)
				cout << "            " << entry->second[i] << endl;
		}
	}
	exit(1);

}
