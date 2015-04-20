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
//      helper functions that convert between standard ASCII PWA2000
//      .evt files and the new ROOT tree format
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <algorithm>
#include <map>

#include <boost/tokenizer.hpp>
#include <boost/progress.hpp>
#include <boost/bimap.hpp>
#include <boost/assign/list_inserter.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TTreePerfStats.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TObjString.h"
#include "TVector3.h"

#include "reportingUtilsRoot.hpp"
#include "conversionUtils.hpp"
#include "particleDataTable.h"
#include "isobarDecayTopology.h"
#include "isobarHelicityAmplitude.h"
#include "evtTreeHelper.h"


using namespace std;
using namespace boost;
using namespace boost::bimaps;


namespace rpwa {


	bool
	checkParticleCharge(const size_t  lineNmb,
	                    const int     id,
	                    const string& name,
	                    const int     chargeToCheck)
	{
		if (name == "unknown") {
			printWarn << "error reading data line " << lineNmb << ": unknown particle" << endl;
			return false;
		} else {
			// check charge
			int charge;
			particleProperties::chargeFromName(name, charge);
			if (chargeToCheck != charge) {
				printWarn << "error reading data line " << lineNmb << ": "
				          << "GEANT particle ID of " << id << " corresponds to particle "
				          << "'" << name << "'. this inconsistent with the charge of "
				          << "'" << chargeToCheck << "' in data file.";
				return false;
			}
		}
		return true;
	}


	double
	getParticleMass(const string& name)
	{
		rpwa::particleDataTable&  pdt  = rpwa::particleDataTable::instance();
		const particleProperties* prop = 0;
		if (pdt.isInTable(name))
			prop = pdt.entry(name);
		else {
			const string n = particleProperties::stripChargeFromName(name);
			if (pdt.isInTable(n))
				prop = pdt.entry(n);
		}
		if (not prop) {
			printWarn << "neither particle '" << name << "' "
			          << "nor '" << particleProperties::stripChargeFromName(name) << "' "
			          << "are in particle data table. using mass 0." << endl;
			return 0;
		}
		return prop->mass();
	}


	void
	parseLeafAndObjNames(const string& cmdLineString,
	                     string&       prodKinPartNamesObjName,
	                     string&       prodKinMomentaLeafName,
	                     string&       decayKinPartNamesObjName,
	                     string&       decayKinMomentaLeafName)
	{
		typedef tokenizer<char_separator<char> > tokenizer;
		char_separator<char> separator(";");
		tokenizer            nameTokens(cmdLineString, separator);
		tokenizer::iterator  nameToken = nameTokens.begin();
		prodKinPartNamesObjName  = *nameToken;
		prodKinMomentaLeafName   = *(++nameToken);
		decayKinPartNamesObjName = *(++nameToken);
		decayKinMomentaLeafName  = *(++nameToken);
		printInfo << "using the following object/leaf names:" << endl
		          << "        production kinematics: "
		          << "particle names = '" << prodKinPartNamesObjName << "', "
		          << "momenta = '" << prodKinMomentaLeafName << "'" << endl
		          << "        decay kinematics:      "
		          << "particle names = '" << decayKinPartNamesObjName << "', "
		          << "momenta = '" << decayKinMomentaLeafName << "'" << endl;
	}


	bool
	getParticleNamesFromRootFile(TFile&         inFile,
	                             TClonesArray*& prodKinPartNames,   // array of particle names to be filled
	                             TClonesArray*& decayKinPartNames,  // array of particle names to be filled
	                             const string&  inTreeName,
	                             const string&  prodKinPartNamesObjName,
	                             const string&  decayKinPartNamesObjName)
	{
		// read particle names from first root file
		inFile.GetObject(prodKinPartNamesObjName.c_str(),  prodKinPartNames );
		inFile.GetObject(decayKinPartNamesObjName.c_str(), decayKinPartNames);
		if (prodKinPartNames and decayKinPartNames) {
			TVector3::Class()->IgnoreTObjectStreamer(true);
			return true;
		}
		// if particle names are not in file try first tree entry for backward compatibility
		TTree* inTree = 0;
		inFile.GetObject(inTreeName.c_str(), inTree);
		if (not prodKinPartNames) {
			TBranch* prodKinPartNamesBr = 0;
			if (inTree->SetBranchAddress(prodKinPartNamesObjName.c_str(),
			                             &prodKinPartNames, &prodKinPartNamesBr) >= 0)
				if (inTree->GetEntries() > 0)
					if (inTree->LoadTree(0) >= 0)
						prodKinPartNamesBr->GetEntry(0);
		}
		if (not decayKinPartNames) {
			TBranch* decayKinPartNamesBr = 0;
			if (inTree->SetBranchAddress(decayKinPartNamesObjName.c_str(),
			                             &decayKinPartNames, &decayKinPartNamesBr) >= 0)
				if (inTree->GetEntries() > 0)
					if (inTree->LoadTree(0) >= 0)
						decayKinPartNamesBr->GetEntry(0);
		}
		if (prodKinPartNames and decayKinPartNames)
			return true;
		else
			return false;
	}


	bool
	openRootEvtFiles(vector<TTree*>&       inTrees,            // array of trees from .root and .evt files
	                 TClonesArray*&        prodKinPartNames,   // array of particle names to be filled
	                 TClonesArray*&        decayKinPartNames,  // array of particle names to be filled
	                 const vector<string>& rootFileNames,      // .root files to be opened
	                 const vector<string>& evtFileNames,       // .evt files to be converted to trees
	                 const string&         inTreeName,
	                 const string&         prodKinPartNamesObjName,
	                 const string&         prodKinMomentaLeafName,
	                 const string&         decayKinPartNamesObjName,
	                 const string&         decayKinMomentaLeafName,
	                 const bool            debug)
	{
		// open root files and build chain
		TChain* inChain = 0;
		if (rootFileNames.size() > 0) {
			inChain = new TChain(inTreeName.c_str());
			for (unsigned int i = 0; i < rootFileNames.size(); ++i) {
				printInfo << "opening ROOT input file '" << rootFileNames[i] << "'" << endl;
				if (inChain->Add(rootFileNames[i].c_str()) < 1)
					printWarn << "no events in ROOT input file '" << rootFileNames[i] << "'" << endl;
			}
			// read particle names from first root file
			TFile* inFile = inChain->GetFile();  // opens first file
			if (not inFile)
				printWarn << "could not open file '" << rootFileNames[0] << "'" << endl;
			else {
				if (not getParticleNamesFromRootFile(*inFile, prodKinPartNames, decayKinPartNames,
				                                     inTreeName,
				                                     prodKinPartNamesObjName, decayKinPartNamesObjName)) {
					printErr << "cannot find production and/or decay kinematics particle names "
					         << "in input file '" << inChain->GetFile()->GetName() << "'" << endl;
					return false;
				}
			}
		}
		if (inChain)
			inTrees.push_back(inChain);

		// convert .evt files to root trees
		for (unsigned int i = 0; i < evtFileNames.size(); ++i) {
			printInfo << "opening .evt input file '" << evtFileNames[i] << "'" << endl;
			ifstream evtFile(evtFileNames[i].c_str());
			if (not evtFile or not evtFile.good()) {
				printWarn << "cannot open .evt input file '" << evtFileNames[i] << "'. skipping." << endl;
				continue;
			}
			printInfo << "converting .evt input file '" << evtFileNames[i] << "' "
			          << "into memory resident tree. this might reduce performance. "
			          << "ROOT input format is recommended." << endl;
			// create tree
			TTree* tree = new TTree(inTreeName.c_str(), inTreeName.c_str());
			if (not tree) {
				printErr << "problems creating tree '" << inTreeName << "'. skipping." << endl;
				continue;
			}
			TClonesArray prodNames ("TObjString");  // names of production kinematics particles read from .evt file
			TClonesArray decayNames("TObjString");  // names of decay kinematics particles read from .evt file
			if (fillTreeFromEvt(evtFile, *tree, prodNames, decayNames, -1,
			                    prodKinMomentaLeafName, decayKinMomentaLeafName, debug))
				inTrees.push_back(tree);
			else {
				printWarn << "problems creating tree from .evt input file '" << evtFileNames[i] << "' "
				          << "skipping." << endl;
			}

			// set particle names
			if (prodKinPartNames) {
				// check consistency
				int j;
				if (prodNames.GetEntriesFast() != prodKinPartNames->GetEntriesFast())
					goto prodKinPartNamesInconsistent;
				for (j = 0; j < prodNames.GetEntriesFast(); ++j) {
					if (not (((TObjString*)prodNames[j])->GetString()
					         == ((TObjString*)(*prodKinPartNames)[j])->GetString()))
						goto prodKinPartNamesInconsistent;
				prodKinPartNamesInconsistent:
					{
						printErr << "production kinematics particle names read from input files "
						         << "are inconsistent:" << endl
						         << "          expected:";
						for (int k = 0; k < prodKinPartNames->GetEntriesFast(); ++k)
							cout << "  " << ((TObjString*)(*prodKinPartNames)[k])->GetString().Data()
							     << "[" << k << "]";
						cout << endl << "          read from '" << evtFileNames[i] << "':";
						for (int k = 0; k < prodNames.GetEntriesFast(); ++k)
							cout << "  " << ((TObjString*)prodNames[k])->GetString().Data() << "[" << k << "]";
						cout << endl << "    check input files." << endl;
						return false;
					}
				}
			} else {
				prodKinPartNames  = new TClonesArray("TObjString");
				*prodKinPartNames = prodNames;
			}
			if (decayKinPartNames) {
				// check consistency
				int j;
				if (decayNames.GetEntriesFast() != decayKinPartNames->GetEntriesFast())
					goto decayKinPartNamesInconsistent;
				for (j = 0; j < decayNames.GetEntriesFast(); ++j) {
					if (not (((TObjString*)decayNames[j])->GetString()
					         == (((TObjString*)(*decayKinPartNames)[j])->GetString())))
						goto decayKinPartNamesInconsistent;
				decayKinPartNamesInconsistent:
					{
						printErr << "decay kinematics particle names read from input files "
						         << "are inconsistent:" << endl
						         << "          expected:";
						for (int k = 0; k < decayKinPartNames->GetEntriesFast(); ++k)
							cout << "  " << ((TObjString*)(*decayKinPartNames)[k])->GetString().Data()
							     << "[" << k << "]";
						cout << endl << "          read from '" << evtFileNames[i] << "':";
						for (int k = 0; k < decayNames.GetEntriesFast(); ++k)
							cout << "  " << ((TObjString*)decayNames[k])->GetString().Data() << "[" << k << "]";
						cout << endl << "    check input files." << endl;
						return false;
					}
				}
			} else {
				decayKinPartNames  = new TClonesArray("TObjString");
				*decayKinPartNames = decayNames;
			}
		}  // loop over .evt files

		// make sure particle names are specified
		bool success = true;
		if (not prodKinPartNames) {
			printErr << "no production kinematics particle names were found in input files." << endl;
			success = false;
		}
		if (not decayKinPartNames) {
			printWarn << "no decay kinematics particle names were found in input files." << endl;
			success = false;
		}

		return success;
	}


	bool
	fillTreeFromEvt(istream&       inEvt,
	                TTree&         outTree,            // tree to be filled
	                TClonesArray&  prodKinPartNames,   // array of particle names to be filled
	                TClonesArray&  decayKinPartNames,  // array of particle names to be filled
	                const long int maxNmbEvents,
	                const string&  prodKinMomentaLeafName,
	                const string&  decayKinMomentaLeafName,
	                const bool     debug,
	                const long int treeCacheSize)
	{
		if (not inEvt or not inEvt.good()) {
			printWarn << "cannot read from input stream" << endl;
			return false;
		}

		// create leaf variables
		TVector3::Class()->IgnoreTObjectStreamer(true);
		TClonesArray* prodKinMomenta  = new TClonesArray("TVector3");
		TClonesArray* decayKinMomenta = new TClonesArray("TVector3");

		// connect leaf variables to tree branches
		const int splitLevel = 99;
		const int bufSize    = 256000;
		outTree.Branch(prodKinMomentaLeafName.c_str(),  "TClonesArray", &prodKinMomenta,  bufSize, splitLevel);
		outTree.Branch(decayKinMomentaLeafName.c_str(), "TClonesArray", &decayKinMomenta, bufSize, splitLevel);

		// loop over events and fill tree
		vector<string> prodNames;   // names of production kinematics particles read from .evt file
		vector<string> decayNames;  // names of decay kinematics particles read from .evt file
		bool           success     = true;
		long int       countEvents = 0;
		long int       countLines  = 0;
		inEvt.seekg(0, ios::end);
		long int fileLength = inEvt.tellg();
		inEvt.seekg(0, ios::beg);
		progress_display* progressIndicator = (not debug) ? new progress_display(fileLength, cout, "") : 0;
		streampos         lastPos           = inEvt.tellg();
		while (inEvt.good()) {
			string line;

			// read number of particles
			int nmbParticles = 0;
			if (getline(inEvt, line)) {
				++countLines;
				stringstream lineStream(line);
				int          n;
				if (lineStream >> n)
					nmbParticles = n;
				else {
					printWarn << "event " << countEvents + 1 << ": error reading number of particles "
					          << "from line " << countLines << ": " << line << endl;
					success = false;
				}
			} else
				break;
			assert(nmbParticles > 0);
			if (debug)
				printDebug << "# of particles = " << nmbParticles << endl;

			// read production kinematics data (beam + fixed target)
			prodNames.clear();
			prodKinMomenta->Clear();
			if (getline(inEvt, line)) {
				++countLines;
				stringstream lineStream(line);
				int          id = 0, charge = 0;
				double       momX = 0, momY = 0, momZ = 0, E = 0;
				if (lineStream >> id >> charge >> momX >> momY >> momZ >> E) {
					const string partName = particleDataTable::particleNameFromGeantId(id);
					prodNames.push_back(partName);
					if (not checkParticleCharge(countLines, id, partName, charge))
						success = false;
					new((*prodKinMomenta)[0]) TVector3(momX, momY, momZ);
				} else {
					printWarn << "event " << countEvents + 1 << ": error reading beam data "
					          << "from line " << countLines << ": " << line << endl;
					success = false;
				}
			} else
				break;

			// check consistency
			const int nmbProdKinPart = prodNames.size();
			assert((nmbProdKinPart > 0) and (nmbProdKinPart == prodKinMomenta->GetEntriesFast()));
			if (debug) {
				printDebug << nmbProdKinPart << " production kinematics particles:" << endl;
				for (int i = 0; i < nmbProdKinPart; ++i)
					cout << "        particle[" << i << "]: "
					     << prodNames[i] << "; " << *((TVector3*)(*prodKinMomenta)[i]) << endl;
			}

			// read decay kinematics data
			decayNames.clear();
			decayNames.resize(nmbParticles - 1, "");
			decayKinMomenta->Clear();
			for (int i = 0; i < nmbParticles - 1; ++i) {
				if (getline(inEvt, line)) {
					++countLines;
					stringstream lineStream(line);
					int          id, charge;
					double       momX, momY, momZ, E;
					if (lineStream >> id >> charge >> momX >> momY >> momZ >> E) {
						const string partName = particleDataTable::particleNameFromGeantId(id);
						decayNames[i] = partName;
						if (not checkParticleCharge(countLines, id, partName, charge))
							success = false;
						new((*decayKinMomenta)[i]) TVector3(momX, momY, momZ);
					} else {
						printWarn << "event " << countEvents + 1 << ": error reading decay kinematics "
						          << "particle[" << i << "] data from line " << countLines << ": " << line << endl;
						success = false;
					}
				} else
					break;
			}

			// check consistency
			const int nmbDecayKinPart = decayNames.size();
			assert((nmbDecayKinPart > 0) and (nmbDecayKinPart == decayKinMomenta->GetEntriesFast()));
			if (debug) {
				printDebug << nmbDecayKinPart << " decay kinematics particles:" << endl;
				for (int i = 0; i < nmbDecayKinPart; ++i)
					cout << "        particle[" << i << "]: "
					     << decayNames[i] << "; " << *((TVector3*) (*decayKinMomenta)[i]) << endl;
			}

			outTree.Fill();
			if (countEvents == 0) {
				for (unsigned int i = 0; i < prodNames.size(); ++i)
					new(prodKinPartNames[i]) TObjString(prodNames[i].c_str());
				for (unsigned int i = 0; i < decayNames.size(); ++i)
					new(decayKinPartNames[i]) TObjString(decayNames[i].c_str());
			} else {
				// check that particle lists do not change
				const int nmbIsPart         = prodNames.size();
				const int nmbIsPartExpected = prodKinPartNames.GetEntriesFast();
				const int nmbFsPart         = decayNames.size();
				const int nmbFsPartExpected = decayKinPartNames.GetEntriesFast ();
				if ((nmbIsPart != nmbIsPartExpected) || (nmbFsPart != nmbFsPartExpected)) {
					printErr << "number of particles changed in event " << countEvents << ". "
					         << "producion kinematics: " << nmbIsPart << " particles "
					         << "(" << nmbIsPartExpected << " expected); "
					         << "decay kinematics: " << nmbFsPart << " particles "
					         << "(" << nmbFsPartExpected << " expected). "
					         << "check .evt file. Aborting..." << endl;
					throw;
				}
				for (int i = 0; i < nmbIsPart; ++i) {
					const string name         = prodNames[i];
					const string nameExpected = ((TObjString*)prodKinPartNames[i])->GetString().Data();
					if (name != nameExpected) {
						printErr << "particle[" << i << "] name mismatch in production kinematics in event "
						         << countEvents << ": '" << name << "', expected '" << nameExpected << "'. "
						         << "check .evt file. Aborting..." << endl;
						throw;
					}
				}
				for (int i = 0; i < nmbFsPart; ++i) {
					const string name         = decayNames[i];
					const string nameExpected = ((TObjString*)decayKinPartNames[i])->GetString().Data();
					if (name != nameExpected) {
						printErr << "particle[" << i << "] name mismatch in decay kinematics in event "
						         << countEvents << ": '" << name << "', expected '" << nameExpected << "'. "
						         << "check .evt file. Aborting..." << endl;
						throw;
					}
				}
			}

			++countEvents;
			if (progressIndicator)
				(*progressIndicator) += inEvt.tellg() - lastPos;
			lastPos = inEvt.tellg();
			if ((maxNmbEvents > 0) and (countEvents >= maxNmbEvents))
				break;
		}

		printInfo << "optimizing tree" << endl;
		//outTree.Print();
		outTree.OptimizeBaskets(treeCacheSize, 1, "d");
		//outTree.Print();

		printInfo << "read " << countLines << " lines from input stream and wrote "
		          << countEvents << " events to tree '" << outTree.GetName() << "' "
		          << "assuming fixed target" << endl;
		return success;
	}


	bool
	writeEvtFromTree(TTree&              inTree,
	                 ostream&            outEvt,
	                 const TClonesArray& prodKinPartNames,
	                 const TClonesArray& decayKinPartNames,
	                 const long int      maxNmbEvents,
	                 const string&       inTreeName,
	                 const string&       prodKinMomentaLeafName,
	                 const string&       decayKinMomentaLeafName,
	                 const bool          debug)
	{
		const long int nmbEventsTree = inTree.GetEntries();
		if (not outEvt) {
			printWarn << "cannot write to output stream" << endl;
			return false;
		}
		string partClassName = prodKinPartNames.GetClass()->GetName();
		if (partClassName != "TObjString") {
			printWarn << "production kinematics particle names are of type '" << partClassName
			          << "' and not TObjString." << endl;
			return false;
		}
		partClassName = decayKinPartNames.GetClass()->GetName();
		if (partClassName != "TObjString") {
			printWarn << "decay kinematics particle names are of type '" << partClassName
			          << "' and not TObjString." << endl;
			return false;
		}

		// create leaf variables
		TClonesArray* prodKinMomenta  = 0;
		TClonesArray* decayKinMomenta = 0;

		// connect leaf variables to tree branches
		inTree.SetBranchAddress(prodKinMomentaLeafName.c_str(),  &prodKinMomenta );
		inTree.SetBranchAddress(decayKinMomentaLeafName.c_str(), &decayKinMomenta);

		// loop over events
		const long int    nmbEvents         = ((maxNmbEvents > 0) ? min(maxNmbEvents, nmbEventsTree)
		                                       : nmbEventsTree);
		progress_display* progressIndicator = (not debug) ? new progress_display(nmbEvents, cout, "") : 0;
		for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
			if (progressIndicator)
				++(*progressIndicator);

			if (inTree.LoadTree(eventIndex) < 0)
				break;
			inTree.GetEntry(eventIndex);

			assert(prodKinMomenta );
			assert(decayKinMomenta);
			const int nmbProdKinPart  = prodKinPartNames.GetEntriesFast ();
			const int nmbDecayKinPart = decayKinPartNames.GetEntriesFast();
			assert(nmbProdKinPart  == prodKinMomenta->GetEntriesFast ());
			assert(nmbDecayKinPart == decayKinMomenta->GetEntriesFast());
			if (nmbProdKinPart < 1) {
				printWarn << "array of production kinematics particles does not have any entries. "
				          << "at least entry for beam (index 0) is required. skipping event." << endl;
				continue;
			}
			if (nmbDecayKinPart < 1) {
				printWarn << "array of decay kinematics particles does not have any entries. "
				          << "at least one final state particle is required. skipping event." << endl;
				continue;
			}

			// write total number of particles
			outEvt << 1 + nmbDecayKinPart << endl;  // PWA2000 supports only beam in production kinematics

			// write production kinematics (beam only)
			if (debug)
				printDebug << "event[" << eventIndex << "]: " << nmbProdKinPart
				           << " production kinematics particles:" << endl;
			{  // only beam
				assert(prodKinPartNames [0]);
				assert((*prodKinMomenta)[0]);
				const string    name = ((TObjString*)prodKinPartNames[0])->GetString().Data();
				const TVector3* mom  = dynamic_cast<TVector3*>((*prodKinMomenta)[0]);
				assert(mom);
				const double mass = getParticleMass(name);
				int id, charge;
				particleDataTable::geantIdAndChargeFromParticleName(name, id, charge);
				outEvt << setprecision(numeric_limits<double>::digits10 + 1)
				       << id << " " << charge << " " << mom->X() << " " << mom->Y() << " " << mom->Z() << " "
				       << sqrt(mass * mass + mom->Mag2()) << endl;
				if (debug) {
					cout << "        particle[" << 0 << "]: " << name << ", id = " << id << ", "
					     << "charge = " << charge << "; " << *mom << endl;
				}
			}

			// write decay kinematics
			if (debug)
				printDebug << "event[" << eventIndex << "]: " << nmbDecayKinPart
				           << " decay kinematics particles:" << endl;
			for (unsigned int i = 0; i < (unsigned int)nmbDecayKinPart; ++i) {
				assert(decayKinPartNames [i]);
				assert((*decayKinMomenta)[i]);
				const string    name = ((TObjString*)decayKinPartNames[i])->GetString().Data();
				const TVector3* mom  = dynamic_cast<TVector3*>((*decayKinMomenta)[i]);
				assert(mom);
				const double mass = getParticleMass(name);
				int id, charge;
				particleDataTable::geantIdAndChargeFromParticleName(name, id, charge);
				outEvt << setprecision(numeric_limits<double>::digits10 + 1)
				       << id << " " << charge << " " << mom->X() << " " << mom->Y() << " " << mom->Z() << " "
				       << sqrt(mass * mass + mom->Mag2()) << endl;
				if (debug) {
					cout << "        particle[" << i << "]: " << name << ", id = " << id << ", "
					     << "charge = " << charge << "; " << *mom << endl;
				}
			}

			if (debug)
				cout << endl;
		}

		printInfo << "wrote " << nmbEvents << " events to output stream" << endl;
		return true;
	}


	bool
	processTree(TTree&                    tree,
	            const TClonesArray&       prodKinPartNames,
	            const TClonesArray&       decayKinPartNames,
	            const isobarAmplitudePtr& amplitude,
	            vector<complex<double> >& ampValues,
	            const long int            maxNmbEvents,
	            const string&             prodKinMomentaLeafName,
	            const string&             decayKinMomentaLeafName,
	            const bool                printProgress,
	            const string&             treePerfStatOutFileName,
	            const long int            treeCacheSize)
	{
		if (not amplitude) {
			printWarn << "null pointer to isobar decay amplitude. cannot process tree." << endl;
			return false;
		}
		// initialize amplitude
		amplitude->init();
		const isobarDecayTopologyPtr& decayTopo = amplitude->decayTopology();

		// create branch pointers and leaf variables
		TBranch*      prodKinMomentaBr  = 0;
		TBranch*      decayKinMomentaBr = 0;
		TClonesArray* prodKinMomenta    = 0;
		TClonesArray* decayKinMomenta   = 0;

		// connect leaf variables to tree branches
		tree.SetBranchAddress(prodKinMomentaLeafName.c_str(),  &prodKinMomenta,  &prodKinMomentaBr );
		tree.SetBranchAddress(decayKinMomentaLeafName.c_str(), &decayKinMomenta, &decayKinMomentaBr);
		tree.SetCacheSize(treeCacheSize);
		tree.AddBranchToCache(prodKinMomentaLeafName.c_str(),  true);
		tree.AddBranchToCache(decayKinMomentaLeafName.c_str(), true);
		tree.StopCacheLearningPhase();
		TTreePerfStats* treePerfStats = 0;
		if (treePerfStatOutFileName != "")
			treePerfStats = new TTreePerfStats("ioPerf", &tree);

		// loop over events
		if (not decayTopo->initKinematicsData(prodKinPartNames, decayKinPartNames)) {
			printWarn << "problems initializing input data. cannot read input data." << endl;
			return false;
		}
		const long int    nmbEventsTree     = tree.GetEntries();
		const long int    nmbEvents         = ((maxNmbEvents > 0) ? min(maxNmbEvents, nmbEventsTree)
		                                       : nmbEventsTree);
		bool              success           = true;
		progress_display* progressIndicator = (printProgress) ? new progress_display(nmbEvents, cout, "") : 0;
		for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
			if (progressIndicator)
				++(*progressIndicator);

			if (tree.LoadTree(eventIndex) < 0)
				break;
			// read only required branches
			prodKinMomentaBr->GetEntry (eventIndex);
			decayKinMomentaBr->GetEntry(eventIndex);

			if (not prodKinMomenta or not decayKinMomenta) {
				printWarn << "at least one of the input data arrays is a null pointer: "
				          << "        production kinematics: " << "momenta = " << prodKinMomenta  << endl
				          << "        decay kinematics:      " << "momenta = " << decayKinMomenta << endl
				          << "skipping event." << endl;
				success = false;
				continue;
			}

			if (decayTopo->readKinematicsData(*prodKinMomenta, *decayKinMomenta))
				ampValues.push_back((*amplitude)());
			else {
				printWarn << "problems reading event[" << eventIndex << "]" << endl;
				success = false;
			}
		}

		if (printProgress)
			tree.PrintCacheStats();
		if (treePerfStats) {
			treePerfStats->SaveAs(treePerfStatOutFileName.c_str());
			delete treePerfStats;
		}
		return success;
	}


}  // namespace rpwa
