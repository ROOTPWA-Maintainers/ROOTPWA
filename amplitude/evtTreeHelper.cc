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

#include <boost/progress.hpp>

#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TObjString.h"
#include "TVector3.h"

#include "pputil.h"
#include "reportingUtilsRoot.hpp"
#include "conversionUtils.hpp"
#include "particleDataTable.h"
#include "isobarDecayTopology.h"
#include "isobarHelicityAmplitude.h"


using namespace std;
using namespace boost;


namespace rpwa {


	string
	particleNameFromGeantId(const int id,
	                        const int charge)
	{
		assert((charge == -1) or (charge == 0) or (charge == +1));
		string name = id2name((Geant_ID)id);
		name = particleProperties::stripChargeFromName(name);
		stringstream n;
		n << name << sign(charge);
		return n.str();
	}


	void
	idAndChargeFromParticleName(const string& name,
	                            int&          id,
	                            int&          charge)
	{
		particleProperties::chargeFromName(name, charge);
		id = name2id(name, charge);
		if (id == g_Unknown)
			id = name2id(particleProperties::stripChargeFromName(name), charge);
		if (id == g_Unknown)
			printWarn << "unknown particle '" << name << "'" << endl;
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


	bool
	fillTreeFromEvt(istream&       inEvt,
	                TTree&         outTree,
	                const long int maxNmbEvents,
	                const string&  prodKinParticlesLeafName,
	                const string&  prodKinMomentaLeafName,
	                const string&  decayKinParticlesLeafName,
	                const string&  decayKinMomentaLeafName,
	                const string&  targetParticleName,
	                const bool     debug)
	{
		if (not inEvt or not inEvt.good()) {
			printWarn << "cannot read from input stream" << endl;
			return false;
		}

		// create leaf variables
		TClonesArray* prodKinParticles  = new TClonesArray("TObjString");
		TClonesArray* prodKinMomenta    = new TClonesArray("TVector3");
		TClonesArray* decayKinParticles = new TClonesArray("TObjString");
		TClonesArray* decayKinMomenta   = new TClonesArray("TVector3");

		// connect leaf variables to tree branches
		const int split   = 0;
		const int bufSize = 256000;
		outTree.Branch(prodKinParticlesLeafName.c_str(),  "TClonesArray", &prodKinParticles,  bufSize, split);
		outTree.Branch(prodKinMomentaLeafName.c_str(),    "TClonesArray", &prodKinMomenta,    bufSize, split);
		outTree.Branch(decayKinParticlesLeafName.c_str(), "TClonesArray", &decayKinParticles, bufSize, split);
		outTree.Branch(decayKinMomentaLeafName.c_str(),   "TClonesArray", &decayKinMomenta,   bufSize, split);

		// loop over events and fill tree
		bool     success     = true;
		long int countEvents = 0;
		long int countLines  = 0;
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
				int n;
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
				printInfo << "# of particles = " << nmbParticles << endl;
    
			// read production kinematics data (beam + fixed target)
			prodKinParticles->Clear();
			prodKinMomenta->Clear  ();
			if (getline(inEvt, line)) {
				++countLines;
				stringstream lineStream(line);
				int    id = 0, charge = 0;
				double momX = 0, momY = 0, momZ = 0, E = 0;
				if (lineStream >> id >> charge >> momX >> momY >> momZ >> E) {
					const string beamParticleName = particleNameFromGeantId(id , charge);
					new((*prodKinParticles)[0]) TObjString(beamParticleName.c_str());
					new((*prodKinMomenta  )[0]) TVector3  (momX, momY, momZ);
				} else {
					printWarn << "event " << countEvents + 1 << ": error reading beam data "
					          << "from line " << countLines << ": " << line << endl;
					success = false;
				}
			} else
				break;
			new((*prodKinParticles)[1]) TObjString(targetParticleName.c_str());
			new((*prodKinMomenta  )[1]) TVector3  (0, 0, 0);
			const int nmbProdKinPart = prodKinParticles->GetEntriesFast();
			assert((nmbProdKinPart > 0) and (nmbProdKinPart == prodKinMomenta->GetEntriesFast()));
			if (debug) {
				printInfo << nmbProdKinPart << " production kinematics particles:" << endl;
				for (int i = 0; i < nmbProdKinPart; ++i)
					cout << "        particle[" << i << "]: "
					     << ((TObjString*)(*prodKinParticles)[i])->GetString() << "; "
					     << *((TVector3*)(*prodKinMomenta)[i]) << endl;
			}

			// read decay kinematics data
			decayKinParticles->Clear();
			decayKinMomenta->Clear  ();
			for (int i = 0; i < nmbParticles - 1; ++i) {
				if (getline(inEvt, line)) {
					++countLines;
					stringstream lineStream(line);
					int    id, charge;
					double momX, momY, momZ, E;
					if (lineStream >> id >> charge >> momX >> momY >> momZ >> E) {
						const string name = particleNameFromGeantId(id , charge);
						new((*decayKinParticles  )[i]) TObjString(name.c_str());
						new((*decayKinMomenta)[i]) TVector3  (momX, momY, momZ);
					} else {
						printWarn << "event " << countEvents + 1 << ": error reading decay kinematics "
						          << "particle[" << i << "] data from line " << countLines << ": " << line << endl;
						success = false;
					}
				} else
					break;
			}
			const int nmbDecayKinPart = decayKinParticles->GetEntriesFast();
			assert((nmbDecayKinPart > 0) and (nmbDecayKinPart == decayKinMomenta->GetEntriesFast()));
			if (debug) {
				printInfo << nmbDecayKinPart << " decay kinematics particles:" << endl;
				for (int i = 0; i < nmbDecayKinPart; ++i)
					cout << "        particle[" << i << "]: "
					     << ((TObjString*)(*decayKinParticles)[i])->GetString() << "; "
					     << *((TVector3*) (*decayKinMomenta)  [i]) << endl;
			}

			outTree.Fill();
			++countEvents;
			if (progressIndicator)
				(*progressIndicator) += inEvt.tellg() - lastPos;
			lastPos = inEvt.tellg();
			if ((maxNmbEvents > 0) and (countEvents >= maxNmbEvents))
				break;
		}

		printInfo << "read " << countLines << " lines from input stream and wrote "
		          << countEvents << " events to tree '" << outTree.GetName() << "' "
		          << "assuming fixed " << targetParticleName << " target" << endl;
		return success;
	}


	bool
	writeEvtFromTree(TChain&        inTree,
	                 ostream&       outEvt,
	                 const long int maxNmbEvents,
	                 const string&  inTreeName,
	                 const string&  prodKinParticlesLeafName,
	                 const string&  prodKinMomentaLeafName,
	                 const string&  decayKinParticlesLeafName,
	                 const string&  decayKinMomentaLeafName,
	                 const bool     debug)
	{
		const long int nmbEventsTree = inTree.GetEntries();
		if (not outEvt) {
			printWarn << "cannot write to output stream" << endl;
			return false;
		}

		// create leaf variables
		TClonesArray* prodKinParticles  = 0;
		TClonesArray* prodKinMomenta    = 0;
		TClonesArray* decayKinParticles = 0;
		TClonesArray* decayKinMomenta   = 0;

		// connect leaf variables to tree branches
		inTree.SetBranchAddress(prodKinParticlesLeafName.c_str(),  &prodKinParticles );
		inTree.SetBranchAddress(prodKinMomentaLeafName.c_str(),    &prodKinMomenta   );
		inTree.SetBranchAddress(decayKinParticlesLeafName.c_str(), &decayKinParticles);
		inTree.SetBranchAddress(decayKinMomentaLeafName.c_str(),   &decayKinMomenta  );
			 
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

			assert(prodKinParticles );
			assert(prodKinMomenta   );
			assert(decayKinParticles);
			assert(decayKinMomenta  );
			const int nmbProdKinPart  = prodKinParticles->GetEntriesFast ();
			const int nmbDecayKinPart = decayKinParticles->GetEntriesFast();
			assert(nmbProdKinPart  == prodKinMomenta->GetEntriesFast ());
			assert(nmbDecayKinPart == decayKinMomenta->GetEntriesFast());
			if (nmbProdKinPart < 1) {
				printWarn << "arrays for production kinematics particles do not have any entries. "
				          << "at least entry for beam (index 0) is required. skipping event." << endl;
				continue;
			}

			// write total number of particles
			outEvt << 1 + nmbDecayKinPart << endl;  // PWA2000 supports only beam in production kinematics

			// write production kinematics (beam only)
			if (debug)
				printInfo << "event[" << eventIndex << "]: " << nmbProdKinPart
				          << " production kinematics particles:" << endl;
			for (int i = 0; i < 1; ++i) {  // only beam
				assert((*prodKinParticles)[i]);
				assert((*prodKinMomenta  )[i]);
				const string   name = ((TObjString*)(*prodKinParticles)[i])->GetString().Data();
				const TVector3 mom  = *((TVector3*)(*prodKinMomenta)[i]);
				const double   mass = getParticleMass(name);
				int id, charge;
				idAndChargeFromParticleName(name, id, charge);
				outEvt << setprecision(numeric_limits<double>::digits10 + 1)
				       << id << " " << charge << " " << mom.X() << " " << mom.Y() << " " << mom.Z() << " "
				       << sqrt(mass * mass + mom.Mag2()) << endl;
				if (debug) {
					cout << "        particle[" << i << "]: " << name << ", id = " << id << ", "
					     << "charge = " << charge << "; " << mom << endl;
				}
			}

			// write decay kinematics
			if (debug)
				printInfo << "event[" << eventIndex << "]: " << nmbDecayKinPart
				          << " decay kinematics particles:" << endl;
			for (int i = 0; i < nmbDecayKinPart; ++i) {
				assert((*decayKinParticles  )[i]);
				assert((*decayKinMomenta)[i]);
				const string   name = ((TObjString*)(*decayKinParticles)[i])->GetString().Data();
				const TVector3 mom  = *((TVector3*)(*decayKinMomenta)[i]);
				const double   mass = getParticleMass(name);
				int id, charge;
				idAndChargeFromParticleName(name, id, charge);
				outEvt << setprecision(numeric_limits<double>::digits10 + 1)
				       << id << " " << charge << " " << mom.X() << " " << mom.Y() << " " << mom.Z() << " "
				       << sqrt(mass * mass + mom.Mag2()) << endl;
				if (debug) {
					cout << "        particle[" << i << "]: " << name << ", id = " << id << ", "
					     << "charge = " << charge << "; " << mom << endl;
				}
			}
    
			if (debug)
				cout << endl;
		}

		printInfo << "wrote " << nmbEvents << " events to output stream" << endl;
		return true;
	}


	bool
	processTree(TTree&                        tree,
	            const isobarAmplitudePtr&     amplitude,
	            vector<complex<double> >&     ampValues,
	            const long int                maxNmbEvents,
	            const string&                 prodKinParticlesLeafName,
	            const string&                 prodKinMomentaLeafName,
	            const string&                 decayKinParticlesLeafName,
	            const string&                 decayKinMomentaLeafName,
	            const bool                    printProgress)
	{
		if (not amplitude) {
			printWarn << "null pointer to isobar decay amplitude. cannot process tree." << endl;
			return false;
		}

		// create branch pointers and leaf variables
		TBranch*      prodKinParticlesBr  = 0;
		TBranch*      prodKinMomentaBr    = 0;
		TBranch*      decayKinParticlesBr = 0;
		TBranch*      decayKinMomentaBr   = 0;
		TClonesArray* prodKinParticles    = 0;
		TClonesArray* prodKinMomenta      = 0;
		TClonesArray* decayKinParticles   = 0;
		TClonesArray* decayKinMomenta     = 0;
	
		// connect leaf variables to tree branches
		tree.SetBranchAddress(prodKinParticlesLeafName.c_str(),  &prodKinParticles,  &prodKinParticlesBr );
		tree.SetBranchAddress(prodKinMomentaLeafName.c_str(),    &prodKinMomenta,    &prodKinMomentaBr   );
		tree.SetBranchAddress(decayKinParticlesLeafName.c_str(), &decayKinParticles, &decayKinParticlesBr);
		tree.SetBranchAddress(decayKinMomentaLeafName.c_str(),   &decayKinMomenta,   &decayKinMomentaBr  );

		// loop over events
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
			prodKinParticlesBr->GetEntry (eventIndex);
			prodKinMomentaBr->GetEntry   (eventIndex);
			decayKinParticlesBr->GetEntry(eventIndex);
			decayKinMomentaBr->GetEntry  (eventIndex);

			if (   not prodKinParticles  or not prodKinMomenta
			    or not decayKinParticles or not decayKinMomenta) {
				printWarn << "at least one of the input data arrays is a null pointer: "
				          << "        production kinematics: particle names = " << prodKinParticles << ", "
				          << "momenta = " << prodKinMomenta << endl
				          << "        decay kinematics:      particle names = " << decayKinParticles << ", "
				          << "momenta = " << decayKinMomenta << endl
				          << "skipping event." << endl;
				success = false;
				continue;
			}

			const isobarDecayTopologyPtr& decayTopo = amplitude->decayTopology();
			if (decayTopo->readData(*prodKinParticles,  *prodKinMomenta,
			                        *decayKinParticles, *decayKinMomenta))
				ampValues.push_back((*amplitude)());
			else {
				printWarn << "problems reading event[" << eventIndex << "]" << endl;
				success = false;
			}
		}
		return success;
	}


}  // namespace rpwa
