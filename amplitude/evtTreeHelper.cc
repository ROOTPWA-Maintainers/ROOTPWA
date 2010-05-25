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

#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TObjString.h"
#include "TVector3.h"

#include "pputil.h"
#include "particleDataTable.h"
#include "utilities.h"


using namespace std;
using namespace rpwa;


string
particleNameFromGeantId(const int id,
			const int charge)
{
  assert((charge == -1) || (charge == 0) || (charge == +1));
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
  id     = name2id(name, charge);
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
  if (!prop) {
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
		const string&  leafNameIsPartNames,
		const string&  leafNameIsPartMomenta,
		const string&  leafNameFsPartNames,
		const string&  leafNameFsPartMomenta,
		const bool     debug)
{
  if (!inEvt || !inEvt.good()) {
    printWarn << "cannot read from input stream" << endl;
    return false;
  }

  // create leaf variables
  TClonesArray* initialStateNames   = new TClonesArray("TObjString", 1);
  TClonesArray* initialStateMomenta = new TClonesArray("TVector3",   1);
  TClonesArray* finalStateNames     = new TClonesArray("TObjString", 0);
  TClonesArray* finalStateMomenta   = new TClonesArray("TVector3",   0);

  // connect leaf variables to tree branches
  const int split   = 0;
  const int bufSize = 256000;
  outTree.Branch(leafNameIsPartNames.c_str(),   "TClonesArray", &initialStateNames,   bufSize, split);
  outTree.Branch(leafNameIsPartMomenta.c_str(), "TClonesArray", &initialStateMomenta, bufSize, split);
  outTree.Branch(leafNameFsPartNames.c_str(),   "TClonesArray", &finalStateNames,     bufSize, split);
  outTree.Branch(leafNameFsPartMomenta.c_str(), "TClonesArray", &finalStateMomenta,   bufSize, split);

  // loop over events and fill tree
  bool     success     = true;
  long int countEvents = 0;
  long int countLines  = 0;
  while (inEvt.good()) {

    // read event
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
    
    // read initial state particle data (beam only)
    initialStateNames->Clear  ();
    initialStateMomenta->Clear();
    if (getline(inEvt, line)) {
      ++countLines;
      stringstream lineStream(line);
      int    id = 0, charge = 0;
      double momX = 0, momY = 0, momZ = 0, E = 0;
      if (lineStream >> id >> charge >> momX >> momY >> momZ >> E) {
      	const string name = particleNameFromGeantId(id , charge);
      	new((*initialStateNames  )[0]) TObjString(name.c_str());
      	new((*initialStateMomenta)[0]) TVector3  (momX, momY, momZ);
      } else {
      	printWarn << "event " << countEvents + 1 << ": error reading beam data "
      		  << "from line " << countLines << ": " << line << endl;
      	success = false;
      }
    } else
      break;
    assert(initialStateNames->GetEntriesFast() == initialStateMomenta->GetEntriesFast());
    const unsigned int nmbIsPart = initialStateNames->GetEntriesFast();
    if (debug) {
      printInfo << nmbIsPart << " initial state particles:" << endl;
      for (unsigned int i = 0; i < nmbIsPart; ++i)
      	cout << "        particle[" << i << "]: "
      	     << ((TObjString*)(*initialStateNames)[i])->GetString() << "; "
      	     << *((TVector3*)(*initialStateMomenta)[i]) << endl;
    }

    // read final state particle data
    finalStateNames->Clear  ();
    finalStateMomenta->Clear();
    for (int i = 0; i < nmbParticles - 1; ++i) {
      if (getline(inEvt, line)) {
    	++countLines;
    	stringstream lineStream(line);
    	int    id, charge;
    	double momX, momY, momZ, E;
    	if (lineStream >> id >> charge >> momX >> momY >> momZ >> E) {
	  const string name = particleNameFromGeantId(id , charge);
    	  new((*finalStateNames  )[i]) TObjString(name.c_str());
    	  new((*finalStateMomenta)[i]) TVector3  (momX, momY, momZ);
    	} else {
    	  printWarn << "event " << countEvents + 1 << ": error reading final state "
		    << "particle[" << i << "] data from line " << countLines << ": " << line << endl;
    	  success = false;
    	}
      } else
	break;
    }
    assert(finalStateNames->GetEntriesFast() == finalStateMomenta->GetEntriesFast());
    const unsigned int nmbFsPart = finalStateNames->GetEntriesFast();
    if (debug) {
      printInfo << nmbFsPart << " final state particles:" << endl;
      for (unsigned int i = 0; i < nmbFsPart; ++i)
      	cout << "        particle[" << i << "]: "
      	     << ((TObjString*)(*finalStateNames)[i])->GetString() << "; "
      	     << *((TVector3*)(*finalStateMomenta)[i]) << endl;
    }

    outTree.Fill();
    ++countEvents;
    if ((maxNmbEvents > 0) && (countEvents >= maxNmbEvents))
      break;
    if (debug)
      cout << endl;
  }

  printInfo << "read " << countLines << " lines from input stream and wrote "
	    << countEvents << " events to tree '" << outTree.GetName() << "'" << endl;
  return success;
}


bool
writeEvtFromTree(TChain&        inTree,
		 ostream&       outEvt,
		 const long int maxNmbEvents          = -1,
		 const string&  pdgTableFileName      = "./particleDataTable.txt",
		 const string&  inTreeName            = "rootPwaEvtTree",
		 const string&  leafNameIsPartNames   = "initialStateNames",
		 const string&  leafNameIsPartMomenta = "initialStateMomenta",
		 const string&  leafNameFsPartNames   = "finalStateNames",
		 const string&  leafNameFsPartMomenta = "finalStateMomenta",
		 const bool     debug                 = false)
{
  rpwa::particleDataTable& pdt = rpwa::particleDataTable::instance();
  pdt.readFile(pdgTableFileName);

  const long int nmbEventsTree = inTree.GetEntries();
  inTree.GetListOfFiles()->ls();
  if (!outEvt) {
    printWarn << "cannot write to output stream" << endl;
    return false;
  }

  // create leaf variables
  TClonesArray* initialStateNames   = 0;
  TClonesArray* initialStateMomenta = 0;
  TClonesArray* finalStateNames     = 0;
  TClonesArray* finalStateMomenta   = 0;

  // connect leaf variables to tree branches
  inTree.SetBranchAddress(leafNameIsPartNames.c_str(),   &initialStateNames  );
  inTree.SetBranchAddress(leafNameIsPartMomenta.c_str(), &initialStateMomenta);
  inTree.SetBranchAddress(leafNameFsPartNames.c_str(),   &finalStateNames    );
  inTree.SetBranchAddress(leafNameFsPartMomenta.c_str(), &finalStateMomenta  );
			 
  // loop over events
  const long int nmbEvents = ((maxNmbEvents > 0) ? min(maxNmbEvents, nmbEventsTree)
			                         : nmbEventsTree);
  for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
    if (!debug)
      progressIndicator(eventIndex, nmbEvents);

    if (inTree.LoadTree(eventIndex) < 0)
      break;
    inTree.GetEntry(eventIndex);

    assert(initialStateNames  );
    assert(initialStateMomenta);
    assert(initialStateNames->GetEntriesFast() == initialStateMomenta->GetEntriesFast());
    const unsigned int nmbIsPart = initialStateNames->GetEntriesFast();

    assert(finalStateNames  );
    assert(finalStateMomenta);
    assert(finalStateNames->GetEntriesFast() == finalStateMomenta->GetEntriesFast());
    const unsigned int nmbFsPart = finalStateNames->GetEntriesFast();

    outEvt << nmbIsPart + nmbFsPart << endl;

    if (debug)
      printInfo << "event[" << eventIndex << "]: " << nmbIsPart
		<< " initial state particles:" << endl;
    for (unsigned int i = 0; i < nmbIsPart; ++i) {
      assert((*initialStateNames  )[i]);
      assert((*initialStateMomenta)[i]);
      const string   name = ((TObjString*)(*initialStateNames)[i])->GetString().Data();
      const TVector3 mom  = *((TVector3*)(*initialStateMomenta)[i]);
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
      printInfo << "event[" << eventIndex << "]: " << nmbFsPart << " final state particles:" << endl;
    for (unsigned int i = 0; i < nmbFsPart; ++i) {
      assert((*finalStateNames  )[i]);
      assert((*finalStateMomenta)[i]);
      const string   name = ((TObjString*)(*finalStateNames)[i])->GetString().Data();
      const TVector3 mom = *((TVector3*)(*finalStateMomenta)[i]);
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
