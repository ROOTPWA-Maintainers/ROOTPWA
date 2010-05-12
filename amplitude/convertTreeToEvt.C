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
// $Rev:: 220                         $: revision of last commit
// $Author:: bgrube                   $: author of last commit
// $Date:: 2010-05-11 16:49:30 +0200 #$: date of last commit
//
// Description:
//      program that converts standrad ASCII PWA2000 .evt files into
//      new ROOT tree format
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <fstream>
#include <string>
#include <cassert>
#include <iomanip>
#include <limits>

#include "TChain.h"
#include "TClonesArray.h"
#include "TObjString.h"
#include "TVector3.h"

#include "pputil.h"
#include "particleDataTable.h"
#include "particleProperties.h"
#include "utilities.h"


using namespace std;
using namespace rpwa;


void
getParticleId(string name,
	      int&   id,
	      int&   charge)
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
convertTreeToEvt(const string&  inFileNamePattern = "testEvents.root",
		 const string&  outFileName       = "testTree.evt",
		 const long int maxNmbEvents      = -1,
		 const string&  inTreeName        = "rootPwaEvtTree",
		 const string&  pdgTableFileName  = "./particleDataTable.txt",
		 const bool     debug             = false)
{
  rpwa::particleDataTable& pdt = rpwa::particleDataTable::instance();
  pdt.readFile(pdgTableFileName);

  // open input file
  TChain chain(inTreeName.c_str());
  if (chain.Add(inFileNamePattern.c_str()) < 1) {
    printWarn << "cannot open input file(s) '" << inFileNamePattern << "'" << endl;
    return false;
  }
  const long int nmbEventsChain = chain.GetEntries();
  chain.GetListOfFiles()->ls();

  // create output file
  ofstream outFile(outFileName.c_str());
  if (!outFile) {
    printWarn << "cannot open output file '" << outFileName << "'" << endl;
    return false;
  }

  // create leaf variables
  // TClonesArray* initialStateNames   = new TClonesArray("TObjString", 0);
  // TClonesArray* initialStateMomenta = new TClonesArray("TVector3",   0);
  // TClonesArray* finalStateNames     = new TClonesArray("TObjString", 0);
  // TClonesArray* finalStateMomenta   = new TClonesArray("TVector3",   0);
  TClonesArray* initialStateNames   = 0;
  TClonesArray* initialStateMomenta = 0;
  TClonesArray* finalStateNames     = 0;
  TClonesArray* finalStateMomenta   = 0;

  // connect leaf variables to tree branches
  // chain.GetBranch("initialStateNames"  )->SetAutoDelete(kFALSE);
  // chain.GetBranch("initialStateMomenta")->SetAutoDelete(kFALSE);
  // chain.GetBranch("finalStateNames"    )->SetAutoDelete(kFALSE);
  // chain.GetBranch("finalStateMomenta"  )->SetAutoDelete(kFALSE);
  chain.SetBranchAddress("initialStateNames",   &initialStateNames  );
  chain.SetBranchAddress("initialStateMomenta", &initialStateMomenta);
  chain.SetBranchAddress("finalStateNames",     &finalStateNames    );
  chain.SetBranchAddress("finalStateMomenta",   &finalStateMomenta  );

  // loop over events
  const long int nmbEvents = ((maxNmbEvents > 0) ? maxNmbEvents : nmbEventsChain);
  for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
    if (!debug)
      progressIndicator(eventIndex, nmbEvents);

    // initialStateNames->Clear  ();
    // initialStateMomenta->Clear();
    // finalStateNames->Clear    ();
    // finalStateMomenta->Clear  ();

    if (chain.LoadTree(eventIndex) < 0)
      break;
    chain.GetEntry(eventIndex);

    assert(initialStateNames  );
    assert(initialStateMomenta);
    assert(initialStateNames->GetSize() == initialStateMomenta->GetSize());
    const unsigned int nmbIsPart = initialStateNames->GetSize();

    assert(finalStateNames  );
    assert(finalStateMomenta);
    assert(finalStateNames->GetSize() == finalStateMomenta->GetSize());
    const unsigned int nmbFsPart = finalStateNames->GetSize();

    outFile << nmbIsPart + nmbFsPart << endl;

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
      getParticleId(name, id, charge);
      outFile << setprecision(numeric_limits<double>::digits10 + 1)
	      << id << " " << charge << " " << mom.X() << " " << mom.Y() << " " << mom.Z() << " "
	      << sqrt(mass * mass + mom.Mag2()) << endl;
      if (debug) {
	cout << "        particle[" << i << "]: " << name << ", id = " << id << ", "
	     << "charge = " << charge << "; ";
	mom.Print();
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
      getParticleId(name, id, charge);
      outFile << setprecision(numeric_limits<double>::digits10 + 1)
	      << id << " " << charge << " " << mom.X() << " " << mom.Y() << " " << mom.Z() << " "
	      << sqrt(mass * mass + mom.Mag2()) << endl;
      if (debug) {
	cout << "        particle[" << i << "]: " << name << ", id = " << id << ", "
	     << "charge = " << charge << "; ";
	mom.Print();
      }
    }
    
    if (debug)
      cout << endl;
  }

  return true;
}
