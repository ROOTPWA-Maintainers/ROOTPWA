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
//      basic test program for mass dependence functions
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <list>

#include "TStopwatch.h"

#include "particleDataTable.h"
#include "massDependence.h"
#include "isobarDecayVertex.h"
#include "nBodyPhaseSpaceGen.h"
#include "reportingUtilsRoot.hpp"

#include "../pwa2000/libpp/massDep.h"
#include "../pwa2000/libpp/particle.h"


extern particleDataTable PDGtable;


using namespace std;
using namespace boost;
using namespace rpwa;


// constucts Lorentz vector of parent particle
TLorentzVector
constructParent(TRandom3&    random,
                const double mass,            // [GeV/c^2]
                const double maxPt  = 0.150,  // [GeV/c]
                const double maxEta = 1)
{
	double pt      = maxPt * random.Rndm();                 // parent transverse momentum in lab frame
	double eta     = 2 * maxEta * random.Rndm() - maxEta;   // parent pseudorapidity in lab frame
	double phi     = 2 * pi * random.Rndm();                // parent azimuth in lab frame
	double tanhEta = tanh(eta);
	TLorentzVector parent;
	parent.SetPx(pt * cos(phi));
	parent.SetPy(pt * sin(phi));
	parent.SetPz(tanhEta * sqrt((mass * mass + pt * pt) / (1 - tanhEta * tanhEta)));
	parent.SetVectM(parent.Vect(), mass);
	return parent;
}


int
main(int argc, char** argv)
{
	printCompilerInfo();
	printSvnVersion();

	rpwa::particleDataTable& pdt = rpwa::particleDataTable::instance();
	pdt.readFile();
	
	massDependence::setDebug(true);

	// define parameters
	const size_t       nmbEvents = 5;
	const unsigned int seed      = 123456789;
	
	massDependencePtr massDep = createRelativisticBreitWigner();
	// massDependencePtr massDep = createPiPiSWaveAuMorganPenningtonKachaev();
	printDebug << *massDep << endl;

	particlePtr daughter1 = createParticle("pi+");
	particlePtr daughter2 = createParticle("pi-");
	particlePtr parent    = createParticle("rho(770)0");
	// particlePtr parent    = createParticle("sigma");

	const unsigned int L = 2;
	const unsigned int S = 0;

	isobarDecayVertexPtr vertex = createIsobarDecayVertex(parent, daughter1, daughter2, L, S, massDep);

	// setup phase-space generator
	nBodyPhaseSpaceGen psGen;
	const double daughterMasses[2] = {daughter1->mass(), daughter2->mass()};
	const double massRange     [2] = {daughterMasses[0] + daughterMasses[1], 3};    // [GeV/c^2]
	psGen.setWeightType(nBodyPhaseSpaceGen::S_U_CHUNG);
	psGen.setKinematicsType(nBodyPhaseSpaceGen::BLOCK);
	psGen.setDecay(2, daughterMasses);
	psGen.setSeed(seed);
	psGen.setMaxWeight(1.01 * psGen.estimateMaxWeight(massRange[1]));

	// setup PWA2000
	PDGtable.initialize("../keyfiles/key5pi/pdgTable.txt");
	::particle::debug();
	::particle pwa2kDaughter1(PDGtable.get("pi"),      +1);
	::particle pwa2kDaughter2(PDGtable.get("pi"),      -1);
	::particle pwa2kParent   (PDGtable.get("rho(770)"), 0);
	pwa2kDaughter1.print();
	pwa2kDaughter2.print();
	cout << endl;
	::decay pwa2kDecay;
	pwa2kDecay.addChild(pwa2kDaughter1);
	pwa2kDecay.addChild(pwa2kDaughter2);
	pwa2kDecay.setL(L);
	pwa2kDecay.setS(S);
	pwa2kDecay.print();
	cout << endl;
	::massDep* pwa2kMassDep = new breitWigner();
	// ::massDep* pwa2kMassDep = new AMP_kach();
	printDebug << "PWA2000: ";
	pwa2kMassDep->print();
	cout << endl;
	pwa2kParent.setDecay(pwa2kDecay);
	pwa2kParent.setMassDep(pwa2kMassDep);
	pwa2kParent.print();
	cout << endl;

	// loop over events
	TRandom3   random(seed);
	TStopwatch timer;
	timer.Reset();
	timer.Start();
	vector<complex<double> > newAmps(nmbEvents, 0);
	vector<complex<double> > oldAmps(nmbEvents, 0);
	size_t                   countEvts = 0;
	while (countEvts < nmbEvents) {
		const TLorentzVector parentLv
			= constructParent(random, massRange[0] + random.Rndm() * (massRange[1] - massRange[0]));
		if (not psGen.generateDecayAccepted(parentLv))
			continue;
		// calculate new amplitude
		daughter1->setMomentum(psGen.daughter(0).Vect());
		daughter2->setMomentum(psGen.daughter(1).Vect());
		vertex->calcParentLzVec();
		printDebug << *vertex << endl;
		const TLorentzVector mom = parent->lzVec();
		printDebug << "parent = " << mom << "; m = " << mom.M() << endl;
		newAmps[countEvts] = vertex->massDepAmplitude();
		// calculate PWA2000 amplitude
		assert(pwa2kParent.Decay()->_children.size() == 2);
		list< ::particle>::iterator child = pwa2kParent.Decay()->_children.begin();
		::particle& d1 = *child;
    ++child;
    ::particle& d2 = *child;
		fourVec pwa2kMom;
		pwa2kMom.set(daughter1->lzVec().E(), daughter1->lzVec().X(),
		             daughter1->lzVec().Y(), daughter1->lzVec().Z());
		d1.set4P(pwa2kMom);
		pwa2kMom.set(daughter2->lzVec().E(), daughter2->lzVec().X(),
		             daughter2->lzVec().Y(), daughter2->lzVec().Z());
		d2.set4P(pwa2kMom);
		pwa2kMom.set(parent->lzVec().E(), parent->lzVec().X(), parent->lzVec().Y(), parent->lzVec().Z());
		pwa2kParent.set4P(pwa2kMom);
		pwa2kParent.print();
		oldAmps[countEvts] = pwa2kMassDep->val(pwa2kParent);
		cout << endl;
		++countEvts;
	}
	timer.Stop();
	printInfo << "calculated mass dependence for " << nmbEvents << " events" << endl
	          << "    this consumed: ";
	timer.Print();


	double maxAbsDeviation = 0;
	double maxRelDeviation = 0;
	for (unsigned int i = 0; i < newAmps.size(); ++i) {
		const complex<double> delta    = oldAmps[i] - newAmps[i];
		const double          absDelta = abs(delta);
		const double          absVal   = abs(newAmps[i]);
		const double          relDelta = absDelta / absVal;
		if (absDelta > abs(maxAbsDeviation))
			maxAbsDeviation = absDelta;
		if ((absVal != 0) and (relDelta  > maxRelDeviation))
			maxRelDeviation = relDelta;
	}
	printInfo << "maximum absolute deviation is " << maxAbsDeviation << "; "
	          << "maximum relative deviation is " << maxRelDeviation << endl;
	
	// printInfo << maxPrecisionDouble(amp) << " vs. PWA2000 " << maxPrecisionDouble(pwa2kAmp)
	//           << "; Delta = " << maxPrecisionDouble(pwa2kAmp - amp) << endl;
}
