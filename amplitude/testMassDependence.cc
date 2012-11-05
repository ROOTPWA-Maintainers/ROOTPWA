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


#include <boost/timer/timer.hpp>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

#include "particleDataTable.h"
#include "massDependence.h"
#include "isobarDecayVertex.h"
#include "nBodyPhaseSpaceGen.h"
#include "reportingUtilsRoot.hpp"

#include "massDep.h"
#include "particle.h"


extern particleDataTable PDGtable;


using namespace std;
using namespace boost;
using namespace boost::timer;
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


void
compareAmplitudes(const size_t          nmbEvents,
                  const unsigned int    seed,
                  massDependencePtr&    massDep,
                  isobarDecayVertexPtr& vertex,
                  ::massDep*            pwa2kMassDep,
                  ::particle&           pwa2kParent,
                  const string&         name    = "",
                  const double          maxMass = 3)
{
	vertex->setMassDependence(massDep);
	pwa2kParent.setMassDep(pwa2kMassDep);
	printInfo << "comparing '" << name << "' dynamic amplitudes" << endl
	          << "        " << *(vertex->massDependence()) << endl
	          << "        ";
	pwa2kMassDep->print();
	cout << endl
	     << "    using decays" << endl
	     << "        " << *vertex << endl;
	pwa2kParent.print();

	// setup phase-space generator
	particlePtr  daughter1         = vertex->daughter1();
	particlePtr  daughter2         = vertex->daughter2();
	const double daughterMasses[2] = {daughter1->mass(), daughter2->mass()};
	const double massRange     [2] = {daughterMasses[0] + daughterMasses[1], maxMass};    // [GeV/c^2]
	nBodyPhaseSpaceGen psGen;
	psGen.setWeightType(nBodyPhaseSpaceGen::S_U_CHUNG);
	psGen.setKinematicsType(nBodyPhaseSpaceGen::BLOCK);
	psGen.setDecay(2, daughterMasses);
	psGen.setSeed(seed);
	psGen.setMaxWeight(1.01 * psGen.estimateMaxWeight(massRange[1]));

	// loop over events
	TRandom3                 random(seed);
	cpu_timer                timer, timerPwa2k;
	vector<complex<double> > newAmps   (nmbEvents, 0);
	vector<complex<double> > oldAmps   (nmbEvents, 0);
	vector<double>           parentMass(nmbEvents, 0);
	size_t                   countEvts = 0;
	bool                     firstCall = true;
	while (countEvts < nmbEvents) {
		const TLorentzVector parentLv
			= constructParent(random, massRange[0] + random.Rndm() * (massRange[1] - massRange[0]));
		if (not psGen.generateDecayAccepted(parentLv))
			continue;
		// calculate new amplitude
		daughter1->setMomentum(psGen.daughter(0).Vect());
		daughter2->setMomentum(psGen.daughter(1).Vect());
		vertex->calcParentLzVec();
		if (firstCall)
			timer.start();
		else
			timer.resume();
		newAmps[countEvts] = vertex->massDepAmplitude();
		timer.stop();
		parentMass[countEvts] = parentLv.M();
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
		pwa2kMom.set(vertex->parent()->lzVec().E(), vertex->parent()->lzVec().X(),
		             vertex->parent()->lzVec().Y(), vertex->parent()->lzVec().Z());
		pwa2kParent.set4P(pwa2kMom);
		if (firstCall)
			timerPwa2k.start();
		else
			timerPwa2k.resume();
		oldAmps[countEvts] = pwa2kMassDep->val(pwa2kParent);
		timerPwa2k.stop();
		++countEvts;
		if (firstCall)
			firstCall = false;
	}
	printInfo << "calculated mass dependence for " << nmbEvents << " events" << endl
	          << "    this consumed: (new code) "
	          << timer.format(3, "[%u sec user] + [%s sec sys] = [%t sec]")
	          << " vs. (PWA2000) " << timerPwa2k.format(3, "[%t sec] = [%u sec user] + [%s sec sys]")
	          << endl;

	// compare amplitude values
	TH2D* hAmp         = new TH2D(("h" + name + "Amp"        ).c_str(), (name + " Amplitude"                 ).c_str(), 1000, -1, 1, 1000, -0.5, 1.5);
	TH1D* hAmpReal     = new TH1D(("h" + name + "AmpReal"    ).c_str(), (name + " Amplitude (Real Part)"     ).c_str(), 1000, -1.2, 1.2);
	TH1D* hAmpImag     = new TH1D(("h" + name + "AmpImag"    ).c_str(), (name + " Amplitude (Imag Part)"     ).c_str(), 1000, -1.2, 1.2);
	TH2D* hInt         = new TH2D(("h" + name + "Int"        ).c_str(), (name + " Intensity"                 ).c_str(), 1000, massRange[0] - 0.2, massRange[1] + 0.2, 1000,   0, 1.2);
	TH2D* hPhase       = new TH2D(("h" + name + "Phase"      ).c_str(), (name + " Phase"                     ).c_str(), 1000, massRange[0] - 0.2, massRange[1] + 0.2, 1000, -200, 200);
	TH1D* hAbsDiffReal = new TH1D(("h" + name + "AbsDiffReal").c_str(), (name + " Absolute Diff. (Real Part)").c_str(), 100000, -1e-12, 1e-12);
	TH1D* hAbsDiffImag = new TH1D(("h" + name + "AbsDiffImag").c_str(), (name + " Absolute Diff. (Imag Part)").c_str(), 100000, -1e-12, 1e-12);
	TH1D* hRelDiffReal = new TH1D(("h" + name + "RelDiffReal").c_str(), (name + " Relative Diff. (Real Part)").c_str(), 100000, -1e-12, 1e-12);
	TH1D* hRelDiffImag = new TH1D(("h" + name + "RelDiffImag").c_str(), (name + " Relative Diff. (Imag Part)").c_str(), 100000, -1e-12, 1e-12);
	double maxAbsDiffReal = 0;
	double maxRelDiffReal = 0;
	double maxAbsDiffImag = 0;
	double maxRelDiffImag = 0;
	for (unsigned int i = 0; i < newAmps.size(); ++i) {
		const complex<double> absDiff = oldAmps[i] - newAmps[i];
		const complex<double> relDiff = complex<double>(absDiff.real() / oldAmps[i].real(),
		                                                absDiff.imag() / oldAmps[i].imag());
		hAmpReal->Fill(newAmps[i].real());
		hAmpImag->Fill(newAmps[i].imag());
		hAmp->Fill(newAmps[i].real(), newAmps[i].imag());
		hInt->Fill(parentMass[i], abs(newAmps[i]));
		hPhase->Fill(parentMass[i], TMath::RadToDeg() * arg(newAmps[i]));
		hAbsDiffReal->Fill(absDiff.real());
		hAbsDiffImag->Fill(absDiff.imag());
		hRelDiffReal->Fill(relDiff.real());
		hRelDiffImag->Fill(relDiff.imag());
		if (abs(absDiff.real()) > maxAbsDiffReal)
			maxAbsDiffReal = abs(absDiff.real());
		if (abs(absDiff.imag()) > maxAbsDiffImag)
			maxAbsDiffImag = abs(absDiff.imag());
		if ((oldAmps[i].real() != 0) and (abs(relDiff.real()) > maxRelDiffReal))
			maxRelDiffReal = abs(relDiff.real());
		if ((oldAmps[i].imag() != 0) and (abs(relDiff.imag()) > maxRelDiffImag))
			maxRelDiffImag = abs(relDiff.imag());
	}
	printSucc << endl
	          << "    real part: maximum absolute deviation is " << maxAbsDiffReal << "; "
	          << "maximum relative deviation is " << maxRelDiffReal << endl
	          << "    imag part: maximum absolute deviation is " << maxAbsDiffImag << "; "
	          << "maximum relative deviation is " << maxRelDiffImag << endl;
}


int
main(int argc, char** argv)
{
	printCompilerInfo();
	printSvnVersion();
	cout << endl;

	// define parameters
	const size_t       nmbEvents = 1000000;
	//const size_t       nmbEvents = 10;
	const unsigned int seed      = 123456789;
	
	rpwa::particleDataTable& pdt = rpwa::particleDataTable::instance();
	pdt.readFile();
	PDGtable.initialize("../keyfiles/key5pi/pdgTable.txt");
	
	// massDependence::setDebug(true);
	// ::particle::debug();

	// define decay
	// particlePtr daughter1  = createParticle("pi+");
	// ::particle  pwa2kDaughter1(PDGtable.get("pi"), +1);
	// particlePtr daughter2  = createParticle("pi-");
	// ::particle  pwa2kDaughter2(PDGtable.get("pi"), -1);
	// const unsigned int L = 2;
	// const unsigned int S = 0;
	// particlePtr parent  = createParticle("rho(770)0");
	// ::particle  pwa2kParent(PDGtable.get("rho(770)"), 0);
	// const unsigned int L = 0;
	// const unsigned int S = 0;
	// particlePtr parent  = createParticle("sigma");
	// ::particle  pwa2kParent(PDGtable.get("sigma"), 0);

	// particlePtr daughter1  = createParticle("rho(770)0");
	// ::particle  pwa2kDaughter1(PDGtable.get("rho(770)"), 0);
	// particlePtr daughter2  = createParticle("sigma0");
	// ::particle  pwa2kDaughter2(PDGtable.get("sigma"), 0);
	particlePtr daughter1  = createParticle("a1(1260)+");
	::particle  pwa2kDaughter1(PDGtable.get("a1(1260)"), +1);
	particlePtr daughter2  = createParticle("pi-");
	::particle  pwa2kDaughter2(PDGtable.get("pi"), -1);
	const unsigned int L = 0;
	const unsigned int S = 2;
	particlePtr parent  = createParticle("rho(1700)0");
	::particle  pwa2kParent(PDGtable.get("rho(1700)"), 0);

	// setup data structures
	isobarDecayVertexPtr vertex = createIsobarDecayVertex(parent, daughter1, daughter2, L, S);
	::decay pwa2kDecay;
	pwa2kDecay.addChild(pwa2kDaughter1);
	pwa2kDecay.addChild(pwa2kDaughter2);
	pwa2kDecay.setL(L);
	pwa2kDecay.setS(S);
	pwa2kParent.setDecay(pwa2kDecay);

	// compare amplitudes
	TFile*       outFile = TFile::Open("testMassDependence.root", "RECREATE");
	const size_t nmbMassDep = 2;
	const string names[nmbMassDep] = {
		// "FLAT",
		"BW",
		// "AMP_M",
		// "AMP_VES",
		// "AMP_KACH",
		"RHO_PRIME"
	};
	massDependencePtr massDep[nmbMassDep] = {
		// createFlatMassDependence(),
		createRelativisticBreitWigner(),
		// createPiPiSWaveAuMorganPenningtonM(),
		// createPiPiSWaveAuMorganPenningtonVes(),
		// createPiPiSWaveAuMorganPenningtonKachaev(),
		createRhoPrimeMassDep()
	};
	::massDep* pwa2kMassDep[nmbMassDep] = {
		// new flat       (),
		new breitWigner(),
		// new AMP_M      (),
		// new AMP_ves    (),
		// new AMP_kach   (),
		new rhoPrime   ()
	};
	for (size_t i = 0; i < nmbMassDep; ++i) {
		compareAmplitudes(nmbEvents, seed, massDep[i], vertex, pwa2kMassDep[i], pwa2kParent, names[i]);	
		cout << endl;
	}
	outFile->Write();
	outFile->Close();
}
