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


#include "massDep.h"
#include "particle.h"

#include "particleDataTable.h"
#include "massDependence.h"
#include "isobarDecayVertex.h"
#include "reportingUtilsRoot.hpp"


extern particleDataTable PDGtable;


using namespace std;
using namespace boost;
using namespace rpwa;


int
main(int argc, char** argv)
{
	printCompilerInfo();
	printSvnVersion();

	rpwa::particleDataTable& pdt = rpwa::particleDataTable::instance();
	pdt.readFile();
	
	massDependence::setDebug(true);

	massDependencePtr massDep = createPiPiSWaveAuMorganPenningtonKachaev();
	printDebug << *massDep << endl;

	particlePtr pi0   = createParticle("pi-", 0);
	particlePtr pi1   = createParticle("pi+", 0);
	particlePtr sigma = createParticle("sigma");

	isobarDecayVertexPtr vertex = createIsobarDecayVertex(sigma, pi0, pi1, 0, 0, massDep);
	pi0->setMomentum(TVector3(-0.0761465106, -0.116917817, 5.89514709));
	pi1->setMomentum(TVector3(-0.0244305532, -0.106013023, 30.6551865));
	vertex->calcParentLzVec();
	printDebug << *vertex << endl;
	TLorentzVector mom = sigma->lzVec();
	printDebug << "sigma = " << mom << "; m = " << mom.M() << endl;
	complex<double> amp = vertex->massDepAmplitude();

	// compare to PWA2000
	PDGtable.initialize("../keyfiles/key5pi/pdgTable.txt");
	AMP_kach pwa2kMassDep;
	printDebug << "PWA2000: ";
	pwa2kMassDep.print();
	cout << endl;
	::particle pwa2kSigma;
	fourVec    pwa2kMom;
	pwa2kSigma.debug();
	pwa2kMom.set(mom.E(), mom.X(), mom.Y(), mom.Z());
	pwa2kSigma.set4P(pwa2kMom);
	pwa2kSigma.print();
	complex<double> pwa2kAmp = pwa2kMassDep.val(pwa2kSigma);
	
	printInfo << maxPrecisionDouble(amp) << " vs. PWA2000 " << maxPrecisionDouble(pwa2kAmp)
	          << "; Delta = " << maxPrecisionDouble(pwa2kAmp - amp) << endl;
}
