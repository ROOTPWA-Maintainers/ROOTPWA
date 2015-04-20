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
//      basic test program for amplitude classes
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <fstream>

#include <boost/progress.hpp>

#include "TVector3.h"
#include "TLorentzRotation.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"

#include "mathUtils.hpp"
#include "reportingUtilsRoot.hpp"
#include "reportingUtilsEnvironment.h"
#include "conversionUtils.hpp"
#include "particleDataTable.h"
#include "diffractiveDissVertex.h"
#include "massDependence.h"
#include "waveDescription.h"
#include "isobarAmplitude.h"
#include "isobarHelicityAmplitude.h"
#include "isobarCanonicalAmplitude.h"
#include "evtTreeHelper.h"


using namespace std;
using namespace boost;
using namespace rpwa;


int
main(int argc, char** argv)
{
	printCompilerInfo();
	printGitHash();

	rpwa::particleDataTable& pdt = rpwa::particleDataTable::instance();
	pdt.readFile();
	TStopwatch timer;

	if (0) {
		isobarDecayVertex::setDebug(true);
		decayTopology::setDebug(true);
		isobarDecayTopology::setDebug(true);
		massDependence::setDebug(true);
		diffractiveDissVertex::setDebug(true);
		isobarAmplitude::setDebug(true);
		isobarHelicityAmplitude::setDebug(true);
		isobarCanonicalAmplitude::setDebug(true);
		waveDescription::setDebug(true);
	}

	if (0) {

		{
			TLorentzVector p(0.5, 0.75, 1, 2);
			TVector3       n = TVector3(0, 0, 1).Cross(p.Vect());
			TRotation R1;
			R1.RotateZ(-n.Phi());
			R1.RotateY(piHalf - n.Theta());
			R1.RotateZ(piHalf);
			n *= R1;
			p *= R1;
			cout << "R1 -> " << n << "    " << p << endl;
			// rotate about yHfAxis so that daughter momentum is along z-axis
			TRotation R2;
			R2.RotateY(-signum(p.X()) * p.Theta());
			p *= R2;
			cout << "R2 -> " << p << endl;
			// boost into daughter RF
			TLorentzRotation L3(-p.BoostVector());
			cout << "L3 -> " << L3 * p << endl;

			R1.Transform(R2);
			TLorentzRotation L(R1);
			L.Boost(-p.BoostVector());
			p = TLorentzVector(0.5, 0.75, 1, 2);
			p *= L;
			cout << "L -> " << p << endl;
		}

		{
			particlePtr X = createParticle("X");
			TLorentzVector p(0.5, 0.75, 1, 2);
			X->setLzVec(p);
			isobarHelicityAmplitude amp;
			TLorentzRotation L = amp.hfTransform(X->lzVec());
			cout << "!!! L -> " << L * p << endl;
		}
	}

	if (0) {
		TLorentzVector beam(1,   0.5,  180, 182);
		TLorentzVector X   (0.5, 0.75, 1,   3);
		isobarHelicityAmplitude amp;
		TLorentzRotation L = amp.gjTransform(beam, X);
		cout << "!!! L -> " << L * X << endl;
	}

	if (0) {
		// define final state particles
		particlePtr pi0 = createParticle("pi-", 0);
		particlePtr pi1 = createParticle("pi+", 0);
		particlePtr pi2 = createParticle("pi-", 1);
		particlePtr pi3 = createParticle("pi+", 1);
		particlePtr pi4 = createParticle("pi-", 2);
		// define isobars
		particlePtr sigma = createParticle("sigma");
		particlePtr a1    = createParticle("a1(1260)+");
		particlePtr f1    = createParticle("f1(1285)");
		// define X-system
		//                                   2I  G  2J  P  C  2M  refl
		particlePtr X = createParticle("X-", 2, -1, 4, +1, 0, 2, +1);
		// define production vertex
		particlePtr              beam     = createParticle("pi-");
		particlePtr              target   = createParticle("p+");
		diffractiveDissVertexPtr prodVert = createDiffractiveDissVertex(beam, target, X);
		// define vertices
		massDependencePtr    massDep = createRelativisticBreitWigner();
		isobarDecayVertexPtr vert0   = createIsobarDecayVertex(X,     pi4, f1,    2, 2);
		isobarDecayVertexPtr vert1   = createIsobarDecayVertex(f1,    pi2, a1,    2, 2, massDep);
		isobarDecayVertexPtr vert2   = createIsobarDecayVertex(a1,    pi3, sigma, 2, 0, massDep);
		isobarDecayVertexPtr vert3   = createIsobarDecayVertex(sigma, pi0, pi1,   0, 0, massDep);
		// set Lorentz vectors
		beam->setLzVec(TLorentzVector(0.104385398, 0.0132061851, 189.987978, 189.988058));
		pi0->setLzVec(TLorentzVector(-0.0761465106, -0.116917817, 5.89514709, 5.89844947));
		pi1->setLzVec(TLorentzVector(-0.0244305532, -0.106013023, 30.6551865, 30.6556973));
		pi2->setLzVec(TLorentzVector(0.000287952441, 0.10263611, 3.95724077, 3.96103114));
		pi3->setLzVec(TLorentzVector(0.0299586212, 0.176440177, 115.703054, 115.703277));
		pi4->setLzVec(TLorentzVector(0.176323963, -0.0985753246, 30.9972271, 30.9981995));
		// build graph
		vector<isobarDecayVertexPtr> decayVertices;
		decayVertices.push_back(vert3);
		decayVertices.push_back(vert1);
		decayVertices.push_back(vert2);
		decayVertices.push_back(vert0);
		vector<particlePtr> fsParticles;
		fsParticles.push_back(pi0);
		fsParticles.push_back(pi1);
		fsParticles.push_back(pi2);
		fsParticles.push_back(pi3);
		fsParticles.push_back(pi4);
		isobarDecayTopologyPtr topo = createIsobarDecayTopology(prodVert, decayVertices, fsParticles);
		// topo->checkTopology();
		// topo->checkConsistency();
		topo->fillKinematicsDataCache();
		isobarHelicityAmplitude amp(topo);
		cout << topo;
		complex<double>         decayAmp = amp.amplitude();
		cout << "!!! decay amplitude = " << decayAmp << endl;
	}

	if (0) {
		// define final state particles
		particlePtr pi0 = createParticle("pi-", 0);
		particlePtr pi1 = createParticle("pi+", 0);
		// define X-system
		//                                   2I  G  2J  P   C  2M
		particlePtr X = createParticle("X0", 2, +1, 2, -1, -1, 0);
		// define production vertex
		particlePtr              beam     = createParticle("pi-");
		particlePtr              target   = createParticle("p+");
		diffractiveDissVertexPtr prodVert = createDiffractiveDissVertex(beam, target, X);
		// define vertices and final-state particles
		isobarDecayVertexPtr         vert0 = createIsobarDecayVertex(X, pi0, pi1, 2, 0);
		vector<isobarDecayVertexPtr> decayVertices;
		decayVertices.push_back(vert0);
		vector<particlePtr> fsParticles;
		fsParticles.push_back(pi0);
		fsParticles.push_back(pi1);
		isobarDecayTopologyPtr topo = createIsobarDecayTopology(prodVert, decayVertices, fsParticles);
		topo->checkTopology();
		topo->checkConsistency();
		isobarAmplitude::setDebug         (true);
		isobarHelicityAmplitude::setDebug (true);
		isobarCanonicalAmplitude::setDebug(true);
		isobarAmplitudePtr amp[2] = {createIsobarHelicityAmplitude (topo),
		                             createIsobarCanonicalAmplitude(topo)};
		for (unsigned int i = 0; i < 2; ++i) {
			amp[i]->enableBoseSymmetrization(false);
			amp[i]->enableReflectivityBasis (false);
			cout << *(amp[i]);
		}
		TChain chain("rootPwaEvtTree");
		chain.Add("../../../massBins/2004/Q3PiData/template.both/1260.1300/1260.1300.root");
		chain.GetListOfFiles()->ls();
		vector<complex<double> > ampValues[2];

		TClonesArray* prodKinPartNames  = 0;
		TClonesArray* decayKinPartNames = 0;
		if (not getParticleNamesFromRootFile(*(chain.GetFile()), prodKinPartNames, decayKinPartNames,
		                                     "rootPwaEvtTree"))
			cout << "cannot find production and/or decay kinematics particle names "
			     << "in input file '" << chain.GetFile()->GetName() << "'." << endl;
		for (unsigned int i = 0; i < 2; ++i)
			processTree(chain, *prodKinPartNames, *decayKinPartNames, amp[i], ampValues[i], 2);
		for (unsigned int i = 0; i < ampValues[0].size(); ++i)
			cout << "amplitude[" << i << "] = " << ampValues[0][i] << " vs. " << ampValues[1][i] << "; "
			     << "ratio = " << ampValues[0][i] / ampValues[1][i] << endl;
	}

	if (1) {
		const long int maxNmbEvents = 10000;
		//const long int maxNmbEvents = 1;

		// const string   newKeyFileName = "test.key";
		// const string   oldKeyFileName = "testAmplitude.key";

		// const string   newKeyFileName = "../../keyfiles/key3pi/SET1_new/1-0-+0+rho770_11_pi-.key";
		// const string   oldKeyFileName = "../../keyfiles/key3pi/SET1/1-0-+0+rho770_11_pi-.key";
		// const string   evtInFileName  = "../../../massBins/2004/Q3PiData/template.both/1260.1300/1260.1300.evt";
		// const string   rootInFileName = "../../../massBins/2004/Q3PiData/template.both/1260.1300/1260.1300.root";

		// waves w/o isospin symmetrization

		// // rel. delta = (1.4218332033726580e-09, 1.4065747228912715e-09)
		// // rms 1.44e-9, 8.34e-10
		// const string newKeyFileName = "test5pi/charly/nosym/1-0-00+f01500=sigma_00_sigma_00_pi-.key";
		// const string oldKeyFileName = "test5pi/sebastian/nosym/1-0-+0+pi-_00_f01500=sigma_0_sigma.key";

		// // rel. delta = (1.3390541006872489e-09, 1.3820297766255656e-09)
		// // rms 1.15e-9, 1.19e-9
		// const string newKeyFileName = "test5pi/charly/nosym/1-1+00+sigma_22_a21320-=rho770_21_pi-.key";
		// const string oldKeyFileName = "test5pi/sebastian/nosym/1-1++0+sigma_22_a21320=pi-_2_rho770.key";

		// // rel. delta = (1.3969833147493598e-09, 1.3631710002893191e-09)
		// // rms 1.71e-9, 1.46e-9
		// const string newKeyFileName = "test5pi/charly/nosym/1-2-00+rho770_02_a21320-=rho770_21_pi-.key";
		// const string oldKeyFileName = "test5pi/sebastian/nosym/1-2-+0+rho770_02_a21320=pi-_2_rho770.key";

		// rel. delta = (4.2414900268409422e-01, -3.7760164904837606e-01)
		// rel. delta = (-1.1363857545500716e-01, 1.1934564489491958e-03)
		// rel. delta = (-2.1179976175700090e-05, 3.5535687503020203e-06)
		// rel. delta = (-3.4768493036075876e-09, 3.5819563781228545e-07) noBose noRefl
		// rel. delta = (-1.3114336807010326e-07, 3.4071135025521628e-07) full sym
		// rel. delta = (-4.2926706747269698e-05, -6.5415181989729538e-06) full sym; orig .evt
		// rel. delta = (-7.1429252866706858e-13, -8.0913804091617701e-12) noBose noRefl; fixed PDG table
		// rel. delta = (1.3880163165229989e-09, 1.3994113453639803e-09) full sym
		// rel. delta = (-4.2794189442389053e-05, -6.8808316461363639e-06) full sym; orig  .evt
		// rms 1.21e-9, 2.92e-10
		// const string newKeyFileName = "test5pi/charly/nosym/1-2-00+sigma_20_pi1800-=sigma_00_pi-.key";
		// const string oldKeyFileName = "test5pi/sebastian/nosym/1-2-+0+sigma_20_pi1800=pi-_0_sigma.key";

		// waves with isospin symmetrization

		// rel. delta = (1.4271813795388615e-09, 1.3774171038156871e-09)
		// rel. delta = (1.4319815006796074e-09, 1.4120124973618052e-09)
		// rms 1.51e-9, 1.91e-9
		// phi = 0, R = 1 ---> 1 / sqrt(2) * (a1 + a2)
		// const string newKeyFileName = "test5pi/charly/sym/1-1+00+rho1700=a11260-=rho770_01_pi-_01_pi+_01_pi-.key";
		// const string oldKeyFileName = "test5pi/sebastian/sym/1-1++0+pi-_01_rho1700=a11269=pi-_0_rho770_0_pi+.key";
		const string newKeyFileName = "test5pi/charly/sym/1-1+00+rho1700=a11260+=rho770_01_pi+_01_pi-_01_pi-.key";
		const string oldKeyFileName = "test5pi/sebastian/sym/1-1++0+pi-_01_rho1700=a11269=pi+_0_rho770_0_pi-.key";

		const string evtInFileName  = "test5pi/1900.1960.genbod.regen.evt";
		const string rootInFileName = "test5pi/1900.1960.genbod.root";
		// const string evtInFileName  = "test5pi/oneEvent.evt";
		// const string rootInFileName = "test5pi/oneEvent.root";

		// decayTopology::setDebug(true);
		// isobarDecayTopology::setDebug(true);
		//massDependence::setDebug(true);
		//isobarAmplitude::setDebug(true);
		//isobarHelicityAmplitude::setDebug(true);

		waveDescription    waveDesc;
		isobarAmplitudePtr amp;
		if (waveDesc.parseKeyFile(newKeyFileName) and waveDesc.constructAmplitude(amp)) {
			isobarDecayTopologyPtr topo = amp->decayTopology();
			printInfo << *amp;
			amp->init();

			// read data from tree
			const string&            inTreeName               = "rootPwaEvtTree";
			const string&            prodKinPartNamesObjName  = "prodKinParticles";
			const string&            prodKinMomentaLeafName   = "prodKinMomenta";
			const string&            decayKinPartNamesObjName = "decayKinParticles";
			const string&            decayKinMomentaLeafName  = "decayKinMomenta";
			vector<complex<double> > myAmps;
			// open input file
			vector<TTree*> inTrees;
			TClonesArray*  prodKinPartNames  = 0;
			TClonesArray*  decayKinPartNames = 0;
			{
				vector<string> rootFileNames(1, rootInFileName);
				vector<string> evtFileNames;
				if (not openRootEvtFiles(inTrees, prodKinPartNames, decayKinPartNames,
				                         rootFileNames, evtFileNames,
				                         inTreeName, prodKinPartNamesObjName, prodKinMomentaLeafName,
				                         decayKinPartNamesObjName, decayKinMomentaLeafName, true)) {
					printErr << "problems opening input files. Aborting..." << endl;
					exit(1);
				}
			}
			const long int nmbEventsChain = inTrees[0]->GetEntries();

			// create branch pointers and leaf variables
			TBranch*      prodKinMomentaBr  = 0;
			TBranch*      decayKinMomentaBr = 0;
			TClonesArray* prodKinMomenta    = 0;
			TClonesArray* decayKinMomenta   = 0;

			// connect leaf variables to tree branches
			inTrees[0]->SetBranchAddress(prodKinMomentaLeafName.c_str(),  &prodKinMomenta,  &prodKinMomentaBr );
			inTrees[0]->SetBranchAddress(decayKinMomentaLeafName.c_str(), &decayKinMomenta, &decayKinMomentaBr);

			// loop over events
			if (not topo->initKinematicsData(*prodKinPartNames, *decayKinPartNames)) {
				printErr << "problems initializing input data. Aborting..." << endl;
				exit(1);
			}
			const long int   nmbEvents = ((maxNmbEvents > 0) ? min(maxNmbEvents, nmbEventsChain)
			                              : nmbEventsChain);
			progress_display progressIndicator(nmbEvents);
			timer.Reset();
			timer.Start();
			for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
				++progressIndicator;

				if (inTrees[0]->LoadTree(eventIndex) < 0)
					break;
				// read only required branches
				prodKinMomentaBr->GetEntry (eventIndex);
				decayKinMomentaBr->GetEntry(eventIndex);

				if (not prodKinMomenta or not decayKinMomenta) {
					printWarn << "at least one data array is null pointer: "
					          << "prodKinMomenta = "    << prodKinMomenta    << ", "
					          << "decayKinMomenta = "   << decayKinMomenta   << ". "
					          << "skipping event." << endl;
					continue;
				}

				if (topo->readKinematicsData(*prodKinMomenta, *decayKinMomenta)) {
					myAmps.push_back((*amp)());
					if ((myAmps.back().real() == 0) or (myAmps.back().imag() == 0))
						printWarn << "event " << eventIndex << ": " << myAmps.back() << endl;
					topo->productionVertex()->productionAmp();
				}
			} // event loop

			timer.Stop();
			printSucc << "read " << myAmps.size() << " events from file(s) "
			          << "'" << rootInFileName << "' and calculated amplitudes" << endl;
			cout << "needed ";
			printInfo << "myAmps[0] = " << maxPrecisionDouble(myAmps[0]) << endl;
			timer.Print();

		}  // parsing of key file successful

	}
}
