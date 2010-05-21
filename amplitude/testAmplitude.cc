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
//      basic test program for amplitude classes
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <fstream>

#include "TVector3.h"
#include "TLorentzRotation.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include "Vec.h"
#include "lorentz.h"
#include "keyfile.h"
#include "event.h"
#include "wave.h"

#include "svnVersion.h"
#include "utilities.h"
#include "particleDataTable.h"
#include "diffractiveDissVertex.h"
#include "massDependence.h"
#include "isobarHelicityAmplitude.h"
#include "keyFileParser.h"


extern particleDataTable PDGtable;
extern wave              gWave;


using namespace std;
using namespace rpwa;


int
main(int argc, char** argv)
{
  printCompilerInfo();
  printSvnVersion();

  rpwa::particleDataTable& pdt = rpwa::particleDataTable::instance();
  pdt.readFile();
  TStopwatch timer;
  
  // isobarDecayVertex::setDebug(true);
  // decayTopology::setDebug(true);
  // isobarDecayTopology::setDebug(true);
  // isobarHelicityAmplitude::setDebug(true);
  // flatMassDependence::setDebug(true);
  // relativisticBreitWigner::setDebug(true);
  // diffractiveDissVertex::setDebug(true);
  keyFileParser::setDebug(true);

  if (0) {
    {
      fourVec  p(2, threeVec(0.5, 0.75, 1));
      threeVec n = threeVec(0, 0, 1) / p.V();
      cout << "before = " << n << "    " << p << endl;
      rotation         R;
      lorentzTransform L1;
      L1.set(R.set(n.phi(), n.theta() - M_PI / 2, -M_PI / 2));
      n *= R;
      p *= L1;
      cout << "L1 -> " << n << "    " << p << endl;
      lorentzTransform L2;
      L2.set(R.set(0, signof(p.x()) * p.theta(), 0));
      p *= L2;
      cout << "L2 -> " << p << endl;
      lorentzTransform L3;
      L3.set(p);
      p *= L3;
      cout << "L3 -> " << p << endl;

      matrix<double> X(4, 4);
      X = L3 * (L2 * L1);
      lorentzTransform L(X);
      p = fourVec(2, threeVec(0.5, 0.75, 1));
      p *= L;
      cout << "L -> " << p << endl;
    }

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
    particlePtr a1    = createParticle("a1(1269)+");
    particlePtr f1    = createParticle("f1(1285)");
    // define X-system
    //                                   2I  G  2J  P   C  2M
    particlePtr X = createParticle("X-", 2, -1, 4, +1, +1, 2);
    // define production vertex
    particlePtr              beam     = createParticle("pi-");
    diffractiveDissVertexPtr prodVert = createDiffractiveDissVertex(beam, X);
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
    isobarDecayTopology topo(prodVert, decayVertices, fsParticles);
    // topo.checkTopology();
    // topo.checkConsistency();
    isobarHelicityAmplitude amp(topo);
    cout << topo;
    complex<double>         decayAmp = amp.amplitude();
    cout << "!!!< decay amplitude = " << decayAmp << endl;

    if (1) {  // compare to PWA2000
      PDGtable.initialize();
      event    ev;
      ifstream eventData("testEvents.evt");
      ev.setIOVersion(1);
      if (!(eventData >> ev).eof()) {
	keyfile key;
	key.open("1-2++1+pi-_11_f11285=pi-_11_a11269=pi+_1_sigma.key");
	complex<double> pwa2000amp;
	key.run(ev, pwa2000amp, true);
	key.rewind();
	key.close();
	cout << "!!! PWA2000 amplitude = " << pwa2000amp << endl;
	cout << "!!! my amplitude = " << decayAmp << " vs. PWA2000 amplitude = " << pwa2000amp << ", "
	     << "delta = " << decayAmp - pwa2000amp << endl;
      }
    }
  }

  if (1) {
    keyFileParser&         parser = keyFileParser::instance();
    isobarDecayTopologyPtr topo;
    if (parser.parse("test.key", topo)) {
      printInfo << *topo;
      topo->checkTopology();
      topo->checkConsistency();
      isobarHelicityAmplitude amp(*topo);
      parser.setAmplitudeOptions(amp);

      // read data from tree
      const string&            inFileNamePattern     = "testEvents.root";
      const long int           maxNmbEvents          = -1;
      const string&            inTreeName            = "rootPwaEvtTree";
      const string&            leafNameIsPartNames   = "initialStateNames";
      const string&            leafNameIsPartMomenta = "initialStateMomenta";
      const string&            leafNameFsPartNames   = "finalStateNames";
      const string&            leafNameFsPartMomenta = "finalStateMomenta";
      vector<complex<double> > myAmps;
      // open input file
      printInfo << "opening input file(s) '" << inFileNamePattern << "'" << endl;
      TChain chain(inTreeName.c_str());
      if (chain.Add(inFileNamePattern.c_str()) < 1) {
	printWarn << "no events in input file(s) '" << inFileNamePattern << "'" << endl;
	return false;
      }
      const long int nmbEventsChain = chain.GetEntries();
      chain.GetListOfFiles()->ls();

      // create branch pointers and leaf variables
      TBranch*      initialStateNamesBr   = 0;
      TBranch*      initialStateMomentaBr = 0;
      TBranch*      finalStateNamesBr     = 0;
      TBranch*      finalStateMomentaBr   = 0;
      TClonesArray* initialStateNames     = 0;
      TClonesArray* initialStateMomenta   = 0;
      TClonesArray* finalStateNames       = 0;
      TClonesArray* finalStateMomenta     = 0;
	
      // connect leaf variables to tree branches
      chain.SetBranchAddress(leafNameIsPartNames.c_str  (), &initialStateNames,   &initialStateNamesBr  );
      chain.SetBranchAddress(leafNameIsPartMomenta.c_str(), &initialStateMomenta, &initialStateMomentaBr);
      chain.SetBranchAddress(leafNameFsPartNames.c_str  (), &finalStateNames,     &finalStateNamesBr    );
      chain.SetBranchAddress(leafNameFsPartMomenta.c_str(), &finalStateMomenta,   &finalStateMomentaBr  );

      // loop over events
      const long int nmbEvents = ((maxNmbEvents > 0) ? min(maxNmbEvents, nmbEventsChain)
				  : nmbEventsChain);
      timer.Reset();
      timer.Start();
      for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
	//for (long int eventIndex = 1; eventIndex < 2; ++eventIndex) {
	progressIndicator(eventIndex, nmbEvents);
	  
	if (chain.LoadTree(eventIndex) < 0)
	  break;
	// read only required branches
	initialStateNamesBr->GetEntry  (eventIndex);
	initialStateMomentaBr->GetEntry(eventIndex);
	finalStateNamesBr->GetEntry    (eventIndex);
	finalStateMomentaBr->GetEntry  (eventIndex);

	if (!initialStateNames || !initialStateMomenta || !finalStateNames || !finalStateMomenta) {
	  printWarn << "at least one data array is null pointer: "
		    << "initialStateNames = "   << initialStateNames   << ", "
		    << "initialStateMomenta = " << initialStateMomenta << ", "
		    << "finalStateNames = "     << finalStateNames     << ", "
		    << "finalStateMomenta = "   << finalStateMomenta   << ". "
		    << "skipping event." << endl;
	  continue;
	}

	if (topo->readData(*initialStateNames, *initialStateMomenta,
			   *finalStateNames,   *finalStateMomenta)) {
	  // topo->printIsParticles(cout);
	  // topo->printFsParticles(cout);
	  // complex<double> A = amp.amplitude();
	  // complex<double> B = amp.amplitude();
	  // topo->revertMomenta();
	  // complex<double> C = amp.amplitude();
	  // cout << "A = " << A << ", B = " << B << ", C = " << C << ", A - C = " << A - C << endl;
	  myAmps.push_back(amp.amplitude());
	  if ((myAmps.back().real() == 0) || (myAmps.back().imag() == 0))
	    printWarn << "event " << eventIndex << ": " << myAmps.back() << endl;
	}
      }
	
      timer.Stop();
      printInfo << "successfully read " << myAmps.size() << " events from file(s) "
		<< "'" << inFileNamePattern << "' and calculated amplitudes" << endl;
      cout << "needed ";
      timer.Print();

      vector<complex<double> > pwa2kAmps;
      if (1) {  // compare to PWA2000
	PDGtable.initialize();
	const string evtFileName = "testTree.evt";
	ifstream     eventData(evtFileName.c_str());
	keyfile      key;
	event        ev;
	key.open("1-2++1+pi-_11_f11285=pi-_11_a11269=pi+_1_sigma.key");
	ev.setIOVersion(1);
	timer.Reset();
	timer.Start();
	while (!(eventData >> ev).eof()) {
	  complex<double> pwa2kamp;
	  key.run(ev, pwa2kamp, true);
	  pwa2kAmps.push_back(pwa2kamp);
	  key.rewind();
	}
	timer.Stop();
	printInfo << "successfully read " << pwa2kAmps.size() << " events from file(s) "
		  << "'" << evtFileName << "' and calculated amplitudes" << endl;
	cout << "needed ";
	timer.Print();
	printInfo << "myAmps[0] = " << myAmps[0] << " vs. pwa2kAmps[0] = " << pwa2kAmps[0] << ", "
		  << "delta = " << myAmps[0] - pwa2kAmps[0] << endl;
	
	if (1) {
	  const string outFileName = "testDiff.root";
	  printInfo << "writing comparison plots to " << outFileName << endl;
	  TFile* f              = TFile::Open(outFileName.c_str(), "RECREATE");
	  TH1D*  hMyAmpsReal    = new TH1D("hMyAmpsReal",    "hMyAmpsReal;Event Number;#Rgothic[Amplitude]",    myAmps.size(),    -0.5, myAmps.size()    - 0.5);
	  TH1D*  hMyAmpsImag    = new TH1D("hMyAmpsImag",    "hMyAmpsImag;Event Number;#Jgothic[Amplitude]",    myAmps.size(),    -0.5, myAmps.size()    - 0.5);
	  TH1D*  hPwa2kAmpsReal = new TH1D("hPwa2kAmpsReal", "hPwa2kAmpsReal;Event Number;#Rgothic[Amplitude]", pwa2kAmps.size(), -0.5, pwa2kAmps.size() - 0.5);
	  TH1D*  hPwa2kAmpsImag = new TH1D("hPwa2kAmpsImag", "hPwa2kAmpsImag;Event Number;#Jgothic[Amplitude]", pwa2kAmps.size(), -0.5, pwa2kAmps.size() - 0.5);
	  TH1D*  hDiffReal      = new TH1D("hDiffReal", "hDiffReal;#Rgothic[Amplitude] Difference;Count", 10000, -3e-9, 3e-9);
	  TH1D*  hDiffImag      = new TH1D("hDiffImag", "hDiffImag;#Jgothic[Amplitude] Difference;Count", 10000, -3e-9, 3e-9);
	  TH2D*  hCorrReal      = new TH2D("hCorrReal", "hCorrReal;#Rgothic[My Amp];#Rgothic[PWA2000 Amp]", 1000, -2, 2, 1000, -2, 2);
	  TH2D*  hCorrImag      = new TH2D("hCorrImag", "hCorrImag;#Jgothic[My Amp];#Jgothic[PWA2000 Amp]", 1000, -2, 2, 1000, -2, 2);
	  for (unsigned int i = 0; i < myAmps.size(); ++i) {
	    hMyAmpsReal->SetBinContent   (i + 1, myAmps[i].real());
	    hMyAmpsImag->SetBinContent   (i + 1, myAmps[i].imag());
	    hPwa2kAmpsReal->SetBinContent(i + 1, pwa2kAmps[i].real());
	    hPwa2kAmpsImag->SetBinContent(i + 1, pwa2kAmps[i].imag());
	    hDiffReal->Fill(pwa2kAmps[i].real() - myAmps[i].real());
	    hDiffImag->Fill(pwa2kAmps[i].imag() - myAmps[i].imag());
	    hCorrReal->Fill(myAmps[i].real(), pwa2kAmps[i].real());
	    hCorrImag->Fill(myAmps[i].imag(), pwa2kAmps[i].imag());
	  }
	  f->Write();
	  f->Close();
	}
      }
    }
  }
}
