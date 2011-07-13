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
//      basic test program for amplitude calculation on CPU and GPU
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


// #include <fstream>

#include <boost/progress.hpp>

// #include "TVector3.h"
// #include "TLorentzRotation.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
// #include "TFile.h"
// #include "TH1.h"
// #include "TH2.h"
// #include "TSystem.h"

// #include "mathUtils.hpp"
// #include "reportingUtilsRoot.hpp"
// #include "conversionUtils.hpp"
#include "particleDataTable.h"
// #include "diffractiveDissVertex.h"
// #include "massDependence.h"
#include "waveDescription.h"
#include "isobarAmplitude.h"
#include "isobarHelicityAmplitude.h"
#include "helampltestInterface.cuh"
// #include "isobarCanonicalAmplitude.h"
// #include "evtTreeHelper.h"


using namespace std;
using namespace boost;
using namespace rpwa;
using namespace cupwa;


int
main(int argc, char** argv)
{
  printCompilerInfo();
  printSvnVersion();
	
  rpwa::particleDataTable& pdt = rpwa::particleDataTable::instance();
  pdt.readFile("../amplitude/particleDataTable.txt");
  TStopwatch timer;
  
  // isobarDecayVertex::setDebug(true);
  // decayTopology::setDebug(true);
  // isobarDecayTopology::setDebug(true);
  massDependence::setDebug(true);
  // diffractiveDissVertex::setDebug(true);
  isobarAmplitude::setDebug(true);
  isobarHelicityAmplitude::setDebug(true);
  // isobarCanonicalAmplitude::setDebug(true);
  // waveDescription::setDebug(true);

  if (1) {
    const long int maxNmbEvents   = 1;
    const string   keyFileName    = "./1-2-+1+f21270_22_pi-.key";
    const string   rootInFileName = "./1260.1300.root";

    waveDescription    waveDesc;
    isobarAmplitudePtr amp;
    if (waveDesc.parseKeyFile(keyFileName) and waveDesc.constructAmplitude(amp)) {
      isobarDecayTopologyPtr topo = amp->decayTopology();
      printInfo << *amp;
			
      // read data from tree
      const string&            inTreeName                = "rootPwaEvtTree";
      const string&            prodKinParticlesLeafName  = "prodKinParticles";
      const string&            prodKinMomentaLeafName    = "prodKinMomenta";
      const string&            decayKinParticlesLeafName = "decayKinParticles";
      const string&            decayKinMomentaLeafName   = "decayKinMomenta";
      vector<complex<double> > cpuAmps;
      // open input file
      printInfo << "opening input file(s) '" << rootInFileName << "'" << endl;
      TChain chain(inTreeName.c_str());
      if (chain.Add(rootInFileName.c_str()) < 1) {
	printWarn << "no events in input file(s) '" << rootInFileName << "'" << endl;
	return false;
      }
      const long int nmbEventsChain = chain.GetEntries();
      chain.GetListOfFiles()->ls();

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
      chain.SetBranchAddress(prodKinParticlesLeafName.c_str(),  &prodKinParticles,  &prodKinParticlesBr );
      chain.SetBranchAddress(prodKinMomentaLeafName.c_str(),    &prodKinMomenta,    &prodKinMomentaBr   );
      chain.SetBranchAddress(decayKinParticlesLeafName.c_str(), &decayKinParticles, &decayKinParticlesBr);
      chain.SetBranchAddress(decayKinMomentaLeafName.c_str(),   &decayKinMomenta,   &decayKinMomentaBr  );

      // loop over events
      const long int   nmbEvents = ((maxNmbEvents > 0) ? min(maxNmbEvents, nmbEventsChain)
				    : nmbEventsChain);
      progress_display progressIndicator(nmbEvents);
      timer.Reset();
      timer.Start();
      for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
	++progressIndicator;
    
	if (chain.LoadTree(eventIndex) < 0)
	  break;
	// read only required branches
	prodKinParticlesBr->GetEntry (eventIndex);
	prodKinMomentaBr->GetEntry   (eventIndex);
	decayKinParticlesBr->GetEntry(eventIndex);
	decayKinMomentaBr->GetEntry  (eventIndex);
	      
	if (   not prodKinParticles  or not prodKinMomenta
	       or not decayKinParticles or not decayKinMomenta) {
	  printWarn << "at least one data array is null pointer: "
		    << "prodKinParticles = "  << prodKinParticles  << ", "
		    << "prodKinMomenta = "    << prodKinMomenta    << ", "
		    << "decayKinParticles = " << decayKinParticles << ", "
		    << "decayKinMomenta = "   << decayKinMomenta   << ". "
		    << "skipping event." << endl;
	  continue;
	}
	      
	if (topo->readData(*prodKinParticles,  *prodKinMomenta,
			   *decayKinParticles, *decayKinMomenta)) {
	  cpuAmps.push_back((*amp)());
	  
	  {
	    vector<isobarDecayVertexPtr> vertices = topo->isobarDecayVertices();
	    double theta1 	= vertices[0]->daughter1()->lzVec().Theta();
	    double phi1		= vertices[0]->daughter1()->lzVec().Phi();
	    double theta2 	= vertices[1]->daughter1()->lzVec().Theta();
	    double phi2		= vertices[1]->daughter1()->lzVec().Phi();
	    int    JX		= vertices[0]->parent()->J();
	    int    J1		= vertices[0]->daughter1()->J();
	    int    J2		= vertices[0]->daughter2()->J();
	    int    S1		= vertices[0]->S();
	    int    S2		= vertices[1]->S();
	    int    L1		= vertices[0]->L();
	    int    L2		= vertices[1]->L();
//	    int    lambda1	= vertices[0]->daughter1()->spinProj();
//	    int    lambda2	= vertices[0]->daughter2()->spinProj();	    
	    int    Lambda	= vertices[0]->parent()->spinProj();		// in my scope MX
//	    int    P		= vertices[0]->parent()->P();
//	    int    refl		= vertices[0]->parent()->reflectivity();
	    double wf2		= vertices[1]->parent()->mass();
	    double wX		= vertices[0]->parent()->lzVec().M();
// // 	    for (unsigned int i = 0; i < vertices.size(); ++i) 
// // 	    {
// // 		cout << "\n" << i << ". parent width: " << maxPrecision(vertices[i]->parent()->width()) << endl;
// // 		cout << "\n" << i << ". daughter1 mass: " << vertices[i]->daughter1()->mass() << endl;
// // 		cout << "\n" << i << ". daughter2 mass: " << vertices[i]->daughter2()->mass() << endl;
// //		cout << "\n" << i << ". : " << vertices[i]->daughter2()->mass() << endl;
// // 		cout << "theta1 = " << theta1 << "     theta2 = " << theta2 << "  L1 = " << L1 << "   L2 = " << L2 << "  S1 = " << S1 << "   S2 = " << S2 << endl;
// // 		cout << "5! = " << cufac(5) << endl;
// // 	    }

// //	    cout << "address of theta1 = " << &theta1 << endl;
	    GPUAmpcall(theta1,phi1,theta2,phi2,wX,wf2,JX,Lambda,J1,L1,S1,J2,L2,S2);
//	    cout << "\n\n\n TEST1: " << test1 << endl;
	  }

	  if ((cpuAmps.back().real() == 0) or (cpuAmps.back().imag() == 0))
	    printWarn << "event " << eventIndex << ": " << cpuAmps.back() << endl;
	}
      } // event loop
      
      timer.Stop();
      printInfo << "successfully read " << cpuAmps.size() << " events from file(s) "
		<< "'" << rootInFileName << "' and calculated amplitudes" << endl;
      cout << "needed ";
      printInfo << "cpuAmps[0] = " << maxPrecisionDouble(cpuAmps[0]) << endl;
      timer.Print();
    }  // parsing of key file successful
    
  }
}
