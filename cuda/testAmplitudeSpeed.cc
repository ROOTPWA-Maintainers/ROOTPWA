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
	////massDependence::setDebug(true);
	// diffractiveDissVertex::setDebug(true);
	////isobarAmplitude::setDebug(true);
	////isobarHelicityAmplitude::setDebug(true);
	// isobarCanonicalAmplitude::setDebug(true);
	// waveDescription::setDebug(true);

	if (1) {
		const long int maxNmbEvents   = 100000;
		const string   keyFileName    = "./1-2-+1+f21270_22_pi-.key";
		//    const string   rootInFileName = "./1260.1300.root";
		const string   rootInFileName = "./allBins.ps.root";

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
      
      
			timespec cpuamp_beg, cpuamp_end; 
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpuamp_beg);
			if (not topo->initKinematicsData(*prodKinParticles, *decayKinParticles)) {
				printErr << "problems initializing input data. aborting." << endl;
				exit(1);
			}
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
	      
				if (topo->readKinematicsData(*prodKinMomenta, *decayKinMomenta)) {  
					cpuAmps.push_back((*amp)());

					if ((cpuAmps.back().real() == 0) or (cpuAmps.back().imag() == 0))
						printWarn << "event " << eventIndex << ": " << cpuAmps.back() << endl;
				}
			} // event loop
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpuamp_end); 
			printf("\n\n%.3f ms required for calculation on host\n", ( ( (cpuamp_end.tv_sec - cpuamp_beg.tv_sec) ) + (cpuamp_end.tv_nsec - cpuamp_beg.tv_nsec) / 1e9 ) * 1e3 ); 
      
			timer.Stop();
			printInfo << "successfully read " << cpuAmps.size() << " events from file(s) "
			          << "'" << rootInFileName << "' and calculated amplitudes" << endl;
			cout << "needed ";
			//      for(int e=0;e<nmbEvents;e++) {
			//      printInfo << "cpuAmps[" << e << "] = " << maxPrecisionDouble(cpuAmps[e]) << endl;
			//      }
			timer.Print();

			cout << "blub0\n" << endl;
      
			//GPU loop
			timespec gpuamp_beg, gpuamp_end; 
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &gpuamp_beg); 
			cout << "blub\n" << nmbEvents << endl << sizeof(vertexdata) << endl;
			vertexdata* event = new vertexdata[nmbEvents];
			//	vertexdata* event;
			DMAccess(event,nmbEvents);
			cout << "blub2" << endl;
			timespec trafo_beg, trafo_end; 
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &trafo_beg);
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

				if (topo->readKinematicsData(*prodKinMomenta, *decayKinMomenta)) {  
	  
					amp->transformDaughters();
	 	  
					vector<isobarDecayVertexPtr> vertices = topo->isobarDecayVertices();
	    	    
					event[eventIndex].theta1 = vertices[0]->daughter1()->lzVec().Theta();
					event[eventIndex].phi1	 = vertices[0]->daughter1()->lzVec().Phi();
					event[eventIndex].theta2 = vertices[1]->daughter1()->lzVec().Theta();
					event[eventIndex].phi2	 = vertices[1]->daughter1()->lzVec().Phi();
					event[eventIndex].JX		 = vertices[0]->parent()->J();
					event[eventIndex].J1		 = vertices[0]->daughter1()->J();
					event[eventIndex].J2		 = vertices[0]->daughter2()->J();
					event[eventIndex].S1		 = vertices[0]->S();
					event[eventIndex].S2		 = vertices[1]->S();
					event[eventIndex].L1		 = vertices[0]->L();
					event[eventIndex].L2		 = vertices[1]->L();  
					event[eventIndex].Lambda = vertices[0]->parent()->spinProj();		// in my scope MX
					event[eventIndex].massf2 = vertices[1]->parent()->mass();
					event[eventIndex].wX		 = vertices[0]->parent()->lzVec().M();
					event[eventIndex].wf2		 = vertices[0]->daughter1()->lzVec().M();
					event[eventIndex].wpi		 = vertices[1]->daughter1()->lzVec().M();
					event[eventIndex].qX		 = vertices[0]->daughter1()->lzVec().Vect().Mag();
					event[eventIndex].qf2		 = vertices[1]->daughter1()->lzVec().Vect().Mag();
	    
				}

			}// GPU event loop
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &trafo_end); 
			printf("\n\n%.3f ms required for assign + trafo daughters\n", ( ( (trafo_end.tv_sec - trafo_beg.tv_sec) ) + (trafo_end.tv_nsec - trafo_beg.tv_nsec) / 1e9 ) * 1e3 ); 

	
			timespec gpuampcall_beg, gpuampcall_end; 
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &gpuampcall_beg); 
			GPUAmpcall2(event,nmbEvents,cpuAmps);
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &gpuampcall_end); 
			printf("\n\n%.3f ms required for call\n", ( ( (gpuampcall_end.tv_sec - gpuampcall_beg.tv_sec) ) + (gpuampcall_end.tv_nsec - gpuampcall_beg.tv_nsec) / 1e9 ) * 1e3 ); 
 
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &gpuamp_end); 
			printf("\n\n%.3f ms required for calculation on device (all together)\n", ( ( (gpuamp_end.tv_sec - gpuamp_beg.tv_sec) ) + (gpuamp_end.tv_nsec - gpuamp_beg.tv_nsec) / 1e9 ) * 1e3 ); 
 
		} // parsing of key file successful
	}  
}

