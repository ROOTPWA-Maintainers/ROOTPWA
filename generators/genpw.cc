///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
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
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////


/** @brief Simple partial wave event generator (homogeneous in m)
 */


#include <complex>
#include <iostream>
//#include <stringstream>

#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include "TPWWeight.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TDiffractivePhaseSpace.h"
#include <event.h>

using namespace std;

extern particleDataTable PDGtable;

void printUsage(char* prog, int errCode=0) {
cerr << "usage:" << endl
       << prog
       << " -n # -t <theta distribution>"
       << "    where:" << endl
       << "        -n #       number of events to generate" << endl
     << "        -t <file>  ROOT file containing histogram with theta distribution "<< endl
       << endl;
  exit(errCode);
}


int main(int argc, char** argv, const int     errCode = 0)
{

  unsigned int nevents=100;
  string theta_file;

  int c;
  while ((c = getopt(argc, argv, "n:t:h")) != -1)
    switch (c) {
    case 'n':
      nevents = atoi(optarg);
      break;
    case 't':
      theta_file = optarg;
      break;
    case 'h':
      printUsage(argv[0]);
      break;
    }

 
  TFile* outfile=TFile::Open("genhisto.root","RECREATE");
  TH1D* hWeights=new TH1D("hWeights","PW Weights",100,0,100);
  TTree* outtree=new TTree("pwevents","pwevents");
  double weight;
  TClonesArray* p=new TClonesArray("TLorentzVector");
  TLorentzVector beam;
  outtree->Branch("weight",&weight,"weight/d");
  outtree->Branch("p",&p);
  outtree->Branch("beam",&beam);

  PDGtable.initialize();
  TPWWeight weighter;
  weighter.addWave("keyfiles/key3pi/SET1/1-1++0+rho770_01_pi-.key",
		   std::complex<double>(15.23,0),0);
  weighter.addWave("keyfiles/key3pi/SET1/1-2++1+rho770_21_pi-.key",
		   std::complex<double>(12.23,3.4),0);
  weighter.loadIntegrals("src/pwafitTest/amplitudes/norm.int");

  TDiffractivePhaseSpace difPS;
  difPS.SetSeed(1236735);
  difPS.SetBeam();
  difPS.SetTarget(-300,0.2,2);
  difPS.SetMassRange(1.2,1.4);			
  TFile* infile=TFile::Open(theta_file.c_str());
  difPS.SetThetaDistribution((TH1*)infile->Get("h1"));
  const double mpi=0.13957018;
  difPS.AddDecayProduct(particleinfo(9,-1,mpi));  
  difPS.AddDecayProduct(particleinfo(8,1,mpi)); 
  difPS.AddDecayProduct(particleinfo(9,-1,mpi));

  double maxweight=-1;
  unsigned int acc=0;

  for(unsigned int i=0;i<nevents;++i)
    {

      p->Delete(); // clear output array

      ofstream str("/tmp/event.evt");
      difPS.event(str);
      str.close();
      
      for(unsigned int ip=0; ip<3;++ip){
	new((*p)[ip]) TLorentzVector(*difPS.GetDecay(ip));
      }

      beam=*difPS.GetBeam();

      // calculate weight
      event e;
      e.setIOVersion(1);
      
      
      ifstream istr("/tmp/event.evt");
      istr >> e;
      istr.close();

      // cerr << e <<endl;
      
      weight=weighter.weight(e);
      if(weight>maxweight)maxweight=weight;

      hWeights->Fill(weight);
      //cerr << i << endl;
      
      outtree->Fill();
      
      

    }

  cerr << "Maxweight: " << maxweight << endl;
  cerr << "Accepted Events: " << acc << endl;
  

  outfile->cd();
  hWeights->Write();
outtree->Write();
  outfile->Close();
  infile->Close();

  return 0;

}

