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
#include "TString.h"
#include "TH1.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TDiffractivePhaseSpace.h"
#include <event.h>
#include "libconfig.h++"

using namespace std;
using namespace libconfig;

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
  //string theta_file;
  string integrals_file;
  string wavelist_file; // format: name Re Im
  string path_to_keyfiles("./");
  string reactionFile;

  int c;
  while ((c = getopt(argc, argv, "n:t:w:k:i:r:h")) != -1)
    switch (c) {
    case 'n':
      nevents = atoi(optarg);
      break;
//     case 't':
//       theta_file = optarg;
//       break;
   case 'w':
      wavelist_file = optarg;
      break;
   case 'i':
      integrals_file = optarg;
      break;
   case 'k':
      path_to_keyfiles = optarg;
      break;
   case 'r':
      reactionFile = optarg;
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

  // read input wavelist and amplitudes
  ifstream wavefile(wavelist_file.c_str());
  while(wavefile.good()){
    TString wavename;
    double RE, IM;
    wavefile >> wavename >> RE >> IM;

    wavename.ReplaceAll(".amp",".key");
    wavename.Prepend(path_to_keyfiles.c_str());

    std::complex<double> amp(RE,IM);
    cout << wavename << " " << amp << endl;
    wavefile.ignore(256,'\n');

    weighter.addWave(wavename.Data(),amp,0);
    
  }



  weighter.loadIntegrals(integrals_file);


  Config reactConf;
  reactConf.readFile(reactionFile.c_str());

  double Mom=reactConf.lookup("beam.momentum");
  double MomSigma=reactConf.lookup("beam.sigma_momentum");
  double DxDz=reactConf.lookup("beam.DxDz");
  double DxDzSigma=reactConf.lookup("beam.sigma_DxDz");
  double DyDz=reactConf.lookup("beam.DyDz");
  double DyDzSigma=reactConf.lookup("beam.sigma_DyDz");

  double targetz=reactConf.lookup("target.pos.z");
  double targetd=reactConf.lookup("target.length");
  double targetr=reactConf.lookup("target.radius");

  double mmin= reactConf.lookup("finalstate.mass_min");
  double mmax= reactConf.lookup("finalstate.mass_max");

  string theta_file= reactConf.lookup("finalstate.theta_file");

  TDiffractivePhaseSpace difPS;
  difPS.SetSeed(1236735);
  difPS.SetBeam(Mom,MomSigma,DxDz,DxDzSigma,DyDz,DyDzSigma);
  difPS.SetTarget(targetz,targetd,targetr);
  difPS.SetMassRange(mmin,mmax);			
  TFile* infile=TFile::Open(theta_file.c_str());
  difPS.SetThetaDistribution((TH1*)infile->Get("h1"));
  

  const Setting& root = reactConf.getRoot();
  const Setting &fspart = root["finalstate"]["particles"];
  int nparticles = fspart.getLength();

  for(int ifs=0;ifs<nparticles;++ifs){
    const Setting &part = fspart[ifs];
    int id;part.lookupValue("g3id",id);
    int q;part.lookupValue("charge",q);
    double m;part.lookupValue("mass",m);
    difPS.AddDecayProduct(particleinfo(id,q,m));
  }
  
  //difPS.AddDecayProduct(particleinfo(8,1,mpi)); 
  //difPS.AddDecayProduct(particleinfo(9,-1,mpi));

  double maxweight=-1;
  unsigned int acc=0;

  for(unsigned int i=0;i<nevents;++i)
    {

      p->Delete(); // clear output array

      ofstream str("/tmp/event.evt");
      difPS.event(str);
      str.close();
      
      for(int ip=0; ip<nparticles;++ip){
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

