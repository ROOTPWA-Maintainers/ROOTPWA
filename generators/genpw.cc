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
#include "TFitBin.h"
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
     << " -n # -o <file> -w <file> -k <path> -i <file> -r <file>"
     << "    where:" << endl
     << "        -n #       number of events to generate" << endl
     << "        -o <file>  ROOT output file"<< endl
     << "        -w <file>  wavelist file (contains production amplitudes)"<< endl 
     << "        -w <file.root>  to use TFitBin tree as input"<< endl 
     << "        -k <path>  path to keyfile directory (all keyfiles have to be there)"<< endl 
     << "        -i <file>  integral file"<< endl 
     << "        -r <file>  reaction config file"<< endl   
     << endl;
 exit(errCode);
}


int main(int argc, char** argv, const int     errCode = 0)
{

  unsigned int nevents=100;
  string output_file("genpw.root");
  string integrals_file;
  string wavelist_file; // format: name Re Im
  string path_to_keyfiles("./");
  string reactionFile;

  int c;
  while ((c = getopt(argc, argv, "n:o:w:k:i:r:h")) != -1)
    switch (c) {
    case 'n':
      nevents = atoi(optarg);
      break;
    case 'o':
      output_file = optarg;
      break;
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

 

  TFile* outfile=TFile::Open(output_file.c_str(),"RECREATE");
  TH1D* hWeights=new TH1D("hWeights","PW Weights",100,0,100);
  TTree* outtree=new TTree("pwevents","pwevents");
  double weight;
  TClonesArray* p=new TClonesArray("TLorentzVector");
  TLorentzVector beam;
  double qbeam;
  std::vector<int> q; // array of charges

  outtree->Branch("weight",&weight,"weight/d");
  outtree->Branch("p",&p);
  outtree->Branch("beam",&beam);
  outtree->Branch("q",&q);
  outtree->Branch("qbeam",&qbeam,"qbeam/i");

  PDGtable.initialize();

  TPWWeight weighter;
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
  double   binCenter    = 500 * (mmin + mmax);

  if(!reactConf.lookupValue("beam.charge",qbeam))qbeam=-1;

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
    int myq;part.lookupValue("charge",myq);
    double m;part.lookupValue("mass",m);
    q.push_back(myq);
    difPS.AddDecayProduct(particleinfo(id,myq,m));
  }



  // check if TFitBin is used as input
  if(wavelist_file.find(".root")!=string::npos){
    cerr << "Using TFitBin as input!" << endl;
    TFile* fitresults=TFile::Open(wavelist_file.c_str(),"READ");
    TFitBin* Bin     = NULL;
    if (!fitresults || fitresults->IsZombie()){
      cerr << "Cannot open start fit results file " << wavelist_file << endl;
      return 1;
    }
    // get tree with start values
    TTree* tree;
    fitresults->GetObject("pwa", tree);
      if (!tree)
        cerr << "Cannot find fitbin tree '"<< "pwa" << "' "<< endl;
      else {
        Bin = new TFitBin();
        tree->SetBranchAddress("fitbin", &Bin);
	// find entry which is closest to mass bin center
        unsigned int iBest = 0;
        double mBest = 0;
        for (unsigned int i = 0; i < tree->GetEntriesFast(); ++i) {
          tree->GetEntry(i);
          if (fabs(binCenter - Bin->mass()) <= fabs(binCenter - mBest)) {
            iBest = i;
            mBest = Bin->mass();
          }
        }  // end loop over TFitBins
	cerr << "Using data from Mass bin m=" << mBest << endl;
        tree->GetEntry(iBest); 
	// write wavelist file for generator
	string tmpname("/tmp/genamps.txt");
	ofstream tmpfile(tmpname.c_str());
	Bin->printAmpsGenPW(tmpfile);
	tmpfile.close();
	wavelist_file=tmpname;
      }
  } // end root file given


  // read input wavelist and amplitudes
  ifstream wavefile(wavelist_file.c_str());
  while(wavefile.good()){
    TString wavename;
    double RE, IM;
    int rank=0;
    int refl=0;
    wavefile >> wavename >> RE >> IM;
    
    if(wavename.Contains("flat") || wavename.Length()<2)continue;
    if(RE==0 && IM==0)continue;
    // check if there is rank information
    if(wavename(0)=='V'){
      // we multiply rank by to to make space for refl+- production vectors
      rank=2*atoi(wavename(1,1).Data());
      // check reflecitivity to sort into correct production vector
      refl = wavename(9)=='+' ? 0 : 1;
      wavename=wavename(3,wavename.Length());
    }

    wavename.ReplaceAll(".amp",".key");
    wavename.Prepend(path_to_keyfiles.c_str());

    std::complex<double> amp(RE,IM);
    cerr << wavename << " " << amp << " r=" << rank/2 
	 << " eps=" << refl << endl;
    wavefile.ignore(256,'\n');

    // production vector index: rank+refl
    weighter.addWave(wavename.Data(),amp,rank+refl);
    
  }



  weighter.loadIntegrals(integrals_file);



  
  double maxweight=-1;
  unsigned int acc=0;

  unsigned int tenpercent=(unsigned int)(nevents/10);

  for(unsigned int i=0;i<nevents;++i)
    {
      if(i>0 && ( i % tenpercent==0) )cerr << "[" << (double)i/(double)nevents*100. << "%]";
					       
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

    } // end event loop

  cerr << endl;
  cerr << "Maxweight: " << maxweight << endl;
  cerr << "Accepted Events: " << acc << endl;
  

  outfile->cd();
  hWeights->Write();
outtree->Write();
  outfile->Close();
  infile->Close();

  return 0;

}

