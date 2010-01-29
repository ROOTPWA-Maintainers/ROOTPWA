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
#include "TProductionAmp.h"
#include "TBWProductionAmp.h"
#include "TDiffractivePhaseSpace.h"
#include <event.h>
#include "libconfig.h++"

using namespace std;
using namespace libconfig;
using namespace rpwa;

extern particleDataTable PDGtable;

void printUsage(char* prog, int errCode=0) {
cerr << "usage:" << endl
     << prog
     << " -n # [-a # -m #] -o <file> -w <file> -k <path> -i <file> -r <file>" << endl
     << "    where:" << endl
     << "        -n #       (max) number of events to generate (default: 100)" << endl
     << "        -a #       (max) number of attempts to do (default: infty)" \
     << endl
     << "        -m #       maxWeight" << endl
     << "        -o <file>  ROOT output file"<< endl
     << "        -w <file>  wavelist file (contains production amplitudes)"<< endl 
     << "        -w <file.root>  to use TFitBin tree as input"<< endl 
     << "        -k <path>  path to keyfile directory (all keyfiles have to be there)"<< endl 
     << "        -i <file>  integral file"<< endl 
     << "        -r <file>  reaction config file"<< endl   
     << endl;
 exit(errCode);
}


int main(int argc, char** argv)
{

  unsigned int nevents=100;
  unsigned int max_attempts=0;
  string output_file("genpw.root");
  string integrals_file;
  bool hasint=false;
  string wavelist_file; // format: name Re Im
  string path_to_keyfiles("./");
  string reactionFile;
  double maxWeight=0;
  int seed=123456;

  int c;
  while ((c = getopt(argc, argv, "n:a:o:w:k:i:r:m:s:h")) != -1)
    switch (c) {
    case 'n':
      nevents = atoi(optarg);
      break;
    case 'a':
      max_attempts = atoi(optarg);
      break;
    case 's':
      seed = atoi(optarg);
      break;
    case 'o':
      output_file = optarg;
      break;
   case 'w':
      wavelist_file = optarg;
      break;
   case 'i':
      integrals_file = optarg;
      hasint=true;
      break;
   case 'k':
      path_to_keyfiles = optarg;
      break;
   case 'r':
      reactionFile = optarg;
      break;
   case 'm':
      maxWeight = atof(optarg);
      break;
    case 'h':
      printUsage(argv[0]);
      break;
    }
 

  gRandom->SetSeed(seed);

  TFile* outfile=TFile::Open(output_file.c_str(),"RECREATE");
  TH1D* hWeights=new TH1D("hWeights","PW Weights",100,0,100);
  TTree* outtree=new TTree("pwevents","pwevents");
  double weight, impweight;
  TClonesArray* p=new TClonesArray("TLorentzVector");
  TLorentzVector beam;
  double qbeam;
  std::vector<int> q; // array of charges

  outtree->Branch("weight",&weight,"weight/d");
  outtree->Branch("impweight",&impweight,"impweight/d");
  outtree->Branch("p",&p);
  outtree->Branch("beam",&beam);
  outtree->Branch("q",&q);
  outtree->Branch("qbeam",&qbeam,"qbeam/I");

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
  double mrecoil=reactConf.lookup("target.mrecoil");

  double mmin= reactConf.lookup("finalstate.mass_min");
  double mmax= reactConf.lookup("finalstate.mass_max");
  double tslope=reactConf.lookup("finalstate.t_slope");
  double binCenter=500 * (mmin + mmax);

 
  if(!reactConf.lookupValue("beam.charge",qbeam))qbeam=-1;

  //string theta_file= reactConf.lookup("finalstate.theta_file");

  TDiffractivePhaseSpace difPS;
  difPS.SetSeed(1236735);
  difPS.SetBeam(Mom,MomSigma,DxDz,DxDzSigma,DyDz,DyDzSigma);
  difPS.SetTarget(targetz,targetd,targetr,mrecoil);
  difPS.SetTPrimeSlope(tslope);
  difPS.SetMassRange(mmin,mmax);			


  double impMass;
  double impWidth;
 
  if(reactConf.lookupValue("importance.mass",impMass) && reactConf.lookupValue("importance.width",impWidth)){
    int act=reactConf.lookup("importance.active");
      if(act==1){
      difPS.SetImportanceBW(impMass,impWidth);
    }
  }
  

  const Setting& root = reactConf.getRoot();
  const Setting& fspart = root["finalstate"]["particles"];
  int nparticles = fspart.getLength();

  for(int ifs=0;ifs<nparticles;++ifs){
    const Setting &part = fspart[ifs];
    int id;part.lookupValue("g3id",id);
    int myq;part.lookupValue("charge",myq);
    double m;part.lookupValue("mass",m);
    q.push_back(myq);
    difPS.AddDecayProduct(particleinfo(id,myq,m));
  }
   // see if we have a resonance in this wave
  const Setting &bws = root["resonances"]["breitwigners"];
  // loop through breitwigners
  int nbw=bws.getLength();
  cerr << "Found " << nbw << " BreitWigners in Config" << endl;
  map<string,TBWProductionAmp*> bwAmps;
  for(int ibw=0;ibw<nbw;++ibw){
    const Setting &bw=bws[ibw];
    string jpcme;bw.lookupValue("jpcme",jpcme);
    double mass;bw.lookupValue("mass",mass);
    double width;bw.lookupValue("width",width);
    double cRe;bw.lookupValue("coupling_Re",cRe);
    double cIm;bw.lookupValue("coupling_Im",cIm);
    std::complex<double> coupl(cRe,cIm);
    cerr << jpcme << " m="<< mass << " w="<<width<<" c="<< coupl << endl;
    bwAmps[jpcme]=new TBWProductionAmp(mass,width,coupl);
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

    TString jpcme=wavename(2,5);

    wavename.ReplaceAll(".amp",".key");
    wavename.Prepend(path_to_keyfiles.c_str());

    std::complex<double> amp(RE,IM);
    cerr << wavename << " " << amp << " r=" << rank/2 
	 << " eps=" << refl << " qn=" << jpcme << endl;
    wavefile.ignore(256,'\n');
    
    TProductionAmp* pamp;
    if(bwAmps[jpcme.Data()]!=NULL){
      pamp=bwAmps[jpcme.Data()];
      cerr << "Using BW for " << jpcme << endl;
      // production vector index: rank+refl
      weighter.addWave(wavename.Data(),pamp,amp,rank+refl);
    }
    else {
      pamp=new TProductionAmp(amp);
      weighter.addWave(wavename.Data(),pamp,std::complex<double>(1,0),rank+refl);
    }   
  }

  if(hasint){
    // read integral files
    ifstream intfile(integrals_file.c_str());
    while(intfile.good()){
      std::string filename;
      double mass;
      intfile >> filename >> mass;
      weighter.loadIntegrals(filename,mass);
      intfile.ignore(256,'\n');
    } // loop over integralfile
  }// endif hasint

  
  double maxweight=-1;
  unsigned int attempts=0;


  unsigned int tenpercent=(unsigned int)(nevents/10);
  unsigned int i=0;
  //difPS.setVerbose(true);
  while(i<nevents && ((max_attempts>0 && attempts<max_attempts) || max_attempts==0))
    {

      ++attempts;
      
					       
      p->Delete(); // clear output array

      ofstream str("/tmp/event.evt");
      difPS.event(str);
      impweight=difPS.impWeight();
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

      if(maxWeight>0){ // do weighting
	//if(weight>maxWeight)maxWeight=weight;
	if(gRandom->Uniform()>weight/maxWeight)continue;
      }
      //cerr << i << endl;
      
      outtree->Fill();
      if(i>0 && ( i % tenpercent==0) )cerr << "[" << (double)i/(double)nevents*100. << "%]";

      ++i;
    } // end event loop

  cerr << endl;
  cerr << "Maxweight: " << maxweight << endl;
  cerr << "Attempts: " << attempts << endl;
  cerr << "Created Events: " << i << endl;
  cerr << "Efficiency: " << (double)i/(double)attempts << endl;
  

  outfile->cd();
  hWeights->Write();
  outtree->Write();
  outfile->Close();
 
  return 0;
}

