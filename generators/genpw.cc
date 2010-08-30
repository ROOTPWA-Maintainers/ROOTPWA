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


#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cstdlib>
#include <string>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include <event.h>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"
//#include "TH1D.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TClonesArray.h"

#include "libconfig.h++"

#include "utilities.h"
#include "TDiffractivePhaseSpace.h"
#include "TPWWeight.h"
#include "TProductionAmp.h"
#include "TBWProductionAmp.h"
#include "TFitBin.h"

#include "TPrimaryVertexGen.h"


using namespace std;
using namespace libconfig;
using namespace rpwa;


extern particleDataTable PDGtable;


void printUsage(char* prog,
		int   errCode = 0)
{
cerr << "usage:" << endl
     << prog
     << " -n # [-a # -m # -M # -B #] -o <file> -w <file> -k <path> -i <file> -r <file>" << endl
     << "    where:" << endl
     << "        -n #       (max) number of events to generate (default: 100)" << endl
     << "        -a #       (max) number of attempts to do (default: infinity)" \
     << endl
     << "        -m #       maxWeight" << endl
     << "        -o <file>  ROOT output file (if not specified, generated automatically)"<< endl
     << "        -w <file>  wavelist file (contains production amplitudes)"<< endl 
     << "        -w <file.root>  to use TFitBin tree as input"<< endl 
     << "        -c <0/1>   if 1 a comgeant eventfile (.fort.26) is written with same naming as the root file (default 1)" << endl
     << "        -k <path>  path to keyfile directory (all keyfiles have to be there)"<< endl 
     << "        -i <file>  integral file"<< endl 
     << "        -r <file>  reaction config file"<< endl
     << "        -M #   lower boundary of mass range in MeV (overwrites values from config file) " << endl
     << "        -B #   width of mass bin in MeV" << endl
     << endl;
 exit(errCode);
}


int main(int argc, char** argv)
{

  unsigned int nevents=100;
  unsigned int max_attempts=0;
  string output_file(""); // either given by option or generated automatically by mass range
  string output_evt("");
  string output_wht("");
  string output_comgeant("");
  string integrals_file;
  bool hasint=false;
  string wavelist_file; // format: name Re Im
  string path_to_keyfiles("./");
  string reactionFile;
  double maxWeight=0;
  int seed=123456;
  int massLower=0;
  int massBinWidth=0;
  bool overwriteMass=false;
  bool writeComGeantout=false;

  int c;
  while ((c = getopt(argc, argv, "n:a:o:w:k:i:r:m:s:M:B:h:c")) != -1)
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
   case 'c':
	  writeComGeantout = true;
	  break;
   case 'M':
      massLower = atoi(optarg);
      overwriteMass=true;
      break;
   case 'B':
      massBinWidth = atoi(optarg);
      overwriteMass=true;
      break;

    case 'h':
      printUsage(argv[0]);
      break;
    }

  PDGtable.initialize();

  TPWWeight weighter;
  Config reactConf;
  reactConf.readFile(reactionFile.c_str());

  // variable that need to get initialized either by input options
  // or the config file
  // will be stored in the tree later
  double weight, impweight;
  TClonesArray* p=new TClonesArray("TLorentzVector");
  TLorentzVector beam;
  TVector3 vertex;
  double tprime(0.);
  double qbeam;
  vector<int> q; // array of charges

  double Mom=reactConf.lookup("beam.momentum");
  double MomSigma=reactConf.lookup("beam.sigma_momentum");
  double BeamPartMass(0.13957018);
  if (reactConf.exists("beam.mass")){
	  BeamPartMass = reactConf.lookup("beam.mass");
  }
  double DxDz=reactConf.lookup("beam.DxDz");
  double DxDzSigma=reactConf.lookup("beam.sigma_DxDz");
  double DyDz=reactConf.lookup("beam.DyDz");
  double DyDzSigma=reactConf.lookup("beam.sigma_DyDz");

  double targetz=reactConf.lookup("target.pos.z");
  double targetd=reactConf.lookup("target.length");
  double targetr=reactConf.lookup("target.radius");
  double mrecoil=reactConf.lookup("target.mrecoil");

  double tprime_min(0.); 
  reactConf.lookupValue("finalstate.t_min", tprime_min);

  double mmin= reactConf.lookup("finalstate.mass_min");
  double mmax= reactConf.lookup("finalstate.mass_max");
  if(overwriteMass){
    mmin=massLower/1000.0;
    mmax=(massLower+massBinWidth)/1000.0;
  }
  // array of tslopes even when only one is existing
  double* tslope = NULL;
  double* inv_m  = NULL;
  int ntslope = 1;
  if (reactConf.lookup("finalstate.t_slope").isArray()){
	  ntslope = reactConf.lookup("finalstate.t_slope").getLength();
	  if (reactConf.lookup("finalstate.inv_m").getLength()!=ntslope){
		  cout << " Error: please check number of t' values and the corresponding invariant masses in the Configuration File! " << endl;
		  return 0;
	  }
	  tslope = new double[ntslope];
	  inv_m  = new double[ntslope];
	  cout << " found array of t' slopes. Reading " << ntslope << " of values ";
	  for (int i = 0; i < ntslope; i++){
		  tslope[i] = reactConf.lookup("finalstate.t_slope")[i];
		  inv_m[i]  = reactConf.lookup("finalstate.inv_m")[i];
		  //cout << inv_m[i] << " " << tslope[i] << endl;
	  }
	  cout << " done. " << endl;
  } else {
	  tslope = new double[1];
	  tslope[0]=reactConf.lookup("finalstate.t_slope");
	  //cout << " tslope set to " << tslope[0];
  }
  double binCenter=500 * (mmin + mmax);

  // check weather to use a primary vertex generator as requested by the config file
  TPrimaryVertexGen* primaryVertexGen(NULL);
  string histfilename_primvertex("");
  if (reactConf.lookupValue("primvertex.histfilename", histfilename_primvertex)){
  	primaryVertexGen = new TPrimaryVertexGen(
  			histfilename_primvertex,
  			BeamPartMass,
  			Mom,
  			MomSigma
  			);
  	if (!primaryVertexGen->Check()){
  		cerr << " Error: histogram filename with beam properties not loaded! " << endl;
  		delete primaryVertexGen;
  		primaryVertexGen = NULL;
  	}
  }

  if(!reactConf.lookupValue("beam.charge",qbeam))qbeam=-1;
  // generate the filename automatically if not specified
  if (output_file == "") {
	  stringstream _filename;
	  _filename << massLower << "." << massLower+massBinWidth << ".genbod.root";
	  cout << " created output filename: " << _filename.str() << endl;
	  output_file = _filename.str();
  }
  output_evt = output_file;
  output_evt.replace(output_evt.find(".root"),5,".evt");
  output_wht = output_file;
  output_wht.replace(output_wht.find(".root"),5,".wht");
  output_comgeant = output_file;
  output_comgeant.erase(output_comgeant.find(".root"),5);
  output_comgeant.append(".fort.26");

  // now create the root file to store the events
  TFile* outfile=TFile::Open(output_file.c_str(),"RECREATE");
  TH1D* hWeights=new TH1D("hWeights","PW Weights",100,0,100);
  TTree* outtree=new TTree("pwevents","pwevents");

  outtree->Branch("weight",&weight,"weight/d");
  outtree->Branch("impweight",&impweight,"impweight/d");
  outtree->Branch("p",&p);
  outtree->Branch("beam",&beam);
  outtree->Branch("vertex", &vertex);
  outtree->Branch("q",&q);
  outtree->Branch("qbeam",&qbeam,"qbeam/I");
  outtree->Branch("tprime", &tprime);

  //string theta_file= reactConf.lookup("finalstate.theta_file");

  TDiffractivePhaseSpace difPS;
  cerr << "Seed=" << seed << endl;
  difPS.SetSeed(seed);
  difPS.SetBeam(Mom,MomSigma,DxDz,DxDzSigma,DyDz,DyDzSigma);
  difPS.SetTarget(targetz,targetd,targetr,mrecoil);
  difPS.SetTPrimeSlope(tslope, inv_m, ntslope);
  difPS.SetMassRange(mmin,mmax);
  difPS.SetPrimaryVertexGen(primaryVertexGen);
  if (tprime_min >= 0)  
	difPS.SettMin(tprime_min);
  else
	cout << " Error: tprime_min must be positive " << endl;


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
    difPS.AddDecayProduct(particleInfo(id,myq,m));
  }

  // see if we have a resonance in this wave
  map<string, TBWProductionAmp*> bwAmps;
  if(reactConf.exists("resonances")){
    const Setting &bws = root["resonances"]["breitwigners"];
    // loop through breitwigners
    int nbw=bws.getLength();
    printInfo << "found " << nbw << " Breit-Wigners in config" << endl;
    
    for(int ibw = 0; ibw < nbw; ++ibw) {
      const Setting &bw = bws[ibw];
      string jpcme;
      double mass, width;
      double cRe, cIm;
      bw.lookupValue("jpcme",       jpcme);
      bw.lookupValue("mass",        mass);
      bw.lookupValue("width",       width);
      bw.lookupValue("coupling_Re", cRe);
      bw.lookupValue("coupling_Im", cIm);
      complex<double> coupl(cRe, cIm);
      cout << "    JPCME = " << jpcme << ", mass = " << mass << " GeV/c^2, "
	   << "width = " << width << " GeV/c^2, coupling = " << coupl << endl;
      bwAmps[jpcme] = new TBWProductionAmp(mass, width, coupl);
    }
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

    complex<double> amp(RE,IM);
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
      weighter.addWave(wavename.Data(),pamp,complex<double>(1,0),rank+refl);
    }   
  }

  if(hasint){
    // read integral files
    ifstream intfile(integrals_file.c_str());
    while(intfile.good()){
      string filename;
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
  ofstream evtout(output_evt.c_str());
  ofstream evtgeant;
  if (writeComGeantout)
	  evtgeant.open(output_comgeant.c_str());
  ofstream evtwht(output_wht.c_str());
  evtwht << setprecision(10);

  while(i<nevents && ((max_attempts>0 && attempts<max_attempts) || max_attempts==0))
    {

      ++attempts;

      
      p->Delete(); // clear output array

      ofstream str("tmpevent.evt");
      if (writeComGeantout)
    	  difPS.event(str, evtgeant);
      else
    	  difPS.event(str);
      impweight=difPS.impWeight();
      str.close();
      
      for(int ip=0; ip<nparticles;++ip){
	new((*p)[ip]) TLorentzVector(*difPS.GetDecay(ip));
      }

      beam=*difPS.GetBeam();
      vertex=*difPS.GetVertex();
      tprime=difPS.Gettprime();

      // calculate weight
      event e;
      e.setIOVersion(1);
      
      
      ifstream istr("tmpevent.evt");
      istr >> e;
      evtout << e;
      evtwht << impweight << endl;
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
      

      if(i>0 && ( i % tenpercent==0) )cerr << "\x1B[2K" << "\x1B[0E" << "[" << (double)i/(double)nevents*100. << "%]";

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
  evtout.close();
  evtwht.close();
  if (writeComGeantout)
	  evtgeant.close();

  delete [] tslope;
  delete [] inv_m;
 
  return 0;
}

