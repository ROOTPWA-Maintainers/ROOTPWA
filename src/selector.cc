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

// Reads in the results of N fits, calculates the best fits and
// creates output file with likelihood values

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include "TString.h"
#include "TChain.h"
#include "TH2D.h"
#include "TFile.h"
#include "fitResult.h"
#include "reportingUtils.hpp"


using namespace std;
using namespace rpwa;


int
main(int argc, char** argv){

	using rpwa::cout;
	using rpwa::cerr;

  if(argc<4){
    cerr<<"Usage: selector nbins nsurvivors preselection fitdir1 fitdir2 fitdir3 ..."<<endl;
    return 1;
  }

  unsigned int nbins=atoi(argv[1]);
  unsigned int nsurv=atoi(argv[2]);
  string pre=argv[3];

  vector<TString> inputdirectories;
  
  map<double,TString> results; // <logli,index>
  unsigned int bestfit=0;
  double bestLogli=0;
  unsigned int counter=0;
  // load list with preselected waves
  ifstream prefile(argv[3]);
  cerr << "Opened Preselectionfile " << pre << endl;
  while(!prefile.eof() && prefile.good()){
    string dir;
    double likeli;
    prefile >> dir >> likeli;
    if(dir=="")continue;
    cerr << "Pre: " << dir << "  " << likeli << endl;
    inputdirectories.push_back(dir);
    results[likeli]=dir;
  }
  counter=inputdirectories.size();
  prefile.close();

  // register new fit directories
  for(int i=4; i<argc; ++i){
    inputdirectories.push_back(argv[i]);
  }

  if(nsurv>inputdirectories.size()){
    cerr << "Not enough fits to create " << nsurv << " survivors." << endl;
  }

  double minevi=1.78E6;
  double maxevi=1.84E6;
  
  //double minevi=712000;
  //double maxevi=700000;
  
  //double minevi=680000;
  //double maxevi=715000;

  unsigned int ngen=100;

  TH2D* hWavesetSize=new TH2D("hWS","Waveset sizes evolution",ngen,-0.5,(double)ngen-.5,100,0,100);
  TH2D* hEvidences=new TH2D("hEvi","Evidence evolution",ngen,-0.5,(double)ngen-.5,1000,minevi,maxevi);
   TH2D* hEviSize=new TH2D("hEviSize","Evidence vs Waveset size",100,0,100,1000,minevi,maxevi);

TH2D* hLogliSize=new TH2D("hLogliSize","LogLikelihood vs Waveset size",100,0,100,1000,minevi,minevi);

  
 
 
  
  //loop over fits and extract quality information
  // we are using the sum of loglikelyhood per event for this
  // start at position counter-1! 
  for(unsigned int j=counter; j<inputdirectories.size(); ++j){
    cerr << "Examining "<<inputdirectories[j]<<endl;
    // get generation
    unsigned int startgen=inputdirectories[j].Index("gen",0)+3;
    unsigned int endgen=inputdirectories[j].Index("/",startgen);
    TString genS=inputdirectories[j](startgen,endgen-startgen);
    cerr << "generation=" << genS << endl;
    double gen=atof(genS.Data());
    fitResult* bin=new fitResult;
    TChain* chain=new TChain("pwa");
  
    TString f=inputdirectories[j];
    f+="/*.result.root";
    cerr << "Loading " << f << endl;
    if(chain->Add(f)==0){
      cerr << "No fitoutput files found." << nbins 
	   << ". Skipping fit!" << endl;
      delete chain;
      delete bin;
      continue;
    }
    //chain->Print();
    cerr << "Set up Chain"<< endl;
    
    chain->SetBranchAddress("fitResult_v2",&bin);
     cerr << "Set Branch"<< endl;
     
    unsigned int n=chain->GetEntries();
    cerr << "Got Entries"<< endl;
    if(n!=nbins){
      cerr << n << " bins in this fit. Expected " << nbins 
	   << ". Skipping fit!" << endl;
      delete chain;
 delete bin;
      continue;
    }
    double sumlogli=0;
    double sumevi=0;
    unsigned int nwaves=0;
     cerr << "Reading data..."<< endl;
     bool allconverged=true;
     for(unsigned int k=0;k<n;++k){ // loop over bins
      chain->GetEntry(k);
      if(!bin->converged()){
	cerr<<"Fit not converged. Will skip this."<< endl;
	allconverged=false;
      }
      if(k==0)nwaves=bin->nmbWaves();
      sumevi+=bin->evidence();
      sumlogli+=-bin->logLikelihood();//+bin->nmbEvents();
    }// end loop over bins

    cerr<<"SumLogli    ="<<setprecision(9)<<sumlogli<<endl;
    cerr<<"SumEvidence ="<<setprecision(9)<<sumevi<<endl;
    if(sumevi==0 || !(sumevi>std::numeric_limits<double>::min() && sumevi < std::numeric_limits<double>::max())){
      cerr<<"Invalid value. Skipping."<< endl;
      delete chain;
      delete bin;
      continue;
    }
    
    // only register result if it all bins converged!
    if(allconverged){
      hWavesetSize->Fill(gen,nwaves);
      hEvidences->Fill(gen,sumevi);
      hEviSize->Fill(nwaves,sumevi);
      hLogliSize->Fill(nwaves,sumlogli);
      
      results[sumevi]=inputdirectories[j];    
      if(sumevi>bestLogli){
	bestLogli=sumevi;
	bestfit=j;
      }
    }
    delete chain;
    delete bin;
  } // end loop over fits
  

  map<double,TString>::iterator mit=results.begin();
  while(mit!=results.end()){
    cerr <<mit->second<<": "<<mit->first<<endl;
    ++mit;
  }

  cerr <<  "Bestfit: "<<inputdirectories[bestfit]
       <<"  with loglikely/event: "<<bestLogli<<endl;

  map<double,TString>::iterator mrit=results.end();
  if(nsurv>results.size())nsurv=results.size();
  for(unsigned int i=0;i<nsurv;++i){
    --mrit;
    cout << mrit->second;
    if(!(mrit->second).Contains("wavelist"))cout<<"/wavelist";
    cout<<" "<<maxPrecision(mrit->first) << endl;
  }
  
  TFile* outfile=TFile::Open("genetic_stats.root","RECREATE");
  hWavesetSize->Write();
  hEvidences->Write();
  hEviSize->Write();
hLogliSize->Write();
  outfile->Close();

  return 0;
}
