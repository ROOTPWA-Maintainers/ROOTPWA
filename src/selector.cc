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
#include <iomanip>
#include "TString.h"
#include "TChain.h"
#include "TFitBin.h"
using namespace std;



int
main(int argc, char** argv){

  if(argc<4){
    cerr<<"Usage: selector nbins nsurvivors fitdir1 fitdir2 fitdir3 ..."<<endl;
    return 1;
  }

  unsigned int nbins=atoi(argv[1]);
  unsigned int nsurv=atoi(argv[2]);

  vector<TString> inputdirectories;
  for(int i=3; i<argc; ++i){
    inputdirectories.push_back(argv[i]);
  }
  
  if(nsurv>inputdirectories.size()){
    cerr << "Not enough fits to create " << nsurv << " survivors." << endl;
  }

  map<double,unsigned int> results; // <logli,index>
  unsigned int bestfit=0;
  double bestLogli=0;
  TFitBin* bin=new TFitBin;
  //loop over fits and extract quality information
  // we are using the sum of loglikelyhood per event for this
  for(unsigned int j=0; j<inputdirectories.size(); ++j){
    cerr << "Examining "<<inputdirectories[j]<<endl;
    TChain* chain=new TChain("pwa");
    TString f=inputdirectories[j];
    f+="/*.root";
    if(chain->Add(f)==0){
      cerr << "No fitoutput files found." << nbins 
	   << ". Skipping fit!" << endl;
      delete chain;
      continue;
    }
    chain->SetBranchAddress("fitbin",&bin);
    unsigned int n=chain->GetEntries();
    if(n!=nbins){
      cerr << n << " bins in this fit. Expected " << nbins 
	   << ". Skipping fit!" << endl;
      delete chain;
      continue;
    }
    double sumlogli=0;
    for(unsigned int k=0;k<n;++k){
      chain->GetEntry(k);
      sumlogli-=bin->logli()/bin->rawEvents();
    }// end loop over bins
    cerr<<"SumLogli="<<setprecision(9)<<sumlogli<<endl;

    results[sumlogli]=j;    
    if(sumlogli>bestLogli){
      bestLogli=sumlogli;
      bestfit=j;
    }
    delete chain;
  } // end loop over fits
  

  map<double,unsigned int>::iterator mit=results.begin();
  while(mit!=results.end()){
    cerr <<inputdirectories[mit->second]<<": "<<mit->first<<endl;
    ++mit;
  }

  cerr <<  "Bestfit: "<<inputdirectories[bestfit]
       <<"  with loglikely/event: "<<bestLogli<<endl;

  map<double,unsigned int>::iterator mrit=results.end();
  if(nsurv>results.size())nsurv=results.size();
  for(unsigned int i=0;i<nsurv;++i){
    --mrit;
    cout << inputdirectories[mrit->second]<<"/wavelist "<<mrit->first << endl;
  }
  
  return 0;
}
