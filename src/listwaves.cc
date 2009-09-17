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
// tool to list the waves contained in a root fitlog
#include "TFile.h"
#include "TTree.h"
#include "TFitBin.h"
#include <string>
#include <iostream>
#include <map>

using namespace std;
char *progname;
int lineno = 1; // global variables needed for lex (not understood)

int main(int argc, char** argv){
  
  progname=argv[0];



  TFile* f=TFile::Open(argv[1]);
  TTree* tree=(TTree*)f->Get("pwa");

  TFitBin* bin=new TFitBin();

  tree->SetBranchAddress("fitbin",&bin);
  
  tree->GetEntry(0);
  int nw=bin->nwaves();
  std::vector<double> intens(nw,0);

  int n=tree->GetEntries();
  int usedbins=0;
  for(int i=0;i<n;++i){
    tree->GetEntry(i);
    double I=bin->intens();
    if(I>1000){
      for(int j=0;j<nw;++j){
	// calculate relative intensity
	intens[j]+=bin->intens(j)/I;
      } // end loop over waves
      ++usedbins;
    }
  }// end loop over mnass bins
 

  // sort the wavenames
  std::map<TString,double> waves;
  std::map<double, TString> waves2;


  for(int j=0;j<nw;++j){/// loop over waves
    // normalize relatve intensities to number of bins!
    waves[bin->wavename(j)]=intens[j]/(double)usedbins;
    waves2[intens[j]/(double)usedbins]=bin->wavename(j);
  }// end loop over waves

  


  // print outcome:

  std::cout<<"Wavename | ID"<<std::endl;
  for(int j=0;j<nw;++j){/// loop over waves
    std::cout<<bin->wavename(j)<<"  |  "<<j<<std::endl;
  }

  std::cout<<"Wavename | Relative Intensity  (sorted by name)"<<std::endl;
  std::map<TString,double>::iterator wit=waves.begin();
  while(wit!=waves.end()){
    std::cout<<wit->first<<"  |  "<<wit->second<<std::endl;
    ++wit;
  }

  std::cout<<"----------------------------------------------------"<<std::endl;
  std::cout<<"Wavename | Relative Intensity  (sorted by intensity)"<<std::endl;
  std::map<double,TString>::iterator wit2=waves2.end();
  while(wit2!=waves2.begin()){
    --wit2;
    std::cout<<wit2->second<<"  |  "<<wit2->first<<std::endl;
  }

}

// dummy function needed since we link to but do not use minuit.o
int mnparm(int, string, double, double, double, double) {
    cerr << "this is impossible" << endl;
    throw "aFit";
    return 0;
}
