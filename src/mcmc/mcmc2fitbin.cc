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
// ROOT Macro
// converts Narkov Chain MonteCarlo output into TFitBins

//#include <fitlog.h>
//#include <integral.h>
#include <iostream>
#include <list>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <complex>
#include <assert.h>
#include "TBranch.h"
#include "TTree.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TString.h"
#include "TComplex.h"
#include "TRandom.h"
#include "TFitBin.h"
#include "TMCMCMeta.h"
#include "TPWALikelihood.h"

using namespace std;


int lineno = 1; // global variables needed for lex (not understood)
char *progname;

int main(int argc, char** argv){

  if(argc<2){
    cout << "Usage: mcmc2fitbin <infile> <outfile>"<<endl;
    return 1;
  }

  // oopen input file
  TFile* infile=TFile::Open(argv[1],"READ");
  TTree* input=(TTree*)infile->Get("MarkovChain");
  TMCMCMeta* meta=(TMCMCMeta*)infile->Get("MetaInfo");
  assert(input!=0);
  assert(meta!=0);

 
  TString ofile_name(argv[2]);
  unsigned int evt=meta->NEvents;
  
  
  double low_mass=meta->low_mass;
  double high_mass=meta->high_mass;
  //double step=meta->stepsize;
  unsigned int rank=2;//ACHTUNG!!!meta->rank; // rank


  int ndim;
  input->SetBranchAddress("N",&ndim);
  input->GetEntry(0);

  double dL[ndim];   // gradient of loglikelihood
  //double dLlast[ndim]; 
  double E;        // likelihood-value
  double H;        // Hamiltonian
  double X[ndim];  // position in phasespace
  //double Xlast[ndim]; // position in PS for last MCMC step
  double R[ndim];  // vector of velocities in PS
  double Rlast[ndim];
  
  
  input->SetBranchAddress("dL",&dL);
  input->SetBranchAddress("X",&X);
  input->SetBranchAddress("L",&E);
  input->SetBranchAddress("H",&H);
  input->SetBranchAddress("R",&R);
  input->SetBranchAddress("Rlast",&Rlast);

  
  TFitBin* result=new TFitBin();
  TFile* outfile=new TFile(ofile_name,"UPDATE");
  // check if tree exists

  TTree* tree=(TTree*)outfile->Get("pwa");
  if(tree==NULL){ // create new tree
    std::cout << "File empty. Creating new results tree!!!" << std::endl;
    tree=new TTree("pwa","pwa");
    tree->Branch("fitbin",&result);
  }
  else tree->SetBranchAddress("fitbin",&result);

 
  double mass=0.5*(low_mass+high_mass);
  
  vector<TString> wavetitles; // without rank
  vector<TString> wavenames; 
  const vector<TString>& parnames=meta->parnames;

  int j=0;
  for(unsigned int i=0;i<parnames.size();++i){
    TString title;
    unsigned int l=parnames[i].Length()-3;
    TString name;
    {if(parnames[i].Contains("V_"))name=parnames[i](0,l+4); else name=parnames[i](0,l);}
    cout << "Name= " << name << endl;
    
    if(std::find(wavenames.begin(),wavenames.end(),name)==wavenames.end()){
      wavenames.push_back(name);
      ++j;
    }
    else continue;
    
    l=wavenames[j-1].Length();
    {if(wavenames[j-1].Contains("V_"))title=wavenames[j-1](2,l-2); else title=wavenames[j-1](3,l-3);}
   
    //cout << "Title=" << title << endl;
    // only fill in title once!!!
    if(std::find(wavetitles.begin(),wavetitles.end(),title)==wavetitles.end()){
      wavetitles.push_back(title);
      cout << title << endl;
    }

  } 

  // Loop over mcmc samples
  int nmc=input->GetEntriesFast();
  for(int mci=0;mci<nmc;++mci){
  
    input->GetEntry(mci);

    vector<TComplex> amps;
    vector<complex<double> > V;
    
    double re,im;
    unsigned int k=0;
    unsigned int nwaves=wavetitles.size()-1; // -1 ... because of flat
    for(unsigned int r=0;r<rank;++r){
      for(unsigned int i=r;i<nwaves;++i){
	if(i<r){re=0;im=0;}
	else if(i==r){re=X[k++];im=0;} // real parameter
	else {
	  re=X[k++];
	  im=X[k++];
	}
	V.push_back(complex<double>(re,im));
      }
    } // end loop over rank
    V.push_back(complex<double>(X[k],0));
    
    // convert to TComplex;
    for(unsigned int i=0;i<V.size();++i)amps.push_back(TComplex(V[i].real(),V[i].imag()));
    
    //cout << " Created amps." << endl;
    
    // error matrix -> No errors for mcmc
    TMatrixD errm(0,0);
    
    
    
    vector<pair<int,int> > indices; // no errm => no indices
    
    
    // get normalization integral
    
        
    int n=wavetitles.size()-1;
    
    TCMatrix integr(n,n);
    integr=meta->Norm;
    
    
    cout << "filling TFitBin" << endl;
    cout << "Number of amps=" << amps.size()<< endl;
    cout << "Number of indices=" << indices.size() << endl;
    cout << "Number of wavenames=" << wavenames.size()<< endl;
    cout << "Number of wavetitles=" << wavetitles.size()<< endl;
    cout << "Dimension of errormatrix=" << errm.GetNrows() << endl;
    cout << "Dimension of integral matrix=" << integr.nrows() << endl;
    
    
    result->fill(amps,
		 indices,
		 wavenames,
		 evt,
		 evt,
		 mass,
		 integr,
		 errm,
		 E,
		 rank);
    
    //result->PrintParameters();
 
 
    TString binname="fitbin";binname+=low_mass;binname+="_";binname+=high_mass;
    binname.ReplaceAll(" ","");
    tree->Fill();

  }// end loop over mcmcsamples
  tree->Write("",TObject::kOverwrite);
  outfile->Close();
  
  return 0;
  }
  
// dummy function needed since we link to but do not use minuit.o
int mnparm(int, string, double, double, double, double) {
    cerr << "this is impossible" << endl;
    throw "aFit";
    return 0;
}
