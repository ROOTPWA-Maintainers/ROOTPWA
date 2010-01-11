#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TLorentzRotation.h"
#include "TH1D.h"
#include "TClonesArray.h"
#include "TROOT.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "NParticleEvent.h"


using namespace std; 


void sampleWeightedEvents(TTree* mctr, string outputfile,unsigned int n,bool doweight=true){
  

  ofstream outfile(outputfile.c_str());


  TTree* tr=mctr;
  double weight=1;
  double maxweight=0; 
  TClonesArray* p=new TClonesArray("TLorentzVector");
  TLorentzVector* beam=NULL;
  int qbeam;
  std::vector<int>* q=NULL; 
  tr->SetBranchAddress("weight",&weight);
  
  // get max weight
  unsigned int nevt=tr->GetEntries();	
  for(unsigned int i=0;i<nevt;++i){
    tr->GetEntry(i);
    if(weight>maxweight)maxweight=weight;
  }// end loop over events
  cout<<"Maxweight="<<maxweight<<endl;

 tr->SetBranchAddress("p",&p);
 tr->SetBranchAddress("beam",&beam);
 tr->SetBranchAddress("qbeam",&qbeam);
 tr->SetBranchAddress("q",&q);
 TVector3 vertex;
 NParticleEvent event(p,q,beam,&qbeam,&vertex);
 
 unsigned int nselected=0;
 unsigned int attempts=0;
 while(nselected<n){
   for(unsigned int i=0;i<nevt;++i){
     tr->GetEntry(i);
     ++attempts;
     if(doweight && gRandom->Uniform()>weight/maxweight)continue;
     
     ++nselected;
     event.refresh();
     event.writeGAMP(outfile);
     if(nselected==n)break;
     
   }// end loop over events
 }
 cout << "Attempts: " << attempts << endl;
 cout << "Efficiency=" << (double)nselected/(double)attempts << endl;
 
 outfile.close();
 
}
