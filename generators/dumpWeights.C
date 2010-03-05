#include <iomanip>
#include <fstream>
#include "TFile.h"
#include "TTree.h"

using namespace std;

void dumpWeights(TString file){

  TFile* infile=TFile::Open(file,"READ");
  file.ReplaceAll(".root",".wht");

  ofstream outfile(file.Data());

  TTree* tr=(TTree*)infile->Get("pwevents");
  double w=1;
  tr->SetBranchAddress("impweight",&w);

  outfile << std::setprecision(10);

  unsigned int n=tr->GetEntries();
  for(unsigned int i=0; i<n; ++i){
    tr->GetEntry(i);
    outfile << w << endl;
  }

  infile->Close();
  outfile.close();


}
