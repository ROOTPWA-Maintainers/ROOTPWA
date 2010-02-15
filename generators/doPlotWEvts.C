void doPlotWEvts(TString evtfile, TString mcfile, 
		 TString outfile, TString mass){

  gROOT->ProcessLine(".L $ROOTPWA/generators/plotWeightedEvts.C+");
  
  TFile* file1=TFile::Open(evtfile,"READ");
  TFile* file2=TFile::Open(mcfile,"READ");
  TTree* data=(TTree*)file1->Get("events");
  TTree* mc=(TTree*)file2->Get("pwevents");

  plotWeightedEvts(mc,data,outfile,mass);

  file1->Close();
  file2->Close();

}
