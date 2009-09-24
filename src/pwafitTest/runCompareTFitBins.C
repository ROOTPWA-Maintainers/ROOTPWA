{
  gSystem->Load("librootpwa.so");

  gROOT->ProcessLine(".x compareTFitBins.C+");
}
