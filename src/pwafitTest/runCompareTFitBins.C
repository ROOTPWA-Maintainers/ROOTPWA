{
  gSystem->Load("libRootPwa.so");

  gROOT->ProcessLine(".x compareTFitBins.C+");
}
