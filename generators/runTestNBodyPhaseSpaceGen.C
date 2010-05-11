void runTestNBodyPhaseSpaceGen()
{
  gSystem->AddIncludePath("-I$ROOTPWA/src");
  gSystem->Load("../build/lib/libRootPwa.so");
  gSystem->Load("../build/lib/libRootPwaGen.so");

  gROOT->ProcessLine(".x testNBodyPhaseSpaceGen.C+");
}
