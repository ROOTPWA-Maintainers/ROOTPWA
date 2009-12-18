void runTestNBodyPhaseSpaceGen()
{
  gSystem->AddIncludePath("-I$ROOTPWA/src");
  gSystem->Load("../build/lib/librootpwa.so");
  gSystem->Load("../build/lib/librootpwagen.so");

  gROOT->ProcessLine(".x testNBodyPhaseSpaceGen.C+");
}
