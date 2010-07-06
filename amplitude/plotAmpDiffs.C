///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
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
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      plots differences between amplitude files
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <string>
#include <set>

#include "TFile.h"
#include "TChain.h"
#include "TTreeFormula.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"

#include "utilities.h"


using namespace std;


bool
plotAmpDiffs(const string&  inFileNamePattern,
             const string&  outFileName  = "ampDiffPlots.root",
             const long int maxNmbEvents = -1,
             const string&  inTreeName   = "ampDiffTree")
{
  // open input file(s)
  printInfo << "opening input file(s) '" << inFileNamePattern << "'" << endl;
  TChain inTree(inTreeName.c_str());
  if (inTree.Add(inFileNamePattern.c_str()) < 1) {
    printWarn << "no events in input file(s) '" << inFileNamePattern << "'" << endl;
    return false;
  }
  inTree.GetListOfFiles()->ls();

  // open output file
  TFile* outFile = 0;
  if (outFileName != "") {
    printInfo << "writing difference plots to '" << outFileName << "'" << endl;
    outFile = TFile::Open(outFileName.c_str(), "RECREATE");
  }

  // connect tree leafs
  TObjString* ampName = 0;
  UInt_t      eventNmb;
  UInt_t      massBinMin, massBinMax;  // [MeV/c^2]
  double_t    valReal[2], valImag[2];
  double_t    absDiffReal, absDiffImag;
  double_t    relDiffReal, relDiffImag;
  inTree.SetBranchAddress("ampName",     &ampName);
  inTree.SetBranchAddress("eventNmb",    &eventNmb);
  inTree.SetBranchAddress("massBinMin",  &massBinMin);
  inTree.SetBranchAddress("massBinMax",  &massBinMax);
  inTree.SetBranchAddress("valReal",     valReal);
  inTree.SetBranchAddress("valImag",     valImag);
  inTree.SetBranchAddress("absDiffReal", &absDiffReal);
  inTree.SetBranchAddress("absDiffImag", &absDiffImag);
  inTree.SetBranchAddress("relDiffReal", &relDiffReal);
  inTree.SetBranchAddress("relDiffImag", &relDiffImag);

  // book histograms
  const long int nmbEventsTree = inTree.GetEntries();
  const long int nmbEvents     = ((maxNmbEvents > 0) ? min(maxNmbEvents, nmbEventsTree)
                                  : nmbEventsTree);
  TH1F* hAmps1Re   = new TH1F("amps1Real",   "amps1Re;Event Number;#Rgothic[Amp 1]", nmbEvents, -0.5, nmbEvents - 0.5);
  TH1F* hAmps1Im   = new TH1F("amps1Imag",   "amps1Im;Event Number;#Jgothic[Amp 1]", nmbEvents, -0.5, nmbEvents - 0.5);
  TH1F* hAmps2Re   = new TH1F("amps2Real",   "amps2Re;Event Number;#Rgothic[Amp 2]", nmbEvents, -0.5, nmbEvents - 0.5);
  TH1F* hAmps2Im   = new TH1F("amps2Imag",   "amps2Im;Event Number;#Jgothic[Amp 2]", nmbEvents, -0.5, nmbEvents - 0.5);
  TH1F* hAbsDiffRe = new TH1F("absDiffReal", "absDiffRe;#Rgothic[Amp 1 - Amp 2];Count", 1000000, -1e-3, 1e-3);
  TH1F* hAbsDiffIm = new TH1F("absDiffImag", "absDiffIm;#Jgothic[Amp 1 - Amp 2];Count", 1000000, -1e-3, 1e-3);
  TH1F* hRelDiffRe = new TH1F("relDiffReal", "relDiffRe;(#Rgothic[Ampl 1] - #Rgothic[Ampl 2]) / #Rgothic[Ampl 1];Count", 1000000, -1e-2, 1e-2);
  TH1F* hRelDiffIm = new TH1F("relDiffImag", "relDiffIm;(#Jgothic[Ampl 1] - #Jgothic[Ampl 2]) / #Jgothic[Ampl 1];Count", 1000000, -1e-2, 1e-2);
  // TH2F* hCorrRe    = new TH2F("corrReal",    "corrRe;#Rgothic[Amp 1];#Rgothic[Amp 2]", 5000, -5, 5, 5000, -5, 5);
  // TH2F* hCorrIm    = new TH2F("corrImag",    "corrIm;#Jgothic[Amp 1];#Jgothic[Amp 2]", 5000, -5, 5, 5000, -5, 5);

  // loop over tree
  set<string> ampNames;
  for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
	  progressIndicator(eventIndex, nmbEvents);
  
    if (inTree.LoadTree(eventIndex) < 0)
      break;
    inTree.GetEntry(eventIndex);

    ampNames.insert(ampName->GetString().Data());

    hAmps1Re->SetBinContent(eventIndex + 1, valReal[0]);
    hAmps1Im->SetBinContent(eventIndex + 1, valImag[0]);
    hAmps2Re->SetBinContent(eventIndex + 1, valReal[1]);
    hAmps2Im->SetBinContent(eventIndex + 1, valImag[1]);
    hAbsDiffRe->Fill(absDiffReal);
    hAbsDiffIm->Fill(absDiffImag);
    hRelDiffRe->Fill(relDiffReal);
    hRelDiffIm->Fill(relDiffImag);
    // hCorrRe->Fill(valReal[0], valReal[1]);
    // hCorrIm->Fill(valImag[0], valImag[1]);
  }
  
  const string leafsToDraw[] = {"absDiffReal",
                                "absDiffImag",
                                "relDiffReal",
                                "relDiffImag"};
  const unsigned int nmbLeafs = sizeof(leafsToDraw) / sizeof(leafsToDraw[0]);
  TH2F* hLeafs[nmbLeafs];
  for (unsigned int i = 0; i < nmbLeafs; ++i) {
	  printInfo << "filling histogram for leaf '" << leafsToDraw[i] << "'" << endl;
	  const double min = inTree.GetMinimum(leafsToDraw[i].c_str());
	  const double max = inTree.GetMaximum(leafsToDraw[i].c_str());
	  hLeafs[i] = new TH2F((leafsToDraw[i] + "Wave").c_str(),
	                       (leafsToDraw[i] + " vs. Wave;Wave;" + leafsToDraw[i]).c_str(),
	                       1, 0, 1, 10000, min - 0.1 * fabs(min), max + 0.1 * fabs(max));
	  hLeafs[i]->SetBit(TH1::kCanRebin);
	  TTreeFormula* leafVal = new TTreeFormula("", leafsToDraw[i].c_str(), &inTree);
	  for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
		  progressIndicator(eventIndex, nmbEvents);
		  if (inTree.LoadTree(eventIndex) < 0)
			  break;
		  inTree.GetEntry(eventIndex);
		  leafVal->UpdateFormulaLeaves();
		  hLeafs[i]->Fill(ampName->GetString().Data(), leafVal->EvalInstance(), 1.);
	  }
	  hLeafs[i]->LabelsDeflate("X");
	  hLeafs[i]->SetLabelSize(0.02, "X");
	  hLeafs[i]->SetDrawOption("COLZ");
  }

  if (outFile) {
	  outFile->Write();
    outFile->Close();
    delete outFile;
    printInfo << "wrote difference plots to '" << outFileName << "'" << endl;
  }
  return true;
}
