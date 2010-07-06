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

#include "TChain.h"
#include "TH1.h"
#include "TH2.h"

#include "utilities.h"


using namespace std;


bool
plotAmpDiffs(const string&  inFileNamePattern,
             const long int maxNmbEvents = -1,
             const string&  inTreeName   = "ampDiffTree")
{
  // open input file
  printInfo << "opening input file(s) '" << inFileNamePattern << "'" << endl;
  TChain inTree(inTreeName.c_str());
  if (inTree.Add(inFileNamePattern.c_str()) < 1) {
    printWarn << "no events in input file(s) '" << inFileNamePattern << "'" << endl;
    return false;
  }
  inTree.GetListOfFiles()->ls();

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
  TH1D* hAmps1Re   = new TH1D("hAmps1Re",   "hAmps1Re;Event Number;#Rgothic[Amp 1]", nmbEvents, -0.5, nmbEvents - 0.5);
  TH1D* hAmps1Im   = new TH1D("hAmps1Im",   "hAmps1Im;Event Number;#Jgothic[Amp 1]", nmbEvents, -0.5, nmbEvents - 0.5);
  TH1D* hAmps2Re   = new TH1D("hAmps2Re",   "hAmps2Re;Event Number;#Rgothic[Amp 2]", nmbEvents, -0.5, nmbEvents - 0.5);
  TH1D* hAmps2Im   = new TH1D("hAmps2Im",   "hAmps2Im;Event Number;#Jgothic[Amp 2]", nmbEvents, -0.5, nmbEvents - 0.5);
  TH1D* hAbsDiffRe = new TH1D("hAbsDiffRe", "hAbsDiffRe;#Rgothic[Amp 1 - Amp 2];Count", 1000000, -1e-3, 1e-3);
  TH1D* hAbsDiffIm = new TH1D("hAbsDiffIm", "hAbsDiffIm;#Jgothic[Amp 1 - Amp 2];Count", 1000000, -1e-3, 1e-3);
  TH1D* hRelDiffRe = new TH1D("hRelDiffRe", "hRelDiffRe;(#Rgothic[Ampl 1] - #Rgothic[Ampl 2]) / #Rgothic[Ampl 1];Count", 1000000, -1e-2, 1e-2);
  TH1D* hRelDiffIm = new TH1D("hRelDiffIm", "hRelDiffIm;(#Jgothic[Ampl 1] - #Jgothic[Ampl 2]) / #Jgothic[Ampl 1];Count", 1000000, -1e-2, 1e-2);
  TH2D* hCorrRe    = new TH2D("hCorrRe",    "hCorrRe;#Rgothic[Amp 1];#Rgothic[Amp 2]", 5000, -5, 5, 5000, -5, 5);
  TH2D* hCorrIm    = new TH2D("hCorrIm",    "hCorrIm;#Jgothic[Amp 1];#Jgothic[Amp 2]", 5000, -5, 5, 5000, -5, 5);

  // loop over tree
  for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
	  progressIndicator(eventIndex, nmbEvents);
  
    if (inTree.LoadTree(eventIndex) < 0)
      break;
    inTree.GetEntry(eventIndex);

    hAmps1Re->SetBinContent(eventIndex + 1, valReal[0]);
    hAmps1Im->SetBinContent(eventIndex + 1, valImag[0]);
    hAmps2Re->SetBinContent(eventIndex + 1, valReal[1]);
    hAmps2Im->SetBinContent(eventIndex + 1, valImag[1]);
    hAbsDiffRe->Fill(absDiffReal);
    hAbsDiffIm->Fill(absDiffImag);
    hRelDiffRe->Fill(relDiffReal);
    hRelDiffIm->Fill(relDiffImag);
    hCorrRe->Fill(valReal[0], valReal[1]);
    hCorrIm->Fill(valImag[0], valImag[1]);
  }

  return true;
}
