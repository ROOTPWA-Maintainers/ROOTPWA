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
//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Creates overview plot with intensities, phase, and coherence
//      for waves A and B from tree.
//      layout:
//              intensity A | phase A - B
//              ------------+----------------
//              intensity B | coherence A - B
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#include <sstream>
#include <string>

#include "TTree.h"
#include "TCanvas.h"

#include "utilities.h"
#include "plotIntensity.h"
#include "plotPhase.h"
#include "plotCoherence.h"


using namespace std;
using namespace rpwa;


// signature with wave indices
void
plot4(TTree*       tree,             // fitResult tree
      const int    waveIndexA,       // index of first wave
      const int    waveIndexB,       // index of second wave
      const double massMin     = 0,  // [GeV/c^2]
      const double massMax     = 0,  // [GeV/c^2]
      const string& branchName = "fitResult_v2")
{
  if (!tree) {
    printErr << "null pointer to tree. exiting." << endl;
    return;
  }

  // select mass range; convert from GeV/c^2 to MeV/c^2
  stringstream selectExpr;
  if ((massMin != 0) || (massMax != 0))
    selectExpr << "(massBinCenter() >= "<< massMin * 1000 << ") && (massBinCenter() <= " << massMax * 1000 << ")";

  stringstream canvName;
  canvName << "4plot_" << waveIndexA << "_" << waveIndexB;
  TCanvas* canv = new TCanvas(canvName.str().c_str(), canvName.str().c_str(), 10, 10, 1000, 800);
  canv->Divide(2, 2);
 
  // wave A intensity
  canv->cd(1);
  plotIntensity(tree, waveIndexA, false, kBlack, false, "", "APZ", 1, 0,
		selectExpr.str(), branchName);

  // wave A - wave B phase angle
  canv->cd(2);
  plotPhase(tree, waveIndexA, waveIndexB, selectExpr.str(), "", "APZ", kBlack, false, branchName);

  // wave B intensity
  canv->cd(3);
  plotIntensity(tree, waveIndexB, false, kBlack, false, "", "APZ", 1, 0,
		selectExpr.str(), branchName);

  // wave A - wave B coherence
  canv->cd(4);
  plotCoherence(tree, waveIndexA, waveIndexB, selectExpr.str(), "", "APZ", kBlack, false, branchName);
}
 

// signature with wave names
void
plot4(TTree*        tree,            // fitResult tree
      const string& waveNameA,       // name of first wave
      const string& waveNameB,       // name of second wave
      const double  massMin    = 0,  // [GeV/c^2]
      const double  massMax    = 0,  // [GeV/c^2]
      const string& branchName = "fitResult_v2")
{
  if (!tree) {
    printErr << "null pointer to tree. exiting." << endl;
    return;
  }
  // get wave indices (assumes same wave set in all trees)
  fitResult* massBin = new fitResult();
  tree->SetBranchAddress(branchName.c_str(), &massBin);
  tree->GetEntry(0);
  const int indexA = massBin->waveIndex(waveNameA);
  const int indexB = massBin->waveIndex(waveNameB);
  if ((indexA >= 0) && (indexB >= 0))
    return plot4(tree, indexA, indexB, massMin, massMax);
  printErr << "cannot find wave(s) in tree '" << tree->GetName() << "'. exiting." << endl;
}
