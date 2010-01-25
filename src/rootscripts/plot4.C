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

#include "plotIntensity.h"
#include "plotPhase.h"
#include "plotCoherence.h"


using namespace std;


void
plot4(TTree*       tree,         // fitResult tree
      const int    waveIndexA,   // index of first wave
      const int    waveIndexB,   // index of second wave
      const double massMin = 1,  // [GeV/c^2]
      const double massMax = 4)  // [GeV/c^2]
{
  // select mass range; convert from GeV/c^2 to MeV/c^2
  stringstream selectExpr;
  selectExpr << "(massBinCenter() >= "<< massMin * 1000 << ") && (massBinCenter() <= " << massMax * 1000 << ")";

  stringstream canvName;
  canvName << "4plot_" << waveIndexA << "_" << waveIndexB;
  TCanvas* canv = new TCanvas(canvName.str().c_str(), canvName.str().c_str(), 10, 10, 1000, 800);
  canv->Divide(2, 2);
 
  // wave A intensity
  canv->cd(1);
  plotIntensity(tree, waveIndexA, selectExpr.str());

  // wave A - wave B phase angle
  canv->cd(2);
  plotPhase(tree, waveIndexA, waveIndexB, selectExpr.str());

  // wave B intensity
  canv->cd(3);
  plotIntensity(tree, waveIndexB, selectExpr.str());

  // wave A - wave B coherence
  canv->cd(4);
  plotCoherence(tree, waveIndexA, waveIndexB, selectExpr.str());
}
 
