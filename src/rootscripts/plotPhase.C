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
//      Draws phase angle of (wave A - wave B) from tree.
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#include <iostream>
#include <sstream>

#include "TAxis.h"
#include "TPad.h"

#include "utilities.h"
#include "TFitResult.h"
#include "plotPhase.h"


using namespace std;


// signature with wave names
TGraphErrors*
plotPhase(TTree*        tree,        // TFitResult tree
	  const string& waveNameA,   // name of first wave
	  const string& waveNameB,   // name of second wave
	  const string& selectExpr,  // TTree::Draw() selection expression
	  const string& graphTitle,  // name and title of graph
	  const char*   drawOption,  // draw option for graph
	  const int     graphColor,  // color of line and marker
	  const bool    saveEps)     // if set, EPS file with name waveId is created
{
  if (!tree) {
    printErr << "NULL pointer to tree. exiting." << endl;
    return 0;
  }

  // call plotPhase with wave indices
  TFitResult* massBin = new TFitResult();
  tree->SetBranchAddress("fitResult", &massBin);
  tree->GetEntry(0);
  const string waveNames[2]   = {waveNameA, waveNameB};
  int          waveIndices[2] = {-1, -1};
  for (int i = 0; i < 2; ++i) {
    for (unsigned int j = 0; j < massBin->nmbWaves(); ++j)
      if (massBin->waveName(j) == waveNames[i]) {
	waveIndices[i] = j;
	break;
      }
    if (waveIndices[i] < 0) {
      printErr << "cannot find wave '" << waveNames[i] << "' "
	       << "in tree '" << tree->GetName() << "'. exiting." << endl;
      return 0;
    }
  }
  return plotPhase(tree, waveIndices[0], waveIndices[1], selectExpr,
		   graphTitle, drawOption, graphColor, saveEps);
}


// signature with wave names
TGraphErrors*
plotPhase(TTree*        tree,        // TFitResult tree
	  const int     waveIndexA,  // index of first wave
	  const int     waveIndexB,  // index of second wave
	  const string& selectExpr,  // TTree::Draw() selection expression
	  const string& graphTitle,  // name and title of graph
	  const char*   drawOption,  // draw option for graph
	  const int     graphColor,  // color of line and marker
	  const bool    saveEps)     // if set, EPS file with name waveId is created
{
  if (!tree) {
    printErr << "NULL pointer to tree. exiting." << endl;
    return 0;
  }
  // get wave names
  TFitResult* massBin = new TFitResult();
  tree->SetBranchAddress("fitResult", &massBin);
  tree->GetEntry(0);
  const string waveNameA = massBin->waveName(waveIndexA).Data();
  const string waveNameB = massBin->waveName(waveIndexB).Data();
  printInfo << "plotting phase between wave '" << waveNameA << "' [" << waveIndexA << "] "
	    << "and wave '" << waveNameB << "' [" << waveIndexB << "]" << endl;
  if (selectExpr != "")
    cout << "    using selection criterion '" << selectExpr << "'" << endl;

  // build and run TTree::Draw() expression
  stringstream drawExpr;
  drawExpr << "phase("     << waveIndexA << "," << waveIndexB << ")"
	   << ":phaseErr(" << waveIndexA << "," << waveIndexB << ")"
	   << ":massBinCenter() >> h" << waveIndexA << "_" << waveIndexB;
  cout << "    running TTree::Draw() expression '" << drawExpr.str() << "' "
       << "on tree '" << tree->GetName() << "', '" << tree->GetTitle() << "'" << endl;
  tree->Draw(drawExpr.str().c_str(), selectExpr.c_str(), "goff");

  // extract data from TTree::Draw() result and build graph
  const int nmbBins = tree->GetSelectedRows();
  vector<double> x(nmbBins), xErr(nmbBins);
  for (int i = 0; i < nmbBins; ++i) {
    x[i]    = tree->GetV3()[i] * 0.001;  // convert mass to GeV
    xErr[i] = 0;
  }
  TGraphErrors* graph = new TGraphErrors(nmbBins,
					 &(*(x.begin())),  // mass
					 tree->GetV1(),    // intensity
					 0,                // mass error
					 tree->GetV2());   // intensity error
  stringstream graphName;
  graphName << "phase_" << waveNameA << "_" << waveNameB;
  {
    stringstream title;
    title << "#Delta#varphi(" << waveNameA << " [" << waveIndexA << "], "
	  << waveNameB << " [" << waveIndexB << "])";
    if (graphTitle != "") {
      graph->SetName (graphTitle.c_str());
      graph->SetTitle(graphTitle.c_str());
    } else {
      graph->SetName(graphName.str().c_str());
      graph->SetTitle(title.str().c_str());
    }
  }

  // beautify and draw graph
  graph->Draw(drawOption);
  graph->SetMarkerStyle(21);
  graph->SetMarkerSize(0.5);
  graph->SetMarkerColor(graphColor);
  graph->SetLineColor(graphColor);
  graph->GetXaxis()->SetTitle("Mass [GeV]");
  graph->GetYaxis()->SetTitle("Phase Angle [deg]");
  graph->SetMinimum(-180);
  graph->SetMaximum(180);
  // TGraphErrors* clone = (TGraphErrors*);
  // clone->SetName(name.str().c_str());
  gPad->Update();

  // create EPS file
  if (saveEps)
    gPad->SaveAs((graphName.str() + ".eps").c_str());

  return graph;
}
