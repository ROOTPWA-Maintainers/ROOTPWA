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
//      Draws coherence of (wave A - wave B) from tree.
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
#include "TLine.h"
#include "TPad.h"

#include "utilities.h"
#include "fitResult.h"
#include "plotCoherence.h"


using namespace std;
using namespace rpwa;


// signature with wave names
TGraphErrors*
plotCoherence(TTree*        tree,        // fitResult tree
	      const string& waveNameA,   // name of first wave
	      const string& waveNameB,   // name of second wave
	      const string& selectExpr,  // TTree::Draw() selection expression
	      const string& graphTitle,  // name and title of graph
	      const char*   drawOption,  // draw option for graph
	      const int     graphColor,  // color of line and marker
	      const bool    saveEps,     // if set, EPS file with name waveId is created
	      const string& branchName)

{
  if (!tree) {
    printErr << "NULL pointer to tree. exiting." << endl;
    return 0;
  }

  // call plotCoherence with wave indices
  fitResult* massBin = new fitResult();
  tree->SetBranchAddress(branchName.c_str(), &massBin);
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
  return plotCoherence(tree, waveIndices[0], waveIndices[1], selectExpr,
		       graphTitle, drawOption, graphColor, saveEps, branchName);
}


// signature with wave indices
TGraphErrors*
plotCoherence(TTree*        tree,        // fitResult tree
	      const int     waveIndexA,  // index of first wave
	      const int     waveIndexB,  // index of second wave
	      const string& selectExpr,  // TTree::Draw() selection expression
	      const string& graphTitle,  // name and title of graph
	      const char*   drawOption,  // draw option for graph
	      const int     graphColor,  // color of line and marker
	      const bool    saveEps,     // if set, EPS file with name waveId is created
	      const string& branchName)
{
  if (!tree) {
    printErr << "NULL pointer to tree. exiting." << endl;
    return 0;
  }
  // get wave names
  fitResult* massBin = new fitResult();
  tree->SetBranchAddress(branchName.c_str(), &massBin);
  tree->GetEntry(0);
  const string waveNameA = massBin->waveName(waveIndexA).Data();
  const string waveNameB = massBin->waveName(waveIndexB).Data();
  printInfo << "plotting coherence of wave '" << waveNameA << "' [" << waveIndexA << "] "
	    << "and wave '" << waveNameB << "' [" << waveIndexB << "]" << endl;
  if (selectExpr != "")
    cout << "    using selection criterion '" << selectExpr << "'" << endl;

  // build and run TTree::Draw() expression
  stringstream drawExpr;
  drawExpr << branchName << ".coherence("  << waveIndexA << "," << waveIndexB << "):"
	   << branchName << ".coherenceErr(" << waveIndexA << "," << waveIndexB << "):"
	   << branchName << ".massBinCenter() >> h" << waveIndexA << "_" << waveIndexB;
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
					 tree->GetV1(),    // coherence
					 0,                // mass error
					 tree->GetV2());   // coherence error
  stringstream graphName;
  graphName << "coherence_" << waveNameA << "_" << waveNameB;
  {
    stringstream title;
    title << "Coherence(" << waveNameA << " [" << waveIndexA << "], "
	  << ", " << waveNameB << " [" << waveIndexB << "])";
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
  graph->GetYaxis()->SetTitle("Coherence");
  graph->SetMinimum(-0.1);
  graph->SetMaximum( 1.1);
  // TGraphErrors* clone = (TGraphErrors*)graph->DrawClone(drawOption);
  // clone->SetName(graphName.str().c_str());
  gPad->Update();
  TLine line;
  line.SetLineStyle(3);
  line.DrawLine(graph->GetXaxis()->GetXmin(), 0, graph->GetXaxis()->GetXmax(), 0);

  // create EPS file
  if (saveEps)
    gPad->SaveAs((graphName.str() + ".eps").c_str());

  return graph;
}
