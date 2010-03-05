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
//      Draws intensity graph for single wave from tree.
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#include <iostream>
#include <sstream>

#include "TGraphErrors.h"
#include "TAxis.h"
#include "TLine.h"
#include "TPad.h"
#include "TMath.h"

#include "plotIntensity.h"


using namespace std;
using namespace rpwa;


// signature with wave index
TMultiGraph*
plotIntensity(const unsigned int nmbTrees,       // number of fitResult trees
	      TTree**            trees,          // array of fitResult trees
	      const int          waveIndex,      // wave index
	      const string&      selectExpr,     // TTree::Draw() selection expression
	      const string&      graphTitle,     // name and title of graph (default is waveId)
	      const char*        drawOption,     // draw option for graph
	      const double       normalization,  // scale factor for intensities
	      const int*         graphColors,    // array of colors for graph line and marker
	      const double       yAxisRangeMax,  // if != 0; range of y-axis is limited to this value
	      const bool         saveEps,        // if set, EPS file with name waveId{
	      const string&      branchName)
{
  for (unsigned int i = 0; i < nmbTrees; ++i)
    if (!trees[i]) {
      printErr << "NULL pointer to tree " << i << ". exiting." << endl;
      return 0;
    }
  // get wave name (assumes same wave set in all trees)
  fitResult* massBin = new fitResult();
  trees[0]->SetBranchAddress(branchName.c_str(), &massBin);
  trees[0]->GetEntry(0);
  const string waveName = massBin->waveName(waveIndex).Data();
  printInfo << "plotting wave intensity for wave '" << waveName << "' [" << waveIndex << "]";
  if (selectExpr != "")
    cout << " using selection criterion '" << selectExpr << "'";
  cout << endl;
  
  // build multiGraph
  TMultiGraph* graph = new TMultiGraph();
  {
    stringstream suffix;
    suffix << " [" << waveIndex << "]";
    if (graphTitle != "") {
      graph->SetName (graphTitle.c_str());
      graph->SetTitle(graphTitle.c_str());
    } else {
      graph->SetName (waveName.c_str());
      graph->SetTitle((waveName + suffix.str()).c_str());
    }
  }
  double maxY = 0;
  for (unsigned int i = 0; i < nmbTrees; ++i) {
    // build and run TTree::Draw() expression
    stringstream drawExpr;
    drawExpr << branchName << ".intensity(\"" << waveName << "\"):"
	     << branchName << ".intensityErr(\"" << waveName << "\"):"
	     << branchName << ".massBinCenter() >> h" << waveName << "_" << i;
    cout << "    running TTree::Draw() expression '" << drawExpr.str() << "' "
	 << "on tree '" << trees[i]->GetName() << "', '" << trees[i]->GetTitle() << "'" << endl;
    trees[i]->Draw(drawExpr.str().c_str(), selectExpr.c_str(), "goff");
      
    // extract data from TTree::Draw() result and build graph
    const int nmbBins = trees[i]->GetSelectedRows();
    vector<double> x(nmbBins), xErr(nmbBins);
    vector<double> y(nmbBins), yErr(nmbBins);
    for (int j = 0; j < nmbBins; ++j) {
      x   [j] = trees[i]->GetV3()[j] * 0.001;  // convert mass to GeV
      xErr[j] = 0;
      y   [j] = trees[i]->GetV1()[j] * normalization;  // scale intensities
      yErr[j] = trees[i]->GetV2()[j] * normalization;  // scale intensity errors
    }
    TGraphErrors* g = new TGraphErrors(nmbBins,
				       &(*(x.begin())),      // mass
				       &(*(y.begin())),      // intensity
				       &(*(xErr.begin())),   // mass error
				       &(*(yErr.begin())));  // intensity error
    {
      stringstream graphName;
      graphName << ((graphTitle == "") ? waveName : graphTitle) << "_" << i;
      g->SetName (graphName.str().c_str());
      g->SetTitle(graphName.str().c_str());
    }
    g->SetMarkerStyle(21);
    g->SetMarkerSize(0.5);
    if (graphColors) {
      g->SetMarkerColor(graphColors[i]);
      g->SetLineColor  (graphColors[i]);
    }
    graph->Add(g);
    
    // compute maximum for y-axis
    for (int j = 0; j < nmbBins; ++j) {
      const double val = y[j] + yErr[j];
      if (maxY < val)
	maxY = val;
    }
  }
  cout << "    maximum intensity for graph " << graph->GetName() << " is " << maxY << endl;
    
  if ((yAxisRangeMax > 0) && (maxY > yAxisRangeMax))
    maxY = yAxisRangeMax;
  graph->SetMinimum(-maxY * 0.1);
  graph->SetMaximum( maxY * 1.1);
  // draw graph
  graph->Draw(drawOption);
  graph->GetXaxis()->SetTitle("Mass [GeV]");
  graph->GetYaxis()->SetTitle("Intensity");
  TLine line;
  line.SetLineStyle(3);
  line.DrawLine(graph->GetXaxis()->GetXmin(), 0, graph->GetXaxis()->GetXmax(), 0);
  gPad->Update();

  // create EPS file
  if (saveEps)
    gPad->SaveAs((waveName + ".eps").c_str());

  return graph;
}
