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

#include "TAxis.h"
#include "TLine.h"
#include "TPad.h"

#include "../utilities.h"
#include "../TFitResult.h"
#include "plotIntensity.h"


using namespace std;


// signature with wave name
TGraphErrors*
plotIntensity(TTree*        tree,           // TFitResult tree
	      const string& waveName,       // wave name
	      const string& selectExpr,     // TTree::Draw() selection expression
	      const string& graphTitle,     // name and title of graph (default is waveId)
	      const char*   drawOption,     // draw option for graph
	      const double  normalization,  // scale factor for intensities
	      const int     graphColor,     // color of line and marker
	      const bool    saveEps)        // if set, EPS file with name of wave is created
{
  if (!tree) {
    cerr << "plotIntensity() error: Null pointer to tree. Exiting." << endl;
    return 0;
  }

  // call plotIntensity with wave index
  TFitResult* massBin = new TFitResult();
  tree->SetBranchAddress("fitResult", &massBin);
  tree->GetEntry(0);
  for (unsigned int i = 0; i < massBin->nmbWaves(); ++i)
    if (massBin->waveName(i) == waveName)
      return plotIntensity(tree, i, selectExpr, graphTitle, drawOption, normalization, graphColor, saveEps);
  cerr << "plotIntensity() error: Cannot find wave '" << waveName << "' in tree '" << tree->GetName() << "'. Exiting." << endl;
  return 0;
}


// signature with wave index
TGraphErrors*
plotIntensity(TTree*        tree,           // TFitResult tree
	      const int     waveIndex,      // wave index
	      const string& selectExpr,     // TTree::Draw() selection expression
	      const string& graphTitle,     // name and title of graph (default is waveId)
	      const char*   drawOption,     // draw option for graph
	      const double  normalization,  // scale factor for intensities
	      const int     graphColor,     // color of line and marker
	      const bool    saveEps)        // if set, EPS file with name of wave is created
{
  if (!tree) {
    cerr << "plotIntensity() error: Null pointer to tree. Exiting." << endl;
    return 0;
  }
  // get wave name
  TFitResult* massBin = new TFitResult();
  tree->SetBranchAddress("fitResult", &massBin);
  tree->GetEntry(0);
  const string waveName = massBin->waveName(waveIndex).Data();
  cout << "Plotting wave intensity for wave '" << waveName << "' [" << waveIndex << "]";
  if (selectExpr != "")
    cout << " using selection criterion '" << selectExpr << "'";
  cout << "." << endl;

  // build and run TTree::Draw() expression
  stringstream drawExpr;
  drawExpr << "intensity(\"" << waveName << "\"):intensityErr(\"" << waveName << "\"):massBinCenter() >> h" << waveName;
  cout << "    Running TTree::Draw() expression '" << drawExpr.str() << "' "
       << "on tree '" << tree->GetName() << "'" << endl;
  tree->Draw(drawExpr.str().c_str(), selectExpr.c_str(), "goff");

  // extract data from TTree::Draw() result and build graph
  const int nmbBins = tree->GetSelectedRows();
  vector<double> x(nmbBins), xErr(nmbBins);
  vector<double> y(nmbBins), yErr(nmbBins);
  for (int i = 0; i < nmbBins; ++i) {
    x   [i] = tree->GetV3()[i] * 0.001;  // convert mass to GeV
    xErr[i] = 0;
    y   [i] = tree->GetV1()[i] * normalization;  // scale intensities
    yErr[i] = tree->GetV2()[i] * normalization;  // scale intensity errors
  }
  TGraphErrors* g = new TGraphErrors(nmbBins,
				     &(*(x.begin())),      // mass
				     &(*(y.begin())),      // intensity
				     &(*(xErr.begin())),   // mass error
				     &(*(yErr.begin())));  // intensity error
  g->SetName (waveName.c_str());
  g->SetTitle(waveName.c_str());
  if (graphTitle != "") {
    g->SetName (graphTitle.c_str());
    g->SetTitle(graphTitle.c_str());
  }

  // compute maximum for y-axis
  double maxY = 0;
  for (int i = 0; i < nmbBins; ++i)
    if(maxY < g->GetY()[i])
      maxY = g->GetY()[i];
  cout << "    Maximum intensity for graph " << g->GetName() << " is " << maxY << "." << endl;
  
  // draw graph
  g->SetMarkerStyle(21);
  g->SetMarkerSize(0.5);
  g->SetMarkerColor(graphColor);
  g->SetLineColor(graphColor);
  g->SetMinimum(-maxY * 0.1);
  g->SetMaximum( maxY * 1.1);
  g->GetXaxis()->SetTitle("Mass [GeV]");
  g->GetYaxis()->SetTitle("Intensity");
  TGraphErrors* clone = (TGraphErrors*)g->DrawClone(drawOption);
  clone->GetYaxis()->SetRangeUser(-maxY * 0.1, maxY * 1.1);
  clone->SetName(waveName.c_str());
  TLine line;
  line.SetLineStyle(3);
  line.DrawLine(clone->GetXaxis()->GetXmin(), 0, clone->GetXaxis()->GetXmax(), 0);

  // create EPS file
  if (saveEps)
    gPad->SaveAs((waveName + ".eps").c_str());

  return g;
}


// signature with wave name
TMultiGraph*
plotIntensity(const unsigned int nmbTrees,       // number of TFitResult trees
	      TTree**            trees,          // array of TFitResult trees
	      const string&      waveName,       // wave name
	      const string&      selectExpr,     // TTree::Draw() selection expression
	      const string&      graphTitle,     // name and title of graph (default is waveId)
	      const char*        drawOption,     // draw option for graph
	      const double       normalization,  // scale factor for intensities
	      const int*         graphColors,    // array of colors for graph line and marker
	      const bool         saveEps)        // if set, EPS file with name waveId{
{
  for (unsigned int i = 0; i < nmbTrees; ++i)
    if (!trees[i]) {
      printErr << "NULL pointer to tree " << i << ". exiting." << endl;
      return 0;
    }

  // call plotIntensity with wave index (assumes same wave set in all trees)
  TFitResult* massBin = new TFitResult();
  trees[0]->SetBranchAddress("fitResult", &massBin);
  trees[0]->GetEntry(0);
  for (unsigned int i = 0; i < massBin->nmbWaves(); ++i)
    if (massBin->waveName(i) == waveName)
      return plotIntensity(nmbTrees, trees, i, selectExpr, graphTitle,
			   drawOption, normalization, graphColors, saveEps);
  printErr << "cannot find wave '" << waveName << "' in tree '"
	   << trees[0]->GetName() << "'. exiting." << endl;
  return 0;
}


// signature with wave index
TMultiGraph*
plotIntensity(const unsigned int nmbTrees,       // number of TFitResult trees
	      TTree**            trees,          // array of TFitResult trees
	      const int          waveIndex,      // wave index
	      const string&      selectExpr,     // TTree::Draw() selection expression
	      const string&      graphTitle,     // name and title of graph (default is waveId)
	      const char*        drawOption,     // draw option for graph
	      const double       normalization,  // scale factor for intensities
	      const int*         graphColors,    // array of colors for graph line and marker
	      const bool         saveEps)        // if set, EPS file with name waveId{
{
  for (unsigned int i = 0; i < nmbTrees; ++i)
    if (!trees[i]) {
      printErr << "NULL pointer to tree " << i << ". exiting." << endl;
      return 0;
    }
  // get wave name (assumes same wave set in all trees)
  TFitResult* massBin = new TFitResult();
  trees[0]->SetBranchAddress("fitResult", &massBin);
  trees[0]->GetEntry(0);
  const string waveName = massBin->waveName(waveIndex).Data();
  printInfo << "plotting wave intensity for wave '" << waveName << "' [" << waveIndex << "]";
  if (selectExpr != "")
    cout << " using selection criterion '" << selectExpr << "'";
  cout << endl;

  
  // build multiGraph
  TMultiGraph* graph = new TMultiGraph();
  graph->SetName (waveName.c_str());
  graph->SetTitle(waveName.c_str());
  if (graphTitle != "") {
    graph->SetName (graphTitle.c_str());
    graph->SetTitle(graphTitle.c_str());
  }
  double       maxY  = 0;
  for (unsigned int i = 0; i < nmbTrees; ++i) {
    // build and run TTree::Draw() expression
    stringstream drawExpr;
    drawExpr << "intensity(\"" << waveName << "\"):intensityErr(\"" << waveName << "\")"
	     << ":massBinCenter() >> h" << waveName << "_" << i;
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
    stringstream graphName;
    graphName << ((graphTitle == "") ? waveName : graphTitle) << "_" << i;
    g->SetName (graphName.str().c_str());
    g->SetTitle(graphName.str().c_str());
    g->SetMarkerStyle(21);
    g->SetMarkerSize(0.5);
    if (graphColors) {
      g->SetMarkerColor(graphColors[i]);
      g->SetLineColor  (graphColors[i]);
    }
    graph->Add(g);
    
    // compute maximum for y-axis
    for (int j = 0; j < nmbBins; ++j)
      if(maxY < (g->GetY()[j] + g->GetEY()[j]))
	maxY = g->GetY()[j] + g->GetEY()[j];
  }
  cout << "    maximum intensity for graph " << graph->GetName() << " is " << maxY << "." << endl;
    
  // draw graph
  graph->Draw(drawOption);
  graph->GetXaxis()->SetTitle("Mass [GeV]");
  graph->GetYaxis()->SetTitle("Intensity");
  graph->GetYaxis()->SetRangeUser(-maxY * 0.1, maxY * 1.1);
  TLine line;
  line.SetLineStyle(3);
  line.DrawLine(graph->GetXaxis()->GetXmin(), 0, graph->GetXaxis()->GetXmax(), 0);

  // create EPS file
  if (saveEps)
    gPad->SaveAs((waveName + ".eps").c_str());

  return graph;
}
