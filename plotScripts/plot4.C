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

#include "reportingUtils.hpp"
#include "plotIntensity.h"
#include "plotPhase.h"
#include "plotCoherence.h"


using namespace std;
using namespace rpwa;


// .............................................................................
// signature with wave indices
void
plot4(const unsigned int nmbTrees,             // number of fitResult trees
      TTree**            trees,                // array of fitResult trees
      const int          waveIndexA,           // index of first wave
      const int          waveIndexB,           // index of second wave
      const bool         saveEps     = false,  // if set, EPS file with name wave ID is created
      const int*         graphColors = NULL,   // array of colors for graph line and marker
      const double*      graphScales = NULL,   // array of scale factors that scale all graphs of one tree
      const bool         drawLegend  = true,   // if set legend is drawn
      const std::string& graphTitle  = "",     // name and title of graph (default is wave IDs)
      const char*        drawOption  = "AP",   // draw option for graph
      const double       massMin     = 0,      // [GeV/c^2]
      const double       massMax     = 0,      // [GeV/c^2]
      const string&      branchName  = "fitResult_v2",
      TCanvas*           canvas      = NULL)   // fill a given canvas instead of creating one (name/title will be changed)
{
	for (unsigned int i = 0; i < nmbTrees; ++i)
		if (!trees[i]) {
			printErr << "null pointer to tree[" << i << "]. Aborting..." << endl;
			return;
		}

	// select mass range; convert from GeV/c^2 to MeV/c^2
	stringstream selectExpr;
	if ((massMin != 0) || (massMax != 0))
		selectExpr << "(massBinCenter() >= "<< massMin * 1000 << ") "
		           << "&& (massBinCenter() <= " << massMax * 1000 << ")";

	stringstream canvName;
	canvName << "4plot_" << waveIndexA << "_" << waveIndexB;
	// create a new canvas if no canvas given to draw into
	if (!canvas){
		canvas = new TCanvas(canvName.str().c_str(), canvName.str().c_str(), 10, 10, 1000, 800);
	} else {
		canvas->SetName(canvName.str().c_str());
		canvas->SetTitle(canvName.str().c_str());
		canvas->SetCanvasSize(1000, 800);
	}
	canvas->Divide(2, 2);

	// wave A intensity
	canvas->cd(1);
	plotIntensity(nmbTrees, trees, waveIndexA, false, graphColors, graphScales, drawLegend,
	              graphTitle, drawOption, 1, 0, selectExpr.str(), branchName);
	// wave A - wave B phase angle
	canvas->cd(2);
	plotPhase(nmbTrees, trees, waveIndexA, waveIndexB, false, graphColors, drawLegend,
	          graphTitle, drawOption, selectExpr.str(), branchName);
	// wave B intensity
	canvas->cd(3);
	plotIntensity(nmbTrees, trees, waveIndexB, false, graphColors, graphScales, drawLegend,
	              graphTitle, drawOption, 1, 0, selectExpr.str(), branchName);
	// wave A - wave B coherence
	canvas->cd(4);
	plotCoherence(nmbTrees, trees, waveIndexA, waveIndexB, false, graphColors, drawLegend,
	              graphTitle, drawOption, selectExpr.str(), branchName);

	// create EPS file
	if (saveEps)
		gPad->SaveAs(((string)canvas->GetName() + ".eps").c_str());
}


void
plot4(TTree*        tree,                 // fitResult tree
      const int     waveIndexA,           // index of first wave
      const int     waveIndexB,           // index of second wave
      const bool    saveEps    = false,   // if set, EPS file with name wave ID is created
      const int     graphColor = kBlack,  // array of colors for graph line and marker
      const bool    drawLegend = false,   // if set legend is drawn
      const string& graphTitle = "",      // name and title of graph (default is wave IDs)
      const char*   drawOption = "AP",    // draw option for graph
      const double  massMin    = 0,       // [GeV/c^2]
      const double  massMax    = 0,       // [GeV/c^2]
      const string& branchName = "fitResult_v2",
      TCanvas*      canvas     = NULL)    // fill a given canvas instead of creating one (name/title will be changed)
{
	plot4(1, &tree, waveIndexA, waveIndexB, saveEps, &graphColor, NULL, drawLegend,
	      graphTitle, drawOption, massMin, massMax, branchName, canvas);
}


// .............................................................................
// signature with wave names
void
plot4(const unsigned int nmbTrees,             // number of fitResult trees
      TTree**            trees,                // array of fitResult trees
      const string&      waveNameA,            // name of first wave
      const string&      waveNameB,            // name of second wave
      const bool         saveEps     = false,  // if set, EPS file with name wave ID is created
      const int*         graphColors = NULL,   // array of colors for graph line and marker
      const double*      graphScales = NULL,   // array of scale factors that scale all graphs of one tree
      const bool         drawLegend  = true,   // if set legend is drawn
      const std::string& graphTitle  = "",     // name and title of graph (default is wave IDs)
      const char*        drawOption  = "AP",   // draw option for graph
      const double       massMin     = 0,      // [GeV/c^2]
      const double       massMax     = 0,      // [GeV/c^2]
      const string&      branchName  = "fitResult_v2")
{
	if (!trees[0]) {
		printErr << "null pointer to tree. exiting." << endl;
		return;
	}
	// get wave indices (assumes same wave set in all trees)
	fitResult* massBin = new fitResult();
	trees[0]->SetBranchAddress(branchName.c_str(), &massBin);
	trees[0]->GetEntry(0);
	const int indexA = massBin->waveIndex(waveNameA);
	const int indexB = massBin->waveIndex(waveNameB);
	if ((indexA >= 0) && (indexB >= 0))
		return plot4(nmbTrees, trees, indexA, indexB, saveEps, graphColors, graphScales, drawLegend,
		             graphTitle, drawOption, massMin, massMax, branchName);
	printErr << "cannot find wave(s) in tree '" << trees[0]->GetName() << "'. exiting." << endl;
}


void
plot4(TTree*        tree,                 // fitResult tree
      const string& waveNameA,            // name of first wave
      const string& waveNameB,            // name of second wave
      const bool    saveEps    = false,   // if set, EPS file with name wave ID is created
      const int     graphColor = kBlack,  // array of colors for graph line and marker
      const bool    drawLegend = false,   // if set legend is drawn
      const string& graphTitle = "",      // name and title of graph (default is wave IDs)
      const char*   drawOption = "AP",    // draw option for graph
      const double  massMin    = 0,       // [GeV/c^2]
      const double  massMax    = 0,       // [GeV/c^2]
      const string& branchName = "fitResult_v2")
{
	plot4(1, &tree, waveNameA, waveNameB, saveEps, &graphColor, NULL, drawLegend,
	      graphTitle, drawOption, massMin, massMax, branchName);
}
