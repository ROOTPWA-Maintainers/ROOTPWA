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
//      Draws coherence of (wave A - wave B) from tree.
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#ifndef PLOTCOHERENCE_HH
#define PLOTCOHERENCE_HH


#include <string>

#include "TTree.h"
#include "TMultiGraph.h"

#include "fitResult.h"


// .............................................................................
// signature with wave indices
TMultiGraph*
plotCoherence(const unsigned int nmbTrees,             // number of fitResult trees
              TTree**            trees,                // array of fitResult trees
              const int          waveIndexA,           // index of first wave
              const int          waveIndexB,           // index of second wave
              const bool         saveEps     = false,  // if set, EPS file with name wave ID is created
              const int*         graphColors = NULL,   // array of colors for graph line and marker
              const bool         drawLegend  = true,   // if set legend is drawn
              const std::string& graphTitle  = "",     // name and title of graph (default is wave IDs)
              const char*        drawOption  = "AP",   // draw option for graph
              const std::string& selectExpr  = "",     // TTree::Draw() selection expression
              const std::string& branchName  = "fitResult_v2");  // fitResult branch name


TMultiGraph*
plotCoherence(TTree*             tree,                 // fitResult tree
              const int          waveIndexA,           // index of first wave
              const int          waveIndexB,           // index of second wave
              const bool         saveEps    = false,   // if set, EPS file with name waveId is created
              const int          graphColor = kBlack,  // color of line and marker
              const bool         drawLegend = false,   // if set legend is drawn
              const std::string& graphTitle = "",      // name and title of graph
              const char*        drawOption = "AP",    // draw option for graph
              const std::string& selectExpr = "",      // TTree::Draw() selection expression
              const std::string& branchName = "fitResult_v2")
{
	return plotCoherence(1, &tree, waveIndexA, waveIndexB, saveEps, &graphColor,
	                     drawLegend, graphTitle, drawOption, selectExpr, branchName);
}


// .............................................................................
// signature with wave names
TMultiGraph*
plotCoherence(const unsigned int nmbTrees,             // number of fitResult trees
              TTree**            trees,                // array of fitResult trees
              const std::string& waveNameA,            // name of first wave
              const std::string& waveNameB,            // name of second wave
              const bool         saveEps     = false,  // if set, EPS file with name wave ID is created
              const int*         graphColors = NULL,   // array of colors for graph line and marker
              const bool         drawLegend  = true,   // if set legend is drawn
              const std::string& graphTitle  = "",     // name and title of graph (default is wave IDs)
              const char*        drawOption  = "AP",   // draw option for graph
              const std::string& selectExpr  = "",     // TTree::Draw() selection expression
              const std::string& branchName  = "fitResult_v2")  // fitResult branch name
{
	if (!trees[0]) {
		printErr << "null pointer to tree. exiting." << std::endl;
		return 0;
	}
	// get wave indices (assumes same wave set in all trees)
	rpwa::fitResult* massBin = new rpwa::fitResult();
	trees[0]->SetBranchAddress(branchName.c_str(), &massBin);
	trees[0]->GetEntry(0);
	const int indexA = massBin->waveIndex(waveNameA);
	const int indexB = massBin->waveIndex(waveNameB);
	if ((indexA >= 0) && (indexB >= 0))
		return plotCoherence(nmbTrees, trees, indexA, indexB, saveEps, graphColors, drawLegend,
		                     graphTitle, drawOption, selectExpr, branchName);
	printErr << "cannot find wave(s) in tree '" << trees[0]->GetName() << "'. exiting." << std::endl;
	return 0;
}


TMultiGraph*
plotCoherence(TTree*             tree,                 // fitResult tree
              const std::string& waveNameA,            // name of first wave
              const std::string& waveNameB,            // name of second wave
              const bool         saveEps    = false,   // if set, EPS file with name waveId is created
              const int          graphColor = kBlack,  // color of line and marker
              const bool         drawLegend = false,   // if set legend is drawn
              const std::string& graphTitle = "",      // name and title of graph
              const char*        drawOption = "AP",    // draw option for graph
              const std::string& selectExpr = "",      // TTree::Draw() selection expression
              const std::string& branchName = "fitResult_v2")
{
	return plotCoherence(1, &tree, waveNameA, waveNameB, saveEps, &graphColor,
	                     drawLegend, graphTitle, drawOption, selectExpr, branchName);
}


#endif  // PLOTCOHERENCE_HH
