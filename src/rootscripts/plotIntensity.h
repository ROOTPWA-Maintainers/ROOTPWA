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


//
// plots wave intensities for a give number of trees on top of each
// other
//


#ifndef PLOTINTENSITY_HH
#define PLOTINTENSITY_HH


#include <string>
#include <vector>

#include "TTree.h"
#include "TMultiGraph.h"

#include "fitResult.h"


// .............................................................................
// signatures with wave index
TMultiGraph*
plotIntensity(const unsigned int nmbTrees,               // number of fitResult trees
              TTree**            trees,                  // array of fitResult trees
              const int          waveIndex,              // wave index
              const bool         saveEps       = false,  // if set, EPS file with name wave ID is created
              const int*         graphColors   = NULL,   // array of colors for graph line and marker
              const double*      graphScales   = NULL,   // array of scales for graphgroups
              const bool         drawLegend    = true,   // if set legend is drawn
              const std::string& graphTitle    = "",     // name and title of graph (default is wave ID)
              const char*        drawOption    = "AP",   // draw option for graph
              const double       normalization = 1,      // scale factor for intensities
              const double       yAxisRangeMax = 0,      // if != 0; range of y-axis is limited to this value
              const std::string& selectExpr    = "",     // TTree::Draw() selection expression
              const std::string& branchName    = "fitResult_v2");  // fitResult branch name


inline
TMultiGraph*
plotIntensity(std::vector<TTree*>&        trees,                  // array of fitResult trees
              const int                   waveIndex,              // wave index
              const bool                  saveEps       = false,  // if set, EPS file with name wave ID is created
              const std::vector<int>&     graphColors   = std::vector<int>(),  // array of colors for graph line and marker
              const std::vector<double>&  graphScales   = std::vector<double>(),  // array of scales for graphgroups
              const bool                  drawLegend    = true,   // if set legend is drawn
              const std::string&          graphTitle    = "",     // name and title of graph (default is wave ID)
              const char*                 drawOption    = "AP",   // draw option for graph
              const double                normalization = 1,      // scale factor for intensities
              const double                yAxisRangeMax = 0,      // if != 0; range of y-axis is limited to this value
              const std::string&          selectExpr    = "",     // TTree::Draw() selection expression
              const std::string&          branchName    = "fitResult_v2")  // fitResult branch name
{
	return plotIntensity(trees.size(), &(*(trees.begin())), waveIndex, saveEps,
	                     &(*(graphColors.begin())), &(*(graphScales.begin())), drawLegend, graphTitle, drawOption,
	                     normalization, yAxisRangeMax, selectExpr, branchName);
}


inline
TMultiGraph*
plotIntensity(TTree*             tree,                    // fitResult tree
              const int          waveIndex,               // wave index
              const bool         saveEps       = false,   // if set, EPS file with name wave ID is created
              const int          graphColor    = kBlack,  // color of line and marker

              const bool         drawLegend    = false,   // if set legend is drawn
              const std::string& graphTitle    = "",      // name and title of graph (default is wave ID)
              const char*        drawOption    = "AP",    // draw option for graph
              const double       normalization = 1,       // scale factor for intensities
              const double       yAxisRangeMax = 0,       // if != 0; range of y-axis is limited to this value
              const std::string& selectExpr    = "",      // TTree::Draw() selection expression
              const std::string& branchName    = "fitResult_v2")  // fitResult branch name
{
	return plotIntensity(1, &tree, waveIndex, saveEps, &graphColor, NULL, drawLegend, graphTitle,
	                     drawOption, normalization,  yAxisRangeMax, selectExpr, branchName);
}


// .............................................................................
// signatures with wave name
// in contrast to other methods, it does not assume
// to have all the same waves in all bins
//inline
TMultiGraph*
plotIntensity(const unsigned int nmbTrees,               // number of fitResult trees
              TTree**            trees,                  // array of fitResult trees
              const std::string& waveName,               // wave name
              const bool         saveEps       = false,  // if set, EPS file with name wave ID is created
              const int*         graphColors   = NULL,   // array of colors for graph line and marker
              const double*      graphScales   = NULL,   // array of scales for graphgroups
              const bool         drawLegend    = true,   // if set legend is drawn
              const std::string& graphTitle    = "",     // name and title of graph (default is wave ID)
              const char*        drawOption    = "AP",   // draw option for graph
              const double       normalization = 1,      // scale factor for intensities
              const double       yAxisRangeMax = 0,      // if != 0; range of y-axis is limited to this value
              const std::string& selectExpr    = "",     // TTree::Draw() selection expression
              const std::string& branchName    = "fitResult_v2")  // fitResult branch name
	;/*{
	   if (!trees[0]) {
	   printErr << "null pointer to tree[" << 0 << "]. aborting." << endl;
	   return 0;
	   }
	   // get wave index (assumes same wave set in all trees)
	   rpwa::fitResult* massBin = new rpwa::fitResult();
	   trees[0]->SetBranchAddress(branchName.c_str(), &massBin);
	   trees[0]->GetEntry(0);
	   const int index = massBin->waveIndex(waveName);
	   if (index >= 0)
	   return plotIntensity(nmbTrees, trees, index, saveEps, graphColors, drawLegend,  graphTitle,
	   drawOption, normalization, yAxisRangeMax, selectExpr, branchName);
	   printErr << "cannot find wave '" << waveName << "' "
	   << "in tree '" << trees[0]->GetName() << "'. aborting." << endl;
	   return 0;
	   }*/


inline
TMultiGraph*
plotIntensity(std::vector<TTree*>&        trees,                  // array of fitResult trees
              const std::string&          waveName,               // wave name
              const bool                  saveEps       = false,  // if set, EPS file with name wave ID is created
              const std::vector<int>&     graphColors   = std::vector<int>(),  // array of colors for graph line and marker
              const std::vector<double>&  graphScales   = std::vector<double>(),  // array of scales for graphgroups
              const bool                  drawLegend    = true,   // if set legend is drawn
              const std::string&          graphTitle    = "",     // name and title of graph (default is wave ID)
              const char*                 drawOption    = "AP",   // draw option for graph
              const double                normalization = 1,      // scale factor for intensities
              const double                yAxisRangeMax = 0,      // if != 0; range of y-axis is limited to this value
              const std::string&          selectExpr    = "",     // TTree::Draw() selection expression
              const std::string&          branchName    = "fitResult_v2")  // fitResult branch name
{
	return plotIntensity(trees.size(), &(*(trees.begin())), waveName, saveEps,
	                     &(*(graphColors.begin())), &(*(graphScales.begin())), drawLegend, graphTitle, drawOption,
	                     normalization, yAxisRangeMax, selectExpr, branchName);
}


inline
TMultiGraph*
plotIntensity(TTree*             tree,                    // fitResult tree
              const std::string& waveName,                // wave name
              const bool         saveEps       = false,   // if set, EPS file with name wave ID is created
              const int          graphColor    = kBlack,  // color of line and marker
              const bool         drawLegend    = true,   // if set legend is drawn
              const std::string& graphTitle    = "",      // name and title of graph (default is wave ID)
              const char*        drawOption    = "AP",    // draw option for graph
              const double       normalization = 1,       // scale factor for intensities
              const double       yAxisRangeMax = 0,       // if != 0; range of y-axis is limited to this value
              const std::string& selectExpr    = "",      // TTree::Draw() selection expression
              const std::string& branchName    = "fitResult_v2")  // fitResult branch name
{
	return plotIntensity(1, &tree, waveName, saveEps, &graphColor, NULL, drawLegend, graphTitle,
	                     drawOption, normalization, yAxisRangeMax, selectExpr, branchName);
}


#endif  //PLOTINTENSITY_HH
