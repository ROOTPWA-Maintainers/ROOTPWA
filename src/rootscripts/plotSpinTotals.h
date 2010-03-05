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
//      Creates summary plots with intensities of spin totals
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#ifndef plotSpinTotals_hh
#define plotSpinTotals_hh


#include <string>
#include <vector>

#include "TVirtualPad.h"
#include "TTree.h"


std::vector<std::pair<string, TVirtualPad*> >
plotSpinTotals(const unsigned int nmbTrees,              // number of fitResult trees
	       TTree**            trees,                 // array of fitResult trees
	       const int*         colors        = 0,     // array of line and marker colors
	       const double       yAxisRangeMax = 0,     // if != 0; range of y-axis is not allowed to be larger than this value
	       const bool         drawLegend    = true,  // if set legend is drawn
	       const std::string& outFileName   = "spinTotals.root",
	       const std::string& branchName    = "fitResult_v2");


std::vector<std::pair<string, TVirtualPad*> >
plotSpinTotals(std::vector<TTree*>&    trees,                 // array of fitResult trees
	       const std::vector<int>& colors,                // array of line and marker colors
	       const double            yAxisRangeMax = 0,     // if != 0; range of y-axis is not allowed to be larger than this value
	       const bool              drawLegend    = true,  // if set legend is drawn
	       const std::string&      outFileName   = "spinTotals.root",
	       const std::string&      branchName    = "fitResult_v2")
{
  return plotSpinTotals(trees.size(), &(*(trees.begin())), &(*(colors.begin())),
			yAxisRangeMax, drawLegend, outFileName, branchName);
}


std::vector<std::pair<std::string, TVirtualPad*> >
plotSpinTotals(TTree*             tree,                    // fitResult tree
	       const int          color         = kBlack,  // color of line and marker
	       const double       yAxisRangeMax = 0,       // if != 0; range of y-axis is not allowed to be larger than this value
	       const bool         drawLegend    = true,    // if set legend is drawn
	       const std::string& outFileName   = "spinTotals.root",
	       const std::string& branchName    = "fitResult_v2")
{
  return plotSpinTotals(1, &tree, &color, yAxisRangeMax, drawLegend, outFileName, branchName);
}


#endif  // plotSpinTotals_hh
