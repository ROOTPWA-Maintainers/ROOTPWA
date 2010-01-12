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


#ifndef plotAllIntensities_hh
#define plotAllIntensities_hh


#include <string>
#include <vector>

#include "TVirtualPad.h"
#include "TTree.h"


// plots all waves sorted by JPC and returns a list of wave names and
// pointers to the pads they were drawn into
std::vector<std::pair<std::string, TVirtualPad*> >
plotAllIntensities(TTree*             tree,                  // TFitResult tree
		   const bool         createPsFile = false,  // if true, plots are written to waves.ps
		   const std::string& outPath      = "./",   // path for output files
		   const bool         mcmc         = false);


std::vector<std::pair<std::string, TVirtualPad*> >
plotAllIntensities(const unsigned int nmbTrees,              // number of TFitResult trees
		   TTree**            trees,                 // array of TFitResult trees
		   const bool         createPsFile = false,  // if true, plots are written to waves.ps
		   const std::string& outPath      = "./",   // path for output files
		   const int*         graphColors  = NULL,   // array of colors for graph line and marker
		   const bool         drawLegend   = true);  // if set legend is drawn


inline
std::vector<std::pair<std::string, TVirtualPad*> >
plotAllIntensities(std::vector<TTree*>&    trees,                 // array of TFitResult trees
		   const bool              createPsFile = false,  // if true, plots are written to waves.ps
		   const std::string&      outPath      = "./",   // path for output files
		   const std::vector<int>& graphColors  = std::vector<int>(),  // array of colors for graph line and marker
		   const bool              drawLegend   = true)   // if set legend is drawn
{
  return plotAllIntensities(trees.size(), &(*(trees.begin())), createPsFile,
			    outPath, &(*(graphColors.begin())), drawLegend);
}


#endif  // plotAllIntensities_hh
