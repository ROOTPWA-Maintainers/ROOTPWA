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

#include "TGraphErrors.h"
#include "TAxis.h"
#include "TPad.h"
#include "TLegend.h"
#include "TList.h"

#include "utilities.h"
#include "plotPhase.h"


using namespace std;
using namespace rpwa;


// signature with wave names
TMultiGraph*
plotPhase(const unsigned int nmbTrees,     // number of fitResult trees
          TTree**            trees,        // array of fitResult trees
          const int          waveIndexA,   // index of first wave
          const int          waveIndexB,   // index of second wave
          const bool         saveEps,      // if set, EPS file with name wave ID is created
          const int*         graphColors,  // array of colors for graph line and marker
          const bool         drawLegend,   // if set legend is drawn
          const string&      graphTitle,   // name and title of graph (default is wave IDs)
          const char*        drawOption,   // draw option for graph
          const string&      selectExpr,   // TTree::Draw() selection expression
          const string&      branchName)   // fitResult branch name
{
	for (unsigned int i = 0; i < nmbTrees; ++i)
		if (!trees[i]) {
			printErr << "null pointer to tree[" << i << "]. aborting." << endl;
			return 0;
		}
	// get wave names (assume same wave set in all trees)
	fitResult* massBin = new fitResult();
	trees[0]->SetBranchAddress(branchName.c_str(), &massBin);
	trees[0]->GetEntry(0);
	const string waveNameA = massBin->waveName(waveIndexA).Data();
	const string waveNameB = massBin->waveName(waveIndexB).Data();
	printInfo << "plotting phase between wave '" << waveNameA << "' [" << waveIndexA << "] "
	          << "and wave '" << waveNameB << "' [" << waveIndexB << "]" << endl;
	if (selectExpr != "")
		cout << "    using selection criterion '" << selectExpr << "'" << endl;

	// create multiGraph
	TMultiGraph* graph = new TMultiGraph();
	{
		if (graphTitle != "") {
			graph->SetName (graphTitle.c_str());
			graph->SetTitle(graphTitle.c_str());
		} else {
			stringstream name;
			name << "phase_" << waveNameA << "_" << waveNameB;
			stringstream title;
			title << "#Delta#varphi(" << waveNameA << " [" << waveIndexA << "], "
			      << waveNameB << " [" << waveIndexB << "])";
			graph->SetName(name.str().c_str());
			graph->SetTitle(title.str().c_str());
		}
	}

	// fill multiGraph
	for (unsigned int i = 0; i < nmbTrees; ++i) {

		// builf and run TTree::Draw() expression
		stringstream drawExpr;
		drawExpr << branchName << ".phase("     << waveIndexA << "," << waveIndexB << "):"
		         << branchName << ".phaseErr(" << waveIndexA << "," << waveIndexB << "):"
		         << branchName << ".massBinCenter() >> h" << waveIndexA << "_" << waveIndexB
		         << "_" << i;
		cout << "    running TTree::Draw() expression '" << drawExpr.str() << "' "
		     << "on tree '" << trees[i]->GetName() << "', '" << trees[i]->GetTitle() << "'" << endl;
		trees[i]->Draw(drawExpr.str().c_str(), selectExpr.c_str(), "goff");

		// extract data from TTree::Draw() result and build graph
		const int nmbBins = trees[i]->GetSelectedRows();
		vector<double> x(nmbBins), xErr(nmbBins);
		for (int j = 0; j < nmbBins; ++j) {
			x   [j] = trees[i]->GetV3()[j] * 0.001;  // convert mass to GeV
			xErr[j] = 0;
		}
		TGraphErrors* g = new TGraphErrors(nmbBins,
		                                   &(*(x.begin())),     // mass
		                                   trees[i]->GetV1(),   // phase
		                                   &(*(xErr.begin())),  // mass error
		                                   trees[i]->GetV2());  // phase error
		
		// beautify graph
		stringstream graphName;
		graphName << graph->GetName() << "_" << i;
		g->SetName (graphName.str().c_str());
		g->SetTitle(graphName.str().c_str());
		g->SetMarkerStyle(21);
		g->SetMarkerSize(0.5);
		if (graphColors) {
			g->SetMarkerColor(graphColors[i]);
			g->SetLineColor  (graphColors[i]);
		}
		graph->Add(g);
	}

	graph->SetMinimum(-240);
	graph->SetMaximum( 240);
	// draw graph
	graph->Draw(drawOption);
	graph->GetXaxis()->SetTitle("Mass [GeV]");
	graph->GetYaxis()->SetTitle("Phase Angle [deg]");
	gPad->Update();

	// add legend
	if (drawLegend && (nmbTrees > 1)) {
		TLegend* legend = new TLegend(0.65,0.80,0.99,0.99);
		legend->SetFillColor(10);
		legend->SetBorderSize(1);
		legend->SetMargin(0.2);
		for (unsigned int i = 0; i < nmbTrees; ++i) {
			TGraph* g = static_cast<TGraph*>(graph->GetListOfGraphs()->At(i));
			legend->AddEntry(g, trees[i]->GetTitle(), "LPE");
		}
		legend->Draw();
	}
  
	// create EPS file
	if (saveEps)
		gPad->SaveAs(((string)graph->GetName() + ".eps").c_str());

	return graph;
}
