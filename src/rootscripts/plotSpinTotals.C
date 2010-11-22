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


#include <iostream>
#include <sstream>
#include <vector>

#include "TFile.h"
#include "TStyle.h"
#include "TString.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLegend.h"

#include "reportingUtils.hpp"
#include "fitResult.h"
#include "plotSpinTotals.h"


using namespace std;
using namespace rpwa;


vector<pair<string, TVirtualPad*> >
plotSpinTotals(const unsigned int nmbTrees,       // number of fitResult trees
               TTree**            trees,          // array of fitResult trees
               const int*         colors,         // array of line and marker colors
               const double       yAxisRangeMax,  // if != 0; range of y-axis is not allowed to be larger than this value
               const bool         drawLegend,     // if set legend is drawn
               const string&      outFileName,
               const string&      branchName)
{
	vector<pair<string, TVirtualPad*> > wavePads;

	for (unsigned int i = 0; i < nmbTrees; ++i)
		if (!trees[i]) {
			printErr << "null pointer to tree " << i << ". exiting." << endl;
			return wavePads;
		}

	TFile* outFile = NULL;
	if (outFileName != "")
		outFile = TFile::Open(outFileName.c_str(), "RECREATE");
	if (!outFile || outFile->IsZombie()) {
		printErr << "could not create file '" << outFileName << "'. exiting." << endl;
		return wavePads;
	}

	const unsigned int nmbPadsPerCanvMin = 6;  // minimum number of pads each canvas is subdivided into
	// define set of spin totals
	const string waves[] = {"logLikelihood",
	                        ".*",  // total intensity
	                        "flat",
	                        "0-\\+0\\+",
	                        "1\\+\\+0\\+",
	                        //"1\\+\\+1\\+",
	                        "2-\\+0\\+",
	                        "2-\\+1\\+",
	                        "2\\+\\+0-",
	                        "2\\+\\+1\\+",
	                        "1-\\+0-",
	                        "1-\\+1\\+",
	                        "3\\+\\+0\\+",
	                        "4\\+\\+1\\+",
	                        "3-\\+1\\+",
	                        "3-\\+1-",
	                        "3-\\+0-",
	                        "0\\+\\+0-"};
	//   const string waves[] = {"logLikelihood",
	//                           ".*",  // total intensity
	//   			  "flat",
	// 			  "0-\\+",
	// 			  "1\\+\\+",
	// 			  "1-\\+",
	// 			  "2\\+\\+",
	// 			  "2-\\+",
	// 			  "3\\+\\+",
	// 			  "4\\+\\+",
	// 			  "4-\\+",
	// 			  "0-\\+0\\+",
	// 			  "1\\+\\+0\\+",
	// 			  "1\\+\\+1\\+",
	// 			  "1\\+\\+1-",
	// 			  "1-\\+0-",
	// 			  "1-\\+1\\+",
	// 			  "1-\\+1-",
	// 			  "2\\+\\+0-",
	// 			  "2\\+\\+1\\+",
	// 			  "2\\+\\+1-",
	// 			  "2-\\+0\\+",
	// 			  "2-\\+1\\+",
	// 			  "2-\\+1-",
	// 			  "3\\+\\+0\\+",
	// 			  "3\\+\\+1\\+",
	// 			  "4\\+\\+1\\+",
	// 			  "4-\\+0\\+",
	// 			  "4-\\+1\\+"};
	// const string waves[] = {"logLikelihood",
	//                         ".*",  // total intensity
	// 			  "flat",
	// 			  "1--",
	// 			  "1\\+-",
	// 			  "2--",
	// 			  "3--",
	// 			  "3\\+-",
	// 			  "1--0\\+",
	// 			  "1--1\\+",
	// 			  "1\\+-1\\+",
	// 			  "2--1\\+",
	// 			  "2--2\\+",
	// 			  "3--0\\+",
	// 			  "3--1\\+",
	// 			  "3--2\\+",
	// 			  "3--3\\+",
	// 			  "3\\+-1\\+",
	// 			  "3\\+-2\\+",
	// 			  "3\\+-3\\+"};
	const unsigned int nmbWaves = sizeof(waves) / sizeof(string);
	printInfo << "plotting spin totals for:" << endl;
	for (unsigned int i = 0; i < nmbWaves; ++i)
		cout << "    " << ((waves[i] != "") ? waves[i] : "total") << endl;

	// set plot style
	gStyle->SetOptStat(0);
	gStyle->SetFillColor(0);
	gStyle->SetPadColor(0);

	// calculate optimum canvas subdivision
	const unsigned int nmbPadsVert    = (int)ceil(sqrt(nmbPadsPerCanvMin));
	const unsigned int nmbPadsHor     = (int)ceil((double)nmbPadsPerCanvMin / (double)nmbPadsVert);
	const unsigned int nmbPadsPerCanv = nmbPadsHor * nmbPadsVert;

	// plot spin totals
	wavePads.resize(nmbWaves, pair<string, TVirtualPad*>("", NULL));
	unsigned int countPad  = 0;
	unsigned int countCanv = 0;
	TCanvas*     canv      = 0;
	string       drawOpt;
	for (unsigned int i = 0; i < nmbWaves; ++i) {
		const TString waveName = unescapeRegExpSpecialChar((waves[i] != ".*") ? waves[i] : "total");
		// create new pad, if necessary
		if (countPad == 0) {
			stringstream canvName;
			canvName << "spinTotals" << countCanv;
			canv = static_cast<TCanvas*>(gROOT->FindObject(canvName.str().c_str()));
			if (!canv) {
				canv = new TCanvas(canvName.str().c_str(), "Spin Totals", 10, 10, 1000, 900);
				canv->Divide(nmbPadsHor, nmbPadsVert);
				drawOpt = "AP";
			} else
				drawOpt = "P SAME";
		}
    
		// build multiGraph
		TMultiGraph* graph = new TMultiGraph();
		double maxY = 0;
		for (unsigned int j = 0; j < nmbTrees; ++j) {
			// build and run TTree::Draw() expression
			string drawExpr;
			if (waves[i] == "logLikelihood")  // draw likelihood
				drawExpr = "-" + branchName + ".logLikelihood() / " + branchName + ".nmbEvents():0:"
					+ branchName + ".massBinCenter()";
			else
				drawExpr = branchName + ".intensity(\"" + waves[i] + "\"):"
					+ branchName + ".intensityErr(\"" + waves[i] + "\"):"
					+ branchName + ".massBinCenter()";
			cout << "    running TTree::Draw() expression '" << drawExpr << "' "
			     << "on tree '" << trees[j]->GetName() << "', '" << trees[j]->GetTitle() << "'" << endl;
			trees[j]->Draw(drawExpr.c_str(), "", "goff");

			// extract data from TTree::Draw() result and build graph
			const int nmbBins = trees[j]->GetSelectedRows();
			vector<double> x(nmbBins), xErr(nmbBins);
			vector<double> y(nmbBins), yErr(nmbBins);
			for (int k = 0; k < nmbBins; ++k) {
				x   [k] = trees[j]->GetV3()[k] * 0.001;  // convert mass to GeV
				xErr[k] = 0;
				y   [k] = trees[j]->GetV1()[k];
				yErr[k] = (waves[i] == "logLikelihood") ? 0 : trees[j]->GetV2()[k];
			}
			TGraphErrors* g = new TGraphErrors(nmbBins,
			                                   &(*(x.begin())),      // mass
			                                   &(*(y.begin())),      // intensity
			                                   &(*(xErr.begin())),   // mass error
			                                   &(*(yErr.begin())));  // intensity error

			{
				stringstream graphName;
				graphName << waveName << "_" << j;
				g->SetName (graphName.str().c_str());
				g->SetTitle(graphName.str().c_str());
			}
			g->SetMarkerStyle(21);
			g->SetMarkerSize(0.5);
			if (colors) {
				g->SetMarkerColor(colors[j]);
				g->SetLineColor  (colors[j]);
			}
			graph->Add(g, "P");

			// compute maximum for y-axis
			for (int k = 0; k < nmbBins; ++k) {
				const double val = y[k] + yErr[k];
				if (maxY < val)
					maxY = val;
			}
		}
		cout << "    maximum intensity for graph " << graph->GetName() << " is " << maxY << endl;

		canv->cd(++countPad);
		graph->SetTitle(waveName);
		graph->SetName(waveName);
		if ((yAxisRangeMax > 0) && (maxY > yAxisRangeMax))
			maxY = yAxisRangeMax;
		graph->SetMinimum(-maxY * 0.1);
		graph->SetMaximum( maxY * 1.1);
		// draw graph
		graph->Draw(drawOpt.c_str());
		graph->GetXaxis()->SetTitle("Mass [GeV]");
		if (waves[i] == "logLikelihood")
			graph->GetYaxis()->SetTitle("ln(Likelihood)");
		else
			graph->GetYaxis()->SetTitle("Intensity");
		TLine* line = new TLine(graph->GetXaxis()->GetXmin(), 0, graph->GetXaxis()->GetXmax(), 0);
		line->SetLineStyle(3);
		graph->GetListOfFunctions()->Add(line);
		// add legend
		if (drawLegend && (nmbTrees > 1)) {
			TLegend* legend = new TLegend(0.65,0.80,0.99,0.99);
			legend->SetFillColor(10);
			legend->SetBorderSize(1);
			legend->SetMargin(0.2);
			for (unsigned int j = 0; j < nmbTrees; ++j) {
				TGraph* g = static_cast<TGraph*>(graph->GetListOfGraphs()->At(j));
				legend->AddEntry(g, trees[j]->GetTitle(), "LPE");
			}
			graph->GetListOfFunctions()->Add(legend);
		}
    
		// memorize pad
		wavePads[i].first  = waves[i];
		wavePads[i].second = gPad;

		if (outFile)
			graph->Write();

		if (countPad == nmbPadsPerCanv) {
			canv->Update();
			countPad = 0;
			++countCanv;
		}
	}

	if (outFile)
		outFile->Close();
  
	return wavePads;
}
