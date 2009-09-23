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

#include "TTree.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TLine.h"


using namespace std;


void
plotSpinTotals(TTree*    tree,  // TFitBin tree
	       const int color = kBlack)
{
  const unsigned int nmbPadsPerCanvMin = 6;  // minimum number of pads each canvas is subdivided into
  // define set of spin totals
  const string waves[] = {"",  // total intensity
			  "flat",
			  "0++0-",
			  "0-+0+",
			  "1++0+",
			  "2-+0+",
			  "2++0-",
			  "2++1+",
			  "1-+0-",
			  "1-+1+",
			  "3++0+",
			  "4++1+",
			  "3-+1+",
			  "3-+1-",
			  "3-+0-"};
  const unsigned int nmbWaves = sizeof(waves) / sizeof(string);
  cout << "Plotting spin totals for:" << endl;
  for (unsigned int i = 0; i < nmbWaves; ++i)
    cout << "    " << ((waves[i] != "") ? waves[i] : "total") << endl;

  // set plot style
  gStyle->SetOptStat(0);
  gStyle->SetMarkerStyle(23);
  gStyle->SetMarkerSize(0.5);
  gStyle->SetMarkerColor(color);
  gStyle->SetLineColor(color);
  gStyle->SetFillColor(0);
  gStyle->SetPadColor(0);

  // calculate optimum canvas subdivision
  const unsigned int nmbPadsVert    = (int)ceil(sqrt(nmbPadsPerCanvMin));
  const unsigned int nmbPadsHor     = (int)ceil((double)nmbPadsPerCanvMin / (double)nmbPadsVert);
  const unsigned int nmbPadsPerCanv = nmbPadsHor * nmbPadsVert;

  // plot spin totals
  unsigned int countPad  = 0;
  unsigned int countCanv = 0;
  TCanvas*     canv      = 0;
  string       drawOpt;
  for (unsigned int i = 0; i < nmbWaves; ++i) {
    // create new pad, if necessary
    if (countPad == 0) {
      stringstream canvName;
      canvName << "spinTotals" << countCanv;
      canv = static_cast<TCanvas*>(gROOT->FindObject(canvName.str().c_str()));
      if (!canv) {
	canv = new TCanvas(canvName.str().c_str(), "Spin Totals", 10, 10, 1000, 900);
	canv->Divide(nmbPadsHor, nmbPadsVert);
	drawOpt = "APZ";
      } else
	drawOpt = "PZ SAME";
    }
    
    // build and run TTree::Draw() expression
    const string drawExpr = "intensity(\"" + waves[i] + "\"):intensityErr(\"" + waves[i] + "\"):massBinCenter()";
    tree->Draw(drawExpr.c_str(), "", "goff");

    // extract data from TTree::Draw() result and build graph
    const int nmbBins = tree->GetSelectedRows();
    vector<double> x(nmbBins), xErr(nmbBins);
    for (int j = 0; j < nmbBins; ++j) {
      x[j]    = tree->GetV3()[j] * 0.001;  // convert mass to GeV
      xErr[j] = 0;
    }
    TGraphErrors* g = new TGraphErrors(nmbBins,
				       &(*(x.begin())),  // mass
				       tree->GetV1(),    // intensity
				       0,                // mass error
				       tree->GetV2());   // intensity error
    // plot graph
    canv->cd(++countPad);
    if (waves[i] != "")
      g->SetTitle(waves[i].c_str());
    else
      g->SetTitle("total");
    g->GetXaxis()->SetTitle("Mass [GeV]");
    g->GetYaxis()->SetTitle("Intensity");
    // compute maximum for y-axis
    double maxY = 0;
    for (int i = 0; i < nmbBins; ++i)
      if(maxY < g->GetY()[i])
	maxY = g->GetY()[i];
    g->SetMinimum(-maxY * 0.1);
    g->SetMaximum( maxY * 1.1);
    g->Draw(drawOpt.c_str());
    TLine line;
    line.SetLineStyle(3);
    line.DrawLine(g->GetXaxis()->GetXmin(), 0, g->GetXaxis()->GetXmax(), 0);

    if (countPad == nmbPadsPerCanv) {
      canv->Update();
      countPad = 0;
      ++countCanv;
    }
  }
}
