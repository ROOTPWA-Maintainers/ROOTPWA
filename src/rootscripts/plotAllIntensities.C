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
//      Draws intensity graphs for all waves.
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#include <iostream>
#include <sstream>
#include <map>
#include <algorithm>
#include <cmath>

#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLatex.h"
#include "TList.h"
#include "TPostScript.h"
#include "TSystem.h"
#include "TLegend.h"

#include "fitResult.h"
#include "plotIntensity.h"
#include "plotAllIntensities.h"


using namespace std;
using namespace rpwa;


string
buildCanvName(const string& jpc,
	      const int     index)
{
  string name;
  {
    stringstream n;
    n << jpc << index;
    name = n.str();
  }
  // remove all spaces
  string::size_type pos = name.find(" ");
  while (pos != string::npos) {
    name.replace(pos, 1, "");
    pos = name.find(" ", pos + 1);
  }
  return name;
}


// predicate for sort
bool
compareIntensities(const pair<string, double>& a,
		   const pair<string, double>& b)
{
  return a.second > b.second;
}


vector<pair<string, TVirtualPad*> >
plotAllIntensities(const unsigned int nmbTrees,       // number of fitResult trees
		   TTree**            trees,          // array of fitResult trees
		   const bool         createPsFile,   // if true, plots are written to waves.ps
		   const string&      outPath,        // path for output files
		   const int*         graphColors,    // array of colors for graph line and marker
		   const bool         drawLegend,     // if set legend is drawn
		   const double       yAxisRangeMax,  // if != 0; range of y-axis is limited to this value
		   const string&      branchName)
{
  const double intensityThr      = 0;            // threshold for total intensity in mass bin
  const int    nmbPadsPerCanvMin = 4;            // minimum number of pads each canvas is subdivided into
  vector<pair<string, TVirtualPad*> > wavePads;  // return value

  for (unsigned int i = 0; i < nmbTrees; ++i)
    if (!trees[i]) {
      printErr << "null pointer to tree " << i << ". exiting." << endl;
      return wavePads;
    }
  // assume that all mass trees have same wave set
  fitResult* massBin = new fitResult();
  trees[0]->SetBranchAddress(branchName.c_str(), &massBin);
  const int nmbMassBins = trees[0]->GetEntries();
  trees[0]->GetEntry(0);
  const int nmbWaves = massBin->nmbWaves();  // assumes that number of waves is the same for all bins
  printInfo << "drawing wave intensities from tree '" << trees[0]->GetName() << "' for "
	    << nmbWaves << " waves in " << nmbMassBins << " mass bins." << endl;

  printInfo << "calculating total wave intensities" << endl;
  int            nmbBinsAboveThr = 0;
  double         totIntensity    = 0;          // total intensity of all waves
  vector<double> totIntensities(nmbWaves, 0);  // total intensity of individual waves
  for (int i = 0; i < nmbMassBins; ++i) {
    trees[0]->GetEntry(i);
    const double binIntensity = massBin->intensity();
    totIntensity += binIntensity;
    if (binIntensity > intensityThr) {
      for (int j = 0; j < nmbWaves; ++j)
	totIntensities[j] += massBin->intensity(j);
      ++nmbBinsAboveThr;
    }
  }

  printInfo << "sorting waves according to their intensity" << endl;
  vector<pair<string, double> > waveIntensities(nmbWaves, pair<string, double>("", 0));
  for (int i = 0; i < nmbWaves; ++i) {
    const double relIntensity = totIntensities[i] / totIntensity;
    const string waveName     = massBin->waveName(i).Data();
    waveIntensities[i] = make_pair(waveName, relIntensity);
  }
  sort(waveIntensities.begin(), waveIntensities.end(), compareIntensities);

  printInfo << "grouping waves w.r.t. to their JPC ..." << endl;
  map<string, int> nmbWavesPerJpc;  // number of waves for each JPC
  for (int i = 0; i < nmbWaves; ++i) {
    const string waveName = massBin->waveName(i).Data();
    const string jpc      = (waveName != "flat") ? string(waveName, 2, 3) : waveName;
    cout << "    wave [" << i << "] of " << nmbWaves << ": " << waveName << ", JPC = " << jpc << endl;
    ++nmbWavesPerJpc[jpc];
  }
  cout << "    ... finished" << endl;

  printInfo << "creating canvases for all JPCs ..." << endl;
  const int nmbPadsHor     = (int)ceil(sqrt(nmbPadsPerCanvMin));
  const int nmbPadsVert    = (int)ceil((double)nmbPadsPerCanvMin / (double)nmbPadsHor);
  const int nmbPadsPerCanv = nmbPadsHor * nmbPadsVert;
  map<string, TCanvas*> canvases;
  for (map<string, int>::const_iterator i = nmbWavesPerJpc.begin(); i != nmbWavesPerJpc.end(); ++i) {
    const string jpc     = i->first;
    const int    nmbCanv = (int)ceil((double)i->second / (double)nmbPadsPerCanv);
    for (int j = 0; j < nmbCanv; ++j) {
      const string name = buildCanvName(jpc, j);
      if (canvases[name] == NULL) {
	canvases[name] = new TCanvas(name.c_str(), name.c_str(), 10, 10, 1000, 800);
	canvases[name]->Divide(nmbPadsHor, nmbPadsVert);
      }
    }
  }

  printInfo << "creating output ROOT file ..." << endl;
  const string outFileName = outPath + "waveIntensities.root";
  TFile*       outFile     = TFile::Open(outFileName.c_str(), "RECREATE");
  if (!outFile) {
    printErr << "could not create output file '" << outFileName << "'. exiting." << endl;
    return wavePads;
  }
  gROOT->cd();

  printInfo << "plotting waves ordered by decreasing intensity ..." << endl;
  wavePads.resize(nmbWaves, pair<string, TVirtualPad*>("", NULL));
  map<string, int> canvJpcCounter; 
  for (int i = 0; i < nmbWaves; ++i) {
    // get canvas
    const string waveName  = waveIntensities[i].first;
    const string jpc       = (waveName != "flat") ? string(waveName, 2, 3) : waveName;
    const int    canvCount = canvJpcCounter[jpc];
    const string canvName  = buildCanvName(jpc, (int)floor((double)canvCount / (double)nmbPadsPerCanv));
    const int    padIndex  = (canvCount % nmbPadsPerCanv) + 1;
    TCanvas*     canv      = canvases[canvName];
    canv->cd(padIndex);
    cout << endl
	 << "    Selecting canvas '" << canvName << "'." << endl
	 << "    wave " << canvCount + 1 << " of " << nmbWavesPerJpc[jpc]
	 << " for JPC = " << jpc << ", intensity rank = " << i + 1 << endl;
    ++canvJpcCounter[jpc];

    // draw intensity graph
    TMultiGraph* graph = plotIntensity(nmbTrees, trees, waveName, "", "", "AP", 1,
				       graphColors, yAxisRangeMax, false, branchName);
    if (!graph)
      continue;

    // write graph to file
    outFile->cd();
    graph->Write();
    gROOT->cd();

    // draw additional info
    const double intensity = floor(waveIntensities[i].second * 1000 + 0.5) / 10;
    stringstream label;
    label << "I = " << intensity << "% (" << i + 1 << ")";
    TLatex* text = new TLatex(0.15, 0.85, label.str().c_str());
    text->SetNDC(true);
    text->Draw();
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
      legend->Draw();
    }

    // memorize pad
    wavePads[i].first  = waveName;
    wavePads[i].second = gPad;
  }

  // also write TH3s created by plotIntensity() to file
  TList* histList = gDirectory->GetList();
  for (unsigned int i = 0; i < nmbTrees; ++i)
    histList->Remove(trees[i]);
  printInfo << "writing the following " << histList->GetEntries() <<" objects to '" << outFile->GetName() << "'" << endl;
  histList->Print();
  outFile->cd();
  histList->Write();

  // cleanup
  outFile->Close();
  if (outFile)
    delete outFile;
  outFile = 0;
  gROOT->cd();

  if (createPsFile) {
    const string psFileName = outPath + "waveIntensities.ps";
    TCanvas      dummyCanv("dummy", "dummy");
    dummyCanv.Print((psFileName + "[").c_str());
    for (map<string, TCanvas*>::iterator i = canvases.begin(); i != canvases.end(); ++i) {
      i->second->Print(psFileName.c_str());
      delete i->second;
      i->second = 0;
    }
    dummyCanv.Print((psFileName + "]").c_str());
    gSystem->Exec(("gv " + psFileName).c_str());
  }

  return wavePads;
}
