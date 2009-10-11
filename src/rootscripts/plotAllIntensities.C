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
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TList.h"
#include "TPostScript.h"
#include "TSystem.h"

#include "../TFitResult.h"

#include "plotIntensity.h"
//#include "plotmcmc.h"


using namespace std;


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
  unsigned int pos = name.find(" ");
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


void
plotAllIntensities(TTree*        tree,                  // TFitResult tree
		   const bool    createPsFile = false,  // if true, plots are written to waves.ps
		   const string& outPath      = "./",   // path for output files
		   const bool    mcmc         = false)
{
  const double intensityThr      = 500;  // threshold for total intensity in mass bin
  const int    nmbPadsPerCanvMin = 4;    // minimum number of pads each canvas is subdivided into

  if (!tree) {
    cerr << "plotAllIntensities() error: Null pointer to tree. Exiting." << endl;
    return;
  }
  TFitResult* massBin = new TFitResult();
  tree->SetBranchAddress("fitbin", &massBin);
  const int nmbMassBins = tree->GetEntries();
  tree->GetEntry(0);
  const int nmbWaves = massBin->nmbWaves();  // assumes that number of waves is the same for all bins
  cout << "Drawing wave intensities from tree '" << tree->GetName() << "' for "
       << nmbWaves << " waves in " << nmbMassBins << " mass bins." << endl;

  // calculate total intensities
  int            nmbBinsAboveThr = 0;
  double         totIntensity    = 0;          // total intensity of all waves
  vector<double> totIntensities(nmbWaves, 0);  // total intensity of individual waves
  for (int i = 0; i < nmbMassBins; ++i) {
    tree->GetEntry(i);
    const double binIntensity = massBin->intensity();
    totIntensity += binIntensity;
    if (binIntensity > intensityThr) {
      for (int j = 0; j < nmbWaves; ++j)
	totIntensities[j] += massBin->intensity(j);
      ++nmbBinsAboveThr;
    }
  }
 
  // sort waves by intensities (normalized to total intensity) in descending order
  vector<pair<string, double> > waveIntensities(nmbWaves, pair<string, double>("", 0));
  for (int i = 0; i < nmbWaves; ++i) {
    const double relIntensity = totIntensities[i] / totIntensity;
    const string waveName     = massBin->waveName(i).Data();
    waveIntensities[i] = make_pair(waveName, relIntensity);
  }
  sort(waveIntensities.begin(), waveIntensities.end(), compareIntensities);

  // group waves according to their JPC
  map<string, int> nmbWavesPerJpc;  // number of waves for each JPC
  for (int i = 0; i < nmbWaves; ++i) {
    const string waveName = massBin->waveName(i).Data();
    const string jpc      = (waveName != "flat") ? string(waveName, 2, 3) : waveName;
    ++nmbWavesPerJpc[jpc];
  }

  // calculate optimum canvas subdivision
  const int nmbPadsHor     = (int)ceil(sqrt(nmbPadsPerCanvMin));
  const int nmbPadsVert    = (int)ceil((double)nmbPadsPerCanvMin / (double)nmbPadsHor);
  const int nmbPadsPerCanv = nmbPadsHor * nmbPadsVert;

  cout << "    Creating canvases for all JPCs ..." << endl;
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

  cout << "    Creating output ROOT file ..." << endl;
  const string outFileName = outPath + "waveIntensities.root";
  TFile*       outFile     = TFile::Open(outFileName.c_str(), "RECREATE");
  if (!outFile) {
    cerr << "plotAllIntensities() error: Could not create output file '" << outFileName << "'. Exiting." << endl;
    return;
  }
  gROOT->cd();

  cout << "    Plotting waves ordered by decreasing intensity ..." << endl;
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
    //TGraph* g = (mcmc) ? plotmcmc(tree, waveName, "", waveName.c_str()) : plotIntensity(tree, waveName, "", waveName);
    TGraph* g = plotIntensity(tree, waveName, "", waveName);
    if (!g)
      continue;
    g->SetName(waveName.c_str());

    // write graph to file
    outFile->cd();
    g->Write();
    gROOT->cd();

    // draw additional info
    const double intensity = floor(waveIntensities[i].second * 1000 + 0.5) / 10;
    stringstream label;
    label << "I = " << intensity << "% (" << i + 1 << ")";
    TLatex* text = new TLatex(0.15, 0.85, label.str().c_str());
    text->SetNDC(true);
    text->Draw();
  }

  // also write TH3s created by plotIntensity() to file
  TList* histList = gDirectory->GetList();
  histList->Remove(tree);
  cout << endl
       << "Writing the following " << histList->GetEntries() <<" objects to '" << outFile->GetName() << "'" << endl;
  histList->Print();
  outFile->cd();
  histList->Write();

  // cleanup
  outFile->Close();
  if (outFile)
    delete outFile;
  outFile = 0;
  gROOT->cd();

  // plot graphs into PS file
  if (createPsFile) {
    // create a postscript file and set the paper size
    const string psFileName = outPath + "waveIntensities.ps";
    const int    psType     = 112;
    TPostScript  ps(psFileName.c_str(), psType);
    //ps.Range(16,24);  // set x,y of printed page
    for (map<string, TCanvas*>::iterator i = canvases.begin(); i != canvases.end(); ++i) {
      ps.NewPage();
      i->second->Update();
      delete i->second;
      i->second = 0;
    }
    ps.Close();
    //const string command = ;
    gSystem->Exec(("gv " + psFileName).c_str());
  }

}
