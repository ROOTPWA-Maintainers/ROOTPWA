//
// plots result from Dima's fit on top of ROOTpwa histograms
//


#include <string>
#include <iostream>
#include <map>
#include <algorithm>

#include "TFile.h"
#include "TPRegexp.h"
#include "TH1.h"
#include "TKey.h"
#include "TObject.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"

#include "plotAllIntensities.h"
#include "plotPhase.h"
#include "TFitResult.h"


using namespace std;


// constructs wave name in Dima's convention from wave name in ROOTpwa convention
string translateWaveName(const string& waveName)
{
  TString dimaName = waveName;

  if (dimaName == "flat")
    return "FLAT";

  // add braces and space
  dimaName.Insert(2, "(");
  dimaName.Insert(6, ")");
  dimaName.Insert(9, " ");

  // translate isobars
  {
    const unsigned int nmbDict                = 5;
    const string       dictionary[nmbDict][2] = {{"sigma",    "eps1"},
						 {"f0980",    "f0"},
						 {"rho770",   "rho"},
						 {"f21270",   "f2"},
						 {"rho31690", "rho3"}};
    for (unsigned int i = 0; i < nmbDict; ++i)
      dimaName.ReplaceAll(dictionary[i][0].c_str(), dictionary[i][1].c_str());
  }
  
  // translate L
  {
    map<char, char> dictionary;
    dictionary['0'] = 'S';
    dictionary['1'] = 'P';
    dictionary['2'] = 'D';
    dictionary['3'] = 'F';
    dictionary['4'] = 'G';
    dictionary['5'] = 'H';
    dictionary['6'] = 'I';
    dictionary['7'] = 'K';
    dictionary['8'] = 'L';

    int  pos = dimaName.First("_");
    char L   = dimaName(pos + 1);
    
    // cut off last part
    dimaName  = dimaName(0, pos);
    dimaName += " pi ";
    dimaName += dictionary[L];
  }

  return dimaName.Data();
}


void
compareIntensitiesWithDima(TTree* tree,  // TFitResult tree
			   TFile* dimaFile)
{
  // get Dima's intensity histograms from file
  TPRegexp          intensityHistNamePattern("^h10\\d{2}$");  // intensities w/o acceptance correction (see histgroups.txt)
  map<string, TH1*> dimaIntensityHists;
  TIterator*        keys = dimaFile->GetListOfKeys()->MakeIterator();
  while (TKey* k = static_cast<TKey*>(keys->Next())) {
    if (!k) {
      cerr << "    error: null pointer to TKey." << endl;
      continue;
    }
    TObject* o = k->ReadObj();
    if (!o) {
      cerr << "    error: cannot read object from TKey '" << k->GetName() << "'." << endl;
      continue;
    }
    const string type = o->IsA()->GetName();
    if (type.find("TH1F") != string::npos) {
      TH1* h              = static_cast<TH1*>(o);
      const TString hName = h->GetName();
      if (hName.Contains(intensityHistNamePattern)) {
	TString hTitle = h->GetTitle();
	hTitle         = hTitle.Strip(TString::kTrailing, ' ');
	const string waveName = hTitle.Data();
	dimaIntensityHists[waveName] = h;
	cout << "        found intensity histogram '"<< hName << "' for wave '" << waveName << "'" << endl;
      }
    }
  }

  // get ROOTpwa intensity histograms
  vector<pair<string, TVirtualPad*> > wavePads = plotAllIntensities(tree, false, "./", false);

  // draw Dima's histograms on top of ROOTpwa's
  cout << endl << "drawing Dima's histograms" << endl;
  for (unsigned int i = 0; i < wavePads.size(); ++i) {
    const string waveName = wavePads[i].first;
    const string dimaName = translateWaveName(waveName);
    cout << "    drawing intensity for wave '" << waveName << "' -> '" << dimaName << "'" << endl;
    
    // draw Dima's plot
    TVirtualPad* pad = wavePads[i].second;
    pad->cd();
    TH1* hDima = dimaIntensityHists[dimaName];
    if (!hDima) {
      cerr << "    no histogram with name '" << dimaName << "'" << endl;
      continue;
    }
    hDima->Scale(3);
    hDima->SetMarkerStyle(21);
    hDima->SetMarkerSize(0.5);
    hDima->SetLineColor(2);
    hDima->SetMarkerColor(2);
    hDima->Draw("SAME");
    pad->Update();
    
    // scale y-axis
    TGraphErrors* graph = static_cast<TGraphErrors*>(pad->GetPrimitive(waveName.c_str()));
    if (graph) {
      const double max = std::max(hDima->GetMaximum(), graph->GetHistogram()->GetMaximum());
      graph->SetMaximum(1.1 * max);
      pad->Modified();
      pad->Update();
    }
  }
}


void
comparePhasesWithDima(TTree* tree,  // TFitResult tree
		      TFile* dimaFile)
{
  // get Dima's phase histograms from file
  TPRegexp          phaseHistNamePattern("^h4\\d{4}$");  // relative phases (see histgroups.txt)
  map<string, TH1*> dimaPhaseHists;
  TIterator*        keys = dimaFile->GetListOfKeys()->MakeIterator();
  while (TKey* k = static_cast<TKey*>(keys->Next())) {
    if (!k) {
      cerr << "    error: null pointer to TKey." << endl;
      continue;
    }
    TObject* o = k->ReadObj();
    if (!o) {
      cerr << "    error: cannot read object from TKey '" << k->GetName() << "'." << endl;
      continue;
    }
    const string type = o->IsA()->GetName();
    if (type.find("TH1F") != string::npos) {
      TH1* h              = static_cast<TH1*>(o);
      const TString hName = h->GetName();
      if (hName.Contains(phaseHistNamePattern)) {
	TString hTitle           = h->GetTitle();
	hTitle                   = hTitle(4, hTitle.Length() - 5);
	int     pos              = hTitle.Index(" - ");
	TString waveNames[2]     = {hTitle(0, pos - 1), hTitle(pos + 3, hTitle.Length() - (pos + 3))};
	waveNames[0]             = waveNames[0].Strip(TString::kTrailing, ' ');
	waveNames[1]             = waveNames[1].Strip(TString::kTrailing, ' ');
	string histName          = (waveNames[0] + "|" + waveNames[1]).Data();
	dimaPhaseHists[histName] = h;
	cout << "        found phase histogram '"<< hName << "' for waves '" << waveNames[0] << "', '" << waveNames[1] << "'" << endl;
      }
    }
  }

  // get ROOTpwa phase histograms
  TFitResult* massBin = new TFitResult();
  tree->SetBranchAddress("fitResult", &massBin);
  tree->GetEntry(0);
  //for (unsigned int i = 0; i < massBin->nmbWaves(); ++i) {
  for (unsigned int i = 0; i < 5; ++i) {
    for (unsigned int j = 0; j < i; ++j) {
      const string waveNames[2] = {massBin->waveName(i).Data(),     massBin->waveName(j).Data()};
      const string dimaNames[2] = {translateWaveName(waveNames[0]), translateWaveName(waveNames[1])};
      cout << "    drawing phase for waves ('" << waveNames[0] << "', '" << waveNames[1] << "')"
	   << " -> ('" << dimaNames[0] << "', '" << dimaNames[1] << "')" << endl;
      
      // get Dima's plot
      const string histName = dimaNames[0] + "|" + dimaNames[1];
      TH1* h = dimaPhaseHists[histName];
      if (!h) {
	cerr << "    no histogram with name '" << histName << "'" << endl;
	continue;
      }

      // draw ROOTpwa phase histograms
      TCanvas*      canv  = new TCanvas();
      TGraphErrors* graph = plotPhase(tree, i, j);
      
      // draw Dima's plot on top of it
      h->SetMarkerStyle(21);
      h->SetMarkerSize(0.5);
      h->SetLineColor(2);
      h->SetMarkerColor(2);
      h->Draw("SAME");
    }
  }

}


void
compareWithDima(TTree*        tree,  // TFitResult tree
		const string& dimaFileName = "/afs/e18.ph.tum.de/user/bgrube/compass/pwa/compassPWA/work/hfit.root")
{
  if (!tree) {
    cerr << "compareWithDima() error: NULL pointer to tree." << endl;
    return;
  }

  cout << "    opening Dima's file '" << dimaFileName << "'" << endl;
  TFile* dimaFile = TFile::Open(dimaFileName.c_str(), "READ");
  if (!dimaFile || dimaFile->IsZombie()) {
    cerr << "    error opening file '" << dimaFileName << "'" << endl;
    return;
  }
  
  //compareIntensitiesWithDima(tree, dimaFile);
  comparePhasesWithDima(tree, dimaFile);
  
}
