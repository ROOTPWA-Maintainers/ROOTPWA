//
// plots result from Dima's fit on top of ROOTPWA histograms
//


#include <string>
#include <iostream>
#include <map>
#include <algorithm>
#include <vector>
#include <utility>

#include "TFile.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TRegexp.h"
#include "TPRegexp.h"
#include "TH1.h"
#include "TLegend.h"
#include "TKey.h"
#include "TObject.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TPostScript.h"
#include "TSystem.h"

#include "reportingUtils.hpp"
#include "fitResult.h"
#include "plotAllIntensities.h"
#include "plotIntensity.h"
#include "plotSpinTotals.h"
#include "plotPhase.h"


using namespace std;
using namespace rpwa;


// reads Dima's histograms from file into map (histogram title, histogram pointer)
map<string, TH1*>
readDimaHists(TFile*   dimaFile,
              TPRegexp histNamePattern)  // histogram name pattern
{
	map<string, TH1*> dimaHists;
	TIterator*        keys = dimaFile->GetListOfKeys()->MakeIterator();
	while (TKey* k = static_cast<TKey*>(keys->Next())) {
		if (!k) {
			printWarn << "NULL pointer to TKey. skipping." << endl;
			continue;
		}
		TObject* o = k->ReadObj();
		if (!o) {
			printWarn << "cannot read object from TKey '" << k->GetName() << "'. skipping." << endl;
			continue;
		}
		const string type = o->IsA()->GetName();
		if (type.find("TH1F") != string::npos) {
			TH1* h              = static_cast<TH1*>(o);
			const TString hName = h->GetName();
			if (hName.Contains(histNamePattern))
				dimaHists[h->GetTitle()] = h;
		}
	}
	return dimaHists;
}


// constructs wave name in Dima's convention from wave name in ROOTPWA convention
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
		const string       dictionary[nmbDict][2] = {//{"sigma",    "eps1"},
			{"sigma",    "f0(1400)"},
			//{"f0980",    "f0"},
			{"f0980",    "f0(980)"},
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
compareIntensityWithDima(TTree* tree, TFile* dimaFile, std::string waveName) {
	printInfo << "comparing intensities with Dima's result for wave " << waveName << endl;

	// get Dima's intensity histograms from file
	map<string, TH1*> dimaIntensityHists;  // intensity histograms w/o acceptance correction (see histgroups.txt)
	{
		map<string, TH1*> dimaHists = readDimaHists(dimaFile, TString("^h10\\d{2}$"));
		// transform map with histogram titles as keys to map with wave names as keys
		for (map<string, TH1*>::iterator i = dimaHists.begin(); i != dimaHists.end(); ++i) {
			cout << "    found intensity histogram '"<< i->first << "' " << flush;
			TString waveName = i->first;
			waveName         = waveName.Strip(TString::kTrailing, ' ');
			dimaIntensityHists[waveName.Data()] = i->second;
			cout << "for wave '" << waveName << "'" << endl;
		}
	}

	TCanvas *pad = new TCanvas(waveName.c_str());
	pad->cd();
	TMultiGraph* rootgraph = plotIntensity(tree, waveName);
	rootgraph->Draw("A*");
	const string dimaName = translateWaveName(waveName);
	cout << "    drawing intensity for wave '" << waveName << "' -> '" << dimaName << "'" << endl;

	// draw Dima's plot
	pad->cd();
	TH1* hDima = dimaIntensityHists[dimaName];
	if (!hDima) {
		printWarn << "cannot find histogram with name '" << dimaName << "'. skipping." << endl;
		return;
	}
	hDima->SetMarkerStyle(21);
	hDima->SetMarkerSize(0.5);
	hDima->SetLineColor(2);
	hDima->SetMarkerColor(2);
	hDima->Draw("SAME");

	TLegend *leg = new TLegend(0.82,0.97,0.97,0.82);
	//leg->SetHeader("The Legend Title");
	leg->AddEntry(rootgraph, "ROOTPWA","flpa");
	leg->AddEntry(hDima,"COMPASSPWA","flpa");
	leg->Draw();

	pad->Update();

	// readjust y-axis range
	TMultiGraph* graph = static_cast<TMultiGraph*>(pad->GetPrimitive(waveName.c_str()));
	if (graph) {
		const double max = std::max(hDima->GetMaximum(), graph->GetHistogram()->GetMaximum());
		graph->SetMaximum(1.1 * max);
		pad->Modified();
		pad->Update();
	}
}

void
compareIntensitiesWithDima(TTree* tree,  // fitResult tree
                           TFile* dimaFile, const bool createPsFile)
{
	printInfo << "comparing intensities with Dima's result." << endl;

	// get Dima's intensity histograms from file
	map<string, TH1*> dimaIntensityHists;  // intensity histograms w/o acceptance correction (see histgroups.txt)
	{
		map<string, TH1*> dimaHists = readDimaHists(dimaFile, TString("^h10\\d{2}$"));
		// transform map with histogram titles as keys to map with wave names as keys
		for (map<string, TH1*>::iterator i = dimaHists.begin(); i != dimaHists.end(); ++i) {
			cout << "    found intensity histogram '"<< i->first << "' " << flush;
			TString waveName = i->first;
			waveName         = waveName.Strip(TString::kTrailing, ' ');
			dimaIntensityHists[waveName.Data()] = i->second;
			cout << "for wave '" << waveName << "'" << endl;
		}
	}

	// get ROOTPWA intensity histograms
	vector<pair<string, TVirtualPad*> > wavePads = plotAllIntensities(tree, false, "./");

	// draw Dima's histograms on top of ROOTPWA's
	cout << endl;
	printInfo << "drawing Dima's histograms" << endl;
	for (unsigned int i = 0; i < wavePads.size(); ++i) {
		const string waveName = wavePads[i].first;
		const string dimaName = translateWaveName(waveName);
		cout << "    drawing intensity for wave '" << waveName << "' -> '" << dimaName << "'" << endl;
    
		// draw Dima's plot
		TLegend *leg = new TLegend(0.82,0.97,0.97,0.82);
		TVirtualPad* pad = wavePads[i].second;
		pad->cd();
		TMultiGraph* list = (TMultiGraph*)pad->GetListOfPrimitives();
		leg->AddEntry(&list[0], "ROOTPWA","flpa");
		TH1* hDima = dimaIntensityHists[dimaName];
		if (!hDima) {
			printWarn << "cannot find histogram with name '" << dimaName << "'. skipping." << endl;
			continue;
		}
		//hDima->Scale(3);
		hDima->SetMarkerStyle(21);
		hDima->SetMarkerSize(0.5);
		hDima->SetLineColor(2);
		hDima->SetMarkerColor(2);
		hDima->Draw("SAME");

		leg->AddEntry(hDima,"COMPASSPWA","flpa");
		leg->Draw();

		pad->Update();
    
		// readjust y-axis range
		TMultiGraph* graph = static_cast<TMultiGraph*>(pad->GetPrimitive(waveName.c_str()));
		if (graph) {
			const double max = std::max(hDima->GetMaximum(), graph->GetHistogram()->GetMaximum());
			graph->SetMaximum(1.1 * max);
			pad->Modified();
			pad->Update();
		}
	}
	if (createPsFile) {
		const string psFileName = "waveIntensities_compareDima.ps";
		TCanvas      dummyCanv("dummy", "dummy");
		dummyCanv.Print((psFileName + "[").c_str());
		vector<TCanvas*> canvases;
		for (unsigned int i = 0; i < wavePads.size(); ++i) {
			for (unsigned int j = 0; j < canvases.size(); j++) {
				if(canvases[j] == wavePads[i].second->GetCanvas()) {
					goto next;
				}
			}
			canvases.push_back(wavePads[i].second->GetCanvas());
			wavePads[i].second->GetCanvas()->Print(psFileName.c_str());
		next: ;
		}
		dummyCanv.Print((psFileName + "]").c_str());
		gSystem->Exec(("gv " + psFileName).c_str());
	}
}


void
compareSpinTotalsWithDima(TTree* tree,  // fitResult tree
                          TFile* dimaFile)
{
	printInfo << "comparing spin totals with Dima's result." << endl;

	// get Dima's spin total histograms from file
	map<string, TH1*> dimaSpinTotalHists;
	{
		map<string, TH1*> dimaHists = readDimaHists(dimaFile, TString("^h4\\d{3}$"));
		// transform map with histogram titles as keys to map with wave names as keys
		// this is still very simplistic and works only for part of the histograms
		for (map<string, TH1*>::iterator i = dimaHists.begin(); i != dimaHists.end(); ++i) {
			TString waveName  = i->first;
			bool    waveMatch = true;
			if (waveName == " Events                                                     ")
				waveName = "";
			else if ((waveName(1, 6) == "J^PC!=") && (waveName(18, 5) == "total")) {
				waveName = waveName(8, 4);
				waveName.ReplaceAll("^", "");
			} else if ((waveName(1, 10) == "J^PC!M[c]=") && (waveName(22, 5) == "total")) {
				waveName = waveName(12, 7);
				waveName.ReplaceAll("^", "");
				waveName.ReplaceAll("!", "");
			} else
				waveMatch = false;
			if (waveMatch) {
				cout << "    found spin total histogram '"<< i->first << "' " << flush;
				dimaSpinTotalHists[waveName.Data()] = i->second;
				cout << "for wave '" << waveName << "'" << endl;
			}
		}
	}

	// get ROOTPWA spinTotal histograms
	vector<pair<string, TVirtualPad*> > wavePads = plotSpinTotals(tree, kBlack, 0, true, "");

	// draw Dima's histograms on top of ROOTPWA's
	cout << endl;
	printInfo<< "drawing Dima's histograms" << endl;
	for (unsigned int i = 0; i < wavePads.size(); ++i) {
		const string waveName = wavePads[i].first;
		cout << "    drawing spinTotal for wave '" << waveName << "'" << endl;
    
		// draw Dima's plot
		TVirtualPad* pad = wavePads[i].second;
		pad->cd();
		TH1* hDima = dimaSpinTotalHists[waveName];
		if (!hDima) {
			printWarn << "cannot find histogram with name '" << waveName << "'. skipping." << endl;
			continue;
		}
		hDima->Scale(3);
		hDima->SetMarkerStyle(21);
		hDima->SetMarkerSize(0.5);
		hDima->SetLineColor(2);
		hDima->SetMarkerColor(2);
		hDima->Draw("SAME");
		pad->Update();
    
		// readjust y-axis range
		TString graphName = "total";
		if (waveName != "") {
			graphName = "g";
			graphName.Append(waveName);
			graphName.ReplaceAll("+", "p");
			graphName.ReplaceAll("-", "m");
		}
		TGraphErrors* graph = static_cast<TGraphErrors*>(pad->GetPrimitive(graphName));
		if (graph) {
			const double max = std::max(hDima->GetMaximum(), graph->GetHistogram()->GetMaximum());
			graph->SetMaximum(1.1 * max);
			pad->Modified();
			pad->Update();
		}
	}
}


void
comparePhasesWithDima(TTree*        tree,  // fitResult tree
                      TFile*        dimaFile,
                      const vector<pair<string, string> >& wavenames_phases,
                      const string& branchName = "fitResult_v2")
{
	printInfo << "comparing phases with Dima's result." << endl;

	// get Dima's phase histograms from file
	map<string, TH1*> dimaPhaseHists;  // relative phases (see histgroups.txt)
	{
		map<string, TH1*> dimaHists = readDimaHists(dimaFile, TString("^h4\\d{4}$"));
		// transform map with histogram titles as keys to map with wave names as keys
		for (map<string, TH1*>::iterator i = dimaHists.begin(); i != dimaHists.end(); ++i) {
			cout << "    found phase histogram '"<< i->first << "' " << flush;
			TString   hTitle       = i->first;
			hTitle                 = hTitle(4, hTitle.Length() - 5);
			const int pos          = hTitle.Index(" - ");
			TString   waveNames[2] = {hTitle(0, pos - 1), hTitle(pos + 3, hTitle.Length() - (pos + 3))};
			waveNames[0]           = waveNames[0].Strip(TString::kTrailing, ' ');
			waveNames[1]           = waveNames[1].Strip(TString::kTrailing, ' ');
			dimaPhaseHists[(waveNames[0] + "|" + waveNames[1]).Data()] = i->second;
			cout << "for waves '" << waveNames[0] << "', '" << waveNames[1] << "'" << endl;
		}
	}

	// get ROOTPWA phase histograms
	fitResult* massBin = new fitResult();
	tree->SetBranchAddress(branchName.c_str(), &massBin);
	tree->GetEntry(0);

	for (unsigned int i = 0; i < wavenames_phases.size(); ++i) {
		// get Dima's plot
		const string histName = translateWaveName(wavenames_phases[i].first) + "|" + translateWaveName(wavenames_phases[i].second);
		TH1* h = dimaPhaseHists[histName];
		if (!h) {
			printWarn << "cannot find histogram with name '" << histName << "'. skipping." << endl;
			continue;
		}

		// draw ROOTPWA phase histograms
		new TCanvas();
		TMultiGraph* rootgraph =  plotPhase(tree, wavenames_phases[i].first, wavenames_phases[i].second);

		// draw Dima's plot on top of it
		h->SetMarkerStyle(21);
		h->SetMarkerSize(0.5);
		h->SetLineColor(2);
		h->SetMarkerColor(2);
		h->Draw("SAME");

		TLegend *leg = new TLegend(0.82,0.97,0.97,0.82);
		//leg->SetHeader("The Legend Title");
		leg->AddEntry(rootgraph, "ROOTPWA","flpa");
		leg->AddEntry(h,"COMPASSPWA","flpa");
		leg->Draw();

	}
}


void
compareWithDima(TTree*        tree,  // fitResult tree
                const string& dimaFileName = "/afs/e18.ph.tum.de/data/compass/bgrube/lambda/pwa/compassPWA/work/hfit.root",
                const string& wavename = "",
                const bool createPsFile =false)
{
	if (!tree) {
		printErr << "NULL pointer to tree. exiting." << endl;
		return;
	}

	printInfo << "opening Dima's file '" << dimaFileName << "'" << endl;
	TFile* dimaFile = TFile::Open(dimaFileName.c_str(), "READ");
	if (!dimaFile || dimaFile->IsZombie()) {
		printErr << "cannot open file '" << dimaFileName << "'. exiting." << endl;
		return;
	}
  
	compareIntensityWithDima(tree, dimaFile, wavename);
	//compareIntensitiesWithDima(tree, dimaFile, createPsFile);
	//compareSpinTotalsWithDima(tree, dimaFile);
	std::vector<std::pair<std::string, std::string> > phases;
	phases.push_back(std::make_pair("1-1++0+rho770_01_pi0.amp", "1-2-+0+f21270_02_pi-.amp"));
	phases.push_back(make_pair("1-2++1+rho770_21_pi0.amp", "1-1++0+rho770_01_pi0.amp"));
	phases.push_back(make_pair("1-2++1+rho770_21_pi0.amp", "1-2-+0+f21270_02_pi-.amp"));
	phases.push_back(make_pair("1-1-+1+rho770_11_pi0.amp", "1-2++1+rho770_21_pi0.amp"));
	phases.push_back(make_pair("1-1-+1+rho770_11_pi0.amp", "1-1++0+rho770_01_pi0.amp"));

	comparePhasesWithDima(tree, dimaFile, phases);

}
