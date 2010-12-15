/*
 * TrpwaWaveSelectFrame.cc
 *
 *  Created on: Aug 26, 2010
 *      Author: Promme
 */

#include "TrpwaPlotAmpsFrame.h"
#include <sstream>
#include <iostream>
#include "TGButtonGroup.h"
#include "TGButton.h"
#include "TGTab.h"
#include "TGLabel.h"
#include "TrpwaCommonTools.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TChain.h"
//#include "TFitBin.h"
#include "plotIntensity.h"
#include "plotPhase.h"
#include "plotCoherence.h"
#include "plotAllIntensities.h"

using namespace std;
using namespace TrpwaCommonTools;

#ifndef TRPWAPLOTAMPSFRAME_CC_
#define TRPWAPLOTAMPSFRAME_CC_

TrpwaPlotAmpsFrame::TrpwaPlotAmpsFrame(
		vector<string> fit_result_paths, // paths containing fit results (all file with the ending .root will be used
		vector<string> fit_result_titles,           // titles of the fits (must be unique)
		vector<string> fit_result_descriptions      // descriptions of the fits
	): TGTransientFrame(gClient->GetRoot(),0,600,600, kVerticalFrame){
	available_fit_results.clear();
	if (fit_result_descriptions.size() == fit_result_paths.size() && fit_result_paths.size() == fit_result_titles.size()){
		for (unsigned int i = 0; i < fit_result_paths.size(); i++){
			//if (i < 5) continue;
			//cout << fit_result_paths[i] << endl;
			if (!Add_Fit_Result(fit_result_paths[i], fit_result_titles[i], fit_result_descriptions[i])){
				cout << " Warning: problems adding fit results from " << fit_result_paths[i] << endl;
			}
		}
		Build();
		fClient->WaitFor(this); // wait till the user closes this window
	} else {
		cout << " Error in TrpwaPlotAmpsFrame::TrpwaPlotAmpsFrame(): wrong input paramters to constructor! " << endl;
		this->CloseWindow();
	}
}

TrpwaPlotAmpsFrame::~TrpwaPlotAmpsFrame(){
	for (Tfilemapit it = available_fit_results.begin(); it != available_fit_results.end(); it++){
		//it->second->Close();
		delete it->second;
	}
	Cleanup();
}

void TrpwaPlotAmpsFrame::Build(){
	// set the root layout for the objects to draw
    gROOT->SetStyle("Plain");
    gStyle->SetOptTitle(kFALSE);

    gStyle->SetOptStat(0);
    //gStyle->SetOptStat(111); // to see histogram statistics in plots
    //gStyle->SetOptFit(111);  // to see fit output in plots
    /*
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetFillColor(0);
    gStyle->SetTitleFillColor(0);*/

    gROOT->ForceStyle();

    // group for root file selection
	TGGroupFrame* frame_rootfile_selections = new TGGroupFrame(this, " select fits ", kHorizontalFrame);

	TGVerticalFrame* subframe_available_fits = new TGVerticalFrame(frame_rootfile_selections);
	subframe_available_fits->AddFrame(
				new TGLabel(subframe_available_fits, "available fit results"));
	// drop down box with available fit results
	box_available_fits = new TGComboBox(subframe_available_fits);
	box_available_fits->Connect("Selected(Int_t)", "TrpwaPlotAmpsFrame", this, "Add_Fit(Int_t)");
	box_available_fits->AddEntry("NONE",-1);
	for (Tfilemapit it = available_fit_results.begin(); it != available_fit_results.end(); it++){
		box_available_fits->AddEntry(it->first.c_str(), (long int) it->second);
	}

	box_available_fits->SetHeight(20);
	subframe_available_fits->AddFrame(box_available_fits, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX,1,1,1,1));

	TGVerticalFrame* subframe_selected_fits = new TGVerticalFrame(frame_rootfile_selections);
	subframe_selected_fits->AddFrame(
				new TGLabel(subframe_selected_fits, "selected fit results"));
	box_selected_fits = new TGComboBox(subframe_selected_fits);
	box_selected_fits->Connect("Selected(Int_t)", "TrpwaPlotAmpsFrame", this, "Remove_Fit(Int_t)");
	box_selected_fits->AddEntry("NONE",-1);
	box_selected_fits->SetHeight(20);
	subframe_selected_fits->AddFrame(box_selected_fits, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX,1,1,1,1));

	frame_rootfile_selections->AddFrame(subframe_available_fits, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX,1,1,1,1));
	frame_rootfile_selections->AddFrame(subframe_selected_fits, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX,1,1,1,1));

	// button to draw all intensities (external script will be called)
	TGTextButton* plot_all_button = new TGTextButton(this, new TGHotString(" draw all waves "));
	plot_all_button->Connect("Clicked()","TrpwaPlotAmpsFrame",this,"Plot_All_selected()");

	// group with buttons to select specific (anchor) waves
	TGGroupFrame* frame_partial_wave_selections = new TGGroupFrame(this, " select partial waves ", kHorizontalFrame);
	TGVerticalFrame* subframe_partial_wave_selections = new TGVerticalFrame(frame_partial_wave_selections);
	subframe_partial_wave_selections->AddFrame(
					new TGLabel(subframe_partial_wave_selections, "available partial waves"));
	// drop down box with available waves given by the selected root files
	box_available_waves = new TGComboBox(subframe_partial_wave_selections);
	box_available_waves->Connect("Selected(Int_t)", "TrpwaPlotAmpsFrame", this, "Select_Wave(Int_t)");
	box_available_waves->AddEntry("NONE",-1);
	box_available_waves->SetHeight(20);
	subframe_partial_wave_selections->AddFrame(box_available_waves, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX,1,1,1,1));

	TGVerticalFrame* subframe_anchor_wave_selections = new TGVerticalFrame(frame_partial_wave_selections);
		subframe_anchor_wave_selections->AddFrame(
						new TGLabel(subframe_anchor_wave_selections, "selected anchor wave"));
	box_available_anchor_waves = new TGComboBox(subframe_anchor_wave_selections);
	box_available_anchor_waves->Connect("Selected(Int_t)", "TrpwaPlotAmpsFrame", this, "Select_Anchor_Wave(Int_t)");
	box_available_anchor_waves->AddEntry("NONE",-1);
	box_available_anchor_waves->SetHeight(20);
	subframe_anchor_wave_selections->AddFrame(box_available_anchor_waves, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX,1,1,1,1));

	frame_partial_wave_selections->AddFrame(subframe_partial_wave_selections, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX,1,1,1,1));
	frame_partial_wave_selections->AddFrame(subframe_anchor_wave_selections, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX,1,1,1,1));

	// button to plot the selection
	TGTextButton* plot_button = new TGTextButton(this, new TGHotString(" draw selected waves "));

	// canvas where the results are drawn to
	TRootEmbeddedCanvas* frame_selected_waves = new TRootEmbeddedCanvas(" canvas_selected_waves ",
			this, 600, 600, kSunkenFrame|kDoubleBorder , 0xffffff);
	canvas_selected_waves = frame_selected_waves->GetCanvas();
	canvas_selected_waves->Divide(2,2);


	this->AddFrame(frame_rootfile_selections, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX,1,1,1,1));
	this->AddFrame(plot_all_button, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX,1,1,1,1));
	this->AddFrame(frame_partial_wave_selections, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX,1,1,1,1));
	this->AddFrame(plot_button, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX,1,1,1,1));
	this->AddFrame(frame_selected_waves);

	SetWindowName("Inspect Results");
	Resize(GetDefaultSize());
	SetWMSize(fWidth, fHeight);
	CenterOnParent();
	MapSubwindows();
	MapWindow();
}

bool TrpwaPlotAmpsFrame::Add_Fit_Result(string fit_result_path, // path containing .root files with fit result trees
			string fit_result_title,			// (unique) title of the fit
			string fit_result_description){		// description of the fit
	bool result(true);
	if (!DirExists(fit_result_path)) return false;
	if (available_fit_results.find(fit_result_title) == available_fit_results.end()){
		// search for .root files in this path
		vector<string> rootfiles;
		if (!(GetDir(fit_result_path, rootfiles, ".root", false) > 0)) result = false;
		TTree* fitresulttree = new TChain(fit_result_title.c_str());
		for (vector<string>::iterator it = rootfiles.begin(); it != rootfiles.end(); it++){
			if (!((TChain*)fitresulttree)->AddFile((fit_result_path+"/"+(*it)).c_str(), -1, "pwa")){
				cout << " Error in TrpwaPlotAmpsFrame::Add_Fit_Result(): Could not add file " << fit_result_path+"/"+(*it) << endl;
			} else {
				//cout << " added " << fit_result_path+"/"+(*it) << endl;
			}
		}
		//cout << " found " << fitresulttree->GetEntries() << " entries " << endl;
		available_fit_results[fit_result_title] = fitresulttree;
		//cout << " tree name is " << fitresulttree->GetName() << endl;
	} else {
		cout << " Error in TrpwaPlotAmpsFrame::Add_Fit_Result(): Fit result with the title " << fit_result_title << " does already exist! " << endl;
		return false;
	}
	return result;
}

// will be called when an available fit file was selected
// an item will be added to the list of selected fit files
void TrpwaPlotAmpsFrame::Add_Fit(int pFitFile){
	if (pFitFile > 0){
		// get the name and add the entry to the list of selected fits
		TTree* selected_fit = (TTree*) pFitFile;
		string fitname = selected_fit->GetName();
		if (selected_fit_results.find(fitname)==selected_fit_results.end()){
			selected_fit_results[fitname]=selected_fit;
			box_selected_fits->AddEntry(fitname.c_str(), pFitFile);
		}
	//selected_fit_results[];
	//cout << pFitFile << endl;
	//cout << " calling add fit " << endl;
	}
}

// will be called when an added fit file was selected
// if the user confirms the selected item will be removed from the list
void TrpwaPlotAmpsFrame::Remove_Fit(int pFitFile){
	cout << pFitFile << endl;
	cout << " calling remove fit " << endl;
}

// a wave will be selected from the list of available waves
void TrpwaPlotAmpsFrame::Select_Wave(int pWave){
	cout << pWave << endl;
	cout << " calling select wave " << endl;
}

// a wave will be selected for the anchor wave (phase determination)
void TrpwaPlotAmpsFrame::Select_Anchor_Wave(int pWave){
	cout << pWave << endl;
	cout << " calling select anchor wave " << endl;
}

vector<string>& TrpwaPlotAmpsFrame::Scan_Fit_Result(TTree* fit_results){
	vector<string>* result = new vector<string>();

	return *result;
}

void TrpwaPlotAmpsFrame::Plot_All_selected(){
	for (Tfilemapit it = selected_fit_results.begin(); it != selected_fit_results.end(); it++){
		TTree* selected_tree = (TTree*)it->second;
		plotAllIntensities(selected_tree, true);
		//box_available_fits->AddEntry(it->first.c_str(), (long int) it->second);
	}

}


#endif /* TRPWAPLOTAMPSFRAME_CC_ */
