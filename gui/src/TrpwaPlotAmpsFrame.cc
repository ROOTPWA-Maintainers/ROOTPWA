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
#include "plotSpinTotals.h"
#include "fitResult.h"
#include "TGTextEntry.h"
#include "TGMsgBox.h"
#include "TAxis.h"
#include "TGFileDialog.h"
#include "TText.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TParameter.h"

using namespace std;
using namespace TrpwaCommonTools;
using namespace rpwa;

#ifndef TRPWAPLOTAMPSFRAME_CC_
#define TRPWAPLOTAMPSFRAME_CC_

TPaveText* label_preliminary;
TPaveText* label_data_info;
TPaveText* label_analysis_info;

void Create_Info(){
	// "preliminary" to plot
	label_preliminary = new TPaveText(0.,0.,1.,1.,"NDC");
	if (1){
		TText* text_preliminary = label_preliminary->AddText("preliminary");
		//TText* text_preliminary = label_preliminary->AddText("ongoing analysis");
		//TText* text_beam_info = label_preliminary->AddText(1,1,"COMPASS 2008 negative hadron beam");
		text_preliminary->SetTextColor(17);
		text_preliminary->SetTextSize(0.10);
		text_preliminary->SetTextAngle(26.15998);
		text_preliminary->SetTextAlign(22);
	}
	label_preliminary->SetFillStyle(4000);
	label_preliminary->SetBorderSize(0);

	label_analysis_info = new TPaveText(0.01,0.95,.99,.99,"NDC");
	label_analysis_info->SetFillColor(0);
	//TText* text_beam_info = label_analysis_info->AddText("COMPASS 2008 negative hadron beam");
	//text_beam_info->SetTextSize(0.05);
	label_analysis_info->SetBorderSize(2);

	//label_data_info = new TPaveText(0.58, 0.78, 0.96, 0.87, "NDC");
	//label_data_info = new TPaveText(0.58, 0.72, 0.96, 0.89, "NDC");
	if (1)
		label_data_info = new TPaveText(0.5, 0.72, 0.89, 0.89, "NDC");
	else
		label_data_info = new TPaveText(0.59, 0.72, 0.89, 0.89, "NDC");
	if (1){
		TText* text_beam_info = label_data_info->AddText("COMPASS 2008 negative hadron beam");
		text_beam_info->SetTextSize(0.03);
		TText* text_channel_info = label_data_info->AddText("K^{-} p #rightarrow K^{-} #pi^{+} #pi^{-} p_{recoil}");
		text_channel_info->SetTextSize(0.04);
		//TText* text_acc_info = label_data_info->AddText("w/o acceptance correction");
		//text_acc_info->SetTextSize(0.03);
		//TText* text_data_info = label_data_info->AddText("W37 = 55% of 2008 data");
		/*
		TText* text_data_info;
		if (0){
			text_data_info = label_data_info->AddText("");//title_MC.c_str());
		} else {
			text_data_info = label_data_info->AddText("");//title_data.c_str());
		}
		text_data_info->SetTextSize(0.03);
		text_data_info->SetTextColor(2);*/
		//text_beam_info->SetNDC();
		//text_beam_info->SetTextAlign(11);
		//label_data_info->SetFillStyle(4000);
		label_data_info->SetFillColor(0);
	} else {
		label_data_info->SetFillStyle(4000);
	}
	label_data_info->SetBorderSize(0);
}

TrpwaPlotAmpsFrame::TrpwaPlotAmpsFrame(
		vector<string> fit_result_paths, // paths containing fit results (all file with the ending .root will be used
		vector<string> fit_result_titles,           // titles of the fits (must be unique)
		vector<string> fit_result_descriptions      // descriptions of the fits
	): TGTransientFrame(gClient->GetRoot(),0,600,600, kVerticalFrame){
	available_fit_results.clear();
	selected_fit_results.clear();
	available_waves.clear();
	available_colors.clear();
	Create_available_colors();
	current_wave = "";
	current_anchor_wave = "";
	current_fit_result = NULL;
	masscutlow  = -1;
	masscuthigh = -1;
	Create_Info();
	for (int ipad = 0; ipad < 5; ipad++){
		plotted_graphs[ipad] = NULL;
		plotted_most_likely_graphs[ipad] = NULL;
	}
	draw_most_likely = true;
	draw_datainfo    = false;
	if (fit_result_descriptions.size() == fit_result_paths.size() && fit_result_paths.size() == fit_result_titles.size()){
		cout << " loading fit results " << endl;
		for (unsigned int i = 0; i < fit_result_paths.size(); i++){
			DrawProgressBar(50,(double)(i+1)/((double)fit_result_paths.size()));
			//if (i < 5) continue;
			//cout << fit_result_paths[i] << endl;
			if (!Add_Fit_Result(fit_result_paths[i], fit_result_titles[i], fit_result_descriptions[i])){
				cout << " Warning: problems adding fit results from " << fit_result_paths[i] << endl;
			}
		}
		cout << " done " << endl;
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
		// store the created graphs
		delete it->second;
	}
	for (Tfitresultmapit it = selected_fit_result_graphs.begin(); it != selected_fit_result_graphs.end(); it++) {
		
		Tdirmapit fit_result_path_it = available_fit_result_paths.find(it->first);
		if (fit_result_path_it != available_fit_result_paths.end()){
			if (it->first == "Current"){
				cout << " graphs of temporary fit results will be not stored to file " << endl;
				continue;
			}
			if (!it->second->Write_to_file(fit_result_path_it->second+"/fit_result_plots.root")){
				cout << " Warning in TrpwaPlotAmpsFrame::~TrpwaPlotAmpsFrame(): could not save graphs to file " << endl;
			}
		}
		delete it->second;
	}
	Cleanup();
}

void TrpwaPlotAmpsFrame::Build(){
	// set the root layout for the objects to draw
    gROOT->SetStyle("Plain");
    // this style is set to fit my phd thesis style:
    if (1){
		gStyle->SetTitleFont(10*13+2,"xyz");
		gStyle->SetTitleSize(0.06, "xyz");
		gStyle->SetTitleOffset(1.3,"y");
		gStyle->SetLabelFont(10*13+2,"xyz");
		gStyle->SetLabelSize(0.06,"xyz");
		gStyle->SetLabelOffset(0.009,"xyz");
		gStyle->SetPadBottomMargin(0.16);
		gStyle->SetPadTopMargin(0.16);
		gStyle->SetPadLeftMargin(0.16);
		gStyle->SetPadRightMargin(0.16);
		gStyle->SetFrameFillColor(0);
	   	gStyle->SetFrameFillStyle(0);
		//gStyle->SetOptTitle(0);
		//gStyle->SetOptStat(0);
    }
    //gStyle->SetOptTitle(kFALSE);

    //gStyle->SetOptStat(0);
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
				new TGLabel(subframe_selected_fits, "list of selected fit results"));
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
	TGTextButton* plot_all_button = new TGTextButton(this, new TGHotString(" draw overview of selected fit results "));
	plot_all_button->Connect("Clicked()","TrpwaPlotAmpsFrame",this,"Plot_All_selected()");
	// Drawing of the spin totals only
	//TGTextButton* plot_spin_totals_button = new TGTextButton(this, new TGHotString(" draw all wave's spin totals "));
	//plot_spin_totals_button->Connect("Clicked()","TrpwaPlotAmpsFrame",this,"Plot_All_selected_Spin_totals()");

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

	button_show_most_likely = new TGCheckButton(this, new TGHotString("plot most likely solution only"));
	// if use normalization not given then it is assumed not to be available
	if (!draw_most_likely){
		button_show_most_likely->SetState(kButtonUp);
	} else {
		button_show_most_likely->SetState(kButtonDown);
	}
	button_show_most_likely->Connect("Clicked()","TrpwaPlotAmpsFrame",this,"Set_show_most_likely()");

	button_draw_datainfo = new TGCheckButton(this, new TGHotString("draw data info"));
	// if use normalization not given then it is assumed not to be available
	if (!draw_datainfo){
		button_draw_datainfo->SetState(kButtonUp);
	} else {
		button_draw_datainfo->SetState(kButtonDown);
	}
	button_draw_datainfo->Connect("Clicked()","TrpwaPlotAmpsFrame",this,"Set_draw_datainfo()");


	// button to plot the selection
	TGTextButton* plot_button = new TGTextButton(this, new TGHotString(" draw selected waves "));
	plot_button->Connect("Clicked()","TrpwaPlotAmpsFrame",this,"Plot_selected_wave()");

	// canvas where the results are drawn to
	TRootEmbeddedCanvas* frame_selected_waves = new TRootEmbeddedCanvas(" canvas_selected_waves ",
			this, 1000, 600, kSunkenFrame|kDoubleBorder , 0xffffff);
	canvas_selected_waves = frame_selected_waves->GetCanvas();
	canvas_selected_waves->Divide(2,2);

	// slider to select a mass range for all 4 plots
	slider_mass_range = new TGDoubleHSlider(this);
	slider_mass_range->Connect("Released()","TrpwaPlotAmpsFrame",this,"Set_Mass_range()");
	slider_mass_range->SetPosition(0.,1.);

	// button to save the plot
	TGTextButton* save_plot_button = new TGTextButton(this, new TGHotString(" save current plot "));
	save_plot_button->Connect("Clicked()","TrpwaPlotAmpsFrame",this,"Save_plot()");

	this->AddFrame(frame_rootfile_selections, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX,1,1,1,1));
	this->AddFrame(plot_all_button, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX,1,1,1,1));
	//this->AddFrame(plot_spin_totals_button, new TGLayoutHints(kLHintsTop | kLHintsLeft |
	//			kLHintsExpandX,1,1,1,1));
	this->AddFrame(frame_partial_wave_selections, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX,1,1,1,1));
	this->AddFrame(button_show_most_likely,new TGLayoutHints(kLHintsCenterX,1,1,1,1));
	this->AddFrame(button_draw_datainfo, new TGLayoutHints(kLHintsCenterX,1,1,1,1));
	this->AddFrame(plot_button, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX,1,1,1,1));
	this->AddFrame(frame_selected_waves);
	this->AddFrame(slider_mass_range, new TGLayoutHints(kLHintsTop | kLHintsLeft |
				kLHintsExpandX,1,1,1,1));
	this->AddFrame(save_plot_button, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX,1,1,1,1));

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
		available_fit_result_paths[fit_result_title] = fit_result_path;
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
			// add also a new fit result to the fitresult container
			vector<int> _one_fitresult;
			_one_fitresult.push_back(selected_fit_results.size()-1);
			TGIconLBEntry* entry = new TGIconLBEntry(box_selected_fits, pFitFile, fitname.c_str(), Get_Icon(_one_fitresult, false));
			box_selected_fits->AddEntry(entry,new TGLayoutHints(kLHintsTop | kLHintsLeft |
									kLHintsExpandX,1,1,1,1));
			//box_selected_fits->AddEntry(fitname.c_str(), pFitFile);
			current_fit_result = selected_fit;
			box_selected_fits->Select(pFitFile, false);

			/*
			TImage *img = TImage::Create();
			img->FillRectangle("#FFFFFF", 0, 0, 100, 100);
			img->BeginPaint();
			img->DrawCircle(50, 50, 35, "#FF0000", -1);
			img->EndPaint();
			img->Scale(16,16); // resize to create a smaller icon
			// now create a TGPicture from our image
			const TGPicture *pic = gClient->GetPicturePool()->GetPicture("pic_name", img->GetPixmap(), img->GetMask());
			TGIconLBEntry* entry = new TGIconLBEntry(box_available_waves, -10, "testentry", pic);//, UInt_t w = 0, Style_t s = 0, UInt_t options = kHorizontalFrame, Pixel_t back = GetWhitePixel())
			box_available_waves->AddEntry(entry,new TGLayoutHints(kLHintsTop | kLHintsLeft |
					kLHintsExpandX,1,1,1,1));
			delete img;*/

			// search now for waves to add to the list of waves
			selected_fit_result_graphs[fitname] = new Tfitresult(fitname, current_fit_result, available_colors[(int)selected_fit_results.size()-1]);
			bool scan_fit_result(true);
			Tdirmapit fit_result_path_it = available_fit_result_paths.find(fitname);
			if (fit_result_path_it != available_fit_result_paths.end() && FileExists(fit_result_path_it->second+"/fit_result_plots.root")){
				if (!selected_fit_result_graphs[fitname]->Read_from_file(fit_result_path_it->second+"/fit_result_plots.root")){
					cout << " Warning in TrpwaPlotAmpsFrame::Add_Fit(): could not read graphs from file " << endl;
				} else {
					scan_fit_result = false;
				}
			}
			vector<string> wavelist;
			if (scan_fit_result){
				wavelist = selected_fit_result_graphs[fitname]->Scan_Fit_Result();
			} else {
				wavelist = selected_fit_result_graphs[fitname]->Get_available_waves();
			}
			for (vector<string>::iterator it = wavelist.begin(); it != wavelist.end(); it++){
				available_waves[*it].push_back((int)selected_fit_results.size()-1);
			}
			// renew the drop down boxes
			box_available_waves->RemoveAll();
			box_available_anchor_waves->RemoveAll();
			current_anchor_wave = "";
			current_wave = "";
			int wavecounter(0);
			for (Twavemapit it = available_waves.begin(); it != available_waves.end(); it++){
				wavecounter++; // index is not used
				TGIconLBEntry* entry = new TGIconLBEntry(box_available_waves, wavecounter, it->first.c_str(), Get_Icon(it->second));
				box_available_waves->AddEntry(entry,new TGLayoutHints(kLHintsTop | kLHintsLeft |
						kLHintsExpandX,1,1,1,1));
				entry = new TGIconLBEntry(box_available_anchor_waves, wavecounter, it->first.c_str(), Get_Icon(it->second));
				box_available_anchor_waves->AddEntry(entry,new TGLayoutHints(kLHintsTop | kLHintsLeft |
				kLHintsExpandX,1,1,1,1));
			}
		}
	//selected_fit_results[];
	//cout << pFitFile << endl;
	//cout << " calling add fit " << endl;
	}
}

// available colors for index markers
//vector<TColor_struct> available_colors;

// creates some available colors for markers in available_colors
void TrpwaPlotAmpsFrame::Create_available_colors(){
	if (available_colors.size() > 0) return; // initialize only once
	TColor_struct color;
	for (int i = 0; i < 50; i++){ // thus colors will repeat 50 times
		/*
		color.rootcolorindex = kBlack;
		color.colorhexcode   = "#000000";
		available_colors.push_back(color);
		color.rootcolorindex = kRed;
		color.colorhexcode   = "#FF0000";
		available_colors.push_back(color);
		color.rootcolorindex = kGreen;
		color.colorhexcode   = "#00FF00";
		available_colors.push_back(color);
		color.rootcolorindex = kBlue;
		color.colorhexcode   = "#0000FF";
		available_colors.push_back(color);*/

		//ULong_t RGB2Pixel(Float_t r, Float_t g, Float_t b)
		//const char * PixelAsHexString(ULong_t pixel)
		// Int_t GetColor(Int_t r, Int_t g, Int_t b)
		for (int i = 1; i < 10; i++){
			if (i == 5 || i == 7) i++; // skip yellow and magenta
			TColor* rootcolor = gROOT->GetColor(i);
			color.rootcolorindex = i;
			color.colorhexcode   = rootcolor->AsHexString();
			available_colors.push_back(color);
		}
	}
}

const TGPicture* TrpwaPlotAmpsFrame::Get_Icon(vector<int>& fit_references, bool divide){
	if ((divide && selected_fit_results.size() == 0) || (!divide && fit_references.size() < 1)){
		cout << " Errot in TrpwaPlotAmpsFrame::Get_Icon(): no fit results are given to refer to! " << endl;
		return NULL;
	}
	TImage *img = TImage::Create();
	if (divide)
		img->FillRectangle("#FFFFFF", 0, 0, 10 * selected_fit_results.size(), 12);
	else
		img->FillRectangle("#FFFFFF", 0, 0, 10 * fit_references.size(), 12);
	img->BeginPaint();
	for (unsigned int i = 0; i < fit_references.size(); i++){
		int posx;
		if (divide)
			posx = 10 * fit_references[i] + 5;
		else
			posx = 10 * i + 5;
		//cout << " drawing " << posx << " " << available_colors[fit_references[i]].colorhexcode << endl;
		img->DrawCircle(posx, 6, 4, available_colors[fit_references[i]].colorhexcode.c_str(), -1);
		//img->Scale(16,16); // resize to create a smaller icon
	}
	img->EndPaint();
	// create an arbitrary name for the new icon in the pool of pictures
	stringstream pic_name;
	static int counter;
	pic_name << "generated_reference_icon_" << counter++;
	const TGPicture *pic = gClient->GetPicturePool()->GetPicture(pic_name.str().c_str(), img->GetPixmap(), img->GetMask());
	delete img;
	return pic;
}


// will be called when an added fit file was selected
// if the user confirms the selected item will be removed from the list
void TrpwaPlotAmpsFrame::Remove_Fit(int pFitFile){
	cout << pFitFile << endl;
	cout << " call to remove a fit is not implemented yet" << endl;
}

// a wave will be selected from the list of available waves
void TrpwaPlotAmpsFrame::Select_Wave(int pWave){
	//cout << pWave << endl;
	//cout << " calling select wave " << endl;
	TGTextLBEntry *filePointer = (TGTextLBEntry *)box_available_waves->GetSelectedEntry();
	string selected_wave = filePointer->GetTitle();
	if (available_waves.find(selected_wave)!=available_waves.end()){
		cout << " selected wave is " << selected_wave << endl;
		current_wave = selected_wave;
	} else {
		current_wave = "";
	}
}

// a wave will be selected for the anchor wave (phase determination)
void TrpwaPlotAmpsFrame::Select_Anchor_Wave(int pWave){
	//string selected_wave = box_available_anchor_waves->GetTextEntry()->GetText();
	TGTextLBEntry *filePointer = (TGTextLBEntry *)box_available_anchor_waves->GetSelectedEntry();
	string selected_wave = filePointer->GetTitle();
	if (available_waves.find(selected_wave)!=available_waves.end()){
		cout << " selected anchor wave is " << selected_wave << endl;
		current_anchor_wave = selected_wave;
	} else {
		current_anchor_wave = "";
	}
}

void TrpwaPlotAmpsFrame::Plot_All_selected(){
	const string branchName = "fitResult_v2";
	static int some_name_counter(0); // due to shitty root implementation for unique names
	stringstream some_name;
	canvas_selected_waves->cd(1); gPad->Clear();
	canvas_selected_waves->cd(2); gPad->Clear();
	canvas_selected_waves->cd(3); gPad->Clear();
	canvas_selected_waves->cd(4); gPad->Clear();
	// to do set the pad values
	TText label;
	canvas_selected_waves->cd(1);
	gPad->Range(0, 0, 1, 1);
	label.DrawText(0.1,0.5, "likelihood distributions");

	canvas_selected_waves->cd(2);
	gPad->Range(0, 0, 1, 1);
	int linecounter(0);
	for (Tfilemapit it = selected_fit_results.begin(); it != selected_fit_results.end(); it++){
		TTree* selected_tree = (TTree*)it->second;
		Tfitresult* fitresult = selected_fit_result_graphs[selected_tree->GetName()];
		label.SetTextColor(fitresult->graphcolor.rootcolorindex);
		label.DrawText(0.1,0.1+linecounter*0.06, fitresult->name.c_str());
		linecounter++;
	}
	label.SetTextColor(kBlack);
	canvas_selected_waves->Print("Fit_overview.ps(");
	canvas_selected_waves->cd(1); gPad->Clear();
	canvas_selected_waves->cd(2); gPad->Clear();
	canvas_selected_waves->cd(3); gPad->Clear();
	canvas_selected_waves->cd(4); gPad->Clear();

	for (int ipad = 0; ipad < 5; ipad++){
		if (plotted_graphs[ipad]){
			plotted_graphs[ipad]->Clear();
			delete plotted_graphs[ipad];
			plotted_graphs[ipad] = NULL;
		}
		if (plotted_most_likely_graphs[ipad]){
			plotted_most_likely_graphs[ipad]->Clear();
			delete plotted_graphs[ipad];
			plotted_most_likely_graphs[ipad]=NULL;
		}
	}
	masscutlow = -1.;
	masscuthigh = -1.;
	slider_mass_range->SetPosition(0.,1.);

	cout << " drawing all intensities into a file " << endl;

	TGraphErrors* _graph;
	vector<TMultiGraph*> _multi_graphs;

	// draw the total intensities
	int ipad(0);
	for (int mostlikely = 0; mostlikely < 2; mostlikely++){
		ipad++;
		canvas_selected_waves->cd(ipad);
		//string last_JPCM = "";
		TMultiGraph* _graphs = NULL;
		for (Tfilemapit it = selected_fit_results.begin(); it != selected_fit_results.end(); it++){
			TTree* selected_tree = (TTree*)it->second;
			Tfitresult* fitresult = selected_fit_result_graphs[selected_tree->GetName()];
			_graph = fitresult->Get_Total_Intensity(mostlikely);
			if (_graph){
				if (!_graphs){
					_graphs = new TMultiGraph();
					some_name.str(""); some_name << "multi_graph_number_" << some_name_counter++;
					_graphs->SetName(some_name.str().c_str());
					_multi_graphs.push_back(_graphs);
				}
				_graphs->Add(_graph);
				_graphs->SetTitle(_graph->GetTitle());
				gPad->Clear();
				_graphs->Draw("APZ");
				_graphs->GetXaxis()->SetTitle(_graph->GetXaxis()->GetTitle());
				_graphs->GetYaxis()->SetTitle(_graph->GetYaxis()->GetTitle());
			} else {
				cout << " no total intensity found in " << selected_tree->GetName()<< endl;
			}
		}
		if (draw_datainfo){
			label_preliminary->Draw();
			label_data_info->Draw();
		}
	}

	// draw the likelihood distributions
	for (int mostlikely = 0; mostlikely < 2; mostlikely++){
		ipad++;
		canvas_selected_waves->cd(ipad);
		//string last_JPCM = "";
		TMultiGraph* _graphs = NULL;
		vector<TGraphErrors*> _temp_container; // to keep the graphs temporary accessible
		for (Tfilemapit it = selected_fit_results.begin(); it != selected_fit_results.end(); it++){
			TTree* selected_tree = (TTree*)it->second;
			Tfitresult* fitresult = selected_fit_result_graphs[selected_tree->GetName()];
			_graph = fitresult->Get_Likelihood(mostlikely);
			if (_graph){
				if (!_graphs){
					_graphs = new TMultiGraph();
					some_name.str(""); some_name << "multi_graph_number_" << some_name_counter++;
					_graphs->SetName(some_name.str().c_str());
					_multi_graphs.push_back(_graphs);
				}
				_temp_container.push_back(_graph);
				_graphs->Add(_graph);
				_graphs->SetTitle(_graph->GetTitle());
				gPad->Clear();
				_graphs->Draw("APZ");
				_graphs->GetXaxis()->SetTitle(_graph->GetXaxis()->GetTitle());
				_graphs->GetYaxis()->SetTitle(_graph->GetYaxis()->GetTitle());
			} else {
				cout << " no likelihood distribution found in " << selected_tree->GetName()<< endl;
			}
		}
		// show only the most likely distributions among several fit results
		// by deleting points with lower likelihoods
		if (mostlikely == 1 && _temp_container.size() > 1){
			cout << " filtering only the most likely log likelihoods among all most likely" << endl;
			for (unsigned int igraph = 0; igraph < _temp_container.size()-1; igraph++){
				TGraphErrors* _graph = _temp_container[igraph];
				_graph->Sort();
				for (unsigned int jgraph = igraph+1; jgraph < _temp_container.size(); jgraph++){
					TGraphErrors* _comp_graph = _temp_container[jgraph];
					_comp_graph->Sort();
					for (int i = 0; i < _graph->GetN(); i++){
						for (int j = 0; j < _comp_graph->GetN(); j++){
							if (_graph->GetX()[i] == _comp_graph->GetX()[j]){ // found to same mass bins
								if (_graph->GetY()[i] < _comp_graph->GetY()[j]){ // _graph is more likely
									_comp_graph->SetPoint(j, _comp_graph->GetX()[j], 100.);
								} else { // _comp_graph is more likely
									_graph->SetPoint(i, _graph->GetX()[i], 100.);
								}
								break; // only one entry per mass bin is expected here
							}
						}
					}
				}
			}
		}
		if (draw_datainfo){
			label_preliminary->Draw();
			label_data_info->Draw();
		}
	}

	canvas_selected_waves->Print("Fit_overview.ps");
	for (ipad = 4; ipad > 0; ipad --){
		canvas_selected_waves->cd(ipad);
		gPad->Clear();
	}

	// draw the evidence distributions
	ipad = 0;
	for (int mostlikely = 0; mostlikely < 2; mostlikely++){
		ipad++;
		canvas_selected_waves->cd(ipad);
		//string last_JPCM = "";
		TMultiGraph* _graphs = NULL;
		vector<TGraphErrors*> _temp_container; // to keep the graphs temporary accessible
		for (Tfilemapit it = selected_fit_results.begin(); it != selected_fit_results.end(); it++){
			TTree* selected_tree = (TTree*)it->second;
			Tfitresult* fitresult = selected_fit_result_graphs[selected_tree->GetName()];
			_graph = fitresult->Get_Evidence(mostlikely);
			if (_graph){
				if (!_graphs){
					_graphs = new TMultiGraph();
					some_name.str(""); some_name << "multi_graph_number_" << some_name_counter++;
					_graphs->SetName(some_name.str().c_str());
					_multi_graphs.push_back(_graphs);
				}
				_temp_container.push_back(_graph);
				_graphs->Add(_graph);
				_graphs->SetTitle(_graph->GetTitle());
				gPad->Clear();
				_graphs->Draw("APZ");
				_graphs->GetXaxis()->SetTitle(_graph->GetXaxis()->GetTitle());
				_graphs->GetYaxis()->SetTitle(_graph->GetYaxis()->GetTitle());
			} else {
				cout << " no evidence distribution found in " << selected_tree->GetName()<< endl;
			}
		}
		// show only the most evident distributions among several fit results
		// by deleting points with lower likelihoods
		if (mostlikely == 1 && _temp_container.size() > 1){
			cout << " filtering only the most evident solutions among all evidences" << endl;
			for (unsigned int igraph = 0; igraph < _temp_container.size()-1; igraph++){
				TGraphErrors* _graph = _temp_container[igraph];
				_graph->Sort();
				for (unsigned int jgraph = igraph+1; jgraph < _temp_container.size(); jgraph++){
					TGraphErrors* _comp_graph = _temp_container[jgraph];
					_comp_graph->Sort();
					for (int i = 0; i < _graph->GetN(); i++){
						for (int j = 0; j < _comp_graph->GetN(); j++){
							if (_graph->GetX()[i] == _comp_graph->GetX()[j]){ // found to same mass bins
								if (_graph->GetY()[i] > _comp_graph->GetY()[j]){ // _graph is more likely
									_comp_graph->SetPoint(j, _comp_graph->GetX()[j], -100.);
								} else { // _comp_graph is more likely
									_graph->SetPoint(i, _graph->GetX()[i], -1.e3);
								}
								break; // only one entry per mass bin is expected here
							}
						}
					}
				}
			}
		}
		if (draw_datainfo){
			label_preliminary->Draw();
			label_data_info->Draw();
		}
	}

	canvas_selected_waves->Print("Fit_overview.ps");
	for (ipad = 4; ipad > 0; ipad --){
		canvas_selected_waves->cd(ipad);
		gPad->Clear();
	}

	// Collect the available spin totals
	map<string, int> _available_jpcs; // list of available JPC
	for (Tfilemapit it = selected_fit_results.begin(); it != selected_fit_results.end(); it++){
		TTree* selected_tree = (TTree*)it->second;
		Tfitresult* fitresult = selected_fit_result_graphs[selected_tree->GetName()];
		vector<string>& result = fitresult->Get_available_Spin_Totals();
		// view the available spin totals
		for (vector<string>::iterator it = result.begin(); it != result.end(); it++){
			_available_jpcs[*it]++;
		}
	}

	// now draw the spin totals
	ipad = 0;
	for (map<string, int>::iterator it = _available_jpcs.begin(); it != _available_jpcs.end(); it++){
		ipad++;
		if (ipad > 4){
			canvas_selected_waves->Print("Fit_overview.ps");
			for (ipad = 4; ipad > 0; ipad --){
				canvas_selected_waves->cd(ipad);
				gPad->Clear();
			}
			ipad = 1;
		}
		canvas_selected_waves->cd(ipad);
		TMultiGraph* _graphs = new TMultiGraph();
		some_name.str(""); some_name << "multi_graph_number_" << some_name_counter++;
		_graphs->SetName(some_name.str().c_str());
		_multi_graphs.push_back(_graphs);
		// collect all spin totals from the available results
		for (Tfilemapit itresult = selected_fit_results.begin(); itresult != selected_fit_results.end(); itresult++){
			TTree* selected_tree = (TTree*)itresult->second;
			Tfitresult* fitresult = selected_fit_result_graphs[selected_tree->GetName()];
			_graph = fitresult->Get_Spin_Total(it->first);
			if(_graph){
				_graphs->Add(_graph);
				_graphs->SetTitle(_graph->GetTitle());
				gPad->Clear();
				_graphs->Draw("APZ");
				_graphs->GetXaxis()->SetTitle(_graph->GetXaxis()->GetTitle());
				_graphs->GetYaxis()->SetTitle(_graph->GetYaxis()->GetTitle());
				gPad->Update();
			}
		}
		if (draw_datainfo){
			label_preliminary->Draw();
			label_data_info->Draw();
		}
	}
	canvas_selected_waves->Print("Fit_overview.ps");

	canvas_selected_waves->cd(1); gPad->Clear();
	canvas_selected_waves->cd(2); gPad->Clear();
	canvas_selected_waves->cd(3); gPad->Clear();
	canvas_selected_waves->cd(4); gPad->Clear();

	//sort(available_waves.begin(), available_waves.end());

	// drawing all intensities
	for (int mostlikely = 0; mostlikely < 2; mostlikely++){
		int ipad(0);
		string last_JPCM = "";
		for (Twavemapit it = available_waves.begin(); it != available_waves.end(); it++){
			string wavename = it->first;

			int J,P,C,M,refl,l,s;
			string iso1, iso2;
			GetJPCMreflISO1lsISO2(wavename,J,P,C,M,refl,iso1,iso2,l,s);
			stringstream _jpcm;// = (wavename.substr(2,4));
			_jpcm << J << P;
			if (last_JPCM != _jpcm.str()){
				if (last_JPCM != ""){
					ipad = 4;
				}
				last_JPCM = _jpcm.str();
			}

			ipad++;
			if (ipad > 4){
				canvas_selected_waves->Print("Fit_overview.ps");
				for (ipad = 4; ipad > 0; ipad --){
					canvas_selected_waves->cd(ipad);
					gPad->Clear();
				}
				ipad = 1;
			}
			canvas_selected_waves->cd(ipad);
			gPad->Clear();
			// for all selected waves
			TMultiGraph* _graphs = NULL;
			for (Tfilemapit it = selected_fit_results.begin(); it != selected_fit_results.end(); it++){
				TTree* selected_tree = (TTree*)it->second;
				Tfitresult* fitresult = selected_fit_result_graphs[selected_tree->GetName()];
				_graph = fitresult->Get_Intensity(wavename, mostlikely);
				if (_graph){
					if (!_graphs){
						_graphs = new TMultiGraph();
						some_name.str(""); some_name << "multi_graph_number_" << some_name_counter++;
						_graphs->SetName(some_name.str().c_str());
						_multi_graphs.push_back(_graphs);
					}
					_graphs->Add(_graph);
					_graphs->SetTitle(_graph->GetTitle());
					gPad->Clear();
					_graphs->Draw("APZ");
					_graphs->GetXaxis()->SetTitle(_graph->GetXaxis()->GetTitle());
					_graphs->GetYaxis()->SetTitle(_graph->GetYaxis()->GetTitle());
				} else {
					cout << " no wave named " << wavename << " in " << selected_tree->GetName() << " found "<< endl;
				}
			}
			if (draw_datainfo){
				label_preliminary->Draw();
				label_data_info->Draw();
			}
			gPad->Update();
		}
		canvas_selected_waves->Print("Fit_overview.ps");
		for (ipad = 4; ipad > 0; ipad --){
			canvas_selected_waves->cd(ipad);
			gPad->Clear();
		}
	}

	canvas_selected_waves->cd(1);
	gPad->Range(0, 0, 1, 1);
	label.DrawText(0.1,0.5, "to do: summary");
	canvas_selected_waves->Print("Fit_overview.ps)");
	gPad->Clear();
	// show the result
	//if (system("acroread Fit_overview.pdf")){
		if (system("evince Fit_overview.ps")){
			if (system("gv Fit_overview.ps")){
				cout << " no reader for Fit_overview.ps found! " << endl;
			}
		}
	//}

	// delete all multi graphs
	for (vector<TMultiGraph*>::iterator it = _multi_graphs.begin(); it != _multi_graphs.end(); it++){
		(*it)->Clear();
		delete (*it);
	}
/*
	for (Tfilemapit it = selected_fit_results.begin(); it != selected_fit_results.end(); it++){
		TTree* selected_tree = (TTree*)it->second;
		Tfitresult* fitresult = selected_fit_result_graphs[selected_tree->GetName()];
		if (!fitresult){
			cout << " TrpwaPlotAmpsFrame::Plot_selected_wave(): Hm, no fit result found! " << endl;
			continue;
		}
		if (current_wave == "") continue;
		ipad++;
		if ()
		canvas_selected_waves->cd(ipad);
		_graph = fitresult->Get_Intensity(current_wave, false);
		if (_graph){
			if (!plotted_graphs[ipad]){
				plotted_graphs[ipad] = new TMultiGraph();
				some_name.str(""); some_name << "multi_graph_number_" << some_name_counter++;
				plotted_graphs[ipad]->SetName(some_name.str().c_str());
			}
			if (masscutlow < 0 || fitresult->waves[current_wave].mass_low < masscutlow)
				 masscutlow = fitresult->waves[current_wave].mass_low;
			if (masscuthigh < 0 || fitresult->waves[current_wave].mass_high > masscuthigh)
				masscuthigh = fitresult->waves[current_wave].mass_high;
			plotted_graphs[ipad]->Add(_graph);
			plotted_graphs[ipad]->SetTitle(_graph->GetTitle());
			gPad->Clear();
			plotted_graphs[ipad]->Draw("APZ");
			plotted_graphs[ipad]->GetXaxis()->SetTitle(_graph->GetXaxis()->GetTitle());
			plotted_graphs[ipad]->GetYaxis()->SetTitle(_graph->GetYaxis()->GetTitle());
		} else {
			cout << " no wave named " << current_wave << " in " << selected_tree->GetName() << " found "<< endl;
		}
		_graph = fitresult->Get_Intensity(current_wave, true);
		if (_graph){
			if (!plotted_most_likely_graphs[ipad]){
				plotted_most_likely_graphs[ipad] = new TMultiGraph();
				some_name.str(""); some_name << "multi_graph_number_" << some_name_counter++;
				plotted_most_likely_graphs[ipad]->SetName(some_name.str().c_str());
			}
			plotted_most_likely_graphs[ipad]->Add(_graph);
			plotted_most_likely_graphs[ipad]->SetTitle(_graph->GetTitle());
			gPad->Clear();
			plotted_most_likely_graphs[ipad]->Draw("APZ");
			plotted_most_likely_graphs[ipad]->GetXaxis()->SetTitle(_graph->GetXaxis()->GetTitle());
			plotted_most_likely_graphs[ipad]->GetYaxis()->SetTitle(_graph->GetYaxis()->GetTitle());
		}
	}*/
	return;
	/*
	for (Tfilemapit it = selected_fit_results.begin(); it != selected_fit_results.end(); it++){
		TTree* selected_tree = (TTree*)it->second;
		cout << " Drawing all results from " << it->first << ". Please be patient... " << endl;
		plotAllIntensities(selected_tree, true);
		cout << " done " << endl;
		//box_available_fits->AddEntry(it->first.c_str(), (long int) it->second);
	}*/
	if (current_fit_result > 0){
		cout << " Drawing all results from " << current_fit_result->GetName() << ". Please be patient... " << endl;
		plotAllIntensities(current_fit_result, true);
		cout << " done " << endl;
	} else {
		cout <<  " no fit selected " << endl;
	}
}

void TrpwaPlotAmpsFrame::Plot_All_selected_Spin_totals(){
	cout << " plotSpinTotals() macro does not work yet properly and may crash! " << endl;
	//return;
	if (current_fit_result > 0){
		cout << " Drawing all spin totals from " << current_fit_result->GetName() << ". Please be patient... " << endl;
		plotSpinTotals(current_fit_result);
		cout << " done " << endl;
	} else {
		cout <<  " no fit selected " << endl;
	}
}

void TrpwaPlotAmpsFrame::Save_plot(){
	TGFileInfo fileinfo;
	const char* filetypes[] = {"pdf","*.pdf","eps","*.eps","gif","*.gif","png","*.png","root","*.root", 0, 0};
	fileinfo.fFileTypes = filetypes;
	stringstream filename;
	filename << current_wave;
	if (current_anchor_wave != "")
		filename << "_vs_" << current_anchor_wave;
	char *cstr = new char [filename.str().size()+1];
	strcpy (cstr, filename.str().c_str());
	fileinfo.fFilename = cstr;
	TGFileDialog* filedialog = new TGFileDialog(gClient->GetRoot(), this, kFDSave, &fileinfo);
	if (filedialog && fileinfo.fFilename) {
		filename.str("");
		filename << fileinfo.fFilename;
		filename << "." << fileinfo.fFileTypes[fileinfo.fFileTypeIdx];
		canvas_selected_waves->Print(filename.str().c_str());
	}
	/*
	int returncode;
	stringstream filename;
	filename << current_wave;
	if (current_anchor_wave != "")
		filename << "_vs_" << current_anchor_wave;
	TGMsgBox* userrespondbox = new TGMsgBox(gClient->GetRoot(), this, "save current output",
			("Do you want to save this plot as\n"+filename.str()+".pdf/.root?").c_str(),
			kMBIconQuestion, (kMBYes | kMBNo), &returncode);
	if (!userrespondbox) cout << " this will be not executed " << endl; // to prevent compiler warnings
	if (returncode == kMBYes){
		canvas_selected_waves->Print((filename.str()+".pdf").c_str());
		canvas_selected_waves->Print((filename.str()+".root").c_str());
	}*/
}

// I could try to copy the fitResult class but I prefer to store
// only the for me relevant values
struct Twavebin_info {
	double massBinCenter;
	double logLikelihood;
	double waveIntensityA;
	double waveIntensityB;
	double waveIntensityErrA;
	double waveIntensityErrB;
	double phase;
	double phaseErr;
	double coherence;
	double coherenceErr;
};

void TrpwaPlotAmpsFrame::Plot_selected_wave(){
	const string branchName = "fitResult_v2";
	static int some_name_counter(0); // due to shitty root implementation for unique names
	stringstream some_name;
	canvas_selected_waves->cd(1); gPad->Clear();
	canvas_selected_waves->cd(2); gPad->Clear();
	canvas_selected_waves->cd(3); gPad->Clear();
	canvas_selected_waves->cd(4); gPad->Clear();
	for (int ipad = 0; ipad < 5; ipad++){
		if (plotted_graphs[ipad]){
			plotted_graphs[ipad]->Clear();
			delete plotted_graphs[ipad];
			plotted_graphs[ipad] = NULL;
		}
		if (plotted_most_likely_graphs[ipad]){
			plotted_most_likely_graphs[ipad]->Clear();
			delete plotted_graphs[ipad];
			plotted_most_likely_graphs[ipad]=NULL;
		}
	}
	masscutlow = -1.;
	masscuthigh = -1.;
	slider_mass_range->SetPosition(0.,1.);

	TGraphErrors* _graph;

	for (Tfilemapit it = selected_fit_results.begin(); it != selected_fit_results.end(); it++){
		TTree* selected_tree = (TTree*)it->second;
		Tfitresult* fitresult = selected_fit_result_graphs[selected_tree->GetName()];
		if (!fitresult){
			cout << " TrpwaPlotAmpsFrame::Plot_selected_wave(): Hm, no fit result found! " << endl;
			continue;
		}
		if (current_wave == "") continue;
		int ipad = 1;
		canvas_selected_waves->cd(ipad);
		_graph = fitresult->Get_Intensity(current_wave, false);
		if (_graph){
			if (!plotted_graphs[ipad]){
				plotted_graphs[ipad] = new TMultiGraph();
				some_name.str(""); some_name << "multi_graph_number_" << some_name_counter++;
				plotted_graphs[ipad]->SetName(some_name.str().c_str());
			}
			if (masscutlow < 0 || fitresult->waves[current_wave].mass_low < masscutlow)
				 masscutlow = fitresult->waves[current_wave].mass_low;
			if (masscuthigh < 0 || fitresult->waves[current_wave].mass_high > masscuthigh)
				masscuthigh = fitresult->waves[current_wave].mass_high;
			plotted_graphs[ipad]->Add(_graph);
			plotted_graphs[ipad]->SetTitle(_graph->GetTitle());
			gPad->Clear();
			plotted_graphs[ipad]->Draw("APZ");
			plotted_graphs[ipad]->GetXaxis()->SetTitle(_graph->GetXaxis()->GetTitle());
			plotted_graphs[ipad]->GetYaxis()->SetTitle(_graph->GetYaxis()->GetTitle());
		} else {
			cout << " no wave named " << current_wave << " in " << selected_tree->GetName() << " found "<< endl;
		}
		_graph = fitresult->Get_Intensity(current_wave, true);
		if (_graph){
			if (!plotted_most_likely_graphs[ipad]){
				plotted_most_likely_graphs[ipad] = new TMultiGraph();
				some_name.str(""); some_name << "multi_graph_number_" << some_name_counter++;
				plotted_most_likely_graphs[ipad]->SetName(some_name.str().c_str());
			}
			plotted_most_likely_graphs[ipad]->Add(_graph);
			plotted_most_likely_graphs[ipad]->SetTitle(_graph->GetTitle());
			gPad->Clear();
			plotted_most_likely_graphs[ipad]->Draw("APZ");
			plotted_most_likely_graphs[ipad]->GetXaxis()->SetTitle(_graph->GetXaxis()->GetTitle());
			plotted_most_likely_graphs[ipad]->GetYaxis()->SetTitle(_graph->GetYaxis()->GetTitle());
		}

		if (current_anchor_wave == "") continue;
		ipad = 3;
		canvas_selected_waves->cd(ipad);
		_graph = fitresult->Get_Intensity(current_anchor_wave, false);
		if (_graph){
			if (!plotted_graphs[ipad]){
				plotted_graphs[ipad] = new TMultiGraph();
				some_name.str(""); some_name << "multi_graph_number_" << some_name_counter++;
				plotted_graphs[ipad]->SetName(some_name.str().c_str());
			}
			//if (masscutlow < 0 || fitresult->waves[current_anchor_wave].mass_low < masscutlow)
			//	 masscutlow = fitresult->waves[current_anchor_wave].mass_low;
			//if (masscuthigh < 0 || fitresult->waves[current_anchor_wave].mass_high > masscuthigh)
			//	masscuthigh = fitresult->waves[current_anchor_wave].mass_high;
			plotted_graphs[ipad]->Add(_graph);
			plotted_graphs[ipad]->SetTitle(_graph->GetTitle());
			gPad->Clear();
			plotted_graphs[ipad]->Draw("APZ");
			plotted_graphs[ipad]->GetXaxis()->SetTitle(_graph->GetXaxis()->GetTitle());
			plotted_graphs[ipad]->GetYaxis()->SetTitle(_graph->GetYaxis()->GetTitle());
		} else {
			cout << " no wave named " << current_anchor_wave << " in " << selected_tree->GetName() << " found "<< endl;
		}
		_graph = fitresult->Get_Intensity(current_anchor_wave, true);
		if (_graph){
			if (!plotted_most_likely_graphs[ipad]){
				plotted_most_likely_graphs[ipad] = new TMultiGraph();
				some_name.str(""); some_name << "multi_graph_number_" << some_name_counter++;
				plotted_most_likely_graphs[ipad]->SetName(some_name.str().c_str());
			}
			plotted_most_likely_graphs[ipad]->Add(_graph);
			plotted_most_likely_graphs[ipad]->SetTitle(_graph->GetTitle());
			gPad->Clear();
			plotted_most_likely_graphs[ipad]->Draw("APZ");
			plotted_most_likely_graphs[ipad]->GetXaxis()->SetTitle(_graph->GetXaxis()->GetTitle());
			plotted_most_likely_graphs[ipad]->GetYaxis()->SetTitle(_graph->GetYaxis()->GetTitle());
		}

		ipad = 2;
		canvas_selected_waves->cd(ipad);
		_graph = fitresult->Get_Phase(current_wave, current_anchor_wave, false);
		if (_graph){
			if (!plotted_graphs[ipad]){
				plotted_graphs[ipad] = new TMultiGraph();
				some_name.str(""); some_name << "multi_graph_number_" << some_name_counter++;
				plotted_graphs[ipad]->SetName(some_name.str().c_str());
			}
			//if (masscutlow < 0 || fitresult->waves[current_wave].mass_low < masscutlow)
			//	 masscutlow = fitresult->waves[current_wave].mass_low;
			//if (masscuthigh < 0 || fitresult->waves[current_wave].mass_high > masscuthigh)
			//	masscuthigh = fitresult->waves[current_wave].mass_high;
			plotted_graphs[ipad]->Add(_graph);
			plotted_graphs[ipad]->SetTitle(_graph->GetTitle());
			gPad->Clear();
			plotted_graphs[ipad]->Draw("APZ");
			plotted_graphs[ipad]->GetXaxis()->SetTitle(_graph->GetXaxis()->GetTitle());
			plotted_graphs[ipad]->GetYaxis()->SetTitle(_graph->GetYaxis()->GetTitle());
		} else {
			cout << " no wave named " << current_wave << " in " << selected_tree->GetName() << " found "<< endl;
		}
		_graph = fitresult->Get_Phase(current_wave, current_anchor_wave, true);
		if (_graph){
			if (!plotted_most_likely_graphs[ipad]){
				plotted_most_likely_graphs[ipad] = new TMultiGraph();
				some_name.str(""); some_name << "multi_graph_number_" << some_name_counter++;
				plotted_most_likely_graphs[ipad]->SetName(some_name.str().c_str());
			}
			plotted_most_likely_graphs[ipad]->Add(_graph);
			plotted_most_likely_graphs[ipad]->SetTitle(_graph->GetTitle());
			gPad->Clear();
			plotted_most_likely_graphs[ipad]->Draw("APZ");
			plotted_most_likely_graphs[ipad]->GetXaxis()->SetTitle(_graph->GetXaxis()->GetTitle());
			plotted_most_likely_graphs[ipad]->GetYaxis()->SetTitle(_graph->GetYaxis()->GetTitle());
		}

		ipad = 4;
		canvas_selected_waves->cd(ipad);
		_graph = fitresult->Get_Coherence(current_wave, current_anchor_wave, false);
		if (_graph){
			if (!plotted_graphs[ipad]){
				plotted_graphs[ipad] = new TMultiGraph();
				some_name.str(""); some_name << "multi_graph_number_" << some_name_counter++;
				plotted_graphs[ipad]->SetName(some_name.str().c_str());
			}
			//if (masscutlow < 0 || fitresult->waves[current_wave].mass_low < masscutlow)
			//	 masscutlow = fitresult->waves[current_wave].mass_low;
			//if (masscuthigh < 0 || fitresult->waves[current_wave].mass_high > masscuthigh)
			//	masscuthigh = fitresult->waves[current_wave].mass_high;
			plotted_graphs[ipad]->Add(_graph);
			plotted_graphs[ipad]->SetTitle(_graph->GetTitle());
			gPad->Clear();
			plotted_graphs[ipad]->Draw("APZ");
			plotted_graphs[ipad]->GetXaxis()->SetTitle(_graph->GetXaxis()->GetTitle());
			plotted_graphs[ipad]->GetYaxis()->SetTitle(_graph->GetYaxis()->GetTitle());
		} else {
			cout << " no wave named " << current_wave << " in " << selected_tree->GetName() << " found "<< endl;
		}
		_graph = fitresult->Get_Coherence(current_wave, current_anchor_wave, true);
		if (_graph){
			if (!plotted_most_likely_graphs[ipad]){
				plotted_most_likely_graphs[ipad] = new TMultiGraph();
				some_name.str(""); some_name << "multi_graph_number_" << some_name_counter++;
				plotted_most_likely_graphs[ipad]->SetName(some_name.str().c_str());
			}
			plotted_most_likely_graphs[ipad]->Add(_graph);
			plotted_most_likely_graphs[ipad]->SetTitle(_graph->GetTitle());
			gPad->Clear();
			plotted_most_likely_graphs[ipad]->Draw("APZ");
			plotted_most_likely_graphs[ipad]->GetXaxis()->SetTitle(_graph->GetXaxis()->GetTitle());
			plotted_most_likely_graphs[ipad]->GetYaxis()->SetTitle(_graph->GetYaxis()->GetTitle());
		}
	}

	Set_Mass_range();
}

void TrpwaPlotAmpsFrame::Set_show_most_likely(){
	draw_most_likely = button_show_most_likely->GetState();
	Set_Mass_range();
}

void TrpwaPlotAmpsFrame::Set_draw_datainfo(){
	draw_datainfo = button_draw_datainfo->GetState();
	Set_Mass_range();
}

void TrpwaPlotAmpsFrame::Set_Mass_range(){
	float mass_min, mass_max;
	slider_mass_range->GetPosition(mass_min, mass_max);
	mass_min = mass_min * (masscuthigh-masscutlow) + masscutlow;
	mass_max = mass_max * (masscuthigh-masscutlow) + masscutlow;
	//mass_min /= 1000.; // drawn in GeV
	//mass_max /= 1000.;
	//cout << " setting from "<< mass_min << " to " << mass_max << endl;
	for (int ipad = 1; ipad < 5; ipad++){
		canvas_selected_waves->cd(ipad);
		gPad->Clear();
		if (!draw_most_likely && plotted_graphs[ipad]){
			plotted_graphs[ipad]->Draw("APZ");
			plotted_graphs[ipad]->GetXaxis()->SetRangeUser(mass_min, mass_max);
			//gPad->Update();
			switch (ipad){
			case 1: // intensity
				plotted_graphs[ipad]->GetYaxis()->SetRangeUser(0.,plotted_graphs[ipad]->GetYaxis()->GetXmax());
				break;
			case 2: // phase
				//plotted_graphs[ipad]->GetYaxis()->SetRangeUser(-600., 600);
				break;
			case 3: // anchor wave
				plotted_graphs[ipad]->GetYaxis()->SetRangeUser(0.,plotted_graphs[ipad]->GetYaxis()->GetXmax());
				break;
			case 4:
				plotted_graphs[ipad]->GetYaxis()->SetRangeUser(0., 1.);
			default:
				break; // do nothing
			}
		}
		if (draw_most_likely && plotted_most_likely_graphs[ipad]){
			//canvas_selected_waves->cd(ipad);
			//gPad->Clear();
			plotted_most_likely_graphs[ipad]->Draw("APZ");
			plotted_most_likely_graphs[ipad]->GetXaxis()->SetRangeUser(mass_min, mass_max);
			switch (ipad){
			case 1: // intensity
				plotted_most_likely_graphs[ipad]->GetYaxis()->SetRangeUser(0.,plotted_most_likely_graphs[ipad]->GetYaxis()->GetXmax());
				break;
			case 2: // phase
				//plotted_most_likely_graphs[ipad]->GetYaxis()->SetRangeUser(-600., 600);
				break;
			case 3: // anchor wave
				plotted_most_likely_graphs[ipad]->GetYaxis()->SetRangeUser(0.,plotted_most_likely_graphs[ipad]->GetYaxis()->GetXmax());
				break;
			case 4:
				plotted_most_likely_graphs[ipad]->GetYaxis()->SetRangeUser(0., 1.);
			default:
				break; // do nothing
			}
			//gPad->Update();
		}
		if (draw_datainfo){
			label_preliminary->Draw();
			label_data_info->Draw();
		}
		gPad->Update();
	}
}

void Tfitresult::Rectify_phase(TGraphErrors* graph_phase){
	if (!graph_phase) return;
	int npoints = graph_phase->GetN();
	double* m = graph_phase->GetX();
	double* phase = graph_phase->GetY();
	double* phase_err = graph_phase->GetEY();
	int istart = 0;
	// first iteration: find the closest point to 0 degree (taking the error into account)
	for (int i = 0; i < npoints; i++){
		if (fabs(phase[i])+phase_err[i] < fabs(phase[istart])+phase_err[istart]){
			istart = i;
		}
	}
	// make the point not jump in phase
	for (int i = istart+1; i < npoints; i++){
		do{
			// define the "direction" to correct phase into
			double phase_diff = phase[i-1]-phase[i];
			int dir = (phase_diff < 0) ? -1 : +1;
			// check that a change in phase of 360° would effect
			// a point to sit nearer to the previous one
			phase_diff = fabs(phase_diff);
			if (phase_diff < 180) break;
			// put the point nearer to the previous one
			phase[i] = phase[i]+dir*360;
			// set also the point of the graph
			graph_phase->SetPoint(i, m[i], phase[i]);
		}while (1); // continue to move points in case of many phase jumps
	}
	// now the other direction
	for (int i = istart-1; i > -1; i--){
		do{
			// define the "direction" to correct phase into
			double phase_diff = phase[i+1]-phase[i];
			int dir = (phase_diff < 0) ? -1 : +1;
			// check that a change in phase of 360° would effect
			// a point to sit nearer to the previous one
			phase_diff = fabs(phase_diff);
			if (phase_diff < 180) break;
			// put the point nearer to the previous one
			phase[i] = phase[i]+dir*360;
			// set also the point of the graph
			graph_phase->SetPoint(i, m[i], phase[i]);
		}while (1); // continue to move points in case of many phase jumps
	}
}

Tfitresult::Tfitresult(
		string unique_name, // needs a unique title
		TTree *loadedfitresult, // needs a pointer to a loaded fit result
		TColor_struct color)
{
	graphcolor 	= color;
	fitresult 	= loadedfitresult;
	name	 	= unique_name;
	// initialize graphs
	all_loglikelihoods = NULL;
	most_likely_likelihoods = NULL;
	all_total_intensities = NULL;
	most_likely_total_intensity = NULL;
	all_evidences = NULL;
	most_likely_evidence = NULL;
	mass_low = -1;
	mass_high = -1;
	//if (fitresult) Scan_Fit_Result();
}



Tfitresult::~Tfitresult(){
	cout << " deleting graphs from " << name << endl;
	delete all_loglikelihoods;
	delete most_likely_likelihoods;
	delete all_total_intensities;
	delete most_likely_total_intensity;
	delete all_evidences;
	delete most_likely_evidence;
	for (Twavegraphmapit it = waves.begin(); it != waves.end(); it++){
		delete it->second.all_intensities;
		delete it->second.most_likely_intensity;
		for (Tphasegraphmapit itphase = it->second.phase.begin(); itphase != it->second.phase.end(); itphase++){
			delete itphase->second.all_coherences;
			delete itphase->second.all_phases;
			delete itphase->second.most_likely_coherence;
			delete itphase->second.most_likely_phase;
		}
	}
	for (Tgraphmapit it = spin_totals.begin(); it != spin_totals.end(); it++){
		delete it->second;
	}
	cout << " done " << endl;
}

TGraphErrors* Tfitresult::Get_Total_Intensity(
					bool most_likely
				){
	TGraphErrors* result = NULL;
	TGraphErrors* _graph = NULL;
	if (most_likely){
		_graph = most_likely_total_intensity;
	} else {
		_graph = all_total_intensities;
	}
	if (_graph){
		result = new TGraphErrors(*_graph); // copy the graph
		stringstream some_name;
		static int some_name_counter(0);
		some_name << "total_intensity_graph_number_" << some_name_counter++;
		result->SetName(some_name.str().c_str());
		result->SetTitle(_graph->GetTitle());
		result->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
		result->GetYaxis()->SetTitle("Total intensity");
	}
	return result;
}

TGraphErrors* Tfitresult::Get_Likelihood(
					bool most_likely
				){
	TGraphErrors* result = NULL;
	TGraphErrors* _graph = NULL;
	if (most_likely){
		_graph = most_likely_likelihoods;
	} else {
		_graph = all_loglikelihoods;
	}
	if (_graph){
		result = new TGraphErrors(*_graph); // copy the graph
		stringstream some_name;
		static int some_name_counter(0);
		some_name << "likelihood_graph_number_" << some_name_counter++;
		result->SetName(some_name.str().c_str());
		result->SetTitle(_graph->GetTitle());
		result->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
		result->GetYaxis()->SetTitle("log likelihood");
	}
	return result;
}

TGraphErrors* Tfitresult::Get_Evidence(
					bool most_likely
				){
	TGraphErrors* result = NULL;
	TGraphErrors* _graph = NULL;
	if (most_likely){
		_graph = most_likely_evidence;
	} else {
		_graph = all_evidences;
	}
	if (_graph){
		result = new TGraphErrors(*_graph); // copy the graph
		stringstream some_name;
		static int some_name_counter(0);
		some_name << "evidence_graph_number_" << some_name_counter++;
		result->SetName(some_name.str().c_str());
		result->SetTitle(_graph->GetTitle());
		result->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
		result->GetYaxis()->SetTitle("evidence");
	}
	return result;
}

TGraphErrors *Tfitresult::Get_Intensity(string wavename, // get an intensity of a given partial wave
bool most_likely){
	Twavegraphmapit it = waves.find(wavename);
	TGraphErrors* result = NULL;
	if (it == waves.end()){
		cout << "Error in Tfitresult::Get_Intensity(): Intensity for " << wavename << " does not exist! " << endl;
		return result;
	} else {
		TGraphErrors* _graph = NULL;
		if (most_likely){
			_graph = it->second.most_likely_intensity;
		} else {
			_graph = it->second.all_intensities;
		}
		if (_graph){
			result = new TGraphErrors(*_graph); // copy the graph
			stringstream some_name;
			static int some_name_counter(0);
			some_name << "intensity_graph_number_" << some_name_counter++;
			result->SetName(some_name.str().c_str());
			result->SetTitle(_graph->GetTitle());
			result->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
			result->GetYaxis()->SetTitle("Intensity");
		}
	}
	return result;
}



TGraphErrors *Tfitresult::Get_Phase(string wavename, // get a phase of a given partial wave
		string anchorwave, // to a corresponding anchor wave
		bool most_likely){
	TGraphErrors* result = NULL;
	if (!Create_Phase_Coherence_graph(anchorwave)) return result;
	Tphasegraphmapit it = waves[wavename].phase.find(anchorwave);
	if (it == waves[wavename].phase.end()){
		cout << "Unexpected Error in Tfitresult::Get_Phase(): aborting!" << endl;
		return result;
	} else {
		TGraphErrors* _graph = NULL;
		if (most_likely){
			_graph = it->second.most_likely_phase;
		} else {
			_graph = it->second.all_phases;
		}
		if (_graph){
			result = new TGraphErrors(*_graph); // copy the graph
			stringstream some_name;
			static int some_name_counter(0);
			some_name << "phase_graph_number_" << some_name_counter++;
			result->SetName(some_name.str().c_str());
			result->SetTitle(_graph->GetTitle());
			result->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
			result->GetYaxis()->SetTitle("Phase Angle [deg]");
		}
	}
	return result;
}



TGraphErrors *Tfitresult::Get_Coherence(
		string wavename, // get a coherence of a given partial wave
		string anchorwave, // to a corresponding anchor wave
		bool most_likely){
	TGraphErrors* result = NULL;
	if (!Create_Phase_Coherence_graph(anchorwave)) return result;
	Tphasegraphmapit it = waves[wavename].phase.find(anchorwave);
	if (it == waves[wavename].phase.end()){
		cout << "Unexpected Error in Tfitresult::Get_Coherence(): aborting!" << endl;
		return result;
	} else {
		TGraphErrors* _graph = NULL;
		if (most_likely){
			_graph = it->second.most_likely_coherence;
		} else {
			_graph = it->second.all_coherences;
		}
		if (_graph){
			result = new TGraphErrors(*_graph); // copy the graph
			stringstream some_name;
			static int some_name_counter(0);
			some_name << "coherence_graph_number_" << some_name_counter++;
			result->SetName(some_name.str().c_str());
			result->SetTitle(_graph->GetTitle());
			result->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
			result->GetYaxis()->SetTitle("Coherence");
		}
	}
	return result;
}

bool Tfitresult::Create_Phase_Coherence_graph(
		string wavename,
		string anchorwave){
	string branchName = "fitResult_v2";
	if (waves.find(wavename) == waves.end()){
		cout << " Error in Tfitresult::Create_Phase_Coherence_graph(): wave " << wavename << " does not exist! " << endl;
		return false;
	}
	if (waves.find(anchorwave) == waves.end()){
		cout << " Error in Tfitresult::Create_Phase_Coherence_graph(): anchor wave " << anchorwave << " does not exist! " << endl;
		return false;
	}
	Tphasegraphmapit it = waves[wavename].phase.find(anchorwave);
	if (it == waves[wavename].phase.end()){
		// create a new graph else do nothing
		Tphasegraph& phasegraph = waves[wavename].phase[anchorwave];
		phasegraph.mass_low = -1;
		phasegraph.mass_high = -1;
		string graphname;
		phasegraph.all_phases 				= new TGraphErrors(0);
		graphname = name+wavename+"vs"+anchorwave+"all_phases";
		phasegraph.all_phases->SetNameTitle(graphname.c_str(), ("#Delta #phi("+wavename+"_vs_"+anchorwave+")").c_str());
		phasegraph.all_phases->SetMarkerColor(graphcolor.rootcolorindex);
		phasegraph.all_phases->SetMarkerStyle(21);
		phasegraph.all_phases->SetMarkerSize(0.5);
		phasegraph.most_likely_phase 		= new TGraphErrors(0);
		graphname = name+wavename+"vs"+anchorwave+"most_likely_phase";
		phasegraph.most_likely_phase->SetNameTitle(graphname.c_str(), ("#Delta #phi("+wavename+"_vs_"+anchorwave+")").c_str());
		phasegraph.most_likely_phase->SetMarkerColor(graphcolor.rootcolorindex);
		phasegraph.most_likely_phase->SetMarkerStyle(21);
		phasegraph.most_likely_phase->SetMarkerSize(0.5);
		phasegraph.all_coherences 			= new TGraphErrors(0);
		graphname = name+wavename+"vs"+anchorwave+"all_coherences";
		phasegraph.all_coherences->SetNameTitle(graphname.c_str(), ("Coherence("+wavename+"_vs_"+anchorwave+")").c_str());
		phasegraph.all_coherences->SetMarkerColor(graphcolor.rootcolorindex);
		phasegraph.all_coherences->SetMarkerStyle(21);
		phasegraph.all_coherences->SetMarkerSize(0.5);
		phasegraph.most_likely_coherence 	= new TGraphErrors(0);
		graphname = name+wavename+"vs"+anchorwave+"most_likely_coherence";
		phasegraph.most_likely_coherence->SetNameTitle(graphname.c_str(), ("Coherence("+wavename+"_vs_"+anchorwave+")").c_str());
		phasegraph.most_likely_coherence->SetMarkerColor(graphcolor.rootcolorindex);
		phasegraph.most_likely_coherence->SetMarkerStyle(21);
		phasegraph.most_likely_coherence->SetMarkerSize(0.5);

		cout << " calculating phases and coherences " << endl;
		map<double, int> massbin_positions;
		map<double, double> smallest_loglikelihood;
		fitResult* massBin = new fitResult();
		fitresult->SetBranchAddress(branchName.c_str(), &massBin);
		int nmassbins = fitresult->GetEntries();
		int different_massbins(0);
		// check all given mass bins
		for (int imassbin = 0; imassbin < nmassbins; imassbin++){
			DrawProgressBar(50, (double)(imassbin+1)/((double)nmassbins));
			fitresult->GetEntry(imassbin);
			double mass = massBin->massBinCenter()/1000.;
			double loglike = massBin->logLikelihood();
			bool found_newmassbin(false);
			if (massbin_positions.find(mass) == massbin_positions.end()){
				found_newmassbin = true;
				massbin_positions[mass] = different_massbins;
				different_massbins++;
			}
			bool found_morelikely(false);
			if (smallest_loglikelihood[mass] > loglike){
				found_morelikely = true;
				smallest_loglikelihood[mass] = loglike;
			}

			int _indexA = -1;
			if (wavename != "")
			_indexA = massBin->waveIndex(wavename);
			int _indexB = -1;
			if (anchorwave != "")
			_indexB = massBin->waveIndex(anchorwave);
			if (_indexA < 0){
				cout << " wave " << wavename << " not found in bin " << mass << endl;
				continue;
			}				
			if (_indexB < 0){
				cout << " anchor wave " << anchorwave << " not found in bin " << mass << endl;
				continue;
			}		
			// both waves were found adjust the valid mass range
			if (phasegraph.mass_low  == -1 || phasegraph.mass_low > mass ) phasegraph.mass_low  = mass;
			if (phasegraph.mass_high == -1 || phasegraph.mass_high < mass) phasegraph.mass_high = mass;

			// retrieve phases and coherences
			double phase        =
				massBin->phase       ((unsigned int) _indexA, (unsigned int) _indexB);
			double phaseErr     =
				massBin->phaseErr    ((unsigned int) _indexA, (unsigned int) _indexB);
			double coherence    =
				massBin->coherence   ((unsigned int) _indexA, (unsigned int) _indexB);
			double coherenceErr =
				massBin->coherenceErr((unsigned int) _indexA, (unsigned int) _indexB);				

			// fill up all phases and coherences
			int _n = phasegraph.all_phases->GetN();
			phasegraph.all_phases->Set(_n+1);
			phasegraph.all_phases->SetPoint(_n, mass, phase);
			phasegraph.all_phases->SetPointError(_n, 0, phaseErr);

			_n = phasegraph.all_coherences->GetN();
			phasegraph.all_coherences->Set(_n+1);
			phasegraph.all_coherences->SetPoint(_n, mass, coherence);
			phasegraph.all_coherences->SetPointError(_n, 0, coherenceErr);
			
			// add a new for the most likely one in case of a new mass bin
			if (found_newmassbin){
				int _n = phasegraph.most_likely_phase->GetN();
				phasegraph.most_likely_phase->Set(_n+1);
				phasegraph.most_likely_phase->SetPoint(_n, mass, phase);
				phasegraph.most_likely_phase->SetPointError(_n, 0, phaseErr);

				_n = phasegraph.most_likely_coherence->GetN();
				phasegraph.most_likely_coherence->Set(_n+1);
				phasegraph.most_likely_coherence->SetPoint(_n, mass, coherence);
				phasegraph.most_likely_coherence->SetPointError(_n, 0, coherenceErr);				
			} else {
				// check if found a more likely one
				if (found_morelikely){
					// search for the corresponding mass bin
					double* _m = phasegraph.most_likely_phase->GetX();
					int 	_n = phasegraph.most_likely_phase->GetN();
					int _pos = -1;
					for (int i = 0; i < _n; i++){
						if (_m[i] == mass){
							_pos = i;
							break;
						}
					} 
					if (_pos < 0){
						cout << " Unexpected error in Tfitresult::Create_Phase_Coherence_graph(): corresponding mass not found! " << endl;
						continue;
					}
					phasegraph.most_likely_phase->SetPoint(_pos, mass, phase);
					phasegraph.most_likely_phase->SetPointError(_pos, 0, phaseErr);

					phasegraph.most_likely_coherence->SetPoint(_pos, mass, coherence);
					phasegraph.most_likely_coherence->SetPointError(_pos, 0, coherenceErr);						
				}
			}
		}
		cout << " done " << endl;
		// now some cosmetics
		phasegraph.most_likely_phase->Sort();
		Rectify_phase(phasegraph.most_likely_phase);
		phasegraph.most_likely_coherence->Sort();

	}
	return true;
}

bool Tfitresult::Create_Phase_Coherence_graph(
		string anchorwave){
	string branchName = "fitResult_v2";
	if (waves.find(anchorwave) == waves.end()){
		cout << " Error in Tfitresult::Create_Phase_Coherence_graph(): anchor wave " << anchorwave << " does not exist! " << endl;
		return false;
	}
	cout << " creating phase and coherence graphs for anchor wave " << anchorwave << endl;
	//cout << " calculating phases and coherences " << endl;
	map<double, int> massbin_positions;
	map<double, double> smallest_loglikelihood;
	fitResult* massBin = new fitResult();
	fitresult->SetBranchAddress(branchName.c_str(), &massBin);
	int nmassbins = fitresult->GetEntries();
	int different_massbins(0);
	map <string, bool> skipphase; // in case a graph exists
	// check all given mass bins
	for (int imassbin = 0; imassbin < nmassbins; imassbin++){
		// after a first iteration searching for existing graphs
		// it is decided if there are some to be filled or not
		if (imassbin == 1 && skipphase.size() == waves.size()){
			bool finish = true;
			for (map <string, bool>::const_iterator itskip = skipphase.begin(); itskip != skipphase.end(); itskip++){
				if (!itskip->second) {
					finish = false; // at least one wave is not created yet
					continue;
				}
			}
			if (finish){
				DrawProgressBar(50, 1.);
				cout << " all phases had been already calculated. " << endl;
				return true;
			}
		}
		DrawProgressBar(50, (double)(imassbin+1)/((double)nmassbins));
		fitresult->GetEntry(imassbin);
		double mass = massBin->massBinCenter()/1000.;
		double loglike = massBin->logLikelihood();
		bool found_newmassbin(false);
		if (massbin_positions.find(mass) == massbin_positions.end()){
			found_newmassbin = true;
			massbin_positions[mass] = different_massbins;
			different_massbins++;
		}
		bool found_morelikely(false);
		if (smallest_loglikelihood[mass] > loglike){
			found_morelikely = true;
			smallest_loglikelihood[mass] = loglike;
		}

		for (Twavegraphmapit itwave = waves.begin(); itwave != waves.end(); itwave++){
			// create a graph if not existent
			Tphasegraphmapit it = itwave->second.phase.find(anchorwave);
			if (it == itwave->second.phase.end()){
				// create a new graph
				Tphasegraph& phasegraph = itwave->second.phase[anchorwave];
				phasegraph.mass_low = -1;
				phasegraph.mass_high = -1;
				string graphname;
				phasegraph.all_phases 				= new TGraphErrors(0);
				graphname = name+itwave->first+"vs"+anchorwave+"all_phases";
				phasegraph.all_phases->SetNameTitle(graphname.c_str(), ("#Delta #phi("+itwave->first+"_vs_"+anchorwave+")").c_str());
				phasegraph.all_phases->SetMarkerColor(graphcolor.rootcolorindex);
				phasegraph.all_phases->SetMarkerStyle(21);
				phasegraph.all_phases->SetMarkerSize(0.5);
				phasegraph.most_likely_phase 		= new TGraphErrors(0);
				graphname = name+itwave->first+"vs"+anchorwave+"most_likely_phase";
				phasegraph.most_likely_phase->SetNameTitle(graphname.c_str(), ("#Delta #phi("+itwave->first+"_vs_"+anchorwave+")").c_str());
				phasegraph.most_likely_phase->SetMarkerColor(graphcolor.rootcolorindex);
				phasegraph.most_likely_phase->SetMarkerStyle(21);
				phasegraph.most_likely_phase->SetMarkerSize(0.5);
				phasegraph.all_coherences 			= new TGraphErrors(0);
				graphname = name+itwave->first+"vs"+anchorwave+"all_coherences";
				phasegraph.all_coherences->SetNameTitle(graphname.c_str(), ("Coherence("+itwave->first+"_vs_"+anchorwave+")").c_str());
				phasegraph.all_coherences->SetMarkerColor(graphcolor.rootcolorindex);
				phasegraph.all_coherences->SetMarkerStyle(21);
				phasegraph.all_coherences->SetMarkerSize(0.5);
				phasegraph.most_likely_coherence 	= new TGraphErrors(0);
				graphname = name+itwave->first+"vs"+anchorwave+"most_likely_coherence";
				phasegraph.most_likely_coherence->SetNameTitle(graphname.c_str(), ("Coherence("+itwave->first+"_vs_"+anchorwave+")").c_str());
				phasegraph.most_likely_coherence->SetMarkerColor(graphcolor.rootcolorindex);
				phasegraph.most_likely_coherence->SetMarkerStyle(21);
				phasegraph.most_likely_coherence->SetMarkerSize(0.5);
			} else {
				// checked if a graph exists already
				if (imassbin == 0)
					skipphase[itwave->first] = true;
				else
					skipphase[itwave->first] = false;
			}
			if (skipphase[itwave->first]) continue;
			Tphasegraph& phasegraph = itwave->second.phase[anchorwave];
			int _indexA = -1;
			//if (wavename != "")
			_indexA = massBin->waveIndex(itwave->first);
			int _indexB = -1;
			if (anchorwave != "")
			_indexB = massBin->waveIndex(anchorwave);
			if (_indexA < 0){
				cout << " wave " << itwave->first << " not found in bin " << mass << endl;
				continue;
			}
			if (_indexB < 0){
				if (itwave == waves.begin())
					cout << " anchor wave " << anchorwave << " not found in bin " << mass << endl;
				continue;
			}
			// both waves were found adjust the valid mass range
			if (phasegraph.mass_low  == -1 || phasegraph.mass_low > mass ) phasegraph.mass_low  = mass;
			if (phasegraph.mass_high == -1 || phasegraph.mass_high < mass) phasegraph.mass_high = mass;

			// retrieve phases and coherences
			double phase        =
				massBin->phase       ((unsigned int) _indexA, (unsigned int) _indexB);
			double phaseErr     =
				massBin->phaseErr    ((unsigned int) _indexA, (unsigned int) _indexB);
			double coherence    =
				massBin->coherence   ((unsigned int) _indexA, (unsigned int) _indexB);
			double coherenceErr =
				massBin->coherenceErr((unsigned int) _indexA, (unsigned int) _indexB);

			// fill up all phases and coherences
			int _n = phasegraph.all_phases->GetN();
			phasegraph.all_phases->Set(_n+1);
			phasegraph.all_phases->SetPoint(_n, mass, phase);
			phasegraph.all_phases->SetPointError(_n, 0, phaseErr);

			_n = phasegraph.all_coherences->GetN();
			phasegraph.all_coherences->Set(_n+1);
			phasegraph.all_coherences->SetPoint(_n, mass, coherence);
			phasegraph.all_coherences->SetPointError(_n, 0, coherenceErr);

			// add a new for the most likely one in case of a new mass bin
			if (found_newmassbin){
				int _n = phasegraph.most_likely_phase->GetN();
				phasegraph.most_likely_phase->Set(_n+1);
				phasegraph.most_likely_phase->SetPoint(_n, mass, phase);
				phasegraph.most_likely_phase->SetPointError(_n, 0, phaseErr);

				_n = phasegraph.most_likely_coherence->GetN();
				phasegraph.most_likely_coherence->Set(_n+1);
				phasegraph.most_likely_coherence->SetPoint(_n, mass, coherence);
				phasegraph.most_likely_coherence->SetPointError(_n, 0, coherenceErr);
			} else {
				// check if found a more likely one
				if (found_morelikely){
					// search for the corresponding mass bin
					double* _m = phasegraph.most_likely_phase->GetX();
					int 	_n = phasegraph.most_likely_phase->GetN();
					int _pos = -1;
					for (int i = 0; i < _n; i++){
						if (_m[i] == mass){
							_pos = i;
							break;
						}
					}
					if (_pos < 0){
						cout << " Unexpected error in Tfitresult::Create_Phase_Coherence_graph(): corresponding mass not found! " << endl;
						continue;
					}
					phasegraph.most_likely_phase->SetPoint(_pos, mass, phase);
					phasegraph.most_likely_phase->SetPointError(_pos, 0, phaseErr);

					phasegraph.most_likely_coherence->SetPoint(_pos, mass, coherence);
					phasegraph.most_likely_coherence->SetPointError(_pos, 0, coherenceErr);
				}
			}
		}

	}
	// do some cosmetics
	for (Twavegraphmapit itwave = waves.begin(); itwave != waves.end(); itwave++){
		Tphasegraphmapit it = itwave->second.phase.find(anchorwave);
		if (it == itwave->second.phase.end()){
			cout << " Uups: should not happen! " << endl;
			return false;
		}
		// now some cosmetics
		it->second.most_likely_phase->Sort();
		Rectify_phase(it->second.most_likely_phase);
		it->second.most_likely_coherence->Sort();
	}

	cout << " done " << endl;
	return true;
}

// most likely spin totals are returned
// the key is coded as JP(C)
TGraphErrors* Tfitresult::Get_Spin_Total(string jpc){
	Tgraphmapit it = spin_totals.find(jpc);
	TGraphErrors* result = NULL;
	if (it == spin_totals.end()){
		cout << "Error in Tfitresult::Get_Spin_Total(): Spin total " << jpc << " does not exist! " << endl;
		return result;
	} else {
		TGraphErrors* _graph = NULL;
		if (0){//most_likely){
			//_graph = it->second.most_likely_intensity;
		} else {
			_graph = it->second;
		}
		if (_graph){
			result = new TGraphErrors(*_graph); // copy the graph
			stringstream some_name;
			static int some_name_counter(0);
			some_name << "spin_total_graph_number_" << some_name_counter++;
			result->SetName(some_name.str().c_str());
			result->SetTitle(_graph->GetTitle());
			result->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
			result->GetYaxis()->SetTitle("Intensity");
		}
	}
	return result;
}

// get a vector with the available spin totals
vector<string>& Tfitresult::Get_available_Spin_Totals(){
	vector<string>* result = new vector<string>();
	for (Tgraphmapit it = spin_totals.begin(); it != spin_totals.end(); it++){
		result->push_back(it->first);
	}
	return *result;
}

vector<string>& Tfitresult::Get_available_waves(){
	vector<string>* result = new vector<string>();
	if (waves.size() == 0){
		cout << " Warning in Tfitresult::Get_available_waves(): No available waves found! Did you scan already?" << endl;
	}
	for (Twavegraphmapit it = waves.begin(); it != waves.end(); it++){
		result->push_back(it->first);
	}
	return *result;
}

vector<string>& Tfitresult::Scan_Fit_Result(string branchName){
	vector<string>* result = new vector<string>();
	if (fitresult <= 0) return *result;
	// in case the stored fit result is requested to be checked for waves
	// and the waves were already loaded, only the available wave list is returned
	if (waves.size() > 0){
		for (Twavegraphmapit it = waves.begin(); it != waves.end(); it++){
			result->push_back(it->first);
		}
		return *result;
	}
	// else scan not only for existing waves but also create already intensity graphs
	cout << " scanning " << fitresult->GetName() << " for waves " << endl;

	fitResult* massBin = new fitResult();
	fitresult->SetBranchAddress(branchName.c_str(), &massBin);
	// store the position of the first mass bin for most likely solution in the graphs
	map<double, int> massbin_positions;
	// store the smallest log likelihood you find here
	map<double, double> smallest_loglikelihood;
	if (all_loglikelihoods || most_likely_likelihoods || all_total_intensities || most_likely_total_intensity){
		cout << " Error in Tfitresult::Scan_Fit_Result(): some graphs are already initialized. Replacing! " << endl;
	}
	string graphname;
	all_loglikelihoods			= new TGraphErrors(0);
	graphname = name+"_all_loglikelihoods";
	all_loglikelihoods->SetNameTitle(graphname.c_str(), "all log likelihoods");
	all_loglikelihoods->SetMarkerColor(graphcolor.rootcolorindex);
	all_loglikelihoods->SetMarkerStyle(21);
	all_loglikelihoods->SetMarkerSize(0.5);
	most_likely_likelihoods 	= new TGraphErrors(0);
	graphname = name+"most_likely_likelihoods";
	most_likely_likelihoods->SetNameTitle(graphname.c_str(), "best log likelihoods");
	most_likely_likelihoods->SetMarkerColor(graphcolor.rootcolorindex);
	most_likely_likelihoods->SetMarkerStyle(21);
	most_likely_likelihoods->SetMarkerSize(0.5);
	all_total_intensities 		= new TGraphErrors(0);
	graphname = name+"all_total_intensities";
	all_total_intensities->SetNameTitle(graphname.c_str(), "all total intensities");
	all_total_intensities->SetMarkerColor(graphcolor.rootcolorindex);
	all_total_intensities->SetMarkerStyle(21);
	all_total_intensities->SetMarkerSize(0.5);
	most_likely_total_intensity = new TGraphErrors(0);
	graphname = name+"most_likely_total_intensity";
	most_likely_total_intensity->SetNameTitle(graphname.c_str(), "most likely total intensity");
	most_likely_total_intensity->SetMarkerColor(graphcolor.rootcolorindex);
	most_likely_total_intensity->SetMarkerStyle(21);
	most_likely_total_intensity->SetMarkerSize(0.5);
	all_evidences = new TGraphErrors(0);
	graphname = name+"all_evidences";
	all_evidences->SetNameTitle(graphname.c_str(), "all evidences");
	all_evidences->SetMarkerColor(graphcolor.rootcolorindex);
	all_evidences->SetMarkerStyle(21);
	all_evidences->SetMarkerSize(0.5);
	most_likely_evidence = new TGraphErrors(0);
	graphname = name+"most_likely_evidence";
	most_likely_evidence->SetNameTitle(graphname.c_str(), "most likely evidence");
	most_likely_evidence->SetMarkerColor(graphcolor.rootcolorindex);
	most_likely_evidence->SetMarkerStyle(21);
	most_likely_evidence->SetMarkerSize(0.5);
	int nbinelements(0);
	int different_massbins(0);
	int nmassbins = fitresult->GetEntries();
	// check all given mass bins
	for (int imassbin = 0; imassbin < nmassbins; imassbin++){
		DrawProgressBar(50, (double)(imassbin+1)/((double)nmassbins));
		fitresult->GetEntry(imassbin);
		//double evidence = massBin->evidence();
		double mass = massBin->massBinCenter()/1000.;
		if (mass_low  == -1 || mass_low > mass ) mass_low  = mass;
		if (mass_high == -1 || mass_high < mass) mass_high = mass;
		double loglike = massBin->logLikelihood();
		double evidence= massBin->evidence();
		int _pos = all_loglikelihoods->GetN();
		all_loglikelihoods->Set(_pos+1);
		all_loglikelihoods->SetPoint(_pos, mass, loglike);
		all_evidences->Set(_pos+1);
		all_evidences->SetPoint(_pos, mass, evidence);
		bool found_newmassbin(false);
		if (massbin_positions.find(mass) == massbin_positions.end()){
			found_newmassbin = true;
			massbin_positions[mass] = different_massbins;
			int _pos = most_likely_likelihoods->GetN();
			// just a check
			if (_pos != different_massbins){
				cout << " Unexpected error in Tfitresult::Scan_Fit_Result() " << endl;
			}
			different_massbins++;
			// add the most likely likelihood
			most_likely_likelihoods->Set(_pos+1); // will be filled later
			most_likely_total_intensity->Set(_pos+1); // will be filled later
			most_likely_evidence->Set(_pos+1);
		}
		bool found_morelikely(false);
		if (smallest_loglikelihood[mass] > loglike){
			found_morelikely = true;
			smallest_loglikelihood[mass] = loglike;
			int _pos = massbin_positions[mass];
			// just a check
			if (_pos < 0 || _pos >= most_likely_likelihoods->GetN()){
				cout << " Unexpected error in Tfitresult::Scan_Fit_Result() " << endl;
				cout << " aborting! " << endl;
				result->clear();
				return *result;
			} else {
				most_likely_likelihoods->SetPoint(_pos, mass, loglike);
				most_likely_evidence->SetPoint(_pos, mass, evidence);
			}
		}
		nbinelements++;
		double totalintensity        = massBin->intensity(".*");//(0);
		double totalintensity_err    = massBin->intensityErr(".*");
		// now go through all available waves in this bin
		for (unsigned iwave = 0; iwave < massBin->waveNames().size(); iwave++){
			// does this wave name exist already?
			string wavename = massBin->waveNames()[iwave];
			double intensity = massBin->intensity(wavename.c_str());
			double intensity_err = massBin->intensityErr(wavename.c_str());
			//totalintensity += intensity;
			//totalintensity_err = sqrt(totalintensity_err*totalintensity_err + intensity_err*intensity_err);
			if (waves.find(wavename) != waves.end()){ // it exists already
				Twavegraph& wavegraph = waves.find(wavename)->second;
				if (wavegraph.mass_low  > mass) wavegraph.mass_low  = mass;
				if (wavegraph.mass_high < mass) wavegraph.mass_high = mass;
				int _pos = wavegraph.all_intensities->GetN();
				wavegraph.all_intensities->Set(_pos+1);
				wavegraph.all_intensities->SetPoint(_pos, mass, intensity);
				wavegraph.all_intensities->SetPointError(_pos, 0, intensity_err);
				if (found_newmassbin){
					_pos = wavegraph.most_likely_intensity->GetN();
					wavegraph.most_likely_intensity->Set(_pos+1);
					wavegraph.most_likely_intensity->SetPoint(_pos, mass, intensity);
					wavegraph.most_likely_intensity->SetPointError(_pos, 0, intensity_err);
				} else {
					// is it a more likely solution?
					if (found_morelikely){
						// search for the corresponding mass bin
						double* _m = wavegraph.most_likely_intensity->GetX();
						int _n = wavegraph.most_likely_intensity->GetN();
						_pos = -1;
						for (int i = 0; i < _n; i++){
							if (_m[i] == mass){
								_pos = i;
								break;
							}
						}
						// check consistency
						if (_pos != -1){
							// everything fine
							wavegraph.most_likely_intensity->SetPoint(_pos, mass, intensity);
							wavegraph.most_likely_intensity->SetPointError(_pos, 0, intensity_err);
						} else {
							cout << " Warning: wave " << wavename << " does not seem to be present in all fit results of " << name << " in mass bin " << mass << endl;
							wavegraph.most_likely_intensity->Set(_n+1);
							wavegraph.most_likely_intensity->SetPoint(_n, mass, intensity);
							wavegraph.most_likely_intensity->SetPointError(_n, 0, intensity_err);
						}
					}
				} // else of if found new mass bin
			} else { // this wave does not exist yet
				Twavegraph& wavegraph = waves[wavename];
				result->push_back(wavename); // store the new wave name
				// create new graphs
				wavegraph.all_intensities = new TGraphErrors(1);
				string _graphname = name+"_"+wavename+"_all_intensities";
				GetJPCMreflISO1lsISO2(wavename, wavegraph.J, wavegraph.P, wavegraph.C, wavegraph.M, wavegraph.e, wavegraph.iso1, wavegraph.iso2, wavegraph.l, wavegraph.s);
				wavegraph.all_intensities->SetName(_graphname.c_str());
				wavegraph.all_intensities->SetTitle(wavename.c_str());
				wavegraph.all_intensities->SetMarkerColor(graphcolor.rootcolorindex);
				wavegraph.all_intensities->SetMarkerStyle(21);
				wavegraph.all_intensities->SetMarkerSize(0.5);
				wavegraph.all_intensities->SetPoint(0, mass, intensity);
				wavegraph.all_intensities->SetPointError(0, 0, intensity_err);
				wavegraph.most_likely_intensity = new TGraphErrors(1);
				_graphname = name+"_"+wavename+"_most_likely_intensity";
				wavegraph.most_likely_intensity->SetName(_graphname.c_str());
				wavegraph.most_likely_intensity->SetTitle(wavename.c_str());
				wavegraph.most_likely_intensity->SetMarkerColor(graphcolor.rootcolorindex);
				wavegraph.most_likely_intensity->SetMarkerStyle(21);
				wavegraph.most_likely_intensity->SetMarkerSize(0.5);
				wavegraph.most_likely_intensity->SetPoint(0, mass, intensity);
				wavegraph.most_likely_intensity->SetPointError(0, 0, intensity_err);
				wavegraph.mass_low  = mass;
				wavegraph.mass_high = mass;
				// phase graphs are created by later request
				// this would take too long to plot each graph against each

				// see if it is a new spin total
				stringstream _jpc;
				if (wavename == "flat") _jpc << "flat";
				else {
					_jpc << wavegraph.J;
					if (wavegraph.P < 0) _jpc << "-"; else _jpc << "+";
					if (wavegraph.C < 0) _jpc << "-"; else _jpc << "+";
				}
				// does this JPC combination exist already?
				if (spin_totals.find(_jpc.str())==spin_totals.end()){
					// create a new graph
					TGraphErrors* &_graph = spin_totals[_jpc.str()];
					_graph = new TGraphErrors(0);//*wavegraph.most_likely_intensity);
					stringstream _graphname;
					static int graphcounter(0);
					_graphname << "total_intensity_graph_"<< _jpc.str() << graphcounter++;
					_graph->SetNameTitle(_graphname.str().c_str(), _jpc.str().c_str());
					_graph->SetMarkerColor(graphcolor.rootcolorindex);
					_graph->SetMarkerStyle(21);
					_graph->SetMarkerSize(0.5);
				}
			}
		}
		// add the total intensities
		if (found_morelikely){
			int _pos = massbin_positions[mass];
			// just a check
			if (_pos < 0 || _pos >= most_likely_total_intensity->GetN()){
				cout << " Unexpected error in Tfitresult::Scan_Fit_Result() " << endl;
				cout << " aborting! " << endl;
				result->clear();
				return *result;
			} else {
				most_likely_total_intensity->SetPoint(_pos, mass, totalintensity);
				most_likely_total_intensity->SetPointError(_pos, 0, totalintensity_err);
				// error not treatet yet
			}
		}
		_pos = all_total_intensities->GetN();
		all_total_intensities->Set(_pos+1);
		all_total_intensities->SetPoint(_pos, mass, totalintensity);
		all_total_intensities->SetPointError(_pos, 0, totalintensity_err);

		// calculate the spin totals
		for (map<string, TGraphErrors*>::iterator it = spin_totals.begin(); it != spin_totals.end(); it++){
			TGraphErrors* &_graph = spin_totals[it->first];
			if (!_graph){
				cout << " Hä? Spin total graph does not exist? should not happen!" << endl;
				continue;
			}
			// replace the + by \+
			stringstream cutstring;
			for (unsigned int i = 0; i < it->first.size(); i++){
				if (it->first[i] == '+')
					cutstring << '\\';
				cutstring << it->first[i];
			}

			/*
			unsigned int spos = 0;
			while (1){
				spos = cutstring.find("+", spos); // start the search from the previous position
				if (spos == string::npos){
					break;
				}
				cutstring.insert(spos,1, '\\');
				spos++; spos++;
				//cout << cutstring << endl;
			}*/
			//cout << cutstring << endl;
			double intensity = massBin->intensity(cutstring.str().c_str());
			double intensityErr = massBin->intensityErr(cutstring.str().c_str());

			// find the corresponding mass bin (if it exists)
			int _pos = -1;
			for (int j = 0; j < _graph->GetN(); j++){
				if (_graph->GetX()[j] == mass){
					_pos = j;
					break;
				}
			}
			if (_pos < 0){
				// create a new point
				_pos = _graph->GetN();
				_graph->Set(_pos+1);
			}
			// fill if no corresponding mass bin found or a more likely solution is given
			if (_pos == _graph->GetN()-1 ||  found_morelikely){
				_graph->SetPoint(_pos, mass, intensity);
				_graph->SetPointError(_pos, 0,  intensityErr);
			}
		}
	}
	//Calculate_Spin_Totals();
	return *result;
}

bool Tfitresult::Write_to_file(string filename){
	cout << " storing available graphs in " << filename << endl;
	bool result = true;
	if (waves.size() == 0){
		cout << " Warning in Tfitresult::Write_to_file: No waves to store " << endl;
		return false;
	}
	TFile file(filename.c_str(), "Recreate");
	if (!file.IsZombie()){
		TTree* tree_wavelist = new TTree("tree_wavelist", "wavelist");
		char* wavename = new char[1000];
		tree_wavelist->Branch("wavename", wavename, "wavename[1000]/B");
		// first iteration: store a list of available waves
		for (Twavegraphmapit it = waves.begin(); it != waves.end(); it++){
			strcpy(wavename,it->first.c_str());
			cout << " Filling " << it->first.c_str() << " as " << wavename << endl;
			tree_wavelist->Fill();
		}
		tree_wavelist->Write();
		// second iteration: store the available waves and their properties
		for (Twavegraphmapit it = waves.begin(); it != waves.end(); it++){
			cout << " storing properties of " << it->first << endl;
			// create folder for each wave
			file.cd();
			file.mkdir(it->first.c_str());
			file.cd(it->first.c_str());
			// the graph it self
			{
				TGraphErrors _graph(*it->second.all_intensities);
				_graph.SetName("all_intensities");
				_graph.Write();
			}
			{
				TGraphErrors _graph(*it->second.most_likely_intensity);
				_graph.SetName("most_likely_intensity");
				_graph.Write();
			}

			TParameter<int> _J("J",it->second.J);
			_J.Write();
			TParameter<int> _P("P",it->second.P);
			_P.Write();
			TParameter<int> _C("C",it->second.C);
			_C.Write();
			TParameter<int> _M("M",it->second.M);
			_M.Write();
			TParameter<int> _e("e",it->second.e);
			_e.Write();
			TParameter<int> _l("l",it->second.l);
			_l.Write();
			TParameter<int> _s("s",it->second.s);
			_s.Write();
			TParameter<int> _mass_low("mass_low",it->second.mass_low);
			_mass_low.Write();
			TParameter<int> _mass_high("mass_high",it->second.mass_high);
			_mass_high.Write();

			TTree* tree_isobars = new TTree("tree_isobars", "isobars");
			char* _isoname = new char[1000];
			tree_isobars->Branch("isoname", _isoname, "isoname[1000]/B");
			strcpy(_isoname,it->second.iso1.c_str());
			tree_isobars->Fill();
			strcpy(_isoname,it->second.iso2.c_str());
			tree_isobars->Fill();
			tree_isobars->Write();
			delete[] _isoname;

			for (Tphasegraphmapit itphase = it->second.phase.begin(); itphase != it->second.phase.end(); itphase++){
				file.cd(it->first.c_str());
				gDirectory->mkdir(itphase->first.c_str());
				file.cd((it->first+"/"+itphase->first).c_str());
				TParameter<int> _mass_low("mass_low",itphase->second.mass_low);
				_mass_low.Write();
				TParameter<int> _mass_high("mass_high",itphase->second.mass_high);
				_mass_high.Write();
				{
					TGraphErrors _graph(*itphase->second.all_coherences);
					_graph.SetName("all_coherences");
					_graph.Write();
				}
				{
					TGraphErrors _graph(*itphase->second.all_phases);
					_graph.SetName("all_phases");
					_graph.Write();
				}
				{
					TGraphErrors _graph(*itphase->second.most_likely_coherence);
					_graph.SetName("most_likely_coherence");
					_graph.Write();
				}
				{
					TGraphErrors _graph(*itphase->second.most_likely_phase);
					_graph.SetName("most_likely_phase");
					_graph.Write();
				}
			}
		}
		// one more step: store global graphs
		file.cd();
		{
			TGraphErrors _graph(*all_loglikelihoods);
			_graph.SetName("all_loglikelihoods");
			_graph.Write();
		}
		{
			TGraphErrors _graph(*most_likely_likelihoods);
			_graph.SetName("most_likely_likelihoods");
			_graph.Write();
		}
		{
			TGraphErrors _graph(*all_total_intensities);
			_graph.SetName("all_total_intensities");
			_graph.Write();
		}
		{
			TGraphErrors _graph(*most_likely_total_intensity);
			_graph.SetName("most_likely_total_intensity");
			_graph.Write();
		}
		{
			TGraphErrors _graph(*all_evidences);
			_graph.SetName("all_evidences");
			_graph.Write();
		}
		{
			TGraphErrors _graph(*most_likely_evidence);
			_graph.SetName("most_likely_evidence");
			_graph.Write();
		}
		TParameter<int> _mass_low("mass_low",mass_low);
		_mass_low.Write();
		TParameter<int> _mass_high("mass_high",mass_high);
		_mass_high.Write();
		// last step: store spin totals
		TTree* tree_spintotallist = new TTree("tree_spintotallist", "spintotals");
		char* spintotalname = new char[1000];
		tree_spintotallist->Branch("spintotalname", spintotalname, "spintotalname[1000]/B");
		for (Tgraphmapit it = spin_totals.begin(); it != spin_totals.end(); it++){
			file.cd();
			strcpy(spintotalname,it->first.c_str());
			cout << " Filling " << it->first.c_str() << " as " << spintotalname << endl;
			tree_spintotallist->Fill();
			file.mkdir(it->first.c_str());
			file.cd(it->first.c_str());
			{
				TGraphErrors _graph(*it->second);
				_graph.SetName("spin_total");
				_graph.Write();
			}
		}
		file.cd();
		tree_spintotallist->Write();
		delete[] spintotalname;

		file.Write();
		file.Close();
		delete[] wavename;
	} else {
		cout << " Error in Tfitresult::Write_to_file: could not write to file: " << filename << endl;
		result = false;
	}
	file.Close();
	return result;
}

bool Tfitresult::Read_from_file(string filename){
	cout << " reading available graphs from " << filename << endl;
	bool result = true;
	if (waves.size() != 0){
		cout << " Warning in Tfitresult::Read_from_file: class is already initialized: cannot read from file! " << endl;
		return false;
	}
	TFile file(filename.c_str(), "Read");
	if (!file.IsZombie()){
		TTree* tree_wavelist = (TTree*) file.Get("tree_wavelist");
		if (!tree_wavelist){
			cout << " Error in Tfitresult::Read_from_file: could not find a wave list in: " << filename << endl;
			return false;
		}
		char* wavename = new char[1000];
		tree_wavelist->SetBranchAddress("wavename", wavename);
		// first: read all available waves and the corresponding graphs and properties
		for (int iwave = 0; iwave < tree_wavelist->GetEntries(); iwave++){
			tree_wavelist->GetEntry(iwave);
			string wavenamestring(wavename);
			Twavegraphmapit it = waves.find(wavename);
			if (it == waves.end()){
				// create a new graph else do nothing
				Twavegraph& wavegraph = waves[wavename];
				cout << " reading properties of " << wavenamestring << endl;
				// read the graph it self
				wavegraph.all_intensities = (TGraphErrors*) file.Get((wavenamestring+"/all_intensities").c_str());
				SetGraphProperties(wavegraph.all_intensities, wavenamestring+"all_intensities");
				wavegraph.most_likely_intensity = (TGraphErrors*) file.Get((wavenamestring+"/most_likely_intensity").c_str());
				SetGraphProperties(wavegraph.most_likely_intensity, wavenamestring+"most_likely_intensity");
				/*				
				if (!wavegraph.all_intensities){
					cout << " Error: could not find " << wavenamestring+"/all_intensities" << endl;
					result = false;
				} else {
					wavegraph.all_intensities->SetName((wavenamestring+"all_intensities").c_str());
					//wavegraph.all_intensities->SetDirectory(0);
				}
				wavegraph.most_likely_intensity = (TGraphErrors*) file.Get((wavenamestring+"/most_likely_intensity").c_str());
				if (!wavegraph.most_likely_intensity){
					cout << " Error: could not find " << wavenamestring+"/most_likely_intensity" << endl;
					result = false;
				} else {
					wavegraph.most_likely_intensity->SetName((wavenamestring+"most_likely_intensity").c_str());
					//wavegraph.most_likely_intensity->SetDirectory(0);
				}
				*/

				TParameter<int>* _J = (TParameter<int>*)file.Get((wavenamestring+"/J").c_str());
				if (_J){
					wavegraph.J = _J->GetVal();
				} else {
					cout << " Error: could not find " << wavenamestring+"/J" << endl;;
					result = false;
				}
				TParameter<int>* _P = (TParameter<int>*)file.Get((wavenamestring+"/P").c_str());
				if (_P){
					wavegraph.P = _P->GetVal();
				} else {
					cout << " Error: could not find " << wavenamestring+"/P" << endl;
					result = false;
				}
				TParameter<int>* _C = (TParameter<int>*)file.Get((wavenamestring+"/C").c_str());
				if (_C){
					wavegraph.C = _C->GetVal();
				} else {
					cout << " Error: could not find " << wavenamestring+"/C" << endl;
					result = false;
				}
				TParameter<int>* _M = (TParameter<int>*)file.Get((wavenamestring+"/M").c_str());
				if (_M){
					wavegraph.M = _M->GetVal();
				} else {
					cout << " Error: could not find " << wavenamestring+"/M" << endl;
					result = false;
				}
				TParameter<int>* _e = (TParameter<int>*)file.Get((wavenamestring+"/e").c_str());
				if (_e){
					wavegraph.e = _e->GetVal();
				} else {
					cout << " Error: could not find " << wavenamestring+"/e" << endl;
					result = false;
				}
				TParameter<int>* _s = (TParameter<int>*)file.Get((wavenamestring+"/s").c_str());
				if (_s){
					wavegraph.s = _s->GetVal();
				} else {
					cout << " Error: could not find " << wavenamestring+"/s" << endl;
					result = false;
				}
				TParameter<int>* _l = (TParameter<int>*)file.Get((wavenamestring+"/l").c_str());
				if (_l){
					wavegraph.l = _l->GetVal();
				} else {
					cout << " Error: could not find " << wavenamestring+"/l" << endl;
					result = false;
				}
				TParameter<int>* _mass_low = (TParameter<int>*)file.Get((wavenamestring+"/mass_low").c_str());
				if (_mass_low){
					wavegraph.mass_low = _mass_low->GetVal();
				} else {
					cout << " Error: could not find " << wavenamestring+"/mass_low" << endl;
					result = false;
				}
				TParameter<int>* _mass_high = (TParameter<int>*)file.Get((wavenamestring+"/mass_high").c_str());
				if (_mass_high){
					wavegraph.mass_high = _mass_high->GetVal();
				} else {
					cout << " Error: could not find " << wavenamestring+"/mass_high" << endl;
					result = false;
				}

				TTree* tree_isobars = (TTree*)file.Get((wavenamestring+"/tree_isobars").c_str());
				if (tree_isobars){
					char* _isoname = new char[1000];
					tree_isobars->SetBranchAddress("isoname", _isoname);
					tree_isobars->GetEntry(0);
					wavegraph.iso1 = _isoname;
					tree_isobars->GetEntry(1);
					wavegraph.iso2 = _isoname;
					delete[] _isoname;
				} else {
					cout << " Error: could not find " << wavenamestring+"/tree_isobars" << endl;
					result = false;
				}
				//testlist[wavenamestring] = wavegraph;
				//cout << " retrieving " << wavename << " as " << wavenamestring << endl;
				//waves[]=wavegraph;
				// finally the phase differences
				for (int iphase = 0; iphase < tree_wavelist->GetEntries(); iphase++){
					tree_wavelist->GetEntry(iphase);
					string phasenamestring(wavename);
					TParameter<int>* _mass_low = (TParameter<int>*)file.Get((wavenamestring+"/"+phasenamestring+"/mass_low").c_str());
					TParameter<int>* _mass_high = (TParameter<int>*)file.Get((wavenamestring+"/"+phasenamestring+"/mass_high").c_str());
					if (!_mass_low || !_mass_high){
						continue;
					} 
					Tphasegraphmapit itphase = wavegraph.phase.find(phasenamestring);
					if (itphase == wavegraph.phase.end()){
						//if (!file.Get(it->first+"/"+phasenamestring+"/all_coherences").c_str())) continue; // no phase graph stored here
						cout << " reading phase vs " << phasenamestring << endl;
						Tphasegraph& phasegraph = wavegraph.phase[phasenamestring];
						phasegraph.mass_low = _mass_low->GetVal();
						phasegraph.mass_high = _mass_high->GetVal();
						phasegraph.all_coherences = (TGraphErrors*) file.Get((wavenamestring+"/"+phasenamestring+"/all_coherences").c_str());
						SetGraphProperties(phasegraph.all_coherences, wavenamestring+"/"+phasenamestring+"/all_coherences");
						phasegraph.all_phases = (TGraphErrors*) file.Get((wavenamestring+"/"+phasenamestring+"/all_phases").c_str());
						SetGraphProperties(phasegraph.all_phases, wavenamestring+"/"+phasenamestring+"/all_phases");
						phasegraph.most_likely_coherence = (TGraphErrors*) file.Get((wavenamestring+"/"+phasenamestring+"/most_likely_coherence").c_str());
						SetGraphProperties(phasegraph.most_likely_coherence, wavenamestring+"/"+phasenamestring+"/most_likely_coherence");
						phasegraph.most_likely_phase = (TGraphErrors*) file.Get((wavenamestring+"/"+phasenamestring+"/most_likely_phase").c_str());
						SetGraphProperties(phasegraph.most_likely_phase, wavenamestring+"/"+phasenamestring+"/most_likely_phase");
					} else {
						cout << " Error: phase graphs " << phasenamestring << " are already existing! " << endl;
						result = false;
					}
				}
			} else {
				cout << " Error: wave " << wavenamestring << " does already exist! " << endl;
				result = false;
			}
		}
		file.cd();
		all_loglikelihoods = (TGraphErrors*) file.Get("all_loglikelihoods");
		SetGraphProperties(all_loglikelihoods, "all_loglikelihoods");
		most_likely_likelihoods = (TGraphErrors*) file.Get("most_likely_likelihoods");
		SetGraphProperties(most_likely_likelihoods, "most_likely_likelihoods");
		all_total_intensities = (TGraphErrors*) file.Get("all_total_intensities");
		SetGraphProperties(all_total_intensities, "all_total_intensities");
		most_likely_total_intensity = (TGraphErrors*) file.Get("most_likely_total_intensity");
		SetGraphProperties(most_likely_total_intensity, "most_likely_total_intensity");
		all_evidences = (TGraphErrors*) file.Get("all_evidences");
		SetGraphProperties(all_evidences, "all_evidences");
		most_likely_evidence = (TGraphErrors*) file.Get("most_likely_evidence");
		SetGraphProperties(most_likely_evidence, "most_likely_evidence");
		TParameter<int>* _mass_low = (TParameter<int>*)file.Get("mass_low");
		if (_mass_low){
			mass_low = _mass_low->GetVal();
		} else {
			cout << " Error: could not find " << "mass_low" << endl;
			result = false;
		}
		TParameter<int>* _mass_high = (TParameter<int>*)file.Get("mass_high");
		if (_mass_high){
			mass_high = _mass_high->GetVal();
		} else {
			cout << " Error: could not find " << "mass_high" << endl;
			result = false;
		}
/*
		if (!all_loglikelihoods){
			cout << " Error: could not find " << "all_loglikelihoods" << endl;
			result = false;
		} else {
			all_loglikelihoods->SetName("all_loglikelihoods");
			//wavegraph.most_likely_intensity->SetDirectory(0);
		}
		most_likely_likelihoods = (TGraphErrors*) file.Get("most_likely_likelihoods");
		if (!most_likely_likelihoods){
			cout << " Error: could not find " << "most_likely_likelihoods" << endl;
			result = false;
		} else {
			most_likely_likelihoods->SetName("most_likely_likelihoods");
			//wavegraph.most_likely_intensity->SetDirectory(0);
		}
		all_total_intensities = (TGraphErrors*) file.Get("all_total_intensities");
		if (!all_total_intensities){
			cout << " Error: could not find " << "all_total_intensities" << endl;
			result = false;
		} else {
			all_total_intensities->SetName("all_total_intensities");
			//wavegraph.most_likely_intensity->SetDirectory(0);
		}
		most_likely_total_intensity = (TGraphErrors*) file.Get("most_likely_total_intensity");
		if (!most_likely_total_intensity){
			cout << " Error: could not find " << "most_likely_total_intensity" << endl;
			result = false;
		} else {
			most_likely_total_intensity->SetName("most_likely_total_intensity");
			//wavegraph.most_likely_intensity->SetDirectory(0);
		}
		all_evidences = (TGraphErrors*) file.Get("all_evidences");
		if (!all_evidences){
			cout << " Error: could not find " << "all_evidences" << endl;
			result = false;
		} else {
			all_evidences->SetName("all_evidences");
			//wavegraph.most_likely_intensity->SetDirectory(0);
		}
		most_likely_evidence = (TGraphErrors*) file.Get("most_likely_evidence");
		if (!all_evidences){
			cout << " Error: could not find " << "most_likely_evidence" << endl;
			result = false;
		} else {
			most_likely_evidence->SetName("most_likely_evidence");
			//wavegraph.most_likely_intensity->SetDirectory(0);
		}*/
		// reading the spin totals
		TTree* tree_spintotallist = (TTree*) file.Get("tree_spintotallist");
		if (tree_spintotallist){
			char* _spintotalname = new char[1000];
			tree_spintotallist->SetBranchAddress("spintotalname", _spintotalname);
			for (int i = 0; i < tree_spintotallist->GetEntries(); i++){
				tree_spintotallist->GetEntry(i);
				string spintotalname(_spintotalname);
				Tgraphmapit it = spin_totals.find(spintotalname);
				if (it == spin_totals.end()){
					
					TGraphErrors* &spintotalgraph = spin_totals[spintotalname];
					spintotalgraph = (TGraphErrors*) file.Get((spintotalname+"/spin_total").c_str());
					SetGraphProperties(spintotalgraph, spintotalname+"spin_total");
					/*					
					if (!spintotalgraph){
						cout << " Error: could not read spin total graph " << spintotalname+"/spin_total" << endl;
						result = false;
					} else {
						spintotalgraph->SetName((spintotalname+"spin_total").c_str());
					}*/
				} else {
					cout << " Error: spin total " << spintotalname << " already exists! " << endl;
				}			
			}
			delete[] _spintotalname;
		} else {
			cout << " Error reading list of spin totals " << endl;
		}

		file.Close();
		cout << " done " << endl;
		delete[] wavename;
	} else {
		cout << " Error in Tfitresult::Read_from_file: could not read from file: " << filename << endl;
		result = false;
	}
	return result;

}

void Tfitresult::SetGraphProperties(TGraph* graph, string name, string title){
	if (!graph){
		cout << " Error: graph " << name << " does not exist " << endl;
		return;	
	}
	if (title != "")
		graph->SetNameTitle(name.c_str(), title.c_str());
	else 
		graph->SetName(name.c_str());
	graph->SetMarkerColor(graphcolor.rootcolorindex);
	graph->SetMarkerStyle(21);
	graph->SetMarkerSize(0.5);
}


#endif /* TRPWAPLOTAMPSFRAME_CC_ */
