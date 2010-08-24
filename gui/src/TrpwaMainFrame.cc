/*
 * TrpwaMainFrame.cc
 *
 *  Created on: Aug 24, 2010
 *      Author: Promme
 */

#include <TApplication.h>
#include <TGClient.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TRootEmbeddedCanvas.h>
#include "TrpwaMainFrame.h"
#include <TGPack.h>
#include <TGButtonGroup.h>
#include <string>
#include <cmath>

static const int nsteps = 11;
static const string step_titles[nsteps] = {
		"      set up workspace        ",
		" fill flat phase space events ",
		"  run MC acceptance analysis  ",
		"   filter RAW data into bins  ",
		"   filter MC data into bins   ",
		"     generate PWA keyfiles    ",
		"   calculate PWA amplitudes   ",
		"  specify amplitudes for fit  ",
		"       fit partial waves      ",
		"         show results         ",
		"            predict           "
};
static const int nbatches = 5;
static const string step_batchnames[nbatches] = {
		"local ",
		"cern ",
		"gridka ",
		"mainz ",
		"munich "
};
// define which batch send job commands are implemented
// [steps][local,cern,gridka,mainz,munich]
static const bool steps_batchimplemented[nsteps][nbatches] = {
		{1,0,0,0,0},
		{1,0,0,0,0},
		{1,0,0,0,0},
		{1,0,0,0,0},
		{1,0,0,0,0},
		{1,0,0,0,0},
		{1,0,0,0,0},
		{1,0,0,0,0},
		{1,0,0,0,0},
		{1,0,0,0,0},
		{1,0,0,0,0},
};

// status of each step [0..1] (0% till 100% finished)
static float step_status[nsteps] = {
		0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
};

TrpwaMainFrame* TrpwaMainFrame::pInstance = NULL;

TrpwaMainFrame::TrpwaMainFrame(const TGWindow *p,UInt_t w,UInt_t h) : TGMainFrame(p,w,h) {
	Build();
}

void TrpwaMainFrame::CloseWindow(){
	// make also the application close when main window is closed
	TGMainFrame::CloseWindow();
	cout << " bye! " << endl;
	gApplication->Terminate(0);
}

void TrpwaMainFrame::Build(){
	// Create canvas widget
	//fEcanvas = new TRootEmbeddedCanvas("Ecanvas",this,200,200);
	//AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX| kLHintsExpandY, 10,10,10,1));
	// Create a horizontal frame widget with buttons
	TGVerticalFrame *steplist = new TGVerticalFrame(this, 800, 100);
	AddFrame(steplist, new TGLayoutHints(kLHintsCenterX,2,2,2,2));

	//TGHorizontalFrame *frame_session_options = new TGHorizontalFrame(steplist,400,20);
	TGGroupFrame *frame_session_options = new TGGroupFrame(steplist, " session ", kHorizontalFrame);
	TGTextButton *button_load = new TGTextButton(frame_session_options," Load ");
	//button_load->Connect("Clicked()","TrpwaMainFrame",this,"DoDraw()");
	frame_session_options->AddFrame(button_load, new TGLayoutHints(kLHintsCenterX,1,1,1,1));
	TGTextButton *button_modify = new TGTextButton(frame_session_options," Modify ");
	//button_modify->Connect("Clicked()","TrpwaMainFrame",this,"DoDraw()");
	frame_session_options->AddFrame(button_modify, new TGLayoutHints(kLHintsCenterX,1,1,1,1));
	TGTextButton *button_new = new TGTextButton(frame_session_options," New ");
	//button_new->Connect("Clicked()","TrpwaMainFrame",this,"DoDraw()");
	frame_session_options->AddFrame(button_new, new TGLayoutHints(kLHintsCenterX,1,1,1,1));
	steplist->AddFrame(frame_session_options, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX,1,1,1,1));

	// setup the options for each step of analysis and connect to the corresponding calls
	for (int istep = 0; istep < nsteps; istep++){
		//TGGroupFrame *frame_session = new TGGroupFrame(steplist, "", kHorizontalFrame);
		TGHorizontalFrame *frame_session = new TGHorizontalFrame(steplist);
		TGVerticalFrame *subframe = new TGVerticalFrame(frame_session, 100, 20);
		frame_session->AddFrame(subframe, new TGLayoutHints(kLHintsTop | kLHintsLeft |
				kLHintsExpandX,1,1,8,1));
		TGTextButton *button = new TGTextButton(subframe, step_titles[istep].c_str());
		subframe->AddFrame(button, new TGLayoutHints(kLHintsTop | kLHintsLeft |
				kLHintsExpandX,1,1,1,1));
		TGHProgressBar* statusbar = new TGHProgressBar(subframe, TGProgressBar::kStandard,300);
		statusbar->SetHeight(button->GetHeight());
		statusbar->SetMin(0.); statusbar->SetMax(1.);
		statusbar->SetPosition(step_status[istep]);
		progressbars.push_back(statusbar); // save for latter status update
		subframe->AddFrame(statusbar, new TGLayoutHints(kLHintsTop | kLHintsLeft |
				kLHintsExpandX,1,1,1,1));
		statusbar->ShowPos(true);
		statusbar->SetBarColor("red");
		TGButtonGroup* buttongroup = new TGButtonGroup(frame_session," batch farm type ",kHorizontalFrame);
		for (int ibatch = 0; ibatch < nbatches; ibatch++){
			//if (steps_batchimplemented[istep][ibatch]){
				TGRadioButton* radiobutton = new TGRadioButton(buttongroup, new TGHotString(step_batchnames[ibatch].c_str()));
				if (ibatch == 0) radiobutton->SetState(kButtonDown);
				if (!steps_batchimplemented[istep][ibatch])radiobutton->SetState(kButtonDisabled);
			//}
		}
		frame_session->AddFrame(buttongroup);
		steplist->AddFrame(frame_session, new TGLayoutHints(kLHintsTop | kLHintsLeft |
				kLHintsExpandX,1,1,1,1));
		//buttongroup->Show();
	}

	TGHorizontalFrame *frame_low = new TGHorizontalFrame(steplist,200,40);
	TGTextButton *button_update = new TGTextButton(frame_low," Update ");
	button_update->Connect("Clicked()","TrpwaMainFrame",this,"CheckStatus()");
	frame_low->AddFrame(button_update, new TGLayoutHints(kLHintsCenterX,5,5,15,10));
	TGTextButton *exit = new TGTextButton(frame_low," Exit ","gApplication->Terminate(0)");
	frame_low->AddFrame(exit, new TGLayoutHints(kLHintsCenterX,5,5,15,10));
	steplist->AddFrame(frame_low, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX,5,5,5,5));

	// Set a name to the main frame
	SetWindowName("rootpwa main window");
	// Map all sub windows of main frame
	MapSubwindows();
	// Initialize the layout algorithm
	Resize(GetDefaultSize());
	// Map main frame
	MapWindow();
}

void TrpwaMainFrame::CheckStatus() {
	cout << " checking the status ... " << endl;
	for (int istep = 0; istep < nsteps; istep++){
		step_status[istep]+=0.1;
		progressbars[istep]->SetPosition(step_status[istep]);
		if (floor(step_status[istep]) >= 1.){
			progressbars[istep]->SetBarColor("green");
		} else {
			progressbars[istep]->SetBarColor("red");
		}
	}
	cout << " done " << endl;
	/*
	// Draws function graphics in randomly choosen interval
	TCanvas *fCanvas = fEcanvas->GetCanvas();
	fCanvas->cd();
	for (int i = 0; i < 1000; i++){
	TF1 f1("f1","sin(x)/x",0,gRandom->Rndm()*10);
	f1.SetFillColor(19);
	f1.SetFillStyle(1);
	f1.SetLineWidth(3);
	f1.Draw();
	fCanvas->Update();
	}*/
}

TrpwaMainFrame::~TrpwaMainFrame() {
	// Clean up used widgets: frames, buttons, layouthints
	Cleanup();
}
