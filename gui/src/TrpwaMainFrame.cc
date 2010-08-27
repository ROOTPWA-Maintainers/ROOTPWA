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
#include "TrpwaSessionManager.h"
#include <TGFileDialog.h>

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
// list with functions to call when button pressed
static const string func_calls[nsteps] = {
		"Dummy()",
		"Dummy()",
		"Dummy()",
		"Dummy()",
		"Dummy()",
		"Dummy()",
		"Dummy()",
		"SelectWaves()",
		"Dummy()",
		"Dummy()",
		"Dummy()"
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
	current_session = NULL;
	steplist = NULL;
	Build();
}

void TrpwaMainFrame::CloseWindow(){
	// make also the application close when main window is closed
	TGMainFrame::CloseWindow();
	cout << " bye! " << endl;
	gApplication->Terminate(0);
}

void TrpwaMainFrame::Build(){
	//TGHorizontalFrame *frame_session_options = new TGHorizontalFrame(steplist,400,20);
	frame_session_options = new TGGroupFrame(this, " session ", kHorizontalFrame);
	TGTextButton *button_load = new TGTextButton(frame_session_options," Load ");
	button_load->Connect("Clicked()","TrpwaMainFrame",this,"LoadSession()");
	frame_session_options->AddFrame(button_load, new TGLayoutHints(kLHintsCenterX,1,1,1,1));
	TGTextButton *button_modify = new TGTextButton(frame_session_options," Modify ");
	//button_modify->Connect("Clicked()","TrpwaMainFrame",this,"DoDraw()");
	frame_session_options->AddFrame(button_modify, new TGLayoutHints(kLHintsCenterX,1,1,1,1));
	TGTextButton *button_new = new TGTextButton(frame_session_options," New ");
	button_new->Connect("Clicked()","TrpwaMainFrame",this,"NewSession()");
	frame_session_options->AddFrame(button_new, new TGLayoutHints(kLHintsCenterX,1,1,1,1));
	TGTextButton *button_save = new TGTextButton(frame_session_options," Save ");
	button_save->Connect("Clicked()","TrpwaMainFrame",this,"SaveSession()");
	frame_session_options->AddFrame(button_save, new TGLayoutHints(kLHintsCenterX,1,1,1,1));
	TGTextButton *button_save_as = new TGTextButton(frame_session_options," Save as ");
	button_save_as->Connect("Clicked()","TrpwaMainFrame",this,"SaveSessionAs()");
	frame_session_options->AddFrame(button_save_as, new TGLayoutHints(kLHintsCenterX,1,1,1,1));
	this->AddFrame(frame_session_options, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX,1,1,1,1));

	steplist = new TGVerticalFrame(this, 800, 100);
	AddFrame(steplist, new TGLayoutHints(kLHintsCenterX,2,2,2,2));

	// setup the options for each step of analysis and connect to the corresponding calls
	for (int istep = 0; istep < nsteps; istep++){
		//TGGroupFrame *frame_session = new TGGroupFrame(steplist, "", kHorizontalFrame);
		TGHorizontalFrame *frame_session = new TGHorizontalFrame(steplist);
		TGVerticalFrame *subframe = new TGVerticalFrame(frame_session, 100, 20);
		frame_session->AddFrame(subframe, new TGLayoutHints(kLHintsTop | kLHintsLeft |
				kLHintsExpandX,1,1,8,1));
		TGTextButton *button = new TGTextButton(subframe, step_titles[istep].c_str());
		button->Connect("Clicked()","TrpwaMainFrame",this,func_calls[istep].c_str());
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

	TGHorizontalFrame *frame_low = new TGHorizontalFrame(this,200,40);
	TGTextButton *button_update = new TGTextButton(frame_low," Update ");
	button_update->Connect("Clicked()","TrpwaMainFrame",this,"CheckStatus()");
	frame_low->AddFrame(button_update, new TGLayoutHints(kLHintsCenterX,5,5,15,10));
	TGTextButton *exit = new TGTextButton(frame_low," Exit ","gApplication->Terminate(0)");
	frame_low->AddFrame(exit, new TGLayoutHints(kLHintsCenterX,5,5,15,10));
	this->AddFrame(frame_low, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX,5,5,5,5));

	// Set a name to the main frame
	SetWindowName("rootpwa main window");
	// Map all sub windows of main frame
	MapSubwindows();
	// Initialize the layout algorithm
	Resize(GetDefaultSize());
	SetWMSize(fWidth, fHeight);
	// Map main frame
	MapWindow();
}

void TrpwaMainFrame::CheckStatus() {
	cout << " checking the status ... " << endl;
	if (current_session){
		/*
				"  specify amplitudes for fit  ",
				"       fit partial waves      ",
				"         show results         ",
				"            predict           "
		*/
		step_status[0]=current_session->Check_binned_data_structure();
		step_status[1]=current_session->Check_flat_phase_space_events();
		step_status[2]=1.;
		step_status[3]=current_session->Check_real_data_events();
		step_status[4]=current_session->Check_MC_data_events();
		step_status[5]=current_session->Check_PWA_keyfiles();

		step_status[6]=(
				current_session->Check_PWA_real_data_amplitudes()   +
				current_session->Check_PWA_MC_acc_data_amplitudes() +
				current_session->Check_PWA_MC_data_amplitudes()
				)/3.;

		step_status[7]=current_session->Check_wave_lists();
		step_status[8]=current_session->Check_fits();

		step_status[9]=1.;
		step_status[10]=1.;
	}

	for (int istep = 0; istep < nsteps; istep++){
		if (!current_session) step_status[istep]=0.;
		//step_status[istep]+=0.1;
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

void TrpwaMainFrame::NewSession(){
	if (current_session){
		SaveSession();
		delete current_session;
	}
	current_session = new TrpwaSessionManager();
	current_session->Set_title(" test session ");
	current_session->Set_Config_File("default.cfg");
	current_session->Set_bin_low(1000);
	current_session->Set_bin_high(2000);
	current_session->Set_n_bins(10);
	current_session->Initialize();

	Update();
}

void TrpwaMainFrame::SaveSession(){
	if (current_session){
		if (!current_session->Save_Session()){
			cout << " Error while saving session occurred! " << endl;
		}
	}
}

void TrpwaMainFrame::SaveSessionAs(){
	if (current_session){
		TGFileInfo fileinfo;
		const char* filetypes[] = {"Session files","*.cfg"};
		fileinfo.fFileTypes = filetypes;
		TGFileDialog* filedialog = new TGFileDialog(gClient->GetRoot(), this, kFDSave, &fileinfo);
		if (fileinfo.fFilename) {
			if (!current_session->Save_Session(fileinfo.fFilename)){
				cout << " Error while saving session occurred! " << endl;
			}
		}
	}
}

void TrpwaMainFrame::LoadSession(){
	TGFileInfo fileinfo;
	const char* filetypes[] = {"Session files","*.cfg"};
	fileinfo.fFileTypes = filetypes;
	TGFileDialog* filedialog = new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fileinfo);
	if (fileinfo.fFilename) {
		if (current_session){
			SaveSession();
			delete current_session;
		}
		cout << " loading " << fileinfo.fFilename << endl;
		current_session = new TrpwaSessionManager();
		string filename = fileinfo.fFilename;
		if (!current_session->Load_Session(filename)){
			cout << " Error while loading session occurred" << endl;
		} else {
			Update();
		}
	}
	// delete filedialog; will be deleted automatically
}

void TrpwaMainFrame::SelectWaves(){
	if (current_session){
	    // set a list of waves together with the waves that are already selected
		TWaveSelections waveselections;
		for (int i = 0; i < current_session->Get_n_bins(); i++){
			TWaveSelection waveselection;
			waveselection.selected_waves  = current_session->GetSelectedWaves(i, waveselection.available_waves, waveselection.bin_low, waveselection.bin_high);
			waveselections.push_back(waveselection);
			/*
			cout << " found " << waveselection.selected_waves.size() << " selected waves " << waveselection.bin_low << "." << waveselection.bin_high << " of available " << waveselection.available_waves.size() << endl;
			for (int i = 0; i < waveselection.selected_waves.size() ; i++){
				cout << waveselection.selected_waves[i] << " " << waveselection.available_waves[i] << endl;
			}*/
		}
		// set the client frames
		frame_wave_select = new TrpwaWaveSelectFrame(waveselections);
		//waveselections;
	}
}

void TrpwaMainFrame::Update(){
	if (current_session){
		frame_session_options->SetTitle(current_session->Get_title().c_str());
		CheckStatus();
		//steplist->Activate(true); ??
	} else {
		//steplist->SetEditDisabled(); ??
	}
}

void TrpwaMainFrame::Dummy(){
	cout << " not implemented yet! " << endl;
}

TrpwaMainFrame::~TrpwaMainFrame() {
	// Clean up used widgets: frames, buttons, layouthints
	Cleanup();
	delete current_session;
}
