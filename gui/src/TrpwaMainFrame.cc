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
#include <TRootEmbeddedCanvas.h>
#include "TrpwaMainFrame.h"
#include <TGPack.h>
#include <TGButtonGroup.h>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include "TrpwaSessionManager.h"
#include "TrpwaJobManager.h"
#include "TrpwaFitOptionsFrame.h"
#include <TGFileDialog.h>
#include <TGMsgBox.h>
#include <cstdlib>

static const int nsteps = 11;
static const string step_titles[nsteps] = {
		"      set up workspace        ",
		" fill flat phase space events ",
		"  run MC acceptance analysis  ",
		"     filter data into bins    ",
//		"           obsolete           ",
		"     generate PWA keyfiles    ",
		"   calculate PWA amplitudes   ",
		"   integrate PWA amplitudes   ",
		"  specify amplitudes for fit  ",
		"       fit partial waves      ",
		"         show results         ",
		"            predict           "
};

// specify if the user may rerun this step when finished
// false if he may rerun it
static const bool step_disabable[nsteps] = {
		true,
		true,
		true,
		true,
		true,
		true,
		true,
		false,
		false,
		false,
		false
};

// list with functions to call when button pressed
static const string func_calls[nsteps] = {
		"SetupWorkspace()",
		"FillFlatPhaseSpaceEvents()",
		"Dummy()",
		"FilterData()",
//		"Dummy()",
		"GenKeys()",
		"CalcAmps()",
		"IntAmps()",
		"SelectWaves()",
		"FitPartialWaves()",
		"ShowFitResults()",
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
//		{1,0,0,0,0},
		{1,0,0,0,0},
		{1,0,0,0,0},
		{1,0,0,0,0},
		{1,0,0,0,0},
		{1,0,0,0,0},
		{1,0,0,0,0},
};

// status of each step [0..1] (0% till 100% finished)
static float step_status[nsteps] = {
		0., 0., 0., 0., 0., 0., 0., 0., 0., 0.//, 0.
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
		stepbuttons.push_back(button); // save to change the status
		subframe->AddFrame(statusbar, new TGLayoutHints(kLHintsTop | kLHintsLeft |
				kLHintsExpandX,1,1,1,1));
		statusbar->ShowPos(true);
		statusbar->SetBarColor("red");
		/*
		TGButtonGroup* buttongroup = new TGButtonGroup(frame_session," batch farm type ",kHorizontalFrame);
		TrpwaJobManager* jobmanager = TrpwaJobManager::Instance();
		for (int ibatch = 0; ibatch < nbatches; ibatch++){
			//if (steps_batchimplemented[istep][ibatch]){
				TGRadioButton* radiobutton = new TGRadioButton(buttongroup, new TGHotString(step_batchnames[ibatch].c_str()));
				if (ibatch == 0) {
					radiobutton->SetState(kButtonDown);
				} else {
					if (jobmanager->GetFarmType() != step_batchnames[ibatch])radiobutton->SetState(kButtonDisabled);
				}
			//}
		}
		frame_session->AddFrame(buttongroup);*/
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
		step_status[3]=(current_session->Check_real_data_events()+current_session->Check_MC_data_events())/2.;
		//step_status[4]=current_session->Check_MC_data_events();
		step_status[4]=current_session->Check_PWA_keyfiles();

		bool checkentries(false);
		// ask the user whether to check the amplitude entries
		// ask if to remove problematic waves if there are some
		int returncode;
		TGMsgBox* userrespondbox = new TGMsgBox(gClient->GetRoot(), this, "check amplitude entries",
				"Do you want to check the amplitude entries as well?\n This may take some time.",
				kMBIconQuestion, (kMBYes | kMBNo), &returncode);
		if (!userrespondbox) cout << " this will be not executed " << endl; // to prevent compiler warnings
		if (returncode == kMBYes){
			checkentries = true;
		}
		step_status[5]=(
				current_session->Check_PWA_real_data_amplitudes(checkentries)   +
				current_session->Check_PWA_MC_acc_data_amplitudes(checkentries) +
				current_session->Check_PWA_MC_data_amplitudes(checkentries)
				)/3.;
		userrespondbox = new TGMsgBox(gClient->GetRoot(), this, "check Integral entries",
				"Do you want to check the Integral entries as well?\n This may take some time.",
				kMBIconQuestion, (kMBYes | kMBNo), &returncode);
		if (!userrespondbox) cout << " this will be not executed " << endl; // to prevent compiler warnings
		if (returncode == kMBYes){
			checkentries = true;
		} else {
			checkentries = false;
		}
		step_status[6]=(
				current_session->Check_PWA_MC_acc_data_integrals(checkentries)  +
				current_session->Check_PWA_MC_data_integrals(checkentries)
				)/2.;
		step_status[7]=current_session->Check_wave_lists();
		step_status[8]=current_session->Check_fits();

		step_status[9]=0.;
		step_status[10]=0.;
	}

	for (int istep = 0; istep < nsteps; istep++){
		if (!current_session) step_status[istep]=0.;
		//step_status[istep]+=0.1;
		progressbars[istep]->SetPosition(step_status[istep]);
		if (floor(step_status[istep]) >= 1.){
			progressbars[istep]->SetBarColor("green");
			if (step_disabable[istep])
				stepbuttons[istep]->SetEnabled(false);

		} else {
			progressbars[istep]->SetBarColor("red");
			stepbuttons[istep]->SetEnabled(true);
		}
	}
	cout << " done " << endl;

	int nproblematicwaves = current_session->Print_problematic_waves();
	if (nproblematicwaves){
		// ask if to remove problematic waves if there are some
		stringstream message;
		message << "There are "<< nproblematicwaves <<" files marked to be problematic.\n";
		message << "(see terminal output for details)\n";
		message << "Do you want to disable these files during further analysis?";
		int returncode;
		TGMsgBox* userrespondbox = new TGMsgBox(gClient->GetRoot(), this, "remove problematic files",
				message.str().c_str(),
				kMBIconQuestion, (kMBYes | kMBNo), &returncode);
		if (!userrespondbox) cout << " this will be not executed " << endl; // to prevent compiler warnings
		if (returncode == kMBYes){
			current_session->Remove_problematic_waves();
		}
	}
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
		const char* filetypes[] = {"Session files","*.cfg", 0, 0};
		fileinfo.fFileTypes = filetypes;
		TGFileDialog* filedialog = new TGFileDialog(gClient->GetRoot(), this, kFDSave, &fileinfo);
		if (filedialog && fileinfo.fFilename) {
			if (!current_session->Save_Session(fileinfo.fFilename)){
				cout << " Error while saving session occurred! " << endl;
			}
		}
	}
}

void TrpwaMainFrame::LoadSession(){
	TGFileInfo fileinfo;
	const char* filetypes[] = {"Session files","*.cfg", 0, 0};
	fileinfo.fFileTypes = filetypes;
	TGFileDialog* filedialog = new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fileinfo);
	if (filedialog && fileinfo.fFilename) {
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
		//cout << " before: " << endl;
		for (int i = 0; i < current_session->Get_n_bins(); i++){
			TWaveSelection waveselection;
			waveselection.selected_waves  = current_session->GetSelectedWaves(i, waveselection.available_waves, waveselection.bin_low, waveselection.bin_high);
			waveselections.push_back(waveselection);
			//cout << " found " << waveselection.selected_waves.size() << " selected waves " << waveselection.bin_low << "." << waveselection.bin_high << " of available " << waveselection.available_waves.size() << endl;
			//for (int i = 0; i < waveselection.selected_waves.size() ; i++){
				//cout << waveselection.selected_waves[i] << " " << waveselection.available_waves[i] << endl;
			//}
		}
		// set the client frames
		frame_wave_select = new TrpwaWaveSelectFrame(waveselections);
		/*cout << " after: " << endl;
		for (int i = 0; i < current_session->Get_n_bins(); i++){
			TWaveSelection& waveselection = waveselections[i];
			cout << " found " << waveselection.selected_waves.size() << " selected waves " << waveselection.bin_low << "." << waveselection.bin_high << " of available " << waveselection.available_waves.size() << endl;
			for (int i = 0; i < waveselection.selected_waves.size() ; i++){
				cout << waveselection.selected_waves[i] << " " << waveselection.available_waves[i] << endl;
			}
		}*/		

		// write now the waves to the session manager
		vector<int> bin_lowedges;
		vector< vector<string> > waves;
		for (unsigned int i = 0; i < waveselections.size(); i++){
			bin_lowedges.push_back(waveselections[i].bin_low);
			vector <string> _waves;
			vector <string>& _available_waves = waveselections[i].available_waves;
			vector <bool>& _selected_waves = waveselections[i].selected_waves;
			for (unsigned int iwave = 0; iwave < _available_waves.size(); iwave++){
				if (_selected_waves[iwave]){
					_waves.push_back(_available_waves[iwave]);
				}
			}
			waves.push_back(_waves);
		}
		if (!current_session->SetSelectedWaves(bin_lowedges,waves)){
			cout << " error setting selected waves " << endl;
		}
		//waveselections;
		current_session->Save_Session();
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

void TrpwaMainFrame::FitPartialWaves(){
	if (current_session){
		// ask the user for the fit settings he wants
		Tfit_options _fit_options;
		_fit_options.niterations = 1;
		_fit_options.rank = 1;
		_fit_options.seed = 12345;
		_fit_options.use_normalization = current_session->Is_Normalization_available();
		_fit_options.fit_consecutive = false;
		TrpwaFitOptionsFrame* frame_fit_options = new TrpwaFitOptionsFrame(_fit_options);
		if (!frame_fit_options) cout << " dummy " << endl;
		cout << current_session->Set_current_fit_title(_fit_options.title) << endl;
		cout << current_session->Set_current_fit_description(_fit_options.description) << endl;

		// calls will move to a separate class depending on the
		// farm type given, but for now implemented here

		// remove old results
		vector<string> fitresultfiles = current_session->GetFitResults();
		for (unsigned int i = 0; i < fitresultfiles.size(); i++){
			stringstream command;
			command << "mv " << fitresultfiles[i] << " " << fitresultfiles[i] << ".previous" << endl;
			cout << system(command.str().c_str()) << endl;
		}
		
		TrpwaJobManager* jobmanager = TrpwaJobManager::Instance();

		// send one job per bin
		if (!_fit_options.fit_consecutive){		
			for (int i = 0; i < current_session->Get_n_bins(); i++){
				string executedir;
				int seed = _fit_options.seed;
				stringstream command;
				for (unsigned int ifit = 0; ifit < _fit_options.niterations; ifit++){
					string fitcommand = current_session->GetFitCommand(i, executedir, _fit_options.use_normalization, _fit_options.rank, _fit_options.seed);
					command << "cd " << executedir << ";\n";
					if (seed > 0) seed++;
					command << fitcommand << ";\n";
				}
				cout << " sending fit job for bin " << i << endl;
				if (!jobmanager->SendJob(command.str(), "fit")){
					cout << " failed!" << endl;
				} else {
					cout << " done " << endl;
				}
			}
		} else {
			string executedir;
			int seed = _fit_options.seed;
			stringstream command;
			string fitcommand;
			string initialfitresult = "";
			string previousfitresult = "";
			for (unsigned int ifit = 0; ifit < _fit_options.niterations; ifit++){
				// start from 1/3 of the available fit range where the number
				// of available events is usually the biggest
				int nbins = current_session->Get_n_bins();
				int ibinstart = nbins/3;
				// perform the first fit with random start values
				fitcommand = current_session->GetFitCommand(ibinstart, executedir, _fit_options.use_normalization, _fit_options.rank, _fit_options.seed, &initialfitresult);
				command << "cd " << executedir << ";\n";	
				command << fitcommand << ";\n";
				// then move to higher values using the previous fit results
				previousfitresult = initialfitresult;
				for (int i = ibinstart+1; i < nbins; i++){
					fitcommand = current_session->GetFitCommand(i, executedir, _fit_options.use_normalization, _fit_options.rank, _fit_options.seed, &previousfitresult);
					command << "cd " << executedir << ";\n";
					command << fitcommand << ";\n";
				}
				// and finaly to lower values
				previousfitresult = initialfitresult;
				for (int i = ibinstart-1; i > -1; i--){
					fitcommand = current_session->GetFitCommand(i, executedir, _fit_options.use_normalization, _fit_options.rank, _fit_options.seed, &previousfitresult);
					command << "cd " << executedir << ";\n";
					command << fitcommand << ";\n";
				}
				cout << " sending fit job for iteration " << ifit << endl;
				if (!jobmanager->SendJob(command.str(), "fit")){
					cout << " failed!" << endl;
				} else {
					cout << " done " << endl;
				}
				if (seed > 0) seed++;
			}
		}
		cout << " saving current constellation " << current_session->Save_Session() << endl;


/*
		// write a script to submit it
		ofstream script("/tmp/_fitpartialwaves.sh");
		if (script.good()){
			script << "# generated script to call wave fitting \n" << endl;
			for (int i = 0; i < current_session->Get_n_bins(); i++){
				string executedir;
				string fitcommand = current_session->GetFitCommand(i, executedir);
				script << "cd " << executedir << endl;
				script << fitcommand << endl;
				script << "cd -" << endl;
			}
		}
		script.close();
		stringstream command;
		command << "source /tmp/_fitpartialwaves.sh";

		cout << system(command.str().c_str()) << endl;
		*/
	}
}

void TrpwaMainFrame::ShowFitResults(){
	if (current_session){
		// ask whether to save the fit results separately
		bool save(false);
		int returncode;
		TGMsgBox* userrespondbox = new TGMsgBox(gClient->GetRoot(), this, "save fit results",
			"Do you want to save the current fit constellation?",
			kMBIconQuestion, (kMBYes | kMBNo), &returncode);
		if (!userrespondbox) cout << " this will be not executed " << endl; // to prevent compiler warnings
		if (returncode == kMBYes){
			save = true;
		}
		vector<string> fit_result_paths;
		vector<string> fit_titles;
		vector<string> fit_descriptions;
		if (save){
			current_session->Save_Fit();
			current_session->Save_Session();
			current_session->Get_List_of_Fits(fit_result_paths, &fit_titles, &fit_descriptions);
		} else {
			current_session->Get_List_of_Fits(fit_result_paths, &fit_titles, &fit_descriptions);
			fit_result_paths.push_back(current_session->Get_fit_results_dir());
			fit_titles.push_back("Current");
			fit_descriptions.push_back("Current (not saved) fit constellation");
		}

		frame_plot_amps = new TrpwaPlotAmpsFrame(fit_result_paths, fit_titles, fit_descriptions);
		return;

		// obsolete part

		// write a script to submit it
		ofstream script("/tmp/_showresults.sh");
		if (script.good()){
			script << "# generated script to show fit results \n" << endl;
			script << "cd ${ROOTPWA}/src/rootscripts" << endl;
			ofstream rootmacro("/tmp/rootmacro.C");
			if (rootmacro.good()){
				rootmacro << "void rootmacro(){" << endl;
				rootmacro << "TChain* tree = new TChain(\"pwa\",\"pwa\");" << endl;
				vector<string> fitresultfiles = current_session->GetFitResults();
				for (unsigned int i = 0; i < fitresultfiles.size(); i++){
					rootmacro << "tree->Add(\"" << fitresultfiles[i] << "\");"<< endl;
				}
				rootmacro << "plotAllIntensities(tree, true);" << endl;
				rootmacro << "};" << endl;
			}
			rootmacro.close();
			script << "root -l -q /tmp/rootmacro.C" << endl;
			//script << "rm /tmp/rootmacro.C" << endl;
			script << "mv ${ROOTPWA}/src/rootscripts/wave*.ps " << current_session->Get_fit_results_dir() << "/" << endl;
			script << "cd -" << endl;
		}
		script.close();
		stringstream command;
		command << "source /tmp/_showresults.sh";

		cout << system(command.str().c_str()) << endl;
	}
	/*
	# visualization must be performed in the rootscript folder containing all needed (logon) scripts
	# run the lines of code
	# hopefully a ps file was created containing all intensities
	mv ${ROOTPWA}/src/rootscripts/waveIntensities.ps ${KPIPI_FIT_DIR}
	# ps2pdf ${KPIPI_FIT_DIR}/waveintensities.ps*/
}

void TrpwaMainFrame::GenKeys() {
	if (current_session){
		cout << " calling the key file generator " << endl;
		string keyfilegen = current_session->Get_key_file_generator_file();
		string keyfilegendir = keyfilegen.substr( 0, keyfilegen.rfind("/")+1 );
		// run the key generator within it's directory
		stringstream command;
		command << "cd "+keyfilegendir << "; ";
		// move all keyfiles there
		command << "mkdir guibackup;";
		command << "mv *.key guibackup/;";
		// create the new key files
		command << "root -l -q "+keyfilegen << "+;";
		// remove all keyfile in the destination directory
		command << "rm " << current_session->Get_key_files_dir() << "/*.key;";
		// move the key files to the destination directory
		command << "mv *.key " << current_session->Get_key_files_dir() << "/ ;";
		command << "cd -;";
		//cout << command.str();
		cout << system(command.str().c_str()) << endl;
		Update();
	}
	//system();
}

void TrpwaMainFrame::CalcAmps(){
	if (current_session){
		cout << " sending jobs for amplitude calculation " << endl;
		vector<string>  evt_real_miss;
		vector<string>  key_real_miss;
		vector<string>& amp_real_miss  = current_session->Get_PWA_real_data_amplitudes(true, &evt_real_miss, &key_real_miss);
		//vector<string>& amp_real_avail = current_session->Get_PWA_real_data_amplitudes(false);
		//cout << endl << " available amplitudes: " << endl;
		/*
		for (unsigned int i = 0; i < amp_real_avail.size(); i++){
			cout << amp_real_avail[i] << endl;
		}
		cout << endl << " missing amplitudes: " << endl;
		for (unsigned int i = 0; i < amp_real_miss.size(); i++){
			cout << amp_real_miss[i] << endl;
		}*/
		vector<string>  evt_mc_miss;
		vector<string>  key_mc_miss;
		vector<string>& amp_mc_miss  = current_session->Get_PWA_MC_data_amplitudes(true, &evt_mc_miss, &key_mc_miss);
		//vector<string>& amp_mc_avail = current_session->Get_PWA_MC_data_amplitudes(false);
		/*
		cout << endl << " available amplitudes: " << endl;
		for (unsigned int i = 0; i < amp_mc_avail.size(); i++){
			cout << amp_mc_avail[i] << endl;
		}
		cout << endl << " missing amplitudes: " << endl;
		for (unsigned int i = 0; i < amp_mc_miss.size(); i++){
			cout << amp_mc_miss[i] << endl;
		}*/
		vector<string>  evt_mc_acc_miss;
		vector<string>  key_mc_acc_miss;
		vector<string>& amp_mc_acc_miss  = current_session->Get_PWA_MC_acc_data_amplitudes(true, &evt_mc_acc_miss, &key_mc_acc_miss);
		//vector<string>& amp_mc_acc_avail = current_session->Get_PWA_MC_acc_data_amplitudes(false);
		/*
		cout << endl << " available amplitudes: " << endl;
		for (unsigned int i = 0; i < amp_mc_acc_avail.size(); i++){
			cout << amp_mc_acc_avail[i] << endl;
		}
		cout << endl << " missing amplitudes: " << endl;
		for (unsigned int i = 0; i < amp_mc_acc_miss.size(); i++){
			cout << amp_mc_acc_miss[i] << " " << evt_mc_acc_miss[i] << endl;
		}*/

		cout << " missing " <<  amp_real_miss.size() << " real data amplitudes " << endl;
		cout << " missing " <<  amp_mc_miss.size() << " real mc data amplitudes " << endl;
		cout << " missing " <<  amp_mc_acc_miss.size() << " real mc acc data amplitudes " << endl;

		cout << " calculating missing amplitudes ... " << endl;

		string pdg_table = current_session->Get_pdg_table();

		TrpwaJobManager* jobmanager = TrpwaJobManager::Instance();

		stringstream batchcommand;

		for (unsigned int i = 0; i < amp_real_miss.size(); i++){
			stringstream command;
			command << "test -s "<< amp_real_miss[i] <<" || cat "<< evt_real_miss[i] <<" | gamp -P "<< pdg_table <<" "<<key_real_miss[i]<<" > "<< amp_real_miss[i]<< " ;" << '\n';
			cout << " calculating " << amp_real_miss[i] << endl;
			batchcommand << command.str();
			//cout << system(command.str().c_str()) << endl;
			if ((i%10) == 0 || i == amp_real_miss.size()-1){
				jobmanager->SendJob(batchcommand.str(), "calc_data_amp");
				batchcommand.str("");
			}
		}

		//jobmanager->SendJob(batchcommand.str(), "calc_data_amp");
		//batchcommand.str("");

		for (unsigned int i = 0; i < amp_mc_miss.size(); i++){
			stringstream command;
			command << "test -s "<< amp_mc_miss[i] <<" || cat "<< evt_mc_miss[i] <<" | gamp -P "<< pdg_table <<" "<<key_mc_miss[i]<<" > "<< amp_mc_miss[i]<< " ;" << '\n';
			cout << " calculating " << amp_mc_miss[i] << endl;
			batchcommand << command.str();
			//cout << system(command.str().c_str()) << endl;
			if ((i%10) == 0 || i == amp_mc_miss.size()-1){
				jobmanager->SendJob(batchcommand.str(), "calc_mc_amp");
				batchcommand.str("");
			}
		}

		//jobmanager->SendJob(batchcommand.str(), "calc_mc_amp");
		//batchcommand.str("");

		for (unsigned int i = 0; i < amp_mc_acc_miss.size(); i++){
			stringstream command;
			command << "test -s "<< amp_mc_acc_miss[i] <<" || cat "<< evt_mc_acc_miss[i] <<" | gamp -P "<< pdg_table <<" "<<key_mc_acc_miss[i]<<" > "<< amp_mc_acc_miss[i]<< " ;" << '\n';
			cout << " calculating " << amp_mc_acc_miss[i] << endl;
			batchcommand << command.str();
			//cout << system(command.str().c_str()) << endl;
			if ((i%10) == 0 || i == amp_mc_acc_miss.size()-1){
				jobmanager->SendJob(batchcommand.str(), "calc_mc_acc_amp");
				batchcommand.str("");
			}
		}

		//jobmanager->SendJob(batchcommand.str(), "calc_mc_acc_amp");
		//batchcommand.str("");

		cout << " done " << endl;
	}
}

void TrpwaMainFrame::IntAmps(){
	if (current_session){
		cout << " sending jobs for amplitude integration " << endl;
		vector< vector<string> > amp_mc_avail;
		vector< vector<string> > amp_mc_acc_avail;
		vector<string>& int_mc_miss = current_session->Get_PWA_MC_data_integrals(true, &amp_mc_avail);
		vector<string>& int_mc_acc_miss = current_session->Get_PWA_MC_acc_data_integrals(true, &amp_mc_acc_avail);

		cout << " available " <<  int_mc_miss.size() << " real mc data integrals " << endl;
		cout << " available " <<  int_mc_acc_miss.size() << " real mc acc data integrals " << endl;

		cout << " integrating available amplitudes ... " << endl;

		TrpwaJobManager* jobmanager = TrpwaJobManager::Instance();

		for (unsigned int i = 0; i < int_mc_miss.size(); i++){
			if (amp_mc_avail[i].size() == 0) continue;
			stringstream command;
			command << "int ";
			cout << " integrating " << int_mc_miss[i] << endl;// << " with " << endl;
			for (unsigned int j = 0; j < amp_mc_avail[i].size(); j++){
				command << amp_mc_avail[i][j] << " ";
				//cout << amp_mc_avail[i][j] << endl;
			}
			command << " > " << int_mc_miss[i];
			//cout << command.str() << endl;
			jobmanager->SendJob(command.str(), "calintegral");
		}

		for (unsigned int i = 0; i < int_mc_acc_miss.size(); i++){
			if (amp_mc_acc_avail[i].size() == 0) continue;
			stringstream command;
			command << "int ";
			cout << " integrating " << int_mc_acc_miss[i] << endl;// << " with " << endl;
			for (unsigned int j = 0; j < amp_mc_acc_avail[i].size(); j++){
				command << amp_mc_acc_avail[i][j] << " ";
				//cout << amp_mc_acc_avail[i][j] << endl;
			}
			command << " > " << int_mc_acc_miss[i];
			//cout << command.str() << endl;
			jobmanager->SendJob(command.str(), "calcintegral");
		}
		cout << " done " << endl;
	}
}

TrpwaMainFrame::~TrpwaMainFrame() {
	// Clean up used widgets: frames, buttons, layouthints
	Cleanup();
	delete current_session;
}

#include "TrpwaEventTreeHandler.h"
void TrpwaMainFrame::FilterData(){
	if (!current_session) {
		cout << " no session loaded! " << endl;
		return;
	}
	cout << " filtering data into bins " << endl;
	TrpwaEventTreeHandler treehandler; // needed to filter data into bins
	if (!treehandler.Set_bin_path(current_session->Get_binned_data_dir())){
		cout << " Error setting up SET directory " << endl;
		return;
	}
	vector<string> data_files = current_session->Get_data_files();
	if (data_files.size() == 0){
		cout << " Warning: no data files to filter specified! " << endl;
	}
	vector<string> mc_data_files = current_session->Get_MC_data_files();
	if (mc_data_files.size() == 0){
		cout << " Warning: no mc data files to filter specified! " << endl;
	}
	data_files.insert(data_files.end(), mc_data_files.begin(), mc_data_files.end());
	//cout << " total number of files: " << data_files.size() << endl;
	if (!treehandler.Add_eventtreefiles(data_files)){
		cout << " Error loading data files " << endl;
		return;
	} else {

		if (!treehandler.Write_Trees_to_BNL_events()){
			cout << " Error filtering events from given trees! " << endl;
		}
		/*
		if (!treehandler.Write_Trees_to_rootpwa_Trees()){
			cout << " Error filtering events from given trees! " << endl;
		}*/
	}
	Update();
}

void TrpwaMainFrame::SetupWorkspace(){
	if (current_session){
		// create the binned data folders if not already there
		current_session->Check_binned_data_structure(true);
		Update();
	}
}

void TrpwaMainFrame::FillFlatPhaseSpaceEvents(){
	if (current_session){
		// send one job per bin
		TrpwaJobManager* jobmanager = TrpwaJobManager::Instance();

		for (int i = 0; i < current_session->Get_n_bins(); i++){
			string executedir;
			string fitcommand = current_session->GetgenpwCommand(i, executedir);
			stringstream command;
			command << "cd " << executedir << ";\n";
			command << fitcommand << ";\n";
			cout << " sending fit job for bin " << i << endl;
			if (!jobmanager->SendJob(command.str(), "genpw")){
				cout << " failed!" << endl;
			} else {
				cout << " done " << endl;
			}
		}
		Update();
	}
}
