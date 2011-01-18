/*
 * TrpwaFitOptionsFrame.cc
 *
 *  Created on: Jan 18, 2011
 *      Author: Promme
 */

#include "TrpwaFitOptionsFrame.h"
#include "TGLabel.h"
#include <iostream>
#include <cmath>

using namespace std;

TrpwaFitOptionsFrame::TrpwaFitOptionsFrame(Tfit_options& fit_options)
	: TGTransientFrame(gClient->GetRoot(),0,400,400, kVerticalFrame){
	_fit_options = &fit_options;
	time_t rawtime;
	time ( &rawtime );
	string stime = ctime (&rawtime);
	if (_fit_options->title == ""){
		cout << stime << endl;
		_fit_options->title = stime;
	}
	Build();

	fClient->WaitFor(this); // wait till the user closes this window
}

TrpwaFitOptionsFrame::~TrpwaFitOptionsFrame(){

}

//unsigned int rank; // rank of fit
//unsigned int niterations; // number of iterations per bin
//bool use_normalization;	// use normalization integral
//int seed; // -1 : random seed else as given
void TrpwaFitOptionsFrame::Build(){
	TGLabel* label_title = new TGLabel(this, "fit title");
	this->AddFrame(label_title, new TGLayoutHints(kLHintsCenterX,1,1,1,1));
	entry_title = new TGTextEntry(this);
	entry_title->SetWidth(200);
	entry_title->SetHeight(20);
	entry_title->Connect("TextChanged(const char *)","TrpwaFitOptionsFrame",this,"CallSetOptions(const char*)");
	this->AddFrame(entry_title);
	entry_title->SetText(_fit_options->title.c_str(), false);

	TGLabel* label_description = new TGLabel(this, "fit description");
	this->AddFrame(label_description, new TGLayoutHints(kLHintsCenterX,1,1,1,1));
	entry_description = new TGTextEntry(this, _fit_options->description.c_str());
	entry_description->SetWidth(200);
	entry_description->SetHeight(20);
	entry_description->Connect("TextChanged(const char *)","TrpwaFitOptionsFrame",this,"CallSetOptions(const char*)");
	this->AddFrame(entry_description);
	entry_description->SetText(_fit_options->description.c_str(), false);


	TGLabel* label_seed = new TGLabel(this, "seed (-1 == random)");
	this->AddFrame(label_seed, new TGLayoutHints(kLHintsCenterX,1,1,1,1));
	entry_seed = new TGNumberEntry(this, (double) _fit_options->seed);
	this->AddFrame(entry_seed, new TGLayoutHints(kLHintsCenterX,1,1,1,1));
	entry_seed->Connect("ValueSet(Long_t)","TrpwaFitOptionsFrame",this,"CallSetOptions(Long_t)");
	//entry_seed->SetLimitValues(-1, 1000000);
	entry_seed->SetLimits(TGNumberEntry::kNELLimitMinMax, -1, 1000000);

	TGLabel* label_nfits = new TGLabel(this, "number of fits per bin");
	this->AddFrame(label_nfits, new TGLayoutHints(kLHintsCenterX,1,1,1,1));
	entry_nfits = new TGNumberEntry(this, (double) _fit_options->niterations);
	this->AddFrame(entry_nfits, new TGLayoutHints(kLHintsCenterX,1,1,1,1));
	entry_nfits->Connect("ValueSet(Long_t)","TrpwaFitOptionsFrame",this,"CallSetOptions(Long_t)");
	//entry_nfits->SetLimitValues(1,100);
	entry_nfits->SetLimits(TGNumberEntry::kNELLimitMinMax, 1, 100);

	TGLabel* label_rank = new TGLabel(this, "rank of fit");
	this->AddFrame(label_rank, new TGLayoutHints(kLHintsCenterX,1,1,1,1));
	entry_rank = new TGNumberEntry(this, (double) _fit_options->rank);
	this->AddFrame(entry_rank, new TGLayoutHints(kLHintsCenterX,1,1,1,1));
	entry_rank->Connect("ValueSet(Long_t)","TrpwaFitOptionsFrame",this,"CallSetOptions(Long_t)");
	//entry_rank->SetLimitValues(1,10);
	entry_rank->SetLimits(TGNumberEntry::kNELLimitMinMax, 1, 10);

	button_use_normalization = new TGCheckButton(this, new TGHotString("use normalization"));
	// if use normalization not given then it is assumed not to be available
	if (!_fit_options->use_normalization){
		button_use_normalization->SetState(kButtonUp);
		button_use_normalization->SetEnabled(kFALSE);
	} else {
		button_use_normalization->SetState(kButtonDown);
	}
	this->AddFrame(button_use_normalization,new TGLayoutHints(kLHintsCenterX,1,1,1,1));
	button_use_normalization->Connect("Clicked()","TrpwaFitOptionsFrame",this,"SetOptions()");

	// Set a name to the main frame
	SetWindowName("fit job options");
	// Map all sub windows of main frame
	MapSubwindows();
	// Initialize the layout algorithm
	Resize(GetDefaultSize());
	SetWMSize(fWidth, fHeight);
	// Map main frame
	MapWindow();
}

#include <sstream>
void TrpwaFitOptionsFrame::SetOptions(){
	_fit_options->niterations = (unsigned int)floor(entry_nfits->GetNumber());
	_fit_options->rank	= (unsigned int)floor(entry_rank->GetNumber());
	_fit_options->seed  = (unsigned int)floor(entry_seed->GetNumber());
	_fit_options->use_normalization = button_use_normalization->GetState();
	_fit_options->title = entry_title->GetText();
	stringstream _description;
	_description << entry_description->GetText();
	_description << " (";
	_description << " nfits: "<< _fit_options->niterations;
	_description << " rank: "<< _fit_options->rank;
	_description << " seed: "<< _fit_options->seed;
	_description << " acc norm: " << _fit_options->use_normalization;
	time_t rawtime;
	time ( &rawtime );
	string stime = ctime (&rawtime);
	_description << " date: " << stime << ")";
	_fit_options->description = _description.str();

	/*
	cout << " niterations are " << _fit_options->niterations << endl;
	cout << " rank is " << _fit_options->rank << endl;
	cout << " seed is " << _fit_options->seed << endl;
	if (_fit_options->use_normalization)
		cout << " normalization is on " << endl;
	else
		cout << " normalization is off " << endl;
	cout << endl;
	*/
}

void TrpwaFitOptionsFrame::CallSetOptions(long value){
	SetOptions();
}

void TrpwaFitOptionsFrame::CallSetOptions(const char *text){
	SetOptions();
}

