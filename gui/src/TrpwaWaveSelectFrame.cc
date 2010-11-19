/*
 * TrpwaWaveSelectFrame.cc
 *
 *  Created on: Aug 26, 2010
 *      Author: Promme
 */

#include "TrpwaWaveSelectFrame.h"
#include <sstream>
#include <iostream>
#include "TGButtonGroup.h"
#include "TGButton.h"

using namespace std;

#ifndef TRPWAWAVESELECTFRAME_CC_
#define TRPWAWAVESELECTFRAME_CC_

TrpwaWaveSelectFrame::TrpwaWaveSelectFrame(TWaveSelections& waveselections)
	: TGTransientFrame(gClient->GetRoot(),0,600,600, kHorizontalFrame){
	_waveselections = &waveselections; // give a reference to be kept after window distruction
	_selected_bin = -1;
	Build();
	BinSelectClick();
	//WaveSelectClick();
	fClient->WaitFor(this); // wait till the user closes this window
}

TrpwaWaveSelectFrame::~TrpwaWaveSelectFrame(){
	Cleanup();
}

void TrpwaWaveSelectFrame::Build(){
	// horizontal frame containing check boxes with all available bins
	TGGroupFrame* frame_bin_selections = new TGGroupFrame(this, " bins ", kHorizontalFrame);

	TGButtonGroup* buttongroup = new TGButtonGroup(frame_bin_selections," individual bins ",kVerticalFrame);
	// add first the button to select all bins
	TGCheckButton* radiobutton = new TGCheckButton(buttongroup, new TGHotString(" all bins "));
	_button_allbins = radiobutton;
	radiobutton->Connect("Clicked()","TrpwaWaveSelectFrame",this,"BinSelectClick()");
	radiobutton->SetState(kButtonDown);
	unsigned int nbins = _waveselections->size();
	for (unsigned int ibin = 0; ibin < nbins; ibin++){
		stringstream _text;
		_text << (*_waveselections)[ibin].bin_low << "-" << (*_waveselections)[ibin].bin_high;
		radiobutton = new TGCheckButton(buttongroup, new TGHotString(_text.str().c_str()));
		_buttons_binselection.push_back(radiobutton);
		radiobutton->Connect("Clicked()","TrpwaWaveSelectFrame",this,"BinSelectClick()");
		//if (ibin == 0) radiobutton->SetState(kButtonDown);
		// create a new radio button group if a certain number of bins is exceeded
		if (((ibin+2)%50)==0){
			frame_bin_selections->AddFrame(buttongroup);
			buttongroup = new TGButtonGroup(frame_bin_selections," individual bins ",kVerticalFrame);
		}
		//if (!steps_batchimplemented[istep][ibatch])radiobutton->SetState(kButtonDisabled);
	}
	frame_bin_selections->AddFrame(buttongroup);
	// add a button to copy from one bin to another
	_copybinbutton = new TGCheckButton(buttongroup, new TGHotString(" copy bin "));
	_copybinbutton->SetState(kButtonUp);
	buttongroup->AddFrame(_copybinbutton);

	this->AddFrame(frame_bin_selections, new TGLayoutHints(kLHintsTop | kLHintsLeft |
				kLHintsExpandX,1,1,1,1));

	// horizontal frame containing check boxes with the waves to select for the bins
	TGGroupFrame* frame_wave_selections = new TGGroupFrame(this, " waves ", kHorizontalFrame);

	buttongroup = new TGButtonGroup(frame_wave_selections," available waves ",kVerticalFrame);
	// add first the button to select all bins
	//TGRadioButton* radiobutton = new TGRadioButton(buttongroup, new TGHotString(" all bins "));
	//radiobutton->SetState(kButtonDown);
	// assuming that all selectable waves are the same for each bin
	unsigned int nwaves(0);
	if (nbins > 0)
		nwaves = (*_waveselections)[0].available_waves.size();
	for (unsigned int iwave = 0; iwave < nwaves; iwave++){
		stringstream _text;
		_text << (*_waveselections)[0].available_waves[iwave];
		TGCheckButton* button = new TGCheckButton(buttongroup, new TGHotString(_text.str().c_str()));
		_buttons_waveselection.push_back(button);
		button->SetState(kButtonUp);
		button->Connect("Clicked()","TrpwaWaveSelectFrame",this,"WaveSelectClick()");
		// create a new radio button group if a certain number of bins is exceeded
		if (((iwave+1)%50)==0){
			frame_wave_selections->AddFrame(buttongroup);
			buttongroup = new TGButtonGroup(frame_wave_selections," available waves ",kVerticalFrame);
		}
	}
	frame_wave_selections->AddFrame(buttongroup);
	this->AddFrame(frame_wave_selections, new TGLayoutHints(kLHintsTop | kLHintsLeft |
					kLHintsExpandX,1,1,1,1));


	SetWindowName("Wave selection");
	Resize(GetDefaultSize());
	SetWMSize(fWidth, fHeight);
	CenterOnParent();
	MapSubwindows();
	MapWindow();
}

// action to be taken on bin selection
void TrpwaWaveSelectFrame::BinSelectClick(){
	//cout << " Click Bin " << endl;

	bool all_bins(false);
	int _last_selected_bin = _selected_bin;
	_selected_bin = -1;
	if (_button_allbins->GetState() == kButtonDown) all_bins = true;
	for (unsigned int ibin = 0; ibin < _buttons_binselection.size(); ibin++){
		// find the new selected bin
		if (_last_selected_bin != (signed) ibin && _buttons_binselection[ibin]->GetState() == kButtonDown ){
			_selected_bin = ibin;
		}
		_buttons_binselection[ibin]->SetState(kButtonUp);
	}
	// found a new bin?
	if (_selected_bin >= 0){
		// copy from last bin?
		if (_copybinbutton->GetState()==kButtonDown && _last_selected_bin >= 0){
			CopyWaveSelection(_last_selected_bin, _selected_bin);
		}
		_buttons_binselection[_selected_bin]->SetState(kButtonDown);
		_button_allbins->SetState(kButtonUp);
	} else {
		// if not then set the last state
		if (!all_bins){
			_selected_bin = _last_selected_bin;
		}
		if (_selected_bin >= 0){
			_buttons_binselection[_selected_bin]->SetState(kButtonDown);
		} else {
			_button_allbins->SetState(kButtonDown);
		}
	}
	_copybinbutton->SetState(kButtonUp);
	UpdateWaveButtons();
}

// action to be taken on wave selection click
void TrpwaWaveSelectFrame::WaveSelectClick(){
	if (_selected_bin > (signed) (*_waveselections).size()){
		cout << " unexpected inconsistency of wave selections size and selected bin! " << endl;
		return;
	}
	unsigned int nwaves = _buttons_waveselection.size();
	if (_selected_bin != -1 && nwaves != (*_waveselections)[_selected_bin].available_waves.size()){
		cout << " Error in TrpwaWaveSelectFrame::UpdateWaveButtons(): number of initialized buttons does not fit the number of available waves for bin " << _selected_bin << endl;
		return;
	}
	// set the status of the wave in the boolean arrays

	for (unsigned iwave = 0; iwave < nwaves; iwave++){
		bool _waveselected(false);
		if (_buttons_waveselection[iwave]->GetState() == kButtonDown){
			_waveselected = true;
		}
		if (_selected_bin >= 0){ // for the selected wave
			(*_waveselections)[_selected_bin].selected_waves[iwave] = _waveselected;
		} else { // for all waves
			for (unsigned int ibin = 0; ibin < _waveselections->size(); ibin++){
				(*_waveselections)[ibin].selected_waves[iwave] = _waveselected;
			}
		}
	}
	// check for consistency in case of all waves
	if (_selected_bin == -1) UpdateWaveButtons();
}

void TrpwaWaveSelectFrame::UpdateWaveButtons(){
	// uncheck all buttons
	for (unsigned int iwave = 0; iwave < _buttons_waveselection.size(); iwave++){
		_buttons_waveselection[iwave]->SetState(kButtonUp);
		_buttons_waveselection[iwave]->SetTextColor(0x000000, false);
	}
	if (_selected_bin > (signed) (*_waveselections).size()){
		cout << " unexpected inconsistency of wave selections size and selected bin! " << endl;
		return;
	}
	if (_selected_bin >= 0){
		// check the size of the wave list first
		unsigned int nwaves = (*_waveselections)[_selected_bin].available_waves.size();
		if (nwaves != _buttons_waveselection.size()){
			cout << " Error in TrpwaWaveSelectFrame::UpdateWaveButtons(): number of initialized buttons does not fit the number of available waves for bin " << _selected_bin << endl;
			return;
		}
		if (nwaves != (*_waveselections)[_selected_bin].selected_waves.size()){
			cout << " unexpected inconsistency of number of waves and the corresponding selection! " << endl;
			return;
		}
		// update the naming corresponding to the list loaded just to be sure
		// check the selection status
		for (unsigned int iwave = 0; iwave < nwaves; iwave++){
			//cout << (*_waveselections)[_selected_bin].available_waves[iwave].c_str() << endl;
			_buttons_waveselection[iwave]->SetText((*_waveselections)[_selected_bin].available_waves[iwave]);
			if ((*_waveselections)[_selected_bin].selected_waves[iwave]){
				_buttons_waveselection[iwave]->SetState(kButtonDown);
			} else {
				_buttons_waveselection[iwave]->SetState(kButtonUp);
			}
		}
	} else {
		// go through all waves lists and check the button if even one wave is set on true
		// that way global disabled waves will stay unchecked
		if ((*_waveselections).size() > 0){
			// assuming that all bins have the same number of available waves
			unsigned int nwaves = (*_waveselections)[0].available_waves.size();
			for (unsigned int iwave = 0; iwave < nwaves; iwave++){
				string lastwave = ""; // check the waves to be at the same place for every bin
				bool keepchecked(false);
				unsigned int  ichecked(0); // count how many times this wave was checked
				for (unsigned int ibin = 0; ibin < _waveselections->size(); ibin++){
					if (nwaves != (*_waveselections)[ibin].available_waves.size()){
						cout << " Error in TrpwaWaveSelectFrame::UpdateWaveButtons(): wave lists in bin ";
						cout << ibin << " differ in respect to the first list. This is not expected yet! " << endl;
						continue;
					}
					vector <string>& _available_waves = (*_waveselections)[ibin].available_waves;
					vector <bool>&   _selected_waves  = (*_waveselections)[ibin].selected_waves;
					if (lastwave == "") lastwave = _available_waves[iwave];
					if ( _available_waves[iwave] != lastwave ){
						cout << " Error in TrpwaWaveSelectFrame::UpdateWaveButtons(): wave lists differ for bins. This is not expected yet! " << endl;
						continue;
					}
					if (_selected_waves[iwave]) {
						keepchecked = true;
						ichecked++;
					}
				}
				if (keepchecked) {
					_buttons_waveselection[iwave]->SetState(kButtonDown);
					// visualize if the wave is not set to this state for all waves
					if (ichecked != 0 && ichecked != _waveselections->size()){
						_buttons_waveselection[iwave]->SetTextColor(0xff0000, false);
					}
				} else {
					_buttons_waveselection[iwave]->SetState(kButtonUp);
				}
			}
		}
	}
}

void TrpwaWaveSelectFrame::CopyWaveSelection(int frombin, int tobin){
	int nbins = (signed)  (*_waveselections).size();
	if (frombin < 0 || frombin > nbins-1) {
		cout << " frombin is out of range: " << frombin << endl;
		return;
	}
	if (tobin < 0 || tobin > nbins-1){
		cout << " tobin is out of range:" << tobin << endl;
		return;
	}
	TWaveSelection& _frombin = (*_waveselections)[frombin];
	TWaveSelection& _tobin   = (*_waveselections)[tobin];
	if (_frombin.selected_waves.size() != _tobin.selected_waves.size()){
		cout << " Error: number of waves does not match! " << endl;
		return;
	}
	int found(0);
	for (unsigned int iwave = 0; iwave < _frombin.selected_waves.size(); iwave++){
		for (unsigned int jwave = 0; jwave < _tobin.selected_waves.size(); jwave++){
			if (_frombin.available_waves[iwave] == _tobin.available_waves[jwave]){
				_tobin.selected_waves[jwave] = _frombin.selected_waves[iwave];
				found++;
				break;
			}
		}
	}
	if ((signed)_frombin.available_waves.size() != found){
		cout << " Error in TrpwaWaveSelectFrame::CopyWaveSelection(): not all waves were found!" << endl;
	}
	/*
	for (unsigned ibin = 0; ibin < _waveselections.size(); ibin++){
		;
	}*/
}


#endif /* TRPWAWAVESELECTFRAME_CC_ */
