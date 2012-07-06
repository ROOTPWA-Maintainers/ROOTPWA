/*
 * TrpwaWaveSelectFrame.h
 *
 *  Created on: Aug 26, 2010
 *      Author: Promme
 */
#include <TGFrame.h>
#include <vector>
#include <map>
#include <string>
#include "TGButton.h"
#include "TGComboBox.h"

using namespace std;

#ifndef TRPWAWAVESELECTFRAME_H_
#define TRPWAWAVESELECTFRAME_H_

// struct to hold the wave selection per bin
struct TWaveSelection {
	vector<string> available_waves; // all amplitudes that are available to be fit
	vector<bool> selected_waves;    // Users selection for this bin
	int bin_low;                    // lower bound for this bin
	int bin_high;                   // upper bound for this bin
};

typedef vector<TWaveSelection> TWaveSelections;

class TrpwaWaveSelectFrame : public TGTransientFrame {
private:

	TWaveSelections* _waveselections; // set and selected waves

	int _selected_bin; // -1 if all bins

	TGButton* _button_allbins; // reference to button with all bins
	TGButton* _copybinbutton;  // reference to button to copy bins
	TGComboBox* _filter_box; // combobox with filters to select
	vector<TGTextButton*> _buttons_binselection; // holds the buttons to all bins

	vector<TGTextButton*> _buttons_waveselection; // holds the buttons to all waves

	vector< vector< TGTextButton* >* > _addresses_for_ids;

	// references to the buttons are stored here by quantum numbers
	map<int, vector< TGTextButton* > > _buttons_waveselection_by_JP;
	map<int, vector< TGTextButton* > > _buttons_waveselection_by_Mrefl;
	map<int, vector< TGTextButton* > > _buttons_waveselection_by_lorb;
	map<string, vector< TGTextButton* > > _buttons_waveselection_by_iso1;
	map<string, vector< TGTextButton* > > _buttons_waveselection_by_iso2;

	TGTextButton* _select_button;
	TGTextButton* _deselect_button;

	// set the buttons according to the wave list in _selected bin
	void UpdateWaveButtons();

	// copy the wave selection from one bin to an other
	void CopyWaveSelection(int frombin, int tobin);

public:

	// provide a selection of waves to be modified by the user here
	TrpwaWaveSelectFrame(TWaveSelections& waveselections);

	// create all buttons etc.
	void Build();

	virtual ~TrpwaWaveSelectFrame();

	// action to be taken on bin selection
	void BinSelectClick();

	// action to be taken on wave selection click
	void WaveSelectClick();

	// filter the waves by the specified setting
	void FilterSelectClick(int selection);

	// select all waves that are selected by the filter_box
	void SelectFiltered();

	//un select all waves that are selected by the filter_box
	void UnselectFiltered();

	// call the root script for class definition
	ClassDef(TrpwaWaveSelectFrame,0);
};

#endif /* TRPWAWAVESELECTFRAME_H_ */
