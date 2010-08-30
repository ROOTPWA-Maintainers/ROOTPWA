/*
 * TrpwaWaveSelectFrame.h
 *
 *  Created on: Aug 26, 2010
 *      Author: Promme
 */
#include <TGFrame.h>
#include <vector>
#include <string>
#include "TGButton.h"

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
	vector<TGTextButton*> _buttons_binselection; // holds the buttons to all bins

	vector<TGTextButton*> _buttons_waveselection; // holds the buttons to all waves

	// set the buttons according to the wave list in _selected bin
	void UpdateWaveButtons();

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

	// call the root script for class definition
	ClassDef(TrpwaWaveSelectFrame,0);
};

#endif /* TRPWAWAVESELECTFRAME_H_ */
