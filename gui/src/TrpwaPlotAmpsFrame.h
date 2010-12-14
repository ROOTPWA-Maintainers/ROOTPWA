/*
 * TrpwaPlotAmpsFrame.h
 *
 *  Created on: Dec 13, 2010
 *      Author: Promme
 */
#include <TGFrame.h>
#include <vector>
#include <map>
#include <string>
#include "TGButton.h"
#include "TGComboBox.h"
#include "TRootEmbeddedCanvas.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

#ifndef TRPWAPLOTAMPSFRAME_H_
#define TRPWAPLOTAMPSFRAME_H_

typedef map<string, TTree*> Tfilemap;
typedef Tfilemap::iterator  Tfilemapit;

class TrpwaPlotAmpsFrame : public TGTransientFrame {
private:
	TCanvas* canvas_selected_waves; // canvas with drawn objects divided into 4 pads

	TGComboBox* box_available_fits; // box with available fits, item index is a pointer to the file
	TGComboBox* box_selected_fits;  // box with selected  fits, item index is a pointer to the file
	TGComboBox* box_available_waves;        // box with available waves to select, item is a pointer to the wave
	TGComboBox* box_available_anchor_waves; // box with available waves to select, item is a pointer to the wave

	Tfilemap available_fit_results; // map with title as key and a pointer to an opened file

	// returns false if no fit results found
	bool Add_Fit_Result(string fit_result_path, // path containing .root files with fit result trees
			string fit_result_title,			// (unique) title of the fit
			string fit_result_description);		// description of the fit

	// searches for all available waves in the fit results
	vector<string>& Scan_Fit_Result(TTree* fit_results);

public:

	// constructor
	TrpwaPlotAmpsFrame(vector<string> fit_result_paths, // paths containing fit results (all file with the ending .root will be used
			vector<string> fit_result_titles,           // titles of the fits (must be unique)
			vector<string> fit_result_descriptions      // descriptions of the fits
			);

	// create all buttons etc.
	void Build();

	virtual ~TrpwaPlotAmpsFrame();

	// will be called when an available fit file was selected
	// an item will be added to the list of selected fit files
	void Add_Fit(int pFitFile);

	// will be called when an added fit file was selected
	// if the user confirms the selected item will be removed from the list
	void Remove_Fit(int pFitFile);

	// a wave will be selected from the list of available waves
	void Select_Wave(int pWave);

	// a wave will be selected for the anchor wave (phase determination)
	void Select_Anchor_Wave(int pWave);

	// call the root script for class definition
	ClassDef(TrpwaPlotAmpsFrame,0);
};

#endif /* TRPWAPLOTAMPSFRAME_H_ */
