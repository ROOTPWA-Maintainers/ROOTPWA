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
#include "TGDoubleSlider.h"
#include "TMultiGraph.h"

using namespace std;

#ifndef TRPWAPLOTAMPSFRAME_H_
#define TRPWAPLOTAMPSFRAME_H_

typedef map<string, TTree*> Tfilemap;
typedef Tfilemap::iterator  Tfilemapit;
typedef map<string, int> Twavemap;
typedef Twavemap::iterator Twavemapit;

class TrpwaPlotAmpsFrame : public TGTransientFrame {
private:
	TCanvas* canvas_selected_waves; // canvas with drawn objects divided into 4 pads
	// keep the plotted graphs
	TMultiGraph* plotted_graphs[5];
	double masscutlow;  // current valid low  range cut of plotted waves
	double masscuthigh; // current valid high range cut of plotted waves

	TGComboBox* box_available_fits; // box with available fits, item index is a pointer to the file
	TGComboBox* box_selected_fits;  // box with selected  fits, item index is a pointer to the file
	TGComboBox* box_available_waves;        // box with available waves to select, item is a pointer to the wave
	TGComboBox* box_available_anchor_waves; // box with available waves to select, item is a pointer to the wave

	TGDoubleHSlider* slider_mass_range; // comman selected mass range for plotting

	Tfilemap available_fit_results; // map with title as key and a pointer to an opened Tree

	Tfilemap selected_fit_results; // map with title as key and a pointer to an opened Tree selected by the user
	TTree* current_fit_result; // pointer to the fit selected in the list of selected fit results

	Twavemap available_waves; // map with the wave name as key and a counter for the number of found waves in the fit results
	string current_wave , current_anchor_wave;

	// returns false if no fit results found
	bool Add_Fit_Result(string fit_result_path, // path containing .root files with fit result trees
			string fit_result_title,			// (unique) title of the fit
			string fit_result_description);		// description of the fit

	// searches for all available waves in the fit results
	vector<string>& Scan_Fit_Result(TTree* fit_results, string branchName = "fitResult_v2");

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

	// plot all intensities in the list of selected fit results
	void Plot_All_selected();

	// plot only the spin totals
	void Plot_All_selected_Spin_totals();

	// plot selected wave and compare (if given) to anchor wave
	void Plot_selected_wave();

	// save the current displayed plot
	void Save_plot();

	// method to be called when changing the mass range with a slider
	void Set_Mass_range();

	// call the root script for class definition
	ClassDef(TrpwaPlotAmpsFrame,0);
};

#endif /* TRPWAPLOTAMPSFRAME_H_ */
