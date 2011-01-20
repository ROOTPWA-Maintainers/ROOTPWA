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
#include <TASImage.h>

using namespace std;

#ifndef TRPWAPLOTAMPSFRAME_H_
#define TRPWAPLOTAMPSFRAME_H_

typedef map<string, TTree*> Tfilemap;
typedef Tfilemap::iterator  Tfilemapit;
typedef map<string, vector< int > > Twavemap;
typedef Twavemap::iterator Twavemapit;

struct TColor_struct{
	int rootcolorindex;
	string colorhexcode;
};

class TrpwaPlotAmpsFrame : public TGTransientFrame {
private:
	TCanvas* canvas_selected_waves; // canvas with drawn objects divided into 4 pads
	// keep the plotted graphs
	TMultiGraph* plotted_graphs[5];
	// keep the most likely graphs, too
	TMultiGraph* plotted_most_likely_graphs[5];
	// draw the most likely graph?
	bool draw_most_likely;
	double masscutlow;  // current valid low  range cut of plotted waves
	double masscuthigh; // current valid high range cut of plotted waves

	TGComboBox* box_available_fits; // box with available fits, item index is a pointer to the file
	TGComboBox* box_selected_fits;  // box with selected  fits, item index is a pointer to the file
	TGComboBox* box_available_waves;        // box with available waves to select, item is a pointer to the wave
	TGComboBox* box_available_anchor_waves; // box with available waves to select, item is a pointer to the wave

	TGDoubleHSlider* slider_mass_range; // common selected mass range for plotting

	TGCheckButton* button_show_most_likely; // set the state to draw the most likely solution

	Tfilemap available_fit_results; // map with title as key and a pointer to an opened Tree

	Tfilemap selected_fit_results; // map with title as key and a pointer to an opened Tree selected by the user
	TTree* current_fit_result; // pointer to the fit selected in the list of selected fit results

	Twavemap available_waves; // map with the wave name as key and a vector with references to the selected_fit_results as indexes starting from 0
	string current_wave , current_anchor_wave;

	// returns false if no fit results found
	bool Add_Fit_Result(string fit_result_path, // path containing .root files with fit result trees
			string fit_result_title,			// (unique) title of the fit
			string fit_result_description);		// description of the fit

	// searches for all available waves in the fit results
	vector<string>& Scan_Fit_Result(TTree* fit_results, string branchName = "fitResult_v2");

	// giving a vector of indexes referring to selected_fit_results
	// an icon is created containing colors corresponding to the fit results
	// if divide == true, empty spaces for missing references (in case the wave is not existing in the specific set)
	// are being created
	// to provide only one icon (for the fit result it self) provide only one entry in the vector
	// with the index of this fit and set divide to false
	const TGPicture* Get_Icon(vector<int>& fit_references, bool divide = true);

	// available colors for index markers
	vector<TColor_struct> available_colors;

	// creates some available colors for markers in available_colors
	void Create_available_colors();

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

	// set the most likely graph to be drawn
	void Set_show_most_likely();

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
