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
#include "TGraphErrors.h"
#include <TASImage.h>

using namespace std;

#ifndef TRPWAPLOTAMPSFRAME_H_
#define TRPWAPLOTAMPSFRAME_H_

typedef map<string, TTree*> Tfilemap;
typedef Tfilemap::iterator  Tfilemapit;
typedef map<string, string> Tdirmap;
typedef Tdirmap::iterator   Tdirmapit;
typedef map<string, vector< int > > Twavemap;
typedef Twavemap::iterator Twavemapit;

struct TColor_struct{
	int rootcolorindex;
	string colorhexcode;
};

// ****************************************************************
// some type definitions to keep already scanned graphs and results
// this part maybe will be moved to a class later

// container to hold a phase difference
struct Tphasegraph {
	TGraphErrors* all_phases;
	TGraphErrors* most_likely_phase;
	TGraphErrors* all_coherences;
	TGraphErrors* most_likely_coherence;
	double mass_low; // valid range
	double mass_high;
};

// to store phase differences vs an other one given by the wave name in the key
typedef map<string, Tphasegraph>  	Tphasegraphmap;
typedef Tphasegraphmap::iterator 	Tphasegraphmapit;

typedef map<string, TGraphErrors*> 	Tgraphmap;
typedef Tgraphmap::iterator			Tgraphmapit;

// container to hold intensities and by request also phase differences vs other waves
struct Twavegraph {
	TGraphErrors* all_intensities;
	TGraphErrors* most_likely_intensity;
	Tphasegraphmap phase; // entries will be created on request only
	double mass_low;
	double mass_high;
	int J; // total spin
	int P; // parity
	int C; // charge
	int M; // spin projection
	int e; // reflectivity
	string iso1; // first isobar
	string iso2; // second isobar
	int l; // angular orbital momentum between both above
	int s; // total spin of both above
};

typedef map<string, Twavegraph> Twavegraphmap;
typedef Twavegraphmap::iterator	Twavegraphmapit;

// multipurpose
//typedef map<double, int> Tdoubleintmap;
//typedef Tdoubleintmap::iterator Tdoubleintmapit;

// container to hold created graphs from a fit result
class Tfitresult{
	// note for all returned Graphs from Getters:
	// you get a copy of an existing one, you may delete it
public:
	Twavegraphmap waves; // all waves available in this fit
	TGraphErrors* all_loglikelihoods;
	TGraphErrors* most_likely_likelihoods;
	TGraphErrors* all_evidences;
	TGraphErrors* most_likely_evidence;
	TGraphErrors* all_total_intensities;
	TGraphErrors* most_likely_total_intensity;
	Tgraphmap spin_totals; // only most likely are stored here
	double mass_low;
	double mass_high;
	TTree* fitresult;
	string name;
	TColor_struct graphcolor;
	// to do: spin totals
	// set all pointers, load intensities, initialize
	Tfitresult(
				string unique_name,	// needs a unique name
				TTree* loadedfitresult, // needs a pointer to a loaded fit result
				TColor_struct color		// common color to be used for all graphs
			);
	// delete all graphs
	virtual ~Tfitresult();

	TGraphErrors* Get_Intensity(
					string wavename, // get an intensity of a given partial wave
					bool most_likely = true // in case of many solutions only the most likely one
				);
	// the graph with a phase is created on request and then stored to
	// access it faster
	TGraphErrors* Get_Phase(
					string wavename,   // get a phase of a given partial wave
					string anchorwave, // to a corresponding anchor wave
					bool most_likely = true // in case of many solutions only the most likely one
				);

	TGraphErrors* Get_Coherence(
					string wavename,   // get a coherence of a given partial wave
					string anchorwave, // to a corresponding anchor wave
					bool most_likely = true // in case of many solutions only the most likely one
				);

	TGraphErrors* Get_Total_Intensity(
					bool most_likely = true // in case of many solutions only the most likely one
				);

	TGraphErrors* Get_Likelihood(
					bool most_likely = true // in case of many solutions only the most likely one
				);

	TGraphErrors* Get_Evidence(
					bool most_likely = true // in case of many solutions only the most likely one
				);

	// most likely spin totals are returned
	// the key is coded as JP(C)
	TGraphErrors* Get_Spin_Total(string jpc = "flat");

	// get a vector with the available spin totals
	vector<string>& Get_available_Spin_Totals();

	// creates phase and coherence for a given pair of waves
	bool Create_Phase_Coherence_graph(
					string wavename,   // between a given partial wave
					string anchorwave // and a corresponding anchor wave
				);

	// this method is called by Get_Phase and Get_Coherence and
	// does not have to be called individually
	// creates phase and coherence for all waves vs the given anchor wave
	bool Create_Phase_Coherence_graph(
					string anchorwave // and a corresponding anchor wave
				);

	// searches for all available waves in the fit results
	// and loads all intensities into graphs
	// if graphs are already loaded only the existing wave list will be returned
	vector<string>& Scan_Fit_Result(
			string branchName = "fitResult_v2"
		);

	// wavelist (if loaded the results from file or already scanned) will be returned
	vector<string>& Get_available_waves();

	// Write all graphs into a root file. If the file exists, it will be overwritten!
	// false if file could not be written
	bool Write_to_file(string filename);

	// Read graphs from a root file
	// false if file does not exist or could not be read
	bool Read_from_file(string filename);

private:
	// find the point next to 0 degree
	// in case of many solutions take the one with
	// the smallest error
	// start from there to correct the phase of consecutive points
	// to be in order and not to jump for example from +180 to -180 degree
	void Rectify_phase(TGraphErrors* graph_phase);

	// sets the graph color, marker style and marker color
	// if title is "" then title provided by the graph is used
	void SetGraphProperties(TGraph* graph, string name, string title = "");

};

typedef map<string, Tfitresult*> Tfitresultmap;
typedef Tfitresultmap::iterator Tfitresultmapit;

// ****************************************************************

class TrpwaPlotAmpsFrame : public TGTransientFrame {
private:
	TCanvas* canvas_selected_waves; // canvas with drawn objects divided into 4 pads
	// keep the plotted graphs
	TMultiGraph* plotted_graphs[5];
	// keep the most likely graphs, too
	TMultiGraph* plotted_most_likely_graphs[5];
	// draw the most likely graph?
	bool draw_most_likely;
	bool draw_datainfo;
	double masscutlow;  // current valid low  range cut of plotted waves
	double masscuthigh; // current valid high range cut of plotted waves

	TGComboBox* box_available_fits; // box with available fits, item index is a pointer to the file
	TGComboBox* box_selected_fits;  // box with selected  fits, item index is a pointer to the file
	TGComboBox* box_available_waves;        // box with available waves to select, item is a pointer to the wave
	TGComboBox* box_available_anchor_waves; // box with available waves to select, item is a pointer to the wave

	TGDoubleHSlider* slider_mass_range; // common selected mass range for plotting

	TGCheckButton* button_show_most_likely; // set the state to draw the most likely solution
	TGCheckButton* 	button_draw_datainfo; // button to overlay data info as "preliminary"

	Tfilemap available_fit_results; // map with title as key and a pointer to an opened Tree
	Tdirmap  available_fit_result_paths; // map with the title as a key and the path to the files containing the fit result trees

	Tfilemap selected_fit_results; // map with title as key and a pointer to an opened Tree selected by the user
	TTree* current_fit_result; // pointer to the fit selected in the list of selected fit results
	Tfitresultmap selected_fit_result_graphs; // map with title as key and the many created graphs with intensities, phases and so on

	Twavemap available_waves; // map with the wave name as key and a vector with references to the selected_fit_results as indexes starting from 0
	string current_wave , current_anchor_wave;

	// returns false if no fit results found
	bool Add_Fit_Result(string fit_result_path, // path containing .root files with fit result trees
			string fit_result_title,			// (unique) title of the fit
			string fit_result_description);		// description of the fit

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

	// plot "preliminary" and the data info
	void Set_draw_datainfo();

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
