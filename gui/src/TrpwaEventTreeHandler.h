/*
 * TrpwaEventTreeHandler.h
 *
 *  Created on: Nov 11, 2010
 *      Author: Promme
 *
 *      Handling trees with events and writing them into a binned structure.
 *      If you provide trees containing the events and the folder structure
 *      to write to, this class will filter the events into this folder structure.
 *
 *      (11.11.10)
 *      - first declarations and implementations
 *      (22.11.10)
 *      - keeping the number of entries per bin and the distribution histograms
 */

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <map>
#include <utility>
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TH1I.h"

using namespace std;

#ifndef TRPWAEVENTTREEHANDLER_H_
#define TRPWAEVENTTREEHANDLER_H_

struct TParticle{
		int geantid;
		int charge;
		double px;
		double py;
		double pz;
		double E;
};

class TrpwaEventTreeHandler {
public:
	// constructor
	TrpwaEventTreeHandler();

	// destructor
	virtual ~TrpwaEventTreeHandler();

	// add a file to the list of tree files
	// expecting a tree of the following structure
	//		TClonesArray* events_lzvecs = NULL; // lz vectors of 0..X-1 particles out, X particle in, X+1 recoil particle
	//		vector <int>* events_charges = NULL;
	//		vector <int>* events_g3pids = NULL;
	//		int		events_id;
	//		int		events_isMC; // 0 not a MC event; 1 MC event; 2 accepted MC event;
	//		tree_events->SetBranchAddress("id"  , &events_id);
	//		tree_events->SetBranchAddress("p"   , &events_lzvecs);
	//		tree_events->SetBranchAddress("q"   , &events_charges);
	//		tree_events->SetBranchAddress("pid" , &events_g3pids);
	//		tree_events->SetBranchAddress("isMC", &events_isMC);
	// with the passed tree name
	bool Add_eventtreefile(string filename, string treename = "events/tree_events_v2");

	// vector of files with the specified format as given above
	bool Add_eventtreefiles(vector<string> filenames, string treename = "events/tree_events_v2");

	// set the directories where to store the events into bins
	// folders must have the naming binlow.binhigh where
	// binlow is included and binhigh is excluded
	// bins must exist
	bool Set_bin_path(string dirname);

	// write the given events to the bin paths
	// event files are written according
	// to the data given by Add_eventtreefiles()
	// with real, mc, mc accepted events
	// attention:
	// - specify the bin baths first
	// - renames previous existing files into
	//   <filename>.previous
	// - if <filename>.previous is existing it will be
	//   removed!
	// - please add all files containing trees with
	//   real data, mc data and mc accepted data first
	//   before calling this method
	bool Write_Trees_to_BNL_events(bool overwrite = false);

	// reset all settings
	void Reset();

private:
	// storage of trees containing the events
	vector <TFile*> _eventtreefiles;
	// and the corresponding tree names
	vector <string> _eventtreenames;
	// binning is stored here: [binlow, binhigh[ MeV
	map <int, int> _bins;
	// the paths to the bins are stored here map key is the binlow
	// the paths are relative to _dir
	map <int, string> _bin_paths;
	// directory containing bins
	string _dir;
	// needed for the progress bar
	double last_percent;
	// store the output distributions as histograms and the values them self
	TH1I* _hist_data;
	TH1I* _hist_data_mc;
	TH1I* _hist_data_mc_acc;
	// first is the bin low and second is the number of entries
	//map <int, int> _nentries_data;
	//map <int, int> _nentries_data_mc;
	//map <int, int> _nentires_data_mc_acc;

	// return the binning for a path
	// in: path where the folder name is <binlow>.<binhigh>
	// returns false if path is not matching the above's structure
	// in this case binlow and binhigh are -1
	bool Path_to_bin(string path, int &binlow, int &binhigh);

	// check whether a file exists
	bool FileExists(string filename);

	// check whether a directory exists
	bool DirExists(string dirname);

	// check the contents of path
	// return objects there as <files>
	// filter for a filter extension given in <filterext>
	// if <filterext> == "dir" only directories will be filtered
	// if rmext then extentions will be removed
	// returns the objects without the path
	int GetDir (string path,
			vector<string> &files, string filterext = "", bool rmext = false);

	// writes one event into stream
	void WriteEventToStream(vector<TParticle>& particles, ofstream& stream);

	// draw a progress bar only when the length changes significantly
	void DrawProgressBar(int len, double percent);

};


#endif /* TRPWAEVENTTREEHANDLER_H_ */
