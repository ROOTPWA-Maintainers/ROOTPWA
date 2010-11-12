/*
 * TrpwaEventTreeHandler.h
 *
 *  Created on: Nov 11, 2010
 *      Author: Promme
 *
 *      Handling trees with events.
 *
 *      (11.11.10)
 *      - first declarations and implementations
 */

#include <iostream>
#include <vector>
#include <map>
#include <TFile.h>
#include <TTree.h>
#include <string>

using namespace std;

#ifndef TRPWAEVENTTREEHANDLER_H_
#define TRPWAEVENTTREEHANDLER_H_

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
	bool Set_bin_paths(string dirnames);

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

private:
	// storage of trees containing the events
	vector <TFile> _eventtreefiles;
	// binning is stored here: [binlow, binhigh[ MeV
	map <int, int> _bins;
	// the paths to the bins are stored here map key is the binlow
	map <int, string> _bin_paths;

	// return the binning for a path
	// in: path where the folder name is <binlow>.<binhigh>
	void Path_to_bin(string path, int &binlow, int &binhigh);

};


#endif /* TRPWAEVENTTREEHANDLER_H_ */
