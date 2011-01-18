/*
 * TrpwaSessionManager.h
 *
 *  Created on: Aug 24, 2010
 *      Author: Promme
 */

//#include <TGWindow.h>
//#include <TGInputDialog.h>
#include <string>
#include <vector>
#include <map>
#include "TrpwaJobManager.h"

using namespace std;

#ifndef TRPWASESSIONMANAGER_H_
#define TRPWASESSIONMANAGER_H_
/*
class TSessionDialogbox: public TGInputDialog {
public:
	TSessionDialogbox(
			const TGWindow* p = 0, const TGWindow* main = 0,
			const char* prompt = 0, const char* defval = 0,
			char* retstr = 0, UInt_t options = kVerticalFrame);
	virtual ~TSessionDialogbox();
	// call the root script for class definition
	ClassDef(TSessionDialogbox,0);
private:
	void Build();
};*/

struct TrpwaBinInfo {
	int bin_low;  // lower bin included (MeV)
	int bin_high; // upper bin excluded (MeV)
	string bin_folder_name; // bin folder name with path to it
	string wave_list_file; // file with specified waves for fit
	vector<string> wave_list; // list of selected waves as specified in wave_list_file
};

typedef map <int, TrpwaBinInfo> TBinMap;

class TrpwaSessionManager{
public:
	TrpwaSessionManager();
	virtual ~TrpwaSessionManager();

	// set the config file containing the
	// settings of this session
	// if not existing a config file will be created
	// please call Save Session after having set the variables
	// true if succeeded
	bool Set_Config_File(string config_file);

	// Get the current config file name
	string Get_Config_File(){return _config_file;};

	// save (as input) the current Session
	// if no file given the config_file that was set is used
	// the (new) _config_file variable is set, too
	// true if succeeded
	bool Save_Session(string config_file = "");

	// load a config file and override all settings
	// true if succeeded
	bool Load_Session(string config_file);

	// set the path to ROOTPWA
	// if no path given ${ROOTPWA} batch variable is taken
	// true if succeeded
	bool Set_ROOTPWA_dir(string path = "");

	// set the filename with path for the PDG table
	// if no path given trying to find one in ${ROOTPWA}'s default destination
	// true if succeeded
	bool Set_pdg_table(string filename = "");

	// set the path to the binned data
	// true if succeeded
	bool Set_binned_data_dir(string path);

	// set the path to the key files
	// true if succeeded
	bool Set_key_files_dir(string path);

	// set the path to the fit results
	// true if succeeded
	// Call Initialize() aver having set all variables!
	bool Set_fit_results_dir(string path);

	// set the lower bin (in MeV)
	// true if valid entry
	// Call Initialize() aver having set all variables!
	bool Set_bin_low(int bin_low);

	// set the upper bin (in MeV)
	// true if valid entry
	// Call Initialize() aver having set all variables!
	bool Set_bin_high(int bin_high);

	// set the number of bins
	// false if (bin_high-bin_low)%n_bins != 0
	// Call Initialize() aver having set all variables!
	bool Set_n_bins(int n_bins);

	// set the config file for the flat phase space generator
	// true if exists
	bool Set_flat_phasespace_config_file(string filename);

	// set the number of events to generate with the
	// flat phase space generator
	bool Set_n_events_flat_phasespace(int n_events);

	// set the filename of the keyfile generator
	// true if exists
	bool Set_key_file_generator_file(string filename);

	// set the title of the session
	void Set_title(string title);

	// set the description of this session
	bool Set_description(string description);

	// Set the current fit title
	bool Set_current_fit_title(string title){
		_fit_title = title;
		return true;
	};

	// Set the description of the current fit
	bool Set_current_fit_description(string description){
		_fit_description = description;
		return true;
	};

	// Call after having set all Variables
	// creates a map of variables for fast access
	// true if succeeded
	bool Initialize();

	// get the path to ROOTPWA
	string Get_ROOTPWA_dir(){return _dir_ROOTPWA;};

	// get the file with the path to the pdg table
	string Get_pdg_table(){return _pdg_table;};

	// Get the path to the binned data
	string Get_binned_data_dir(){return _dir_binned_data;};

	// Get the path to the key files
	string Get_key_files_dir(){return _dir_key_files;};

	// Get the path to the fit results
	string Get_fit_results_dir(){return _dir_fit_results;};

	// Get the lower bin (in MeV)
	int Get_bin_low(){return _bin_low;};

	// Get the upper bin (in MeV)
	int Get_bin_high(){return _bin_high;};

	// Get the number of bins
	int Get_n_bins(){return _n_bins;};

	// Get the config file for the flat phase space generator
	string Get_flat_phasespace_config_file(){return _file_flat_phasespace_config;};

	// Get the number of events to generate with the
	// flat phase space generator
	int Get_n_events_flat_phasespace(){return _n_events_flat_phasespace;};

	// get the filename of the keyfile generator
	string Get_key_file_generator_file(){return _file_keyfile_generator;};

	// set the title of the session
	string Get_title(){return _title;};

	// Get the description of this session
	string Get_description(){return _description;};

	// Get the title of the current fit
	// note: Might be empty if not set1
	string Get_current_fit_title(){return _fit_title;};

	// Get the description of the current fit
	string Get_current_fit_description(){return _fit_description;};

	// returns the status [0-1] of the folder structure
	// (counting folders)
	// if create then missing folders will be created
	// base folder _dir_binned_data must be set to an
	// existing folder
	float Check_binned_data_structure(bool create = false);

	// returns the status [0-1] of flat phase space events
	// (counting .genbod.evt files) + additional checks
	float Check_flat_phase_space_events();

	// returns the status [0-1] of real data events in the folders
	// (counting .evt files) + additional checks
	float Check_real_data_events();

	// returns the status [0-1] of MC data events in the folders
	// (counting .acc.evt files) + additional checks
	float Check_MC_data_events();

	// returns the status [0-1] of generated keyfiles
	// (searching for .key files)
	float Check_PWA_keyfiles();

	// returns the status [0-1] of calculated amplitudes
	// of real data
	// (comparing number of .amp files with .key files in the real data folder)
	// Check_PWA_keyfiles is called
	float Check_PWA_real_data_amplitudes(bool checkentries = true);

	// Get a list of (missing/available) calculated amplitudes with full path
	// if corresponding_eventfiles is given as an empty vector
	// the filenames of the corresponding eventfiles are written out
	// if corresponding_keyfiles is given as an empty vector
	// the keyfiles with path are written out
	vector<string>& Get_PWA_real_data_amplitudes(bool missing = true,
			vector<string>* corresponding_eventfiles = NULL,
			vector<string>* corresponding_keyfiles = NULL);

	// returns the status [0-1] of calculated amplitudes
	// of flat phase space data
	// (comparing number of .amp files with .key files
	// in the flat phase space data folder)
	// Check_PWA_keyfiles is called
	float Check_PWA_MC_data_amplitudes(bool checkentries = true);

	// Get a list of (missing/available) calculated amplitudes with full path
	// if corresponding_eventfiles is given as an empty vector
	// the filenames of the corresponding eventfiles are written out
	// if corresponding_keyfiles is given as an empty vector
	// the keyfiles with path are written out
	vector<string>& Get_PWA_MC_data_amplitudes(bool missing = true,
			vector<string>* corresponding_eventfiles = NULL,
			vector<string>* corresponding_keyfiles = NULL);

	// returns the status [0-1] of calculated amplitudes
	// of accepted flat phase space data
	// (comparing number of .amp files with .key files
	// in the accpeted events data folder)
	// Check_PWA_keyfiles is called
	float Check_PWA_MC_acc_data_amplitudes(bool checkentries = true);

	// Get a list of (missing/available) calculated amplitudes with full path
	// if corresponding_eventfiles is given as an empty vector
	// the filenames of the corresponding eventfiles are written out
	// if corresponding_keyfiles is given as an empty vector
	// the keyfiles with path are written out
	vector<string>& Get_PWA_MC_acc_data_amplitudes(bool missing = true,
			vector<string>* corresponding_eventfiles = NULL,
			vector<string>* corresponding_keyfiles = NULL);

	// Get a list of (missing/available) calculated amplitudes with full path
	// folder = AMPS | PSPAMPS | ACCAMPS
	// if corresponding_eventfiles is given as an empty vector
	// the filenames of the corresponding eventfiles are written out
	// if corresponding_keyfiles is given as an empty vector
	// the keyfiles with path are written out
	vector<string>& Get_PWA_data_amplitudes(string folder, bool missing = true,
			vector<string>* corresponding_eventfiles = NULL,
			vector<string>* corresponding_keyfiles = NULL);

	// check the entries of an amplitude file
	// some files may have (0,0) entries that lead to 0 entries in the integral files
	// ampfile = filename with full path to the ampfile
	// checkentries = whether to check the entries or only if the file exists
	// nlines  = number of lines to be read and tested
	// return true if no 0 entries found
	bool Is_valid_amplitude(string ampfile, bool checkentries = true, int nlines = 10);

	// Get a list of (missing/available) calculated integrals of
	// the available amplitudes with full path
	// folder = PSPAMPS | ACCAMPS
	// if corresponding_amplitudefiles is given as an empty vector
	// the amplitude files with path are written out
	vector<string>& Get_PWA_data_integrals(string folder, bool missing = true,
			vector < vector<string> >* corresponding_amplitudefiles = NULL, bool checkentries = true);

	// same as above for one bin
	string Get_PWA_data_integral(int bin, string folder, bool missing = true,
			vector<string>* corresponding_amplitudefiles = NULL);

	// returns the status [0-1] of calculated integrals
	// of accepted flat phase space data
	// (searching for .int files)
	float Check_PWA_MC_acc_data_integrals(bool checkentries = true);

	// Get a list of (missing/available) calculated integrals of
	// the available MC accepted amplitudes with full path
	// if corresponding_amplitudefiles is given as an empty vector
	// the amplitude files with path are written out
	vector<string>& Get_PWA_MC_acc_data_integrals(bool missing = true,
			vector < vector<string> >* corresponding_amplitudefiles = NULL);

	// returns the status [0-1] of calculated integrals
	// of flat phase space data
	// (searching for .int files)
	float Check_PWA_MC_data_integrals(bool checkentries = true);

	// Get a list of (missing/available) calculated integrals of
	// the available MC amplitudes with full path
	// if corresponding_amplitudefiles is given as an empty vector
	// the amplitude files with path are written out
	vector<string>& Get_PWA_MC_data_integrals(bool missing = true,
			vector < vector<string> >* corresponding_amplitudefiles = NULL);

	// Check a normalization integral for consistency
	// current implementation checks the dimension and the amplitudes listed in it
	// to do: check the number of events, and maybe the entries themselves
	// writes detailed errors to cout
	// all filenames have to be given with full path
	bool Is_valid_norm_integral(const string norm_file, const vector<string>& corresponding_amplitudefiles, bool checkentries = true);

	// returns the status [0-1] of wave lists
	// (searches for wave lists in the bins)
	float Check_wave_lists();

	// returns the status [0-1] of the fits of the bins
	// (searches and counts fit result files)
	float Check_fits();

	// return a list of the wavenames without an extension in the ith bin
	vector<string> &GetSelectedWaves(int ibin);

	// return a list of selections to the list of "allwaves" in the ith bin
	vector<bool>   &GetSelectedWaves(int ibin, vector<string>& allwaves);

	// return a list of selections to the list of "allwaves" in the ith bin and in addition the bin edges
	vector<bool>   &GetSelectedWaves(int ibin, vector<string>& allwaves, int& bin_low, int& bin_high);

	// return a list of available Waves in the ith bin
	// Check_PWA_keyfiles() is called
	vector<string> &GetAvailableWaves(int ibin);

	// return the lower bound of the ith bin
	int	GetBinLow(int ibin);

	// return the upper bound of the ith bin
	int GetBinHigh(int ibin);

	// return the wave list of selected waves for a bin
	vector<string> &GetWaveList(int ibin);

	// return the wave list of selected waves for a bin with bin low bound and bin high bound
	vector<string> &GetWaveList(int ibin, int& bin_low, int& bin_high);

	// return the absolute path the the file with waves
	// bin_low and bin_high are the lower and upper edge of the bin
	// the fit result file name with path that should be used is also returned
	string GetWaveListFile(int ibin, int& bin_low, int& bin_high, string& fitresultfile);

	// return (if specified) data files containing the PWA vectors
	vector<string> Get_data_files(){return _data_files;};

	// return (if specified) data files containing the MC PWA vectors
	vector<string> Get_MC_data_files(){return _mc_data_files;};

	// return the command to fit a certain bin
	// executedir is the directory to execute this command in
	string GetFitCommand(int ibin, string& executedir,
			bool normalize = true, // use acceptance normalization integral if available
			unsigned int rank = 1, // set the rank of the fit
			int seed = 12345);// set the seed for the random number generator (-1 use a randomly generated number)

	// return the command to generate flat phase space events
	// into the specified bin with the number of events specified
	// in the config file
	string GetgenpwCommand(int ibin, string& executedir);

	// retrieve a list of available fit results with path to it
	vector<string>& GetFitResults();

	// set the selected Waves
	// bin_lowedes are the low edges of the bins to be set
	// selections are the corresponding waves for each bin
	bool SetSelectedWaves(const vector<int>& bin_lowedes,const vector< vector<string> >& waves);

	// save the wave lists
	bool SaveSelectedWaves();

	// sort by JPC iso1 iso2 M reflectivity
	void SortWaves(vector<string>& wavelist);

	/*
	TrpwaSessionManager& operator=(const TrpwaSessionManager& copysource) const{
		//Set_n_bins(copysource.Get_n_bins());
		// etc. to do!
	    //return this;
	};*/

	// show the entries of _keyfiles_blacklist in the std::cout
	// returns the total number of problematic waves
	int Print_problematic_waves();

	// problematic amplitude files and integral files are stored in the
	// black list that will be used to disable these files for further analysis
	// The corresponding key files if causing this troubles should be
	// fixed by the user him self
	bool Remove_problematic_waves();

	// save the current fit constellation to a specified sub folder
	// if the folder name is not given a sub directory in the fit folder
	// will be created
	// returns true if succeeded
	// folder will be put to the list of existing fit results
	// if title or description not given
	// current_title and current_description will be used
	// if title then still not set then title will be generated
	bool Save_Fit(string folder = "", string title = "", string description = "");

	// load a Fit from the specified folder
	// the existing fit will be overwritten!
	bool Load_Fit(string folder);

	void Get_List_of_Fits(vector<string>& folders, // list of folders with fits
			vector<string>* titles = NULL, 		  // optional list of fit titles
			vector<string>* descriptions = NULL); // optional list of fit descriptions

	// true if .int files are available
	// if ibin < 0 then checking for all bins
	bool Is_Normalization_available(int ibin = -1);

private:
	string _config_file; // filename with path to the config file of this session
	string _dir_ROOTPWA; // path to ROOTPWA (determined by ${ROOTPWA})
	string _pdg_table; // filename with path to the PDG table needed for amplitude calculations
	vector<string> _data_files; // files containing the analysis data trees
	vector<string> _mc_data_files; // same as above for mc
	string _dir_binned_data; // path to the bin folders (containing data, amplitudes and the integrals)
	string _dir_key_files;   // path to the key files
	string _dir_fit_results; // path to the fit results
	int _bin_low;  // lower bin in MeV included
	int _bin_high; // upper bin in MeV excluded
	int _n_bins;   // number of bins
	string _file_flat_phasespace_config; // filename with path to the flat phase space config file
	int _n_events_flat_phasespace; // number of events per bin for integration
	string _file_keyfile_generator; // filename with path to the key file generator
	string _title; // the title of the session
	string _description; // users description of this session

	string _fit_title; // name of current fit
	string _fit_description; // description of current fit

	vector<string> _fit_folders; // absolute paths of fit folders
	vector<string> _fit_titles ; // titles of the fits
	vector<string> _fit_descriptions; // descriptions of the fits

	TBinMap _bins; // map with settings to each bin

	vector<string> _keyfiles; // key files without the extension determined by accessing the keyfile folder
	int _n_keyfiles; // will be determined by accessing the key file folder

	map<string, vector<string> > _keyfiles_blacklist; // a map of key files containing files with reported problems

	TrpwaJobManager* jobManager; // to send commands performing analysis on different farm types

	// true if both lists are equal
	bool AreListsEqual(const vector<string>& list1, const vector<string>& list2);

	// gives the number of entries in list1 found in list2
	int CompareLists(const vector<string>& list1, const vector<string>& list2);

	// returns the entries of list1 (not if missing = true) found in list2
	// does not check for dubletts
	vector<string>& CompareLists(const vector<string>& list1, const vector<string>& list2, bool missing = true);

	// read the wave list given
	// put every entry with the ending .amp without the character# in it
	// into the list that is returned
	// if the wave list is not existing wavelist will be written
	// in this case Check_PWA_keyfiles() is called to get a list of expected keyfiles
	vector<string> &ReadWaveList(string filename);

	// simple stats bar
	void DrawProgressBar(int len, double percent);

	// you may provide a key file or amplitude file with path
	// the key will be returned
	string Get_key(string file);
};


#endif /* TRPWASESSIONMANAGER_H_ */
