/*
 * TrpwaSessionManager.cc
 *
 *  Created on: Aug 24, 2010
 *      Author: Promme
 *
 *      Manage option files, read paths, do some checks
 */

#include "TrpwaSessionManager.h"
#include <cstdlib>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include "libconfig.h++"
#include <fstream>
#include <sstream>

using namespace std;
using namespace libconfig;
//#include "TGInputDialog.h"
/*
TSessionDialogbox::TSessionDialogbox(const TGWindow* p,
		const TGWindow* main, const char* prompt, const char* defval,
		char* retstr, UInt_t options) :
		TGInputDialog(p, main, prompt, defval, retstr, options){
	// do nothing;
	Build();
}

void TSessionDialogbox::Build(){
	return;
}*/

TrpwaSessionManager::TrpwaSessionManager(){
	_bin_high = 0;
	_bin_low  = 0;
	_n_bins   = 0;
	_n_events_flat_phasespace = 0;
	_config_file = "";
	Set_ROOTPWA_dir();
}

TrpwaSessionManager::~TrpwaSessionManager(){
	;
}

bool TrpwaSessionManager::Set_Config_File(string config_file){
	bool result(false);

	_config_file = config_file;

	result = true;

	return result;
}

// save (as input) the current Session
// the (new) _config_file variable is set, too
// true if succeeded
bool TrpwaSessionManager::Save_Session(string config_file){
	bool result(false);

	// write the current configuration down
	string filename = config_file;
	if (filename == "") filename = _config_file;
	ofstream file(filename.c_str());
	if (file.good()){
		cout << " saving current session as " << filename << endl;
		file << " # configuration file for rootpwa session manager " << endl;

		file << " title = \"" << _title << "\";" << endl;
		file << " description = \"" << _description << "\";" << endl;
		file << " dir_ROOTPWA = \"" << _dir_ROOTPWA << "\";" << endl;
		file << " dir_binned_data = \"" << _dir_binned_data << "\";" << endl;
		file << " dir_key_files = \"" << _dir_key_files << "\";" << endl;
		file << " dir_fit_results = \"" << _dir_fit_results << "\";" << endl;
		file << " bin_low = " << _bin_low << ";" << endl;
		file << " bin_high = " << _bin_high << ";" << endl;
		file << " n_bins = " << _n_bins << ";" << endl;
		file << " file_flat_phasespace_config = \"" << _file_flat_phasespace_config << "\";" << endl;
		file << " n_events_flat_phasespace = " << _n_events_flat_phasespace << ";" << endl;
		file << " file_keyfile_generator = \"" << _file_keyfile_generator << "\";" << endl;

		//map <int, TrpwaBinInfo> _bins; // map with settings to each bin
		result = true;
	}
	file.close();
	return result;
}

// load a config file and override all settings
// true if succeeded
bool TrpwaSessionManager::Load_Session(string config_file){
	bool result(true);
	cout << " Loading session from " << config_file << endl;
	cout << " Warning: no checks performed yet! " << endl;

	Config _config;
	_config.readFile(config_file.c_str());
	_config_file = config_file;
	if (!_config.lookupValue("title", _title)){result = false;}
	if (!_config.lookupValue("description", _description)){result = false;}
	if (!_config.lookupValue("dir_ROOTPWA", _dir_ROOTPWA)){result = false;}
	if (!_config.lookupValue("dir_binned_data", _dir_binned_data)){result = false;}
	if (!_config.lookupValue("dir_key_files", _dir_key_files)){result = false;}
	if (!_config.lookupValue("dir_fit_results", _dir_fit_results)){result = false;}
	if (!_config.lookupValue("bin_low", _bin_low)){result = false;}
	if (!_config.lookupValue("bin_high", _bin_high)){result = false;}
	if (!_config.lookupValue("n_bins", _n_bins)){result = false;}
	if (!_config.lookupValue("file_flat_phasespace_config", _file_flat_phasespace_config)){result = false;}
	if (!_config.lookupValue("n_events_flat_phasespace", _n_events_flat_phasespace)){result = false;}
	if (!_config.lookupValue("file_keyfile_generator", _file_keyfile_generator)){result = false;}
	Initialize();

	return result;
}

// set the path to ROOTPWA
// if no path given ${ROOTPWA} batch variable is taken
// true if succeeded
bool TrpwaSessionManager::Set_ROOTPWA_dir(string path){
	bool result(false);
    string _path = path;
	// use environment variable in case no path was given
	if (_path == ""){
		cout << " searching for ${ROOTPWA} " << endl;
		_path = getenv("ROOTPWA");
		cout << " using ${ROOTPWA} = " << _path << " as ROOTPWA directory! " << endl;
		//i=system ("dir");
	}
	struct stat st;
	if(stat(_path.c_str(),&st) != 0){
		cout << " Error in TrpwaSessionManager::Set_ROOTPWA_dir(): path " << _path << " does not exist!" << endl;
	} else {
		_dir_ROOTPWA = _path;
		result = true;
	}
	return result;
}

// set the path to the binned data
// folder will be created if does not exist
// true if succeeded
bool TrpwaSessionManager::Set_binned_data_dir(string path){
	bool result(false);

	return result;
}

// set the path to the key files
// folder will be created if does not exist
// true if succeeded
bool TrpwaSessionManager::Set_key_files_dir(string path){
	bool result(false);

	return result;
}

// set the path to the fit results
// folder will be created if does not exist
// true if succeeded
bool TrpwaSessionManager::Set_fit_results_dir(string path){
	bool result(false);

	return result;
}

// set the lower bin (in MeV)
// true if valid entry
bool TrpwaSessionManager::Set_bin_low(int bin_low){
	bool result(false);

	return result;
}

// set the upper bin (in MeV)
// true if valid entry
bool TrpwaSessionManager::Set_bin_high(int bin_high){
	bool result(false);

	return result;
}

// set the number of bins
// false if (bin_high-bin_low)%n_bins != 0
bool TrpwaSessionManager::Set_n_bins(int n_bins){
	bool result(false);

	return result;
}

// set the config file for the flat phase space generator
// true if exists
bool TrpwaSessionManager::Set_flat_phasespace_config_file(string filename){
	bool result(false);

	return result;
}

// set the number of events to generate with the
// flat phase space generator
bool TrpwaSessionManager::Set_n_events_flat_phasespace(int n_events){
	bool result(false);

	return result;
}

// set the filename of the keyfile generator
// true if exists
bool TrpwaSessionManager::Set_key_file_generator_file(string filename){
	bool result(false);

	return result;
}

void TrpwaSessionManager::Set_title(string title){
	_title = title;
}

// set the description of this session
bool TrpwaSessionManager::Set_description(string description){
	bool result(false);

	return result;
}

bool TrpwaSessionManager::Initialize(){
	bool result(false);
	if ( _n_bins <= 0 ||  _bin_high-_bin_low <= 0 || (_bin_high-_bin_low)%_n_bins != 0){
		cout << " Error in TrpwaSessionManager::Initialize(): please check bin settings! " << endl;
		return result;
	}
	int step = (_bin_high-_bin_low) / _n_bins;
	for (int ibin = 0; ibin < _n_bins; ibin++){
		int lowedge = _bin_low+step*ibin;
		int highedge= _bin_low+step*(ibin+1);
		TrpwaBinInfo _bin;
		stringstream _folder_name;
		_folder_name << lowedge << "." << highedge;
		_bin.bin_folder_name = _folder_name.str();
		_bin.wave_list_file = _dir_fit_results + "/wavelist." + _bin.bin_folder_name;
		_bin.bin_high = highedge;
		_bin.bin_low  = lowedge;
		_bins[lowedge] = _bin;
		cout << " initialized bin " << _bin.wave_list_file << endl;
	}
	result = true;

	return result;
}

// returns the status [0-1] of the folder structure
// (counting folders)
float TrpwaSessionManager::Check_binned_data_structure(){
	float result(0.);

	result = 0.1;

	return result;
}

// returns the status [0-1] of flat phase space events
// (counting .genbod.evt files) + additional checks
float TrpwaSessionManager::Check_flat_phase_space_events(){
	float result(0.);

	result = 0.2;

	return result;
}

// returns the status [0-1] of real data events in the folders
// (counting .evt files) + additional checks
float TrpwaSessionManager::Check_real_data_events(){
	float result(0.);

	result = 0.3;

	return result;
}

// returns the status [0-1] of MC data events in the folders
// (counting .acc.evt files) + additional checks
float TrpwaSessionManager::Check_MC_data_events(){
	float result(0.);

	result = 0.4;

	return result;
}

// returns the status [0-1] of generated keyfiles
// (searching for .key files)
float TrpwaSessionManager::Check_PWA_keyfiles(){
	float result(0.);

	result = 0.5;

	return result;
}

// returns the status [0-1] of calculated amplitudes
// of real data
// (comparing number of .amp files with .key files in the real data folder)
float TrpwaSessionManager::Check_PWA_real_data_amplitudes(){
	float result(0.);

	result = 0.6;

	return result;
}

// returns the status [0-1] of calculated amplitudes
// of flat phase space data
// (comparing number of .amp files with .key files
// in the flat phase space data folder)
float TrpwaSessionManager::Check_PWA_MC_data_amplitudes(){
	float result(0.);

	result = 0.7;

	return result;
}

// returns the status [0-1] of calculated amplitudes
// of accepted flat phase space data
// (comparing number of .amp files with .key files
// in the accpeted events data folder)
float TrpwaSessionManager::Check_PWA_MC_acc_data_amplitudes(){
	float result(0.);

	result = 0.8;

	return result;
}

// returns the status [0-1] of wave lists
// (searches for wave lists in the bins)
float TrpwaSessionManager::Check_wave_lists(){
	float result(0.);

	result = 0.9;

	return result;
}

// returns the status [0-1] of the fits of the bins
// (searches and counts fit result files)
float TrpwaSessionManager::Check_fits(){
	float result(0.);

	result = 1.;

	return result;
}
