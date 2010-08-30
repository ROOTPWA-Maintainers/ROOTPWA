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
#include <stdlib.h>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include "libconfig.h++"
#include <fstream>
#include <sstream>
#include <algorithm>

#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

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

		if ( _bins.size() > 0 ){
			if (_bins.size() != (unsigned) _n_bins){
				cout << " Error in TrpwaSessionManager::Save_Session(): number of stored bins does not match the bin number stored! " << endl;
				return result;
			}
			// write also the individual bins down
			file << " bins_lowedge = [ ";
			for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
				if (it != _bins.begin()) file << " , ";
				file << it->second.bin_low;
			}
			file << " ]; " << endl;

			file << " bins_highedge = [ ";
			for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
				if (it != _bins.begin()) file << " , ";
				file << it->second.bin_high;
			}
			file << " ]; " << endl;

			file << " bins_wavelist = [ ";
			for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
				if (it != _bins.begin()) file << " , ";
				file << "\"" << it->second.wave_list_file << "\"";
			}
			file << " ]; " << endl;

			file << " bins_foldername = [ ";
			for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
				if (it != _bins.begin()) file << " , ";
				file << "\"" << it->second.bin_folder_name << "\"";
			}
			file << " ]; " << endl;
		}


		//map <int, TrpwaBinInfo> _bins; // map with settings to each bin
		result = true;
	}
	file.close();

	if (!SaveSelectedWaves()) result = false;
	
	return result;
}

// load a config file and override all settings
// true if succeeded
bool TrpwaSessionManager::Load_Session(string config_file){
	bool result(true);
	if (!FileExists(config_file)){
		cout << " Error in TrpwaSessionManager::Load_Session(): " << config_file << " does not exist! " << endl;
		return false;
	}

	cout << " Loading session from " << config_file << endl;
	Config _config;
	_config.readFile(config_file.c_str());
	_config_file = config_file;
	int _val_int;
	string _val_string;

	if (!_config.lookupValue("title", _val_string)){
		result = false;
	} else {
		Set_title(_val_string);
	}

	if (!_config.lookupValue("description", _val_string)){
		result = false;
	} else {
		if (!Set_description(_val_string)) result = false;
	}

	if (!_config.lookupValue("dir_ROOTPWA", _val_string)){
		result = false;
	} else {
		if (!Set_ROOTPWA_dir(_val_string)) result = false;
	}

	if (!_config.lookupValue("dir_binned_data", _val_string)){
		result = false;
	} else {
		if (!Set_binned_data_dir(_val_string)) result = false;
	}

	if (!_config.lookupValue("dir_key_files", _val_string)){
		result = false;
	} else {
		if (!Set_key_files_dir(_val_string)) result = false;
	}

	if (!_config.lookupValue("dir_fit_results", _val_string)){
		result = false;
	} else {
		if (!Set_fit_results_dir(_val_string)) result = false;
	}

	if (!_config.lookupValue("bin_low", _val_int)){
		result = false;
	} else {
		if (!Set_bin_low(_val_int)) result = false;
	}

	if (!_config.lookupValue("bin_high", _val_int)){
		result = false;
	} else {
		if (!Set_bin_high(_val_int)) result = false;
	}

	if (!_config.lookupValue("n_bins", _val_int)){
		result = false;
	} else {
		if (!Set_n_bins(_val_int)) result = false;
	}

	if (!_config.lookupValue("file_flat_phasespace_config", _val_string)){
		result = false;
	} else {
		if (!Set_flat_phasespace_config_file(_val_string)) result = false;
	}

	if (!_config.lookupValue("n_events_flat_phasespace", _val_int)){
		result = false;
	} else {
		if (!Set_n_events_flat_phasespace(_val_int)) result = false;
	}

	if (!_config.lookupValue("file_keyfile_generator", _val_string)){
		result = false;
	} else {
		if (!Set_key_file_generator_file(_val_string)) result = false;
	}

	if (    _config.exists("bins_lowedge")&&
			_config.exists("bins_highedge")&&
			_config.exists("bins_wavelist")&&
			_config.exists("bins_foldername")
			){
		Setting& _setting_lowedges    = _config.lookup("bins_lowedge");
		Setting& _setting_highedges   = _config.lookup("bins_highedge");
		Setting& _setting_wavelists   = _config.lookup("bins_wavelist");
		Setting& _setting_foldernames = _config.lookup("bins_foldername");
		if (!(  _setting_lowedges.isArray() &&
				_setting_highedges.isArray() &&
				_setting_wavelists.isArray() &&
				_setting_foldernames.isArray())){
			cout << " Error in TrpwaSessionManager::Load_Session(): check bin input arrays in config file! " << endl;
			return false;
		}
		if (!(  _setting_lowedges.getLength() == _n_bins &&
				_setting_highedges.getLength() == _n_bins &&
				_setting_wavelists.getLength() == _n_bins &&
				_setting_foldernames.getLength() == _n_bins)){
			cout << " Error in TrpwaSessionManager::Load_Session(): Number of bins does not match the number of entries! " << endl;
			cout << "                                               Check bin input arrays in the config file! " << endl;
			return false;
		}
		if (_bins.size() != 0){
			cout << " Error in TrpwaSessionManager::Load_Session(): bin structure already exists! Replacing." << endl;
			_bins.clear();
		}
		for (int i = 0; i < _n_bins; i++){
			TrpwaBinInfo bin;
			bin.bin_folder_name = (const char*) _setting_foldernames[i];
			bin.wave_list_file  = (const char*) _setting_wavelists  [i];
			bin.bin_low         = _setting_lowedges   [i];
			bin.bin_high        = _setting_highedges  [i];
			bin.wave_list	    = ReadWaveList(_dir_fit_results +"/"+ bin.wave_list_file);
			_bins[bin.bin_low]  = bin;
		}
		if (_bins.size() != (unsigned) _n_bins){
			cout << " Error in TrpwaSessionManager::Load_Session(): number of bins that differ is not correct! " << endl;
			return false;
		}
	} else {
		cout << " No valid bins found: creating new bin structure! " << endl;
		Initialize();
	}
	//Initialize();

	return result;
}

// set the path to ROOTPWA
// if no path given ${ROOTPWA} batch variable is taken
// true if succeeded
bool TrpwaSessionManager::Set_ROOTPWA_dir(string path){
	bool result(false);
    string _path = path;
	// use environment variable in case no or wrong path was given
	if (_path == "" || !DirExists(path)){
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

bool TrpwaSessionManager::Set_binned_data_dir(string path){
	bool result(false);

	if (DirExists(path)){
		_dir_binned_data = path;
		result = true;
	}

	return result;
}

bool TrpwaSessionManager::Set_key_files_dir(string path){
	bool result(false);

	if (DirExists(path)){
		_dir_key_files = path;
		result = true;
	}

	return result;
}

bool TrpwaSessionManager::Set_fit_results_dir(string path){
	bool result(false);

	if (DirExists(path)){
		_dir_fit_results = path;
		result = true;
	}

	return result;
}

bool TrpwaSessionManager::Set_bin_low(int bin_low){
	bool result(false);

	if (bin_low > 0){
		_bin_low = bin_low;
		result = true;
	}

	return result;
}

bool TrpwaSessionManager::Set_bin_high(int bin_high){
	bool result(false);

	if (bin_high > 0){
		_bin_high = bin_high;
		result = true;
	}

	return result;
}

// set the number of bins
// false if (bin_high-bin_low)%n_bins != 0
bool TrpwaSessionManager::Set_n_bins(int n_bins){
	bool result(false);

	if (n_bins > 0) {
		_n_bins = n_bins;
		result = true;
	}

	return result;
}

// set the config file for the flat phase space generator
// true if exists
bool TrpwaSessionManager::Set_flat_phasespace_config_file(string filename){
	bool result(false);

	if (FileExists(filename)){
		_file_flat_phasespace_config = filename;
		result = true;
	}

	return result;
}

// set the number of events to generate with the
// flat phase space generator
bool TrpwaSessionManager::Set_n_events_flat_phasespace(int n_events){
	bool result(false);

	if (n_events > 0){
		_n_events_flat_phasespace = n_events;
		result = true;
	}

	return result;
}

// set the filename of the keyfile generator
// true if exists
bool TrpwaSessionManager::Set_key_file_generator_file(string filename){
	bool result(false);

	if (FileExists(filename)){
		_file_keyfile_generator = filename;
		result = true;
	}

	return result;
}

void TrpwaSessionManager::Set_title(string title){
	_title = title;
}

// set the description of this session
bool TrpwaSessionManager::Set_description(string description){
	bool result(true);

	_description = description;

	return result;
}

bool TrpwaSessionManager::Initialize(){
	bool result(false);
	if (_bins.size() != 0){
		cout << " Error in TrpwaSessionManager::Initialize(): A bin structure already exists! " << endl;
		return result;
	}
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
		_bin.wave_list_file = /*_dir_fit_results+"/" +*/ "wavelist." + _bin.bin_folder_name;
		_bin.bin_high = highedge;
		_bin.bin_low  = lowedge;
		_bin.wave_list = ReadWaveList(_dir_fit_results+"/"+_bin.wave_list_file);
		_bins[lowedge] = _bin;
		cout << " initialized bin " << _bin.wave_list_file << endl;
	}
	result = true;

	return result;
}

// returns the status [0-1] of the folder structure
// (counting folders)
float TrpwaSessionManager::Check_binned_data_structure(){
	cout << " checking binned data structure in " << _dir_binned_data << endl;
	float result(0.);

	if ( _n_bins <= 0 ) return result;
	int counter(0);
	for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
		bool complete(true);
		string _bindir = _dir_binned_data + "/" + it->second.bin_folder_name;
		if (!DirExists(_bindir)){
			complete = false;
		}
		if (!DirExists(_bindir + "/AMPS")){
			complete = false;
		}
		if (!DirExists(_bindir + "/PSPAMPS")){
			complete = false;
		}
		if (!DirExists(_bindir + "/ACCAMPS")){
			complete = false;
		}
		if (!DirExists(_bindir + "/MC")){
			complete = false;
		}
		if (complete) counter++;
	}
	cout << " found " << counter << " of " << _n_bins << " expected folders" << endl;
	result = (double) counter / ((double) _n_bins);

	return result;
}

// returns the status [0-1] of flat phase space events
// (counting .genbod.evt files) + additional checks
float TrpwaSessionManager::Check_flat_phase_space_events(){
	cout << " searching flat phase space events " << endl;
	float result(0.);

	if ( _n_bins <= 0 ) return result;
	int counter(0);
	for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
		if (FileExists(_dir_binned_data + "/" + it->second.bin_folder_name + "/" + it->second.bin_folder_name + ".genbod.evt")){
			counter++;
		}
	}
	cout << " found " << counter << " of " << _n_bins << " expected .genbod.evt files" << endl;
	result = (double) counter / ((double) _n_bins);

	return result;
}

// returns the status [0-1] of real data events in the folders
// (counting .evt files) + additional checks
float TrpwaSessionManager::Check_real_data_events(){
	cout << " searching real data events " << endl;
	float result(0.);

	if ( _n_bins <= 0 ) return result;
	int counter(0);
	for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
		if (FileExists(_dir_binned_data + "/" + it->second.bin_folder_name + "/" + it->second.bin_folder_name + ".evt")){
			counter++;
		}
	}
	cout << " found " << counter << " of " << _n_bins << " expected .evt files" << endl;
	result = (double) counter / ((double) _n_bins);

	return result;
}

// returns the status [0-1] of MC data events in the folders
// (counting .acc.evt files) + additional checks
float TrpwaSessionManager::Check_MC_data_events(){
	cout << " searching MC acceptance events " << endl;
	float result(0.);

	if ( _n_bins <= 0 ) return result;
	int counter(0);
	for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
		if (FileExists(_dir_binned_data + "/" + it->second.bin_folder_name + "/" + it->second.bin_folder_name + ".acc.evt")){
			counter++;
		}
	}
	cout << " found " << counter << " of " << _n_bins << " expected .acc.evt files" << endl;
	result = (double) counter / ((double) _n_bins);

	return result;
}

// returns the status [0-1] of generated keyfiles
// (searching for .key files)
float TrpwaSessionManager::Check_PWA_keyfiles(){
	cout << " searching for keyfiles in " << _dir_key_files << endl;
	float result(0.);

	_n_keyfiles = GetDir(_dir_key_files, _keyfiles, ".key", true);
	SortWaves(_keyfiles);

	if (_n_keyfiles > 0) result = 1.;

	cout << " found " << _n_keyfiles << " .key files " << endl;

	return result;
}

// returns the status [0-1] of calculated amplitudes
// of real data
// (comparing number of .amp files with .key files in the real data folder)
float TrpwaSessionManager::Check_PWA_real_data_amplitudes(){
	Check_PWA_keyfiles();
	cout << " searching for amplitude files in AMPS folders " << endl;
	float result(0.);

	if ( _n_bins <= 0 ) return result;
	int counter(0);
	for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
		string _dir = _dir_binned_data + "/" + it->second.bin_folder_name + "/" + "AMPS/";
		vector<string> _amps;
		GetDir(_dir, _amps, ".amp", true);
		if (AreListsEqual(_keyfiles, _amps)) {
			counter+=_amps.size();
		}
		//if (FileExists(_dir_binned_data + "/" + it->second.bin_folder_name + "/" + "AMPS/")){

		//}
	}
	cout << " found " << counter << " of " << _n_bins*_keyfiles.size() << " expected .amp files" << endl;
	result = (double) counter / ((double) _n_bins*_keyfiles.size());

	return result;
}

// returns the status [0-1] of calculated amplitudes
// of flat phase space data
// (comparing number of .amp files with .key files
// in the flat phase space data folder)
float TrpwaSessionManager::Check_PWA_MC_data_amplitudes(){
	Check_PWA_keyfiles();
	cout << " searching for amplitude files in PSPAMPS folders " << endl;
	float result(0.);

	if ( _n_bins <= 0 ) return result;
	int counter(0);
	for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
		string _dir = _dir_binned_data + "/" + it->second.bin_folder_name + "/" + "PSPAMPS/";
		vector<string> _amps;
		GetDir(_dir, _amps, ".amp", true);
		if (AreListsEqual(_keyfiles, _amps)) {
			counter+=_amps.size();
		}
		//if (FileExists(_dir_binned_data + "/" + it->second.bin_folder_name + "/" + "AMPS/")){

		//}
	}
	cout << " found " << counter << " of " << _n_bins*_keyfiles.size() << " expected .amp files" << endl;
	result = (double) counter / ((double) _n_bins*_keyfiles.size());

	return result;
}

// returns the status [0-1] of calculated amplitudes
// of accepted flat phase space data
// (comparing number of .amp files with .key files
// in the accpeted events data folder)
float TrpwaSessionManager::Check_PWA_MC_acc_data_amplitudes(){
	Check_PWA_keyfiles();
	cout << " searching for amplitude files in ACCAMPS folders " << endl;
	float result(0.);

	if ( _n_bins <= 0 ) return result;
	int counter(0);
	for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
		string _dir = _dir_binned_data + "/" + it->second.bin_folder_name + "/" + "ACCAMPS/";
		vector<string> _amps;
		GetDir(_dir, _amps, ".amp", true);
		if (AreListsEqual(_keyfiles, _amps)) {
			counter+=_amps.size();
		}
		//if (FileExists(_dir_binned_data + "/" + it->second.bin_folder_name + "/" + "AMPS/")){

		//}
	}
	cout << " found " << counter << " of " << _n_bins*_keyfiles.size() << " expected .amp files" << endl;
	result = (double) counter / ((double) _n_bins*_keyfiles.size());

	return result;
}

// returns the status [0-1] of wave lists
// (searches for wave lists in the bins)
float TrpwaSessionManager::Check_wave_lists(){
	cout << " searching for wave lists in " << _dir_fit_results << endl;
	float result(0.);

	if ( _n_bins <= 0 ) return result;
	int counter(0);
	for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
		if (FileExists(_dir_fit_results + "/" + it->second.wave_list_file)){
			counter++;
		}
	}
	cout << " found " << counter << " of " << _n_bins << " expected wave list files" << endl;
	result = (double) counter / ((double) _n_bins);

	return result;
}

// returns the status [0-1] of the fits of the bins
// (searches and counts fit result files)
float TrpwaSessionManager::Check_fits(){
	cout << " searching for fit results in " << _dir_fit_results << endl;
	float result(0.);

	if ( _n_bins <= 0 ) return result;
	int counter(0);
	for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
		if (FileExists(_dir_fit_results + "/" + "fit_result_" + it->second.wave_list_file + ".root")){
			counter++;
		}
	}
	cout << " found " << counter << " of " << _n_bins << " expected fit result files" << endl;
	result = (double) counter / ((double) _n_bins);

	return result;
}

// return a list of the wavenames without an extension in the ith bin
vector<string> &TrpwaSessionManager::GetSelectedWaves(int ibin){
	vector<string>* result = new vector<string>();
	vector<string> _allwaves;
	vector<bool> _selected = GetSelectedWaves(ibin, _allwaves);
	for (unsigned int i = 0; i < _allwaves.size(); i++){
		if (_selected[i]) result->push_back(_allwaves[i]);
	}
	return *result;
}

	// return a list of selections to the list of "allwaves" in the ith bin
vector<bool>   &TrpwaSessionManager::GetSelectedWaves(int ibin, vector<string>& allwaves){
	int _bin_low;
	int _bin_high;
	return GetSelectedWaves(ibin, allwaves, _bin_low, _bin_high);
}

vector<bool>   &TrpwaSessionManager::GetSelectedWaves(int ibin, vector<string>& allwaves, int& bin_low, int& bin_high){
	Check_PWA_keyfiles(); // to (re)initialize the _kefiles variable
	vector<bool>* result = new vector<bool>();
	if (_keyfiles.size() == 0){
		cout << " Error in TrpwaSessionManager::GetSelectedWaves(): No Keyfiles available! " << endl;
		return *result;
	}
	allwaves = _keyfiles;
	if (_bins.size() > 0){
		int bincounter(0);
		for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
			if (bincounter == ibin){
				bin_low = (*it).second.bin_low;
				bin_high= (*it).second.bin_high;
				vector<string> _selected_waves = (*it).second.wave_list;
				// find the corresponding entries in the list of available keys
				for (unsigned int i = 0; i < allwaves.size(); i++){
					bool found(false);
					for (vector<string>::iterator j = _selected_waves.begin(); j != _selected_waves.end(); j++){
						// the key was found in the mother list
						// write it out and delete it from the list to see afterwards if all keys were also available
						if ((*j) == allwaves[i]) {
							result->push_back(true);
							found = true;
							_selected_waves.erase(j);
							break;
						}
					}
					if (!found) result->push_back(false);
				}
				// the list should be empty now, check this
				// if not we have selected waves that are not available in the key list -> inconsistency
				if (_selected_waves.size() != 0){
					cout << " Error in TrpwaSessionManager::GetSelectedWaves(): there are " << _selected_waves.size() << " amps in ";
					cout << (*it).second.wave_list_file << " specified that are not available in the list of keys! "  << endl;
					result->clear();
					allwaves.clear();
					return *result;
				}
				if (allwaves.size()!=result->size()){
					cout << " unexpected error in TrpwaSessionManager::GetSelectedWaves(), please inform the coder! " << endl;
					result->clear();
					allwaves.clear();
					return *result;
				}
				break; // the correct bin was found
			}
			bincounter++;
		}
	} else {
		if (_n_bins > 0) cout << " Error in TrpwaSessionManager::GetSelectedWaves(): Bin structure not initialized yet! " << endl;
	}
	return *result;
}

	// return a list of available Waves in the ith bin
	// Check_PWA_keyfiles() is called
vector<string> &TrpwaSessionManager::GetAvailableWaves(int ibin){
	vector<string>* result = new vector<string>();
	Check_PWA_keyfiles();
	(*result) = _keyfiles;
	return *result;
}

	// return the lower bound of the ith bin
int	TrpwaSessionManager::GetBinLow(int ibin){
	cout << " TrpwaSessionManager::GetBinLow() not implemented yet " << endl;
	return 0;
}

	// return the upper bound of the ith bin
int TrpwaSessionManager::GetBinHigh(int ibin){
	cout << " TrpwaSessionManager::GetBinHigh() not implemented yet " << endl;
	return 0;
}

vector<string>& TrpwaSessionManager::GetWaveList(int ibin, int& bin_low, int& bin_high){
	vector<string>* result = new vector<string>();
	int bincounter(0);
	for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
		if (bincounter == ibin){
			bin_low = (*it).second.bin_low;
			bin_high= (*it).second.bin_high;
			(*result) = (*it).second.wave_list; // copy the wave list
			break;
		}
		bincounter++;
	}
	return (*result);
}

string TrpwaSessionManager::GetWaveListFile(int ibin, int& bin_low, int& bin_high, string& fitresultfile){
	string result = "";
	int bincounter(0);
	for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
		if (bincounter == ibin){
			bin_low = (*it).second.bin_low;
			bin_high= (*it).second.bin_high;
			result = _dir_fit_results + "/" + (*it).second.wave_list_file;
			fitresultfile = _dir_fit_results + "/" + "fit_result_" + it->second.wave_list_file + ".root";
			break;
		}
		bincounter++;
	}
	return result;
}

string TrpwaSessionManager::GetFitCommand(int ibin, string& executedir){
	string result = "";
	string wavelistfile = "";
	string normalizationfile = "";
	string fitresultfile = "";
	int bincounter(0);
	int bin_low = -1;
	int bin_high = -1;
	for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
		if (bincounter == ibin){
			bin_low = (*it).second.bin_low;
			bin_high= (*it).second.bin_high;
			wavelistfile = _dir_fit_results + "/" + (*it).second.wave_list_file;
			fitresultfile = _dir_fit_results + "/" + "fit_result_" + it->second.wave_list_file + ".root";
			executedir = _dir_binned_data + "/" + it->second.bin_folder_name + "/" + "AMPS/";
			normalizationfile = _dir_binned_data + "/" + it->second.bin_folder_name + "/" + "PSPAMPS/" + "norm.int";
			break;
		}
		bincounter++;
	}
	if (!FileExists(wavelistfile)){
		cout << " Error in TrpwaSessionManager::GetFitCommand(): wave list file " << wavelistfile << " does not exist! " << endl;
		return result;
	}
	if (FileExists(fitresultfile)){
		cout << " Warning in TrpwaSessionManager::GetFitCommand(): fit result " << fitresultfile << " exist already! " << endl;
	}
	if (!FileExists(normalizationfile)){
		cout << " Error in TrpwaSessionManager::GetFitCommand(): normalization integral file " << normalizationfile << " does not exist! " << endl;
		return result;
	}
	stringstream _result;
	_result <<  " pwafit -q -w " + wavelistfile + " -o " + fitresultfile + " -r 1 " + " -l " << bin_low << " -N -u " << bin_high << " -n " + normalizationfile;
	result = _result.str();
	return result;
}

vector<string>& TrpwaSessionManager::GetWaveList(int ibin){
	int bin_low;
	int bin_high;
	GetWaveList(ibin, bin_low, bin_high);
}

// check whether a file exists
bool TrpwaSessionManager::FileExists(string filename){
	  ifstream ifile(filename.c_str());
	  return ifile;
}

// check whether a directory exists
bool TrpwaSessionManager::DirExists(string dirname){
	struct stat st;
	if(stat(dirname.c_str(),&st) == 0)
		return true;
	else
		return false;
}

int TrpwaSessionManager::GetDir (string path,
		vector<string> &files, string filterext, bool rmext){
    DIR *dp;
    struct dirent *dirp;
    files.clear();
    if((dp  = opendir(path.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << path << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
    	string _filename = dirp->d_name;
    	if ((_filename.size() <= filterext.size())||
    		(filterext != "" && _filename.compare(_filename.size()-filterext.size(), filterext.size(), filterext) != 0)) continue;
    	if (rmext) _filename.erase(_filename.size()-filterext.size(), filterext.size());
    	files.push_back(_filename);
    }
    closedir(dp);
    return (signed) files.size();
}

bool TrpwaSessionManager::AreListsEqual(const vector<string>& list1, const vector<string>& list2){
	if (list1.size() != list2.size()) return false;
	// prepare a list for ereasing
	vector<string> _templist = list2;
	// compare every entry of one list with every entry of the other and remove
	// same entries till nothing is left or all elements had been checked
	for (unsigned int i = 0; i < list1.size(); i++){
		bool found(false);
		for (vector<string>::iterator j = _templist.begin(); j != _templist.end(); ++j){
			if (list1[i] == *j){
				_templist.erase(j);
				found = true;
				break;
			}
		}
		// was one element found? else no chance to continue
		if (!found) break;
	}
	// only if all elements were found list2 will be empty
	if (_templist.size() == 0) return true; else return false;
}

vector<string> &TrpwaSessionManager::ReadWaveList(string filename){
	cout << " reading wave list " << filename << endl;
	vector<string>* result = new vector<string>();

	// Check whether the file exists, if not create it if possible
	if (!FileExists(filename)){
		Check_PWA_keyfiles();
		cout << " wave list " << filename << " does not exist: creating a new one with ";
		cout << _keyfiles.size() << " key files " << endl;
		ofstream _file(filename.c_str());
		if (_file.good()){
			for (unsigned int i = 0; i < _keyfiles.size(); i++){
				_file << _keyfiles[i] << ".amp" << endl;
			}

		} else {
			cout << " Error in TrpwaSessionManager::ReadWaveList(): Could not create a new wave list -> wave list ";
			cout << filename << " is missing! " << endl;
			return *result;
		}
		_file.close();
	}
	// read the written wave list
	ifstream _file(filename.c_str());
	if (_file.good()){
		char* line = new char[1024];
		while(1){
			_file.getline(line, 1024);
			if (!_file.good()) break;
			string _wave = line;
			// remove empty spaces
			for (unsigned int i=0;i<_wave.length();i++)
				if (_wave[i]==' ') {
					_wave.erase(i,1);
					i--;
				}
			// skip waves with # in it
			if (_wave.find('#')!=string::npos){
				cout << " omitting line: " << _wave << endl;
				continue;
			}
			// take only entries with .amp ending
			if (_wave.size() < 4 || _wave.compare(_wave.size()-4, 4, ".amp") != 0){
				cout << " omitting line with no .amp ending " << endl;
				continue;
			}
			// remove the .amp ending
			_wave.erase(_wave.size()-4, 4);
			result->push_back(_wave);
		}
	} else {
		cout << " Error in TrpwaSessionManager::ReadWaveList(): Could not read wave list "<< filename;
	}
	cout << " found " << result->size() << " waves " << endl;

	return *result;
}

bool TrpwaSessionManager::SetSelectedWaves(const vector<int>& bin_lowedes,const vector< vector<string> >& waves){
	bool result(true);	
	// perform some checks
	if (bin_lowedes.size() != waves.size()){
		return false;
	}
	// find the corresponding waves and set the status
	for (unsigned int i = 0; i < waves.size(); i++){
		// check if the requested bin is existing
		if (_bins.find(bin_lowedes[i]) == _bins.end()){
			cout << " Error in TrpwaSessionManager::SetSelectedWaves(): bin key ";
			cout << bin_lowedes[i] << " does not exist " << endl;
			continue;
		}
		vector < string >& _wavelist   = _bins[bin_lowedes[i]].wave_list;
		_wavelist.clear();
		for (unsigned int iwave = 0 ; iwave < waves[i].size(); iwave++){
			// search for the matching string in the key files and fill the empty list
			bool found(false);
			for (unsigned int jwave = 0; jwave < _keyfiles.size(); jwave++){
				if (_keyfiles[jwave] == waves[i][iwave]){
					found = true;
					// add the wave if requested
					_wavelist.push_back(_keyfiles[jwave]);
					break;
				}
			}
			if (!found){
				cout << " Error in TrpwaSessionManager::SetSelectedWaves(): no matching wave for ";
				cout << waves[i][iwave] << " in bin " << bin_lowedes[i] << " in the keyfiles found ";
				result = false;
			}
		}
	}
	//if (! SaveSelectedWaves()) return false;
	return result;
}

bool TrpwaSessionManager::SaveSelectedWaves(){
	bool result(true);

	// write for every bin the keyfiles out
	// add an # in front of the wave key if it is not specified as an selected wave

	for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
		vector < string > _wavelist   = it->second.wave_list; // copy the wave list
		string filename = _dir_fit_results + "/" + it->second.wave_list_file;		
		ofstream _file(filename.c_str());
		if (_file.good()){
			cout << " bin " << it->first << " nwaves " << _wavelist.size() << endl;		
			for (unsigned int iwave = 0; iwave < _keyfiles.size(); iwave++){
				// search for the wave in the selected waves to know whether it occurs
				bool found(false);
				for (vector<string>::iterator jwave = _wavelist.begin(); jwave !=_wavelist.end(); jwave++){
					if ((*jwave) == _keyfiles[iwave]){
						found = true;
						// erase this element. At the end no element should be left
						_wavelist.erase(jwave);
						break;
					}					
				}
				if (!found) _file << "# ";
				_file << _keyfiles[iwave] << ".amp" << endl;
			} 
		} else {
				cout << " Error in TrpwaSessionManager::SaveSelectedWaves: could not write wave list " << filename << endl;
		}
		_file.close();
		if (_wavelist.size() > 0){
			cout << " Error in TrpwaSessionManager::SaveSelectedWaves: wave list of bin " << it->first << " contains unsaved keys " << endl;
			for (unsigned int iwave = 0; iwave < _wavelist.size(); iwave++){
				cout << _wavelist[iwave] << endl;			
			}
			cout << endl;
		} else {
			result = true;
		}
	}

	return result;
}

// sort by JPC iso1 iso2 M reflectivity
void TrpwaSessionManager::SortWaves(vector<string>& wavelist){
	//int J,P,C,M,refl,l,s;
	//string iso1,iso2;
	sort(wavelist.begin(), wavelist.end());
	//cout << " sorted " << wavelist.size() << " waves " << endl;
}

// get the corresponding variables to the coded wavename
void TrpwaSessionManager::GetJPCMreflISO1lsISO2(string wavename, int& J, int& P, int& C, int& M, int& refl, string& iso1, string& iso2, int& l, int& s){
	cout << " TrpwaSessionManager::GetJPCMreflISO1lsISO2() not implemented yet " << endl;
	for (unsigned int i = 0; i < wavename.size(); i++){
		cout << wavename[i] << " ";
	}
	cout << endl;
}

