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
			file << " bins_lowedge = [ ";
			// write also the individual bins down
			for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
				if (it != _bins.begin()) file << " , ";
				file << it->second.bin_low;
			}
			file << " ]; " << endl;

			file << " bins_highedge = [ ";
			// write also the individual bins down
			for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
				if (it != _bins.begin()) file << " , ";
				file << it->second.bin_high;
			}
			file << " ]; " << endl;

			file << " bins_wavelist = [ ";
			// write also the individual bins down
			for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
				if (it != _bins.begin()) file << " , ";
				file << "\"" << it->second.wave_list_file << "\"";
			}
			file << " ]; " << endl;

			file << " bins_foldername = [ ";
			// write also the individual bins down
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

