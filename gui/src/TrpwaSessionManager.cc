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
#include <complex>

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
	jobManager = TrpwaJobManager::Instance();
	Set_ROOTPWA_dir();
}

TrpwaSessionManager::~TrpwaSessionManager(){

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
		file << " pdg_table = \"" << _pdg_table << "\";" << endl;
		file << " dir_binned_data = \"" << _dir_binned_data << "\";" << endl;

		file << " file_data = [";
		for (unsigned int i = 0; i < _data_files.size(); i++){
			file << "\"" << _data_files[i] << "\"";
			if (i != _data_files.size()-1) file << " ,\n ";
		}
		file << "]; " << endl;

		file << " file_mc_data = [";
		for (unsigned int i = 0; i < _mc_data_files.size(); i++){
			file << "\"" << _mc_data_files[i] << "\"";
			if (i != _mc_data_files.size()-1) file << " ,\n ";
		}
		file << "]; " << endl;

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

	if (!_config.lookupValue("pdg_table", _val_string)){
		if (!Set_pdg_table("")) result = false;
	} else {
		if (!Set_pdg_table("")) result = false;
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

	_data_files.clear();
	if ( _config.exists("file_data")){
		Setting& _setting_file_data    = _config.lookup("file_data");
		if (!_setting_file_data.isArray()){
			cout << " Error in TrpwaSessionManager::Load_Session(): check files with data in config file! " << endl;
			return false;
		}
		int nfiles = _setting_file_data.getLength();
		if (0 == nfiles){
			cout << " Warning in TrpwaSessionManager::Load_Session(): no data files specified " << endl;
		}
		for (int i = 0; i < nfiles; i++){
			string _file = (const char*) _setting_file_data[i];
			if (FileExists(_file)){
				_data_files.push_back(_file);
			} else {
				cout << " Error in TrpwaSessionManager::Load_Session(): File does not exist: " << _file << endl;
				result = false;
			}
		}
	}

	_mc_data_files.clear();
	if ( _config.exists("file_mc_data")){
		Setting& _setting_file_data    = _config.lookup("file_mc_data");
		if (!_setting_file_data.isArray()){
			cout << " Error in TrpwaSessionManager::Load_Session(): check files with mc data in config file! " << endl;
			return false;
		}
		int nfiles = _setting_file_data.getLength();
		if (0 == nfiles){
			cout << " Warning in TrpwaSessionManager::Load_Session(): no mc data files specified " << endl;
		}
		for (int i = 0; i < nfiles; i++){
			string _file = (const char*) _setting_file_data[i];
			if (FileExists(_file)){
				_mc_data_files.push_back(_file);
			} else {
				cout << " Error in TrpwaSessionManager::Load_Session(): File does not exist: " << _file << endl;
				result = false;
			}
		}
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

bool TrpwaSessionManager::Set_pdg_table(string filename){
	bool result(false);
	string _filename = filename;
	if (_filename == "" || !FileExists(_filename)){
		cout << " searching for pdg table in ${ROOTPWA}: " << _dir_ROOTPWA << endl;
		_filename = _dir_ROOTPWA + "/keygen/pdgTable.txt";
		if (FileExists(_filename)){
			cout << " using " << _filename << " as PDG table input " << endl;
			result = true;
		} else {
			cout << " Error in TrpwaSessionManager::Set_pdg_table: no valid PDG table set! " << endl;
			_filename = "";
		}
	}
	_pdg_table = _filename;
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
float TrpwaSessionManager::Check_binned_data_structure(bool create){
	cout << " checking binned data structure in " << _dir_binned_data << endl;
	if (create)
		cout << " creating missing folders " << endl;
	float result(0.);

	if (!DirExists(_dir_binned_data)){
		cout << " Error: Bin folder " << _dir_binned_data << " does not exist!" << endl;
		return 0.;
	}
	if ( _n_bins <= 0 ) return result;
	int counter(0);
	for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
		bool complete(true);
		string _bindir = _dir_binned_data + "/" + it->second.bin_folder_name;
		string _checkdir = _bindir;
		if (!DirExists(_checkdir)){
			if (create){
				if((mkdir(_checkdir.c_str(),00777))==-1) {
					cout << " Error: could not create " << _checkdir << endl;
					return 0.;
				}
			} else {
				complete = false;
			}
		}
		_checkdir = _bindir + "/AMPS";
		if (!DirExists(_checkdir)){
			if (create){
				if((mkdir(_checkdir.c_str(),00777))==-1) {
					cout << " Error: could not create " << _checkdir << endl;
					return 0.;
				}
			} else {
				complete = false;
			}
		}
		_checkdir = _bindir + "/PSPAMPS";
		if (!DirExists(_checkdir)){
			if (create){
				if((mkdir(_checkdir.c_str(),00777))==-1) {
					cout << " Error: could not create " << _checkdir << endl;
					return 0.;
				}
			} else {
				complete = false;
			}
		}
		_checkdir = _bindir + "/ACCAMPS";
		if (!DirExists(_checkdir)){
			if (create){
				if((mkdir(_checkdir.c_str(),00777))==-1) {
					cout << " Error: could not create " << _checkdir << endl;
					return 0.;
				}
			} else {
				complete = false;
			}
		}
		_checkdir = _bindir + "/MC";
		if (!DirExists(_checkdir)){
			if (create){
				if((mkdir(_checkdir.c_str(),00777))==-1) {
					cout << " Error: could not create " << _checkdir << endl;
					return 0.;
				}
			} else {
				complete = false;
			}
		}
		if (complete) counter++;
	}
	if ((int) counter != (int) _n_bins)
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
	if ((int) counter != (int) _n_bins)
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
	if ((int) counter != (int) _n_bins)
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
	if ((int) counter != (int) _n_bins)
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

	//cout << " found " << _n_keyfiles << " .key files " << endl;

	return result;
}

string TrpwaSessionManager::Get_PWA_data_integral(int bin, string folder, bool missing,
	vector<string>* corresponding_amplitudefiles){
	string result = "";
	if (corresponding_amplitudefiles) corresponding_amplitudefiles->clear();	
	if (bin < 0 || bin >= _n_bins){
			return result;
	}
	int ibin(0);
	//cout << " bin " << bin << endl;
	for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); it++){
		if (ibin != bin){
			ibin++;			
			continue;
		}
		string _dir = _dir_binned_data + "/" + it->second.bin_folder_name + "/" + folder + "/";
		//cout << " searching in " << _dir << endl;
		vector<string> _ints;
		GetDir(_dir, _ints, ".int", false);
		// append the directory path and fill it to the result
		bool found(false);
		for (vector<string>::iterator _it = _ints.begin(); _it != _ints.end(); _it++){
			//if (!missing) // accept all files with the ending .int as not missing input
			//	result->push_back(_dir+(*_it));
			if ((*_it)=="norm.int"){ // requirer norm.int
				found = true;
			}
		}
		// write out norm.int if search for it was not successful
		if ((missing && !found) || (!missing && found)){
			result = _dir+"norm.int";
			if (corresponding_amplitudefiles){
				// retrieve the available amplitude files
				GetDir(_dir, (*corresponding_amplitudefiles), ".amp", false);
				// add the full path to it
				for (vector<string>::iterator _itamps = corresponding_amplitudefiles->begin(); _itamps != corresponding_amplitudefiles->end(); _itamps++){
					(*_itamps) = _dir+(*_itamps);
				}
			}
		}
		break;
	}	
	//if (result=="") cout << " found " << endl; else cout << " not found " << endl;
	return result;
}

vector<string>& TrpwaSessionManager::Get_PWA_data_integrals(string folder, bool missing,
		vector < vector<string> >* corresponding_amplitudefiles){
	vector<string>* result = new vector<string>();
	// create a vector for corresponding_amplitudefiles if not given
	// since it is needed internally for consistency checks
	bool delete_corresponding_ampsfiles(false);
	if(!corresponding_amplitudefiles){
		corresponding_amplitudefiles = new vector < vector<string> >();
		delete_corresponding_ampsfiles = true;
	} else {
		corresponding_amplitudefiles->clear();
	}
	cout << " searching for integral files in "<< folder <<" folders " << endl;
	if ( _n_bins <= 0 ) return *result;
	for (int ibin = 0; ibin < _n_bins; ibin++){
		vector<string> _corresponding_amplitudefiles;	
		string integralfile = Get_PWA_data_integral(ibin, folder, missing, &_corresponding_amplitudefiles);
		if (integralfile != ""){
			result->push_back(integralfile);
			corresponding_amplitudefiles->push_back(_corresponding_amplitudefiles);
		}
	}
/*
	for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
		string _dir = _dir_binned_data + "/" + it->second.bin_folder_name + "/" + folder + "/";
		vector<string> _ints;
		GetDir(_dir, _ints, ".int", false);
		// append the directory path and fill it to the result
		bool found(false);
		bool push_back_amplitudefiles(false);
		for (vector<string>::iterator _it = _ints.begin(); _it != _ints.end(); _it++){
			//if (!missing) // accept all files with the ending .int as not missing input
			//	result->push_back(_dir+(*_it));
			if ((*_it)=="norm.int"){ // requirer norm.int
				found = true;
				result->push_back(_dir+(*_it));
				push_back_amplitudefiles = true;
			}
		}
		// write out norm.int if search for it was not successful
		if (missing && !found){
			result->push_back(_dir+"norm.int");
			push_back_amplitudefiles = true;
		}
		if (push_back_amplitudefiles){
			// retrieve the available amplitude files
			vector<string> _amps;
			GetDir(_dir, _amps, ".amp", false);
			// add the full path to it
			for (vector<string>::iterator _itamps = _amps.begin(); _itamps != _amps.end(); _itamps++){
				(*_itamps) = _dir+(*_itamps);
			}
			// store it to the output
			corresponding_amplitudefiles->push_back(_amps);
		}
	}
*/
	// check now the consistency of this result
	//cout << " integrals before check " << result->size() << " missing " << missing << endl;
	if (corresponding_amplitudefiles->size() != result->size()){
		cout << " Error in TrpwaSessionManager::Get_PWA_data_integrals: consistency of result not ensured! " << endl;
		cout << corresponding_amplitudefiles->size() << " does not fit " << result->size() << endl;
		corresponding_amplitudefiles->clear();
		result->clear();
	} else {
		vector< vector<string> >::iterator it_pair = corresponding_amplitudefiles->begin();
		for (vector<string>::iterator it = result->begin(); it != result->end(); ++it){
			if (!missing&&!Is_valid_norm_integral((*it), (*it_pair))){
				result->erase(it);
				corresponding_amplitudefiles->erase(it_pair);
				it--;
				it_pair--;
			}
			it_pair++;
		}
	}
	//cout << " integrals after check " << result->size() << " missing " << missing << endl;
	if (corresponding_amplitudefiles->size() != result->size()){ cout << " noooooooo ! " << endl;}
	if (delete_corresponding_ampsfiles) {
		delete corresponding_amplitudefiles;
		corresponding_amplitudefiles = NULL;
	}
	return *result;
}

float TrpwaSessionManager::Check_PWA_MC_acc_data_integrals(){
	float result(0.);
	result = Get_PWA_data_integrals("ACCAMPS",false).size();
	if ((int) result != (int) _n_bins)
		cout << " found " << (int) result << " of " << _n_bins << " expected .int files" << endl;
	result = (double) result / ((double) _n_bins);
	return result;
}

vector<string>& TrpwaSessionManager::Get_PWA_MC_acc_data_integrals(bool missing,
		vector < vector<string> >* corresponding_amplitudefiles){
	return Get_PWA_data_integrals("ACCAMPS", missing, corresponding_amplitudefiles);;
}

float TrpwaSessionManager::Check_PWA_MC_data_integrals(){
	float result(0.);
	result = Get_PWA_data_integrals("PSPAMPS",false).size();
	if ((int) result != (int) _n_bins)
		cout << " found " << (int) result << " of " << _n_bins << " expected .int files" << endl;
	result = (double) result / ((double) _n_bins);
	return result;
}

vector<string>& TrpwaSessionManager::Get_PWA_MC_data_integrals(bool missing,
		vector < vector<string> >* corresponding_amplitudefiles){
	return Get_PWA_data_integrals("PSPAMPS", missing, corresponding_amplitudefiles);
}

vector<string>& TrpwaSessionManager::Get_PWA_data_amplitudes(string folder, bool missing,
		vector<string>* corresponding_eventfiles, vector<string>* corresponding_keyfiles){
	Check_PWA_keyfiles();
	vector<string>* result = new vector<string>(0);
	cout << " searching for amplitude files in "<< folder <<" folders " << endl;
	if ( _n_bins <= 0 ) return *result;
	for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
		string _dir = _dir_binned_data + "/" + it->second.bin_folder_name + "/" + folder + "/";
		vector<string> _amps;
		GetDir(_dir, _amps, ".amp", true);
		// retrieve the requested lists
		vector<string>& _temp_result = CompareLists(_keyfiles, _amps, missing);
		// append the directory path and fill it to the result
		for (vector<string>::iterator _it = _temp_result.begin(); _it != _temp_result.end(); _it++){
			result->push_back(_dir+(*_it)+".amp");
			if (corresponding_eventfiles){
				if (folder == "AMPS")
					corresponding_eventfiles->push_back(_dir_binned_data + "/" + it->second.bin_folder_name + "/" + it->second.bin_folder_name + ".evt");
				else
				if (folder == "PSPAMPS")
					corresponding_eventfiles->push_back(_dir_binned_data + "/" + it->second.bin_folder_name + "/" + it->second.bin_folder_name + ".genbod.evt");
				else
				if (folder == "ACCAMPS")
					corresponding_eventfiles->push_back(_dir_binned_data + "/" + it->second.bin_folder_name + "/" + it->second.bin_folder_name + ".acc.evt");
				else
					corresponding_eventfiles->push_back("folder_"+folder+"_not_known");
			}
			if (corresponding_keyfiles)
				corresponding_keyfiles->push_back(_dir_key_files +"/"+(*_it)+".key");
		}
	}
	return *result;
}

float TrpwaSessionManager::Check_PWA_real_data_amplitudes(){
	float result(0.);
	vector<string> amps = Get_PWA_data_amplitudes("AMPS",false);
	int i(0);
	for (vector<string>::iterator it = amps.begin(); it != amps.end(); it++){
		if (Is_valid_amplitude((*it))) result +=1. ;
		i++;
		DrawProgressBar(50, (i+1)/((double)amps.size()));
	}
	//result = amps.size();
	if ((int) result != (int) _n_bins* (int) _keyfiles.size())
		cout << " found " << (int) result << " of " << _n_bins*_keyfiles.size() << " expected .amp files" << endl;
	if (_n_bins*_keyfiles.size() == 0.){
		cout << " no amplitudes found! " << endl;
		result = 0.;
	} else {
		result = (double) result / ((double) _n_bins*_keyfiles.size());
	}
	return result;
	/*
	Check_PWA_keyfiles();
	cout << " searching for amplitude files in AMPS folders " << endl;


	if ( _n_bins <= 0 ) return result;
	int counter(0);
	for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
		string _dir = _dir_binned_data + "/" + it->second.bin_folder_name + "/" + "AMPS/";
		vector<string> _amps;
		GetDir(_dir, _amps, ".amp", true);
		if (AreListsEqual(_keyfiles, _amps)) {
			counter+=_amps.size();
		} else {
			counter += CompareLists(_amps, _keyfiles);
		}
		//if (FileExists(_dir_binned_data + "/" + it->second.bin_folder_name + "/" + "AMPS/")){

		//}
	}

	cout << " found " << counter << " of " << _n_bins*_keyfiles.size() << " expected .amp files" << endl;
	result = (double) counter / ((double) _n_bins*_keyfiles.size());
	return result;
	*/
}

vector<string>& TrpwaSessionManager::Get_PWA_real_data_amplitudes(bool missing,
		vector<string>* corresponding_eventfiles, vector<string>* corresponding_keyfiles){
	return Get_PWA_data_amplitudes("AMPS", missing, corresponding_eventfiles, corresponding_keyfiles);
}

vector<string>& TrpwaSessionManager::CompareLists(const vector<string>& list1, const vector<string>& list2, bool missing){
	vector<string>* result = new vector<string>();
	for (vector<string>::const_iterator i = list1.begin(); i != list1.end(); ++i){
		bool found(false);
		for (vector<string>::const_iterator j = list2.begin(); j != list2.end(); ++j){
			if (*i == *j){
				if (!missing) result->push_back(*i);
				found = true;
				break;
			}
		}
		if (!found && missing) result->push_back(*i);
	}
	return *result;
}


// returns the status [0-1] of calculated amplitudes
// of flat phase space data
// (comparing number of .amp files with .key files
// in the flat phase space data folder)
float TrpwaSessionManager::Check_PWA_MC_data_amplitudes(){
	float result(0.);
	vector<string> amps = Get_PWA_data_amplitudes("PSPAMPS",false);
	int i(0);
	for (vector<string>::iterator it = amps.begin(); it != amps.end(); it++){
		if (Is_valid_amplitude((*it))) result +=1. ;
		i++;
		DrawProgressBar(50, (i+1)/((double)amps.size()));
	}
	//result = Get_PWA_data_amplitudes("PSPAMPS",false).size();
	if ((int) result != (int) _n_bins* (int) _keyfiles.size())
		cout << " found " << (int) result << " of " << _n_bins*_keyfiles.size() << " expected .amp files" << endl;
	if (_n_bins*_keyfiles.size() == 0.){
		cout << " no amplitudes found! " << endl;
		result = 0.;
	} else {
		result = (double) result / ((double) _n_bins*_keyfiles.size());
	}
	return result;
	/*
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
		} else {
			counter += CompareLists(_amps, _keyfiles);
		}
		//if (FileExists(_dir_binned_data + "/" + it->second.bin_folder_name + "/" + "AMPS/")){

		//}
	}
	cout << " found " << counter << " of " << _n_bins*_keyfiles.size() << " expected .amp files" << endl;
	result = (double) counter / ((double) _n_bins*_keyfiles.size());

	return result;
	*/
}

vector<string>& TrpwaSessionManager::Get_PWA_MC_data_amplitudes(bool missing,
		vector<string>* corresponding_eventfiles, vector<string>* corresponding_keyfiles){
	return Get_PWA_data_amplitudes("PSPAMPS", missing, corresponding_eventfiles, corresponding_keyfiles);
}

// returns the status [0-1] of calculated amplitudes
// of accepted flat phase space data
// (comparing number of .amp files with .key files
// in the accpeted events data folder)
float TrpwaSessionManager::Check_PWA_MC_acc_data_amplitudes(){
	float result(0.);
	vector<string> amps = Get_PWA_data_amplitudes("ACCAMPS",false);
	int i(0);
	for (vector<string>::iterator it = amps.begin(); it != amps.end(); it++){
		if (Is_valid_amplitude((*it))) result +=1. ;
		i++;
		DrawProgressBar(50, (i+1)/((double)amps.size()));
	}
	//result = Get_PWA_data_amplitudes("ACCAMPS",false).size();
	if ((int) result != (int) _n_bins* (int) _keyfiles.size())
		cout << " found " << (int) result << " of " << _n_bins*_keyfiles.size() << " expected .amp files" << endl;
	if (_n_bins*_keyfiles.size() == 0.){
		cout << " no amplitudes found! " << endl;
		result = 0.;
	} else {
		result = (double) result / ((double) _n_bins*_keyfiles.size());
	}
	return result;
	/*
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
		} else {
			counter += CompareLists(_amps, _keyfiles);
		}
		//if (FileExists(_dir_binned_data + "/" + it->second.bin_folder_name + "/" + "AMPS/")){

		//}
	}
	cout << " found " << counter << " of " << _n_bins*_keyfiles.size() << " expected .amp files" << endl;
	result = (double) counter / ((double) _n_bins*_keyfiles.size());

	return result;
	*/
}

vector<string>& TrpwaSessionManager::Get_PWA_MC_acc_data_amplitudes(bool missing,
		vector<string>* corresponding_eventfiles, vector<string>* corresponding_keyfiles){
	return Get_PWA_data_amplitudes("ACCAMPS", missing, corresponding_eventfiles, corresponding_keyfiles);
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
	if ((int) counter != (int) _n_bins)
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
	if ((int) counter != (int) _n_bins)
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
	string acceptancefile = "";
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
			acceptancefile = _dir_binned_data + "/" + it->second.bin_folder_name + "/" + "ACCAMPS/" + "norm.int";
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
	if (!FileExists(acceptancefile)){
		cout << " Warning in TrpwaSessionManager::GetFitCommand(): acceptance integral file " << acceptancefile << " does not exist! " << endl;
		cout << " Supplying fit command without normalization! " << endl;
		//_result << " pwafit -q -w " + wavelistfile + " -o " + fitresultfile + " -r 1 " + " -l " << bin_low << " -N -u " << bin_high << " -n " + normalizationfile;
		_result <<  "pwafit -q -w " + wavelistfile + " -o " + fitresultfile + " -r 1 " + " -l " << bin_low << " -A " << _n_events_flat_phasespace << " -a "<< normalizationfile <<" -u " << bin_high << " -N -n " + normalizationfile;
		result = _result.str();
		return result;
	} else {
		//_result <<  " pwafit -q -w " + wavelistfile + " -o " + fitresultfile + " -r 1 " + " -l " << bin_low << " -N -u " << bin_high << " -n " + normalizationfile;
		_result <<  "pwafit -q -w " + wavelistfile + " -o " + fitresultfile + " -r 1 " + " -l " << bin_low << " -A " << _n_events_flat_phasespace << " -a "<< acceptancefile <<" -u " << bin_high << " -N -n " + normalizationfile;
		result = _result.str();
		return result;
	}
}

string TrpwaSessionManager::GetgenpwCommand(int ibin, string& executedir){
	// genpw -n ${KPIPI_NPHASESPACE_MC_EVENTS} -M ${BINLOW} -B ${BINWIDTH} -c -r ${KPIPI_MC_CONFIG_FILE_FLAT}
	string result = "";
	int bincounter(0);
	int bin_low = -1;
	int bin_high = -1;
	for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
		if (bincounter == ibin){
			bin_low = (*it).second.bin_low;
			bin_high= (*it).second.bin_high;
			executedir = _dir_binned_data + "/" + it->second.bin_folder_name;
			break;
		}
		bincounter++;
	}
	stringstream _result;
	_result <<  "genpw -n " << _n_events_flat_phasespace << " -M " << bin_low << " -B " << bin_high-bin_low << " -c -r " << _file_flat_phasespace_config;
	result = _result.str();
	return result;
}

vector<string>& TrpwaSessionManager::GetFitResults(){
	vector<string>* result = new vector<string>();
	for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
		string fitresultfile;
		fitresultfile = _dir_fit_results + "/" + "fit_result_" + it->second.wave_list_file + ".root";
		if (FileExists(fitresultfile)) (*result).push_back(fitresultfile);
		else cout << " Warning in TrpwaSessionManager::GetFitResults(): fit result " << fitresultfile << " does not exist. Skipping " << endl;
	}
	return *result;
}

vector<string>& TrpwaSessionManager::GetWaveList(int ibin){
	int bin_low;
	int bin_high;
	return GetWaveList(ibin, bin_low, bin_high);
}

// check whether a file exists
bool TrpwaSessionManager::FileExists(string filename){
	  ifstream ifile(filename.c_str());
	  return ifile;
}

// check whether a directory exists
bool TrpwaSessionManager::DirExists(string dirname){
	struct stat st;
	if(stat(dirname.c_str(),&st) == 0 && S_ISDIR(st.st_mode))
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

int TrpwaSessionManager::CompareLists(const vector<string>& list1, const vector<string>& list2){
	int result(0);
	vector<string> _templist = list2;
	for (unsigned int i = 0; i < list1.size(); i++){
		bool found(false);
		for (vector<string>::iterator j = _templist.begin(); j != _templist.end(); ++j){
			if (list1[i] == *j){
				result++;
				_templist.erase(j);
				found = true;
				break;
			}
		}
	}
	return result;
}

bool TrpwaSessionManager::Is_valid_norm_integral(const string norm_file, const vector<string> & corresponding_amplitudefiles){
	bool result = true;
	if (!FileExists(norm_file)){
		cout << " Error in TrpwaSessionManager::Is_valid_norm_integral: " << norm_file << " does not exist " << endl;
		return false;
	}

	ifstream integralfile(norm_file.c_str());
	int namps(0);
	int nevents(0);
	int dimensionx(0);
	int dimensiony(0);

	integralfile >> namps;
	integralfile >> nevents;
	integralfile >> dimensionx >> dimensiony;

	//cout << norm_file << ": " << endl;
	//cout << " namps " << namps << " nevents " << nevents << endl;

	if (namps != (int) corresponding_amplitudefiles.size()){
		cout << " Inonsystency of " << norm_file << " dimension is wrong " << endl;
		cout << " expected " << corresponding_amplitudefiles.size() << " amplitudes but got " << namps << endl;
		integralfile.close();
		return false;
	}
	if (namps != dimensionx || namps != dimensiony){
		cout << " Warning in TrpwaSessionManager::Is_valid_norm_integral: " << norm_file << " dimensions do not fit the number of amplitudes." << endl;
	}

	// search the matrix elements
	complex<double> matrix_element;

	// store broken col/row entries
	vector<int> broken_entries;

	// count 0 entries in cols and rows
	vector<int> zero_entries_x(dimensionx);
	vector<int> zero_entries_y(dimensiony);

	// read the matrix
	for (int i = 0; i < dimensionx; i++){
		for (int j = 0; j < dimensiony; j++){
			integralfile >> matrix_element;
			if (!integralfile.good()){
				cout << " Inconsistency of " << norm_file << " found: " << endl;
				cout << " matrix is incomplete! " << endl;
				return false;
			}
			// no 0 entries in the diagonal allowed!
			if (matrix_element.real() == 0. && matrix_element.imag() == 0.){
				if (i == j){
					//cout << " Inconsistency of " << norm_file << " found: " << endl;
					//cout << " 0 entry in diagonal element! row " << i << " col " << j << " " << matrix_element << endl;
					result = false;
					broken_entries.push_back(i); // store the broken entry for identification afterwards
				}
				zero_entries_x[i]++;
				zero_entries_y[i]++;
			}
		}
	}

	// search for amplitudes in the file
	vector<string> integrated_ampslist;

	char* line = new char[100000];
	while(1) {
		integralfile.getline(line, 100000);
		if (!integralfile.good()) break;
		// search for amplitudes
		string amplitudename(line);
		if (amplitudename.find(".amp")!=string::npos){
			int pos_amp = amplitudename.find(".amp");
			// get the numbers behind the .amp name
			string snumber;
			if (pos_amp+5 < (int) amplitudename.size()){
				snumber = amplitudename.substr(pos_amp+5, amplitudename.size()-pos_amp);
			} else {
				cout << " Not a valid amplitude position! " << endl;
				result = false;
				snumber = "-1";
			}
			int number = atoi(snumber.c_str());
			// remove the .amp ending
			amplitudename.erase(pos_amp, amplitudename.size()-pos_amp);
			// find the slash, if not available -> take the first position of the string
			int slashpos = amplitudename.rfind('/');
			if (slashpos != (int) string::npos){
				amplitudename.erase(0, slashpos+1);
			}
			//cout << " filtered: " << amplitudename << " pos " << snumber << endl;
			integrated_ampslist.push_back(amplitudename);
			// give amplitudes corresponding to entries with problems
			for (vector<int>::iterator it = broken_entries.begin(); it != broken_entries.end(); it++){
				if (number == (*it)){
					// Suppress output for amplitudes with known problems
					if (_keyfiles_blacklist.find(amplitudename) == _keyfiles_blacklist.end()){
						cout << " wave with broken entries: " << amplitudename << endl;
					}
					_keyfiles_blacklist[amplitudename].push_back(norm_file);
					break;
				}
			}
			if (zero_entries_x[number] == dimensionx){
				if (_keyfiles_blacklist.find(amplitudename) == _keyfiles_blacklist.end()){
					cout << " all x entries are 0! for " << amplitudename << endl;
				}
				_keyfiles_blacklist[amplitudename].push_back(norm_file);
				result = false;
			}
			if (zero_entries_y[number] == dimensiony){
				if (_keyfiles_blacklist.find(amplitudename) == _keyfiles_blacklist.end()){
					cout << " all y entries are 0! for " << amplitudename << endl;
				}
				_keyfiles_blacklist[amplitudename].push_back(norm_file);
				result = false;
			}
		}
	}

	// create a list of amplitude files for comparison
	vector<string> ampslist;
	for (vector<string>::const_iterator it = corresponding_amplitudefiles.begin();
			it != corresponding_amplitudefiles.end() ; it++){
		string amplitudename = (*it);
		// find the slash, if not available -> take the first position of the string
		int slashpos = amplitudename.rfind('/');
		if (slashpos != (int) string::npos){
			amplitudename.erase(0, slashpos+1);
		}
		// find the .amp ending and erase it
		if (amplitudename.find(".amp")!=string::npos){
			int pos_amp = amplitudename.find(".amp");
			amplitudename.erase(pos_amp, amplitudename.size()-pos_amp);
		}
		//cout << " filtered in amps list: " << amplitudename << endl;
		ampslist.push_back(amplitudename);
	}

	if (ampslist.size() == integrated_ampslist.size() && AreListsEqual(ampslist, integrated_ampslist)){
		//result = true;
	} else {
		result = false;
		cout << " Inconsistency of " << norm_file << " processed amplitudes do not match! " << endl;
	}
	integralfile.close();

	return result;
}

bool TrpwaSessionManager::Is_valid_amplitude(string ampfile, int nlines)
{
	bool result = true;
	ifstream file(ampfile.c_str());
	if (!file) {
		cout << " Error: " << ampfile << " does not exist!" << endl;
		return false;
	}
	int i(0);
	while (1)
	{
		complex<double> value;
		//file >> value;
		file.read((char*)&value, sizeof(complex<double>));
		//cout << value << endl;
		if (!file.good()){
			cout << " Error reading " << ampfile << " at line " << i << endl;
			result = false;
			break;
		}
		if (i > 10) break;
		if (value.real() == 0. && value.imag() == 0.) {
			// display once the error for a 0 entry
			if (result){
				string key = Get_key(ampfile);
				// Suppress output for amplitudes with known problems
				if (_keyfiles_blacklist.find(key) == _keyfiles_blacklist.end()){
					cout << " Amplitudefile " << ampfile << " has a (0, 0) entry within the first " << nlines << endl;
				}
				_keyfiles_blacklist[key].push_back(ampfile);
			}
			result = false;
		}
		i++;
	}
	//cout << " read " << i << " lines " << endl;
	file.close();
	return result;
}

vector<string> &TrpwaSessionManager::ReadWaveList(string filename){
	//cout << " reading wave list " << filename << endl;
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
				//cout << " omitting line: " << _wave << endl;
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
	//cout << " found " << result->size() << " waves " << endl;

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

	// write for every bin the key files out
	// add an # in front of the wave key if it is not specified as an selected wave

	for( TBinMap::const_iterator it = _bins.begin(); it != _bins.end(); ++it){
		vector < string > _wavelist   = it->second.wave_list; // copy the wave list
		string filename = _dir_fit_results + "/" + it->second.wave_list_file;		
		ofstream _file(filename.c_str());
		if (_file.good()){
			//cout << " bin " << it->first << " nwaves " << _wavelist.size() << endl;
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

string TrpwaSessionManager::Get_key(string file){
	// remove the extension if available
	if (file.find(".amp")!=string::npos){
		int pos_amp = file.find(".amp");
		// remove the .amp ending
		file.erase(pos_amp, file.size()-pos_amp);
	}
	// remove the extension if available
	if (file.find(".key")!=string::npos){
		int pos_key = file.find(".key");
		// remove the .amp ending
		file.erase(pos_key, file.size()-pos_key);
	}
	// remove the path if given
	int slashpos = file.rfind('/');
	if (slashpos != (int) string::npos){
		file.erase(0, slashpos+1);
	}
	return file;
}

void TrpwaSessionManager::Print_problematic_waves(){
	cout << endl;
	cout << " ********************************************************* " << endl;
	cout << " *                                                       *"  << endl;
	cout << " *                       wave report                     * " << endl;
	cout << " *                                                       *"  << endl;
	cout << " ********************************************************* " << endl;
	cout << endl;
	for (vector<string>::iterator it = _keyfiles.begin(); it != _keyfiles.end(); it++){
		cout << " " << *it << " \t ";
		if (_keyfiles_blacklist.find(*it) != _keyfiles_blacklist.end()){
			cout << " gave " << _keyfiles_blacklist.find(*it)->second.size() << " times errors " << endl;
		} else {
			cout << " has no problems " << endl;
		}
	}
	cout << endl;
	cout << " ********************************************************* " << endl;
	cout << " *                                                       *"  << endl;
	cout << " *               list of problematic waves               * " << endl;
	cout << " *                                                       *"  << endl;
	cout << " ********************************************************* " << endl;
	for (map<string, vector<string> >::iterator it = _keyfiles_blacklist.begin(); it != _keyfiles_blacklist.end(); it++){
		cout << " " << it->first << endl;
	}
	cout << endl;
}

// get the corresponding variables to the coded wavename
void TrpwaSessionManager::GetJPCMreflISO1lsISO2(string wavename, int& J, int& P, int& C, int& M, int& refl, string& iso1, string& iso2, int& l, int& s){
	cout << " TrpwaSessionManager::GetJPCMreflISO1lsISO2() not implemented yet " << endl;
	for (unsigned int i = 0; i < wavename.size(); i++){
		cout << wavename[i] << " ";
	}
	cout << endl;
}

void TrpwaSessionManager::DrawProgressBar(int len, double percent) {
	  static double last_percent(0.);
	  if ((int)(last_percent*100) == (int)(percent*100)) return;
	  //cout << " drawing " << endl;
	  cout << "\x1B[2K"; // Erase the entire current line.
	  cout << "\x1B[0E"; // Move to the beginning of the current line.
	  string progress;
	  for (int i = 0; i < len; ++i) {
		if (i < static_cast<int>(len * percent)) {
		  progress += "=";
		} else {
		  progress += " ";
		}
	  }
	  cout << "[" << progress << "] " << (static_cast<int>(100 * percent)) << "%";
	  flush(cout); // Required.
	  last_percent = percent;
}
