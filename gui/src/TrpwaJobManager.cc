/*
 * TrpwaJobManager.cc
 *
 *  Created on: Oct 12, 2010
 *      Author: Promme
 *
  *      (12.10.10)
 *      - Implementation of singleton class
 */

#include "TrpwaJobManager.h"
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iostream>

TrpwaJobManager* TrpwaJobManager::_pinstance = NULL;

TrpwaJobManager::TrpwaJobManager(){
	_available_farmtype = arrfarmtypes[0]; // set to local
	cout << system("echo $(pwd)/rootpwa_gui_scripts > /tmp/_pwd") << endl;
	ifstream _filestream("/tmp/_pwd");
	_filestream >> _temp_space;
	_filestream.close();
	cout << system(("mkdir "+_temp_space).c_str()) << endl;
	cout << " using " << _temp_space << " as farm script submission " << endl;
	CheckFarmType();
}

TrpwaJobManager::~TrpwaJobManager(){
	if (_pinstance){
		delete _pinstance;
	}
}

TrpwaJobManager* TrpwaJobManager::Instance(){
	if (!_pinstance){
		_pinstance = new TrpwaJobManager();
	}
	return _pinstance;
}

void TrpwaJobManager::CheckFarmType(){
	cout << " Is Mainz blaster available? ";
	flush(cout);
	{ 	// is it a Mainz Blaster login node?
		// search for "kph" after executing qnodes request
		stringstream command;
		command << "qnodes | grep kph > /tmp/_blaster_test;";
		cout << system(command.str().c_str()) << endl;
		ifstream testfile("/tmp/_blaster_test");
		if (testfile){
			while (1) {
				char* line = new char[1024];
				testfile.getline(line, 1024);
				if (!testfile.good()) break;
				string sline(line);
				if (sline.find("kph")!=string::npos){
					_available_farmtype = arrfarmtypes[3];
				}
			}
		}
	}
	if (_available_farmtype == arrfarmtypes[3]) cout << " yes! " << endl; else cout << " no " << endl;
	cout << " Is lxbatch farm available? ";
	flush(cout);
	{
		// is it there an lxbatch farm available ?
		// search for "atlas" when executing "bqueues" command
		// not tested yet!
		stringstream command;
		command << "bqueues | grep atlas > /tmp/_lxplus_test;";
		cout << system(command.str().c_str()) << endl;
		ifstream testfile("/tmp/_lxplus_test");
		if (testfile){
			while (1) {
				char* line = new char[1024];
				testfile.getline(line, 1024);
				if (!testfile.good()) break;
				string sline(line);
				if (sline.find("atlas")!=string::npos){
					_available_farmtype = arrfarmtypes[1];
				}
			}
		}
	}
	if (_available_farmtype == arrfarmtypes[1]) cout << " yes! " << endl; else cout << " no " << endl;
	cout << " Is gridka farm available? ";
	flush(cout);
	{
		// is it there a gridka login node ?
		// search for "atlas" when executing "bstat" command
		// not tested yet!
		stringstream command;
		command << "qstat | grep atlas > /tmp/_gridka_test;";
		cout << system(command.str().c_str()) << endl;
		ifstream testfile("/tmp/_gridka_test");
		if (testfile){
			while (1) {
				char* line = new char[1024];
				testfile.getline(line, 1024);
				if (!testfile.good()) break;
				string sline(line);
				if (sline.find("atlas")!=string::npos){
					_available_farmtype = arrfarmtypes[2];
				}
			}
		}
	}
	if (_available_farmtype == arrfarmtypes[2]) cout << " yes! " << endl; else cout << " no " << endl;
	//cout << " Is lxbatch farm available? ";
	//flush(cout);

	cout << " test for munich cluster not implemented yet! " << endl;
	/*
	{
		// is it there a gridka login node ?
		// search for "atlas" when executing "bstat" command
		// not tested yet!
		stringstream command;
		command << "bqstat | grep atlas > /tmp/_gridka_test;";
		cout << system(command.str().c_str()) << endl;
		ifstream testfile("/tmp/_gridka_test");
		if (testfile){
			while (1) {
				char* line = new char[1024];
				testfile.getline(line, 1024);
				if (!testfile.good()) break;
				string sline(line);
				if (sline.find("atlas")!=string::npos){
					_available_farmtype = arrfarmtypes[2];
				}
			}
		}
	}*/
	//_available_farmtype = arrfarmtypes[3]; // test purposes
	cout << " farmtype is set to " << _available_farmtype << endl;
}

bool TrpwaJobManager::SendJob(string command, string jobname, int duration){
  if (command == "") return false;
	// make a unique label by time and a counter
	time_t rawtime;
	time ( &rawtime );
	string stime = ctime (&rawtime);
	// replace empty spaces in time string by -
	for (unsigned int i = 0; i < stime.size(); i++){
		if (stime[i] == ' ' || stime[i] == '\n'){
			stime[i]='_';
		}
		//cout << stime[i];
	}
	//cout << endl;
	static int jobcounter(0);
	jobcounter++;
	stringstream batch_script_name;
	batch_script_name << _temp_space << "/" << jobname << "_" << stime << jobcounter << ".sh";
	ofstream batch_script(batch_script_name.str().c_str());
	if (!batch_script) {
		cout << " Error in TrpwaJobManager::SendJob(): could not write script file " <<  batch_script_name.str()<< endl;
		return false;
	}
	if (_available_farmtype == arrfarmtypes[3]){  // case mainz blaster
		batch_script << "#!/bin/bash" << endl;
		batch_script << "#" << endl;
		batch_script << "#PBS -N "<< jobname << endl;
		batch_script << "#PBS -j oe" << endl;
		batch_script << "#PBS -o "<< batch_script_name.str() <<".out" << endl;
		batch_script << "#PBS -V" <<endl;
		batch_script << "#PBS -l nodes=1:x86_64,walltime=12:00:00" << endl;
		batch_script << endl;
		batch_script << "export PATH=$PBS_O_PATH" << endl;
		batch_script << "cd $PBS_O_WORKDIR" << endl;
	}
	if (_available_farmtype == arrfarmtypes[2]){  // case gridka farm
		batch_script << "#!/bin/bash" << endl;
		// execute the login script
		batch_script << "source "<< getenv("HOME") <<"/.bash_profile" << endl;
	}
	batch_script << command << endl;
	batch_script << " echo \"removing " << batch_script_name.str() << "\"" << endl;
	batch_script << "rm " << batch_script_name.str();
	batch_script.close();

	// send the jobs
	if (_available_farmtype == arrfarmtypes[0]){  // case local
		if (system(("source "+batch_script_name.str()).c_str()) == 0) return true; else return false;
	}
	if (_available_farmtype == arrfarmtypes[2]){  // case gridka farm
		if (system(("qsub -q e-long -l max_test_run=1 "+batch_script_name.str()).c_str()) == 0) return true; else return false;
	}
	if (_available_farmtype == arrfarmtypes[3]){  // case mainz blaster
		if (system(("qsub "+batch_script_name.str()).c_str()) == 0) return true; else return false;
	}
	return false;
}
