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
	cout << " farmtype is set to " << _available_farmtype << endl;
}
