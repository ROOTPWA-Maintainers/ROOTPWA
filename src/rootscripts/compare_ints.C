// jasinski@kph.uni-mainz.de, promme@web.de
// created: 03.11.10
// comparing to integral files with each other
// needed since the matrix layout is following the input of amplitudes to "int"
// call with root for example like:
// root -l compare_ints.C+\(\"1400.1420/ACCAMPS/normwo0.int\"\,\"../SET_20MeV/1400.1420/ACCAMPS/norm.int\"\)
// with all slashes!
// copy rights: do what ever you want with it as long as you don't blame me!
#include <complex>
#include <iostream>
#include <string>
#include <fstream>
#include <map>

using namespace std;


void compare_ints(string file1, string file2){
	ifstream filestream[2];
	filestream[0].open(file1.c_str());
	filestream[1].open(file2.c_str());
	if (!filestream[0] || !filestream[1]){
		cout << " Error, could not load files " << endl;
		return;	
	}

	int nentries[2];
	int nwaves[2];
	int dimensionx[2];
	int dimensiony[2];

	// header must be the same
	for (int ifile = 0; ifile < 2; ifile++){
		filestream[ifile] >> nwaves[ifile];
		filestream[ifile] >> nentries[ifile];
		filestream[ifile] >> dimensionx[ifile] >> dimensiony[ifile];
	}

	if (!filestream[0].good() || !filestream[1].good()){
		cout << " Error, could not read header of files " << endl;
		return;	
	}
	
	if (nwaves[0] != nwaves[1] || 
		nentries[0] != nentries[1] || 
		dimensionx[0] != dimensionx[1] || 
		dimensiony[0] != dimensiony[1]) 
	{
		cout << " Error, the header differ! " << endl;
		return;
	}
	if (dimensionx[0] > 100 || dimensiony[0] > 100){
		cout << " cannot process more than 100 x 100 entries, sorry " << endl;
		return;	
	}
	complex<double> entries[2][100][100];

	cout << " reading matrix of size " << dimensionx[0] << " X " << dimensiony[0] << endl;

	// now read all entries into memory
	for (int ix = 0; ix < dimensionx[0] ; ix++){
		for (int iy = 0; iy < dimensiony[0]; iy++){
			for (int ifile = 0; ifile < 2; ifile++){
				filestream[ifile] >> entries[ifile][ix][iy];
				if (!filestream[ifile].good()){
					cout << " Error reading matrix, aborting " << endl;
					return;
				}
			}	
		}	
	}
	
	// the number of waves
	for (int ifile = 0; ifile < 2; ifile++){
		filestream[ifile] >> nwaves[ifile];
		if (!filestream[ifile].good()){
			cout << " Error reading number of waves, aborting " << endl;
			return;
		}
	}
	
	if (nwaves[0] != nwaves[1]){
		cout << " Error: number of waves differs, aborting " << endl;
		return;
	}
	cout << " reading " << nwaves[0] << " waves" << endl;

	map<string,int> waves[2];
	// finaly read the corresponding waves
	for (int ifile = 0; ifile < 2; ifile++){
		for (int iwave = 0; iwave < nwaves[ifile]; iwave++){
			string key;
			int position;
			filestream[ifile] >> key >> position;
			if (!filestream[ifile].good()){
				cout << " Error reading waves, aborting " << endl;
				return;
			}
			// get rid of an eventual path
			int slashpos = key.rfind('/');
			if (slashpos != (int) string::npos){
				key.erase(0, slashpos+1);
			}
			waves[ifile][key] = position;
		}
	}

	// check mapping by size
	if (waves[0].size() != waves[1].size()){
		cout << " Error: waves differ " << endl;
		return;	
	}
	if (waves[0].size() != (unsigned) nwaves[0]){
		cout << " Error: number of different waves is not correct, aborting " << endl;
		return;	
	}

	int ncomperisons(0);
	// compare all entries according to the key
	for ( map<string,int>::iterator itx = waves[0].begin(); itx != waves[0].end(); itx++){
		for ( map<string,int>::iterator ity = waves[0].begin(); ity != waves[0].end(); ity++){
			if (entries[0][itx->second][ity->second] != entries[1][waves[1][itx->first]][waves[1][ity->first]]){
				cout << " Error: wave " << itx->first << " differs " << endl;
				return;
			}
			ncomperisons++;
		}	
	}
	// just a check since everything until now went too smooth!
	if (ncomperisons != nwaves[0]*nwaves[0]){
		cout << " Unexpected error: number of comparisons is too low! " << endl;
	}

	filestream[0].close();
	filestream[1].close();

	cout << " files seem to be equivalent " << endl;
}
