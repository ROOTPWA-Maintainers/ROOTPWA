#include <iostream>
#include <string>
#include <stdio.h>

#include "libconfig.h++"

#include "libConfigUtils.hpp"

#include "TFhh.h"


bool ReadoutTestKeyfile();


using namespace std;
using namespace libconfig;
using namespace rpwa;


int main(int narg, char* carg[]) {
  
	if (not ReadoutTestKeyfile()) {
		cout << "ReadoutTestKexfile failed" << endl;
	}

	
  if (narg < 4) {
    cout << endl
	 << "This program requires 3 input strings for the mother and "<<endl
	 << "the 2 decay particles, each of the form Jp," << endl
	 << "where J is the spin of the particle and p = +/- its parity" 
	 << endl
	 << endl
	 << "options that may follow the three Jp terms:" << endl
	 << "-H     result output also in header file format" << endl
	 << endl;
    
    return 0;
  }
  
  int opt=0;
  for (int oi=4; oi<narg; oi++) {
    int nchar = sizeof(carg[oi])/sizeof(char);
    if (nchar>1 && carg[oi][1]=='H') {
      cout << "H option length:" << nchar << endl;
      opt=2; 
    }
  }
  
  int  jmother;
  char pmother; int pm;
  int  jdecay1;
  char pdecay1; int p1;
  int  jdecay2;
  char pdecay2; int p2;
  
  sscanf(carg[1], "%d%c", &jmother, &pmother);
  cout << "Mother particle: " << jmother << pmother << endl;
  if (pmother=='+') pm= 1;
  else              pm=-1;
  
  sscanf(carg[2], "%1d%c", &jdecay1, &pdecay1);
  cout << "1. decay particle: " << jdecay1 << pdecay1 << endl;
  if (pdecay1=='+') p1= 1;
  else              p1=-1;
  
  sscanf(carg[3], "%1d%c", &jdecay2, &pdecay2);
  cout << "2. decay particle: " << jdecay2 << pdecay2 << endl;
  if (pdecay2=='+') p2= 1;
  else              p2=-1;
  
  cout << jmother << "," << pm << "," 
       << jdecay1 << "," << p1 << "," 
       << jdecay2 << "," << p2 << endl;
  
  TJSS jss(jmother, pm, jdecay1, p1, jdecay2, p2, opt);
  
}

bool ReadoutTestKeyfile() {
	const string KeyfileName = "../relampl/test.key";
	Config key;
	
	bool debug = false;
	
	if(not parseLibConfigFile(KeyfileName, key, debug)) {
		printWarn << "problems reading keyfile" << endl;
		return false;
	}
	
	const Setting& rootKey = key.getRoot();
	const Setting* waveKey = findLibConfigGroup(rootKey, "wave");
	const Setting* XWaveQn = findLibConfigGroup(*waveKey, "XQuantumNumbers");
	int J_parent;
	int P_parent;
	(*XWaveQn).lookupValue("J",J_parent);
	(*XWaveQn).lookupValue("P",P_parent);
	
	cout << "J: " << J_parent << endl;
	cout << "P: " << P_parent << endl;
	cout << endl;
	
	const Setting& XDecay = key.lookup("wave.XDecay");
	const Setting& isobar = key.lookup("wave.XDecay.isobars.[0]");
	
	const char *iso_name = isobar["name"];

	const char *iso_daughter1 = isobar["fsParticles"][0]["name"];
	const char *iso_daughter2 = isobar["fsParticles"][1]["name"];
	
	int iso_L = isobar["L"];
	int iso_S = isobar["S"];
	
	int XDecay_L = XDecay["L"];
	int XDecay_S = XDecay["S"];
	const char *particle = XDecay["fsParticles"][0]["name"];
	
	cout << "isobar name: " << iso_name << endl;
	cout << "isobar 1. daughter name: " << iso_daughter1 << endl;
	cout << "isobar 2. daughter name: " << iso_daughter2 << endl;
	cout << "isobar L: " << iso_L << endl;
	cout << "isobar S: " << iso_S << endl;
	cout << endl;
	cout << "XDecay L: " << XDecay_L << endl;
	cout << "XDecay S: " << XDecay_S << endl;
	cout << "XDecay 2nd particle: " << particle << endl;
	

	return true;
}
