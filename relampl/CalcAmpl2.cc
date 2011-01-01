#include <iostream>
#include <string>
#include <stdio.h>

#include "libconfig.h++"

#include "libConfigUtils.hpp"

#include "TFhh.h"

using namespace std;
using namespace libconfig;
using namespace rpwa;

// declaring functions:
bool ReadoutTestKeyfile();
bool buildDecayTopologyRelAmpl(const Setting* parent, int p_J, int p_P);
bool RelAmplDummy(int p_J, int p_P, int d1_J, int d1_P, int d2_J, int d2_P, int L, int S);
int lookupJDummy(const char* name);
int lookupPDummy(const char* name);


// main function:
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
	int iso_L;
	int iso_S;
	int XDecay_L;
	int XDecay_S;
	const char* iso_name;
	const char* iso_daughter1;
	const char* iso_daughter2;
	const char* particle;

	(*XWaveQn).lookupValue("J",J_parent);
	(*XWaveQn).lookupValue("P",P_parent);

	const Setting* XDecay = findLibConfigGroup(*waveKey,"XDecay");
	const Setting* isobar = findLibConfigList(*XDecay,"isobars");
	const Setting* isobar0 = findLibConfigGroup(*isobar,0);

	iso_name = (*isobar0)["name"];
	iso_daughter1 = (*isobar0)["fsParticles"][0]["name"];
	iso_daughter2 = (*isobar0)["fsParticles"][1]["name"];
	
	iso_L = (*isobar0)["L"];
	iso_S = (*isobar0)["S"];
	
	XDecay_L = (*XDecay)["L"];
	XDecay_S = (*XDecay)["S"];
	particle = (*XDecay)["fsParticles"][0]["name"];
	
//	cout << "J: " << J_parent << endl;
//	cout << "P: " << P_parent << endl;
//	cout << endl;	
//	cout << "isobar name: " << iso_name << endl;
//	cout << "isobar 1. daughter name: " << iso_daughter1 << endl;
//	cout << "isobar 2. daughter name: " << iso_daughter2 << endl;
//	cout << "isobar L: " << iso_L << endl;
//	cout << "isobar S: " << iso_S << endl;
//	cout << endl;
//	cout << "XDecay L: " << XDecay_L << endl;
//	cout << "XDecay S: " << XDecay_S << endl;
//	cout << "XDecay 2nd particle: " << particle << endl;
	
	if (not buildDecayTopologyRelAmpl(XDecay,J_parent,P_parent)) {
		return false;
	}
	
	return true;
}

bool buildDecayTopologyRelAmpl(const Setting* parent, int p_J = 0, int p_P = 0) {
	
	int d1_J = 0, d1_P = 0, d2_J = 0, d2_P = 0, L = 0, S = 0;
//	const char* name;
//	const char* fsParticle1;
//	const char* fsParticle2;
//	const char* fsParticle3;
	const Setting* isobar = findLibConfigList(*parent,"isobars",false);
	
	if (not not isobar) {
		const Setting* isobar0 = findLibConfigGroup(*isobar,0);
		if (not not isobar0) {
			if (not buildDecayTopologyRelAmpl(isobar0)) return false;
		}
		
		const Setting* isobar1 = findLibConfigGroup(*isobar,1);
		if (not not isobar1){
			if (not buildDecayTopologyRelAmpl(isobar1)) return false;
		}
	}
	
	try { 
		const char* name = (*parent)["name"];
		cout << "name: " << name << endl;
		p_J = lookupJDummy(name);
		p_P = lookupPDummy(name);
	} catch (const SettingNotFoundException& NotFound) {
		
	}
	
	L = (*parent)["L"];
	S = (*parent)["S"];
	const Setting* fsParticles = findLibConfigList(*parent,"fsParticles",false);
	if (not not fsParticles) {
		const Setting* fsParticles0 = findLibConfigGroup(*fsParticles,0);
		if (not not fsParticles0) {
			const char* fsParticle1 = (*fsParticles0)["name"];
			cout << "fsParticle1: " << fsParticle1 << endl;
			d1_J = lookupJDummy(fsParticle1);
			d1_P = lookupPDummy(fsParticle1);
		}
		
		const Setting* fsParticles1 = findLibConfigGroup(*fsParticles,1);
		if (not not fsParticles1) {
			const char* fsParticle2 = (*fsParticles1)["name"];
			cout << "fsParticle2: " << fsParticle2 << endl;
			d2_J = lookupJDummy(fsParticle2);
			d2_P = lookupPDummy(fsParticle2);
		} else {
			const char* Particle2 = (*parent)["isobars"][0]["name"];
			d2_J = lookupJDummy(Particle2);
			d2_P = lookupPDummy(Particle2);
		}
		
//		const Setting* fsParticles2 = findLibConfigGroup(*fsParticles,2);
//		if (not not fsParticles2) {
//			const char* fsParticle3 = (*fsParticles2)["name"];
//			cout << "fsParticle3: " << fsParticle3 << endl;
//		}
	}
	
//	cout << "name: " << name << endl;
	cout << "L: " << L << endl;
	cout << "S: " << S << endl;
	
	if (not RelAmplDummy(p_J,p_P,d1_J,d1_P,d2_J,d2_P,L,S)) return false;
	return true;
}

bool RelAmplDummy (int p_J, int p_P, int d1_J, int d1_P, int d2_J, int d2_P, int L, int S) {
	// This dummy doesn't calculate the amplitudes, it just returns constants
	// This function calcuates the relativistic amplitudes and writes them into the key file
//	cout << "RelAmplDummy running" << endl;
	return true;
}

int lookupJDummy(const char* name) {
//	cout << "lookupJDummy running" << endl;
	return 0;
}

int lookupPDummy(const char* name) {
//	cout << "lookupPDummy running" << endl;
	return 0;
}
