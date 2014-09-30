#include <iostream>
#include <string>
#include <stdio.h>

#include "libconfig.h++"

#include "libConfigUtils.hpp"
#include "particleDataTable.h"
#include "particle.h"

#include "TFhh.h"

using namespace std;
using namespace libconfig;
using namespace rpwa;

// declaring functions:
int JansRelAmpl(int narg, char* carg[]);
bool ReadoutTestKeyfile(const char* filename);
bool buildDecayTopologyRelAmpl(const Setting* parent, int p_J, int p_P);
bool RelAmplDummy(int p_J, int p_P, int d1_J, int d1_P, int d2_J, int d2_P, int L, int S, const Setting* parent);
int lookupJ(const char* name);
int lookupP(const char* name);

// main function:
int main(){
	if (not ReadoutTestKeyfile("../relativisticAmpCorrections/test.key")) {
		cout << "ReadoutTestKeyfile failed" << endl;
	}

	return 0;
}

int JansRelAmpl(int narg, char* carg[]) {

	if (not ReadoutTestKeyfile("../relativisticAmpCorrections/test.key")) {
		cout << "ReadoutTestKeyfile failed" << endl;
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

	return 0;
}

bool ReadoutTestKeyfile(const char* filename) {
	const string KeyfileName = filename;
	Config key;

	bool debug = false;

	if(not parseLibConfigFile(KeyfileName, key, debug)) {
		printWarn << "problems reading keyfile" << endl;
		return false;
	}

	const Setting& rootKey      = key.getRoot();
	const Setting* decayVertKey = findLibConfigGroup(rootKey,       "decayVertex");
	const Setting* XWaveQn      = findLibConfigGroup(*decayVertKey, "XQuantumNumbers");
	const Setting* XDecay       = findLibConfigGroup(*decayVertKey, "XDecay");

	int J_parent = (*XWaveQn)["J"];
	int P_parent = (*XWaveQn)["P"];

	particleDataTable& pdt = particleDataTable::instance();
	pdt.readFile("../particleData/particleDataTable.txt");


	if (not buildDecayTopologyRelAmpl(XDecay,J_parent/2,P_parent)) {
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
		const Setting* isobar0 = findLibConfigGroup(*isobar,0,false);
		if (not not isobar0) {
			if (not buildDecayTopologyRelAmpl(isobar0)) return false;
		}

		const Setting* isobar1 = findLibConfigGroup(*isobar,1,false);
		if (not not isobar1){
			if (not buildDecayTopologyRelAmpl(isobar1)) return false;
		}
	}

	try {
		const char* name = (*parent)["name"];
		//		cout << "name: " << name << endl;
		p_J = lookupJ(name);
		p_P = lookupP(name);
	} catch (SettingNotFoundException& NotFound) {

	}

	L = (*parent)["L"];
	S = (*parent)["S"];
	const Setting* fsParticles = findLibConfigList(*parent,"fsParticles",false);
	if (not not fsParticles) {
		const Setting* fsParticles0 = findLibConfigGroup(*fsParticles,0,false);
		if (not not fsParticles0) {
			const char* fsParticle1 = (*fsParticles0)["name"];
			//			cout << "fsParticle1: " << fsParticle1 << endl;
			d1_J = lookupJ(fsParticle1);
			d1_P = lookupP(fsParticle1);
		}

		const Setting* fsParticles1 = findLibConfigGroup(*fsParticles,1,false);
		if (not not fsParticles1) {
			const char* fsParticle2 = (*fsParticles1)["name"];
			//			cout << "fsParticle2: " << fsParticle2 << endl;
			d2_J = lookupJ(fsParticle2);
			d2_P = lookupP(fsParticle2);
		} else {
			const char* Particle2 = (*parent)["isobars"][0]["name"];
			d2_J = lookupJ(Particle2);
			d2_P = lookupP(Particle2);
		}

		//		const Setting* fsParticles2 = findLibConfigGroup(*fsParticles,2);
		//		if (not not fsParticles2) {
		//			const char* fsParticle3 = (*fsParticles2)["name"];
		//			cout << "fsParticle3: " << fsParticle3 << endl;
		//		}
	}

	//	cout << "name: " << name << endl;
	//	cout << "L: " << L << endl;
	//	cout << "S: " << S << endl;

	if (not RelAmplDummy(p_J,p_P,d1_J,d1_P,d2_J,d2_P,L,S,parent)) return false;
	return true;
}

bool RelAmplDummy (int p_J, int p_P, int d1_J, int d1_P, int d2_J, int d2_P, int L, int S, const Setting* parent) {
	// This dummy doesn't calculate the amplitudes, it just returns constants
	// This function calcuates the relativistic amplitudes and writes them into the key file
	//	cout << "RelAmplDummy running" << endl;
	cout << endl;
	cout << "Parent J: "<< p_J << endl;
	cout << "Parent P: "<< p_P << endl;
	cout << "Daughter 1 J: "<< d1_J << endl;
	cout << "Daughter 1 P: "<< d1_P << endl;
	cout << "Daughter 2 J: "<< d2_J << endl;
	cout << "Daughter 2 P: "<< d2_P << endl;
	cout << "L: " << L << endl;
	cout << "S: " << S << endl;

	//	int opt = 0;
	//	TJSS jss(p_J, p_P, d1_J, d1_P, d2_J, d2_P, opt);

	const char* filename = "../relativisticAmpCorrections/test.key";
	Config key;
	parseLibConfigFile(filename, key, false);
	const string path = (*parent).getPath();
	//	cout << "Path: " << path << endl;

	Setting& parent0 = key.lookup(path);

	if (parent0.exists("Relampl")) {
		parent0.remove("Relampl");
	}

	parent0.add("Relampl", Setting::TypeString) = path;


	key.writeFile(filename);

	return true;
}

int lookupJ(const char* name) {
	//	cout << "lookupJDummy running" << endl;
	//	particleDataTable& pdt = particleDataTable::instance();
	//	pdt.readFile("../particleData/particleDataTable.txt");
	particleProperties partProp;
	if (not partProp.fillFromDataTable(name)) {
		printErr << "canot find particle '" << name << "' in particle data table. aborting." << endl;
		throw;
	}
	return (partProp.J()/2);
}

int lookupP(const char* name) {
	//	cout << "lookupPDummy running" << endl;
	//	particleDataTable& pdt = particleDataTable::instance();
	//	pdt.readFile("../particleData/particleDataTable.txt");
	particleProperties partProp;
	if (not partProp.fillFromDataTable(name)) {
		printErr << "canot find particle '" << name << "' in particle data table. aborting." << endl;
		throw;
	}
	return partProp.P();
}
