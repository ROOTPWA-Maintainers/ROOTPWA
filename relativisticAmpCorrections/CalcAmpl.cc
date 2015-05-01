
#include <iostream>
#include <string>
#include <stdio.h>

#include "TFhh.h"
#include "TJSS.h"

#include "primeNumbers.h"
#include "reportingUtils.hpp"

using namespace std;

int main(int narg, char* carg[]) {

	if (narg < 4) {
		cout << endl
		     << "This program requires 3 input strings for the mother and " << endl
		     << "the 2 decay particles, each of the form Jp," << endl
		     << "where J is the spin of the particle and p = +/- its parity" << endl;

		return 0;
	}
	string fileName = "primeNumberCache.root";
	if(narg == 5) {
		fileName = carg[4];
	}
	if(not rpwa::primeNumbers::instance().readCacheFile(fileName)) {
		printErr << "could not read prime number cache file. Aborting..." << endl;
		return 1;
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

	TJSS jss(jmother, pm, jdecay1, p1, jdecay2, p2);
	jss.CalcAmpl();


}

