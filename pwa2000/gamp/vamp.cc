#line 20 "../vamp.nw"
#include <complex>
#include <iostream>
#include <cstdlib>
#include <unistd.h>

using namespace std;

#line 81 "../vamp.nw"
void printUsage(char* prog) {
	cerr << "usage:" << endl;
	cerr << "  " << prog << " [-h] [-m maxEvents] [-B] < inputfile" << endl;
	cerr << "-B\tconvert ascii to binary mode" << endl;

}

#line 29 "../vamp.nw"
int main(int argc, char** argv) {

	int maxEvents = 0;
	int nRead = 0;
	int binaryMode = 0;

	
#line 63 "../vamp.nw"
    extern char* optarg;
    int c;
    while ( (c = getopt(argc,argv, "m:h:B")) != -1 )
        switch(c) {
            case 'm':
                maxEvents = atoi(optarg);
                break;
		case 'B':
			binaryMode = 1;
			break;
            case 'h':
            case '?':
                printUsage(argv[0]);
				exit(0);
	       
        }


#line 37 "../vamp.nw"
	complex<double> a;

	if (binaryMode) {
	while ((cin >> a) && (maxEvents?nRead<maxEvents:1) ) {
		nRead++;
		cout.write((char*) &a,sizeof(a));
	}
	}
	else {

	while( (cin.read((char*) &a,sizeof(a))) && (maxEvents?nRead<maxEvents:1) ) {
		nRead++;
		cout << a << endl;
	}
}

	return 0;

}

