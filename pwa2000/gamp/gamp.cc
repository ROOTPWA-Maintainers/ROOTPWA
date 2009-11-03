#line 13 "../gamp.nw"
#include <cstdlib>
#include <unistd.h>
#include <pputil.h>
#include <keyfile.h>

extern int keydebug;
extern particleDataTable PDGtable;
void printUsage(char* pname);

int main(int argc, char** argv) {
	char *pdgFile = (char *) NULL;
	char* pname = argv[0];
	extern char *optarg;
	extern int optind;
	int c;
	keyfile keyf;
	event e;
	int io_ver = 1;
	int printPDG = 0;
	int dumpPDG = 0;

	
#line 74 "../gamp.nw"
	while ( (c = getopt(argc,argv, "htP:vi:D")) != -1 )
		switch(c) {
		case 'P':
			pdgFile = optarg;
			break;
		case 'v':
			printPDG = 1;
			break;
		case 'D':
			dumpPDG = 1;
			break;
		case 'i':
			io_ver = atoi(optarg);
			break;
		case 't':
			keydebug = 1;
			break;
		case 'h':
		case '?':
			printUsage(pname);
		}

#line 35 "../gamp.nw"
        if(pdgFile)
		PDGtable.initialize(pdgFile);
        else
	        PDGtable.initialize();

	if (dumpPDG) {
		PDGtable.dump();
		cout << flush;
	}
	if (printPDG) {
		PDGtable.print();
		cout << flush;
	}	
		
 
	if (!argv[optind]) {
		cerr << pname << ": ERROR: no keyfile specified" << endl;
		printUsage(pname);
	}
	keyf.open(argv[optind]);	

	
#line 102 "../gamp.nw"
	e.setIOVersion(io_ver);
	while(!(cin>>e).eof()) {
		keyf.run(e);
		keyf.rewind();
	}



#line 58 "../gamp.nw"
	return 0;

}

#line 110 "../gamp.nw"
void printUsage(char* pname) {
	cerr << "usage: " << pname << " [-t] keyfile [-P pdgFile] [-v] [-D] < datafile" << endl;
	cerr << "\t-t:\tparser trace" << endl;
	cerr << "\t-P pdgFile:\tread PDG table from pdgFile" << endl;
	cerr << "\t-v:\tprint pdgTable" << endl;
	cerr << "\t-i io_ver:\tset event io version to io_ver (1 or 2)" << endl;
	cerr << "\t-D:\tdump pdgTable" << endl;
	cerr << "\tkeyfile: decay amplitude specification" << endl;

	exit(1);
}

