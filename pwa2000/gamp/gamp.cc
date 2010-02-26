#include <cstdlib>
#include <unistd.h>

#include "event.h"
#include "keyfile.h"


using namespace std;


extern int               keydebug;
extern particleDataTable PDGtable;


void printUsage(char* pname);


int main(int argc, char** argv)
{
  char* pdgFile = (char*) NULL;
  char* pname   = argv[0];
  int   io_ver   = 1;
  int   printPDG = 0;
  int   dumpPDG  = 0;
	
  int          c;
  extern char* optarg;
  extern int   optind;
  while ((c = getopt(argc, argv, "htP:vi:D")) != -1)
    switch (c) {
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

  if (pdgFile)
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
  keyfile keyf;
  keyf.open(argv[optind]);	

  event e;
  e.setIOVersion(io_ver);
  while (!(cin >> e).eof()) {
    keyf.run(e);
    keyf.rewind();
  }

  return 0;
}


void printUsage(char* pname)
{
  cerr << "usage: " << pname << " [-t] keyfile [-P pdgFile] [-v] [-D] < datafile" << endl
       << "\t-t:\tparser trace" << endl
       << "\t-P pdgFile:\tread PDG table from pdgFile" << endl
       << "\t-v:\tprint pdgTable" << endl
       << "\t-i io_ver:\tset event io version to io_ver (1 or 2)" << endl
       << "\t-D:\tdump pdgTable" << endl
       << "\tkeyfile: decay amplitude specification" << endl;
  exit(1);
}
