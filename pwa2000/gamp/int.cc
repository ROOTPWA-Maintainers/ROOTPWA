#include <iostream>
#include <string>
#include <list>
#include <cstdlib>
#include <unistd.h>

#include "integral.h"


using namespace std;


void
printUsage(const char* prog)
{
  cerr << "usage:"  << endl
       << "  creating     : " << prog << " [-d] [-m max] [-r n] [-w weightfile] files" << endl
       << "  adding       : " << prog << " [-d] [-m max] [-r n] -a intfile"            << endl
       << "  display evts : " << prog << " [-q] -i intfile"                            << endl
       << "  renormalizing: " << prog << " [-r n] -i intfile"                          << endl
       << "  help         : " << prog << " [-h]"                                       << endl
       << "where:" << endl
       << "    max       : maximum number of events"                 << endl
       << "    n         : number of events to renormalize to"       << endl
       << "    weightfile: file containing weights for de-weighting" << endl
       << "                (values will be divided by weights!)"     << endl
       << "    intfile   : integral file to read"                    << endl
       << "                (for adding to or renormalizing)"         << endl
       << endl;
  exit(0);
}


int main(int argc, char** argv)
{

	//int    debug          = 0;
  int    display_events = 0;
  int    maxevents      = 0;
  int    renorm         = 0;
  string oldint;
  int    add            = 1;
  string weightfilename;
  extern int   optind;
  extern char* optarg;
  int c;
  if (argc == 1)
    printUsage(argv[0]);
  while ((c = getopt(argc, argv, "dm:r:a:i:h:w:q")) != -1)
    switch(c) {
    case 'd':
	    //debug = 1;    
      break;
    case 'q':
      display_events = 1;
      break;
    case 'm':
      maxevents = atoi(optarg);    
      break;
    case 'r':
      renorm = atoi(optarg);    
      break;
    case 'a':
      oldint = optarg;    
      break;
    case 'i':
      add = 0;
      oldint = optarg;    
      break;
    case 'w':
      weightfilename = optarg;    
      break;
    case 'h':
    case '?':
      printUsage(argv[0]);
    }

  integral ni;
  if (oldint.size() != 0) {
    ifstream oldfile(oldint.c_str());
    ni.scan(oldfile);
  } else
    ni.files(argv + optind);

  if (weightfilename.size() != 0)
    ni.weightfile(weightfilename);

  if (display_events)
    ni.print_events();
  else {
    if (add) {
      ni.max(maxevents);
      try {
        ni.integrate();
      } catch (const char* m) {
        cerr << m << endl;
        return 0;
      }
    }
    if (renorm)
      ni.renormalize(renorm);
    ni.print();
  }
  return 0;
}
