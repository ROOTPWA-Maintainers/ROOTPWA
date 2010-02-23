#include <complex>
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <cstdlib>
#include <unistd.h>


using namespace std;


void
printUsage(char* prog)
{
  cerr << "usage:" << endl
       << "  " << prog << " [-h] [-m maxEvents] [-B] < inputfile" << endl
       << "-B\tconvert ascii to binary mode" << endl;

}


int
main(int argc, char** argv)
{
  int maxEvents  = 0;
  int nRead      = 0;
  int binaryMode = 0;
	
  extern char* optarg;
  int c;
  while ((c = getopt(argc, argv, "m:h:B")) != -1)
    switch (c) {
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

  complex<double>    a;
  const unsigned int size = sizeof(a);
  if (binaryMode) {
    while ((cin >> a) && (maxEvents ? (nRead < maxEvents) : 1)) {
      nRead++;
      cout.write((char*) &a, size);
    }
  } else {
    while ((cin.read((char*) &a, size)) && (maxEvents ? (nRead < maxEvents) : 1)) {
      nRead++;
      const unsigned int nmbDigits = numeric_limits<double>::digits10 + 1;
      ostringstream s;
      s.precision(nmbDigits);
      s.setf(ios_base::scientific, ios_base::floatfield);
      s << a;
      cout << s.str() << endl;
    }
  }
  
  return 0;
}
