
#include <iostream>

#include <boost/progress.hpp>

#include <TFile.h>
#include <TTree.h>

#include <reportingUtils.hpp>

#include "primeNumbers.h"


using namespace std;


void usage(const string& progName,
           const int     errCode = 0)
{
	cerr << "generate a prime number cache file" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " -o outputFileName [ -l cacheSize ]"<< endl
	     << "    where:" << endl
	     << "        -o file  the name of the output file to generate (default: primeNumberCache.root)" << endl
	     << "        -l the number of prime numbers to generate (default: 100000)" << endl
	     << "        -h print help " << endl
	     << endl;
	exit(errCode);
}


int main(int argc, char** argv)
{
	std::string fileName  = "primeNumberCache.root";
	size_t      cacheSize = 100000;
	extern char* optarg;
	int c;
	while ((c = getopt(argc, argv, "o:l:h")) != -1)
	{
		switch (c) {
		case 'o':
			fileName = optarg;
			break;
		case 'l':
			cacheSize = atoi(optarg);
			break;
		case 'h':
			usage(argv[0]);
			break;
		}
	}

	TFile* cacheFile = TFile::Open(fileName.c_str(), "NEW");
	if(not cacheFile) {
		printErr << "could not create file '" << fileName << "'. Aborting..." << endl;
		return 1;
	}
	TTree* tree = new TTree(rpwa::primeNumbers::TREE_NAME.c_str(), rpwa::primeNumbers::TREE_NAME.c_str());
	rpwa::primeNumbers::entryType entry = 0;
	tree->Branch(rpwa::primeNumbers::BRANCH_NAME.c_str(), &entry, "entry/l");
	entry = 2;
	tree->Fill();
	entry = 3;
	size_t nmbFoundPrimes = 1;
	boost::progress_display progressIndicator(cacheSize, cout, "");
	while(nmbFoundPrimes < cacheSize) {
		if(rpwa::primeNumbers::isPrimeNumber(entry)) {
			tree->Fill();
			++nmbFoundPrimes;
			++progressIndicator;
		}
		entry += 2;
	}
	cout << endl;
	tree->Write();
	cacheFile->Close();
	return 0;
}
