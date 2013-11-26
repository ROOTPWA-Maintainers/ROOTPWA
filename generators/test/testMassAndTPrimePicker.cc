

#include<cstdlib>
#include<iostream>
#include<unistd.h>

#include <boost/progress.hpp>

#include<TFile.h>
#include<TH1D.h>
#include<TH2D.h>

#include "generatorPickerFunctions.h"
#include "randomNumberGenerator.h"
#include "reportingUtils.hpp"


using namespace libconfig;
using namespace rpwa;
using namespace std;


void printUsage(char* prog, int errCode = 0)
{
	cerr << "usage:" << endl
	     << prog
	     << " -n # [ -s #] -o <file> -p <file> -r <file>" << endl
	     << "    where:" << endl
	     << "        -n #       (max) number of events to generate (default: 100)" << endl
	     << "        -o <file>  ROOT output file" << endl
	     << "        -p         path to particle data table file (default: ./particleDataTable.txt)" << endl
	     << "        -r <file>  reaction config file" << endl
	     << "        -s #   set seed " << endl
	     << endl;
	exit(errCode);
}


int main(int argc, char** argv)
{

	unsigned int nEvents = 100;
	string outputFileName = "";
	string pdgFileName = "./particleDataTable.txt";
	string reactionFile;
	int seed = 123456;
	bool seedSet = false;

	int c;
	while ((c = getopt(argc, argv, "n:o:p:r:s")) != -1) {
		switch (c) {
			case 'n':
				nEvents = atoi(optarg);
				break;
			case 's':
				seed = atoi(optarg);
				seedSet = true;
				break;
			case 'o':
				outputFileName = optarg;
				break;
			case 'p':
				pdgFileName = optarg;
				break;
			case 'r':
				reactionFile = optarg;
				break;

			case 'h':
				printUsage(argv[0]);
				break;
			default:
				printUsage(argv[0]);
				break;
		}
	}

	if(not seedSet) {
		printInfo << "Setting random seed to " << seed << endl;
	}
	randomNumberGenerator::instance()->setSeed(seed);

	TFile* rootfile = TFile::Open(outputFileName.c_str(), "RECREATE");
	TH1D* masses = new TH1D("masses", "masses", 1000, 0, 7);
	TH1D* tPrimeH = new TH1D("tPrime", "tPrime", 1000, 0, 5);
	TH2D* massAgTPrime = new TH2D("5_pi_mass_ag_tp", "5_pi_mass_ag_tp", 1000, 1, 5., 1000, 0, 1.1);
	massAndTPrimePicker* picker = new polynomialMassAndTPrimeSlopePicker();
	Config reactConf2;
	reactConf2.readFile(reactionFile.c_str());
	const Setting& rootS = reactConf2.getRoot();
	picker->init(rootS["t_and_m_dependence"]["settings"]);

	boost::progress_display* progressIndicator = new boost::progress_display(nEvents, cout, "");

	for(unsigned int i = 0; i < nEvents; ++i) {

		double mass, tPrime;
		(*picker)(mass, tPrime);
		masses->Fill(mass);
		tPrimeH->Fill(tPrime);
		massAgTPrime->Fill(mass, tPrime);
		++(*progressIndicator);

	}

	rootfile->Write();
	rootfile->Close();

}
