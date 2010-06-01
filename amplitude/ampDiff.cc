///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      compares two amplitude files and calculates differences
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <iostream>
#include <fstream>
#include <unistd.h>
#include <vector>
#include <complex>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include "svnVersion.h"
#include "utilities.h"


using namespace std;


void
usage(const string& progName,
      const int     errCode = 0)
{
  cerr << "usage:" << endl
       << progName
       << " [-o out file -v -h] .amp file 1 .amp file 2" << endl
       << "    where:" << endl
       << "        -o file    path to output ROOT file (default: none)" << endl
       << "        -v         verbose; print debug output (default: false)" << endl
       << "        -h         print help" << endl
       << endl;
  exit(errCode);
}


int
main(int    argc,
     char** argv)
{
  printCompilerInfo();
  printSvnVersion();

  // parse command line options
  const string progName        = argv[0];
  string       outFileName     = "";
  bool         debug           = false;
  string       ampFileNames[2] = {"", ""};
  extern char* optarg;
  extern int   optind;
  int          c;
  while ((c = getopt(argc, argv, "o:vh")) != -1)
    switch (c) {
    case 'o':
      outFileName = optarg;
      break;
    case 'v':
      debug = true;
      break;
    case 'h':
    default:
      usage(progName);
    }

  // get amplitude file names
  if (argc - optind != 2) {
    printErr << "you have to specify two amplitude files. aborting." << endl;;
    usage(progName, 1);
  }
  ampFileNames[0] = argv[optind++];
  ampFileNames[1] = argv[optind++];

  // open amplitude files
  printInfo << "comparing amplitude files '" << ampFileNames[0] << "' and "
	    << "'" << ampFileNames[1] << endl;
  ifstream ampFiles[2];
  for (unsigned int i = 0; i < 2; ++i) {
    ampFiles[i].open(ampFileNames[i].c_str());
    if (!ampFiles[i]) {
      printErr << "cannot open amplitude file '" << ampFileNames[i] << "'. aborting." << endl;
      exit(1);
    }
  }

  // read amplitudes into memory
  vector<complex<double> > amps[2];
  for (unsigned int i = 0; i < 2; ++i) {
    complex<double> amp;
    while (ampFiles[i].read((char*) &amp, sizeof(complex<double>)))
      amps[i].push_back(amp);
  }  
  if (amps[0].size() != amps[1].size())
    printWarn << "the two amplitude files have different number of amplitudes "
	      << "(" << amps[0].size() << " vs. " << amps[1].size() << ")." << endl;
  const unsigned int nmbAmps = min(amps[0].size(), amps[1].size());

  // open output file and book histograms
  TFile* outFile = 0;
  if (outFileName != "") {
    printInfo << "writing difference plots to '" << outFileName << "'" << endl;
    outFile = TFile::Open(outFileName.c_str(), "RECREATE");
  }
  TH1D* hAmps1Re   = new TH1D("hAmps1Re",   "hAmps1Re;Event Number;#Rgothic[Amp 1]", nmbAmps, -0.5, nmbAmps - 0.5);
  TH1D* hAmps1Im   = new TH1D("hAmps1Im",   "hAmps1Im;Event Number;#Jgothic[Amp 1]", nmbAmps, -0.5, nmbAmps - 0.5);
  TH1D* hAmps2Re   = new TH1D("hAmps2Re",   "hAmps2Re;Event Number;#Rgothic[Amp 2]", nmbAmps, -0.5, nmbAmps - 0.5);
  TH1D* hAmps2Im   = new TH1D("hAmps2Im",   "hAmps2Im;Event Number;#Jgothic[Amp 2]", nmbAmps, -0.5, nmbAmps - 0.5);
  TH1D* hAbsDiffRe = new TH1D("hAbsDiffRe", "hAbsDiffRe;#Rgothic[Amp 1 - Amp 2];Count", 100000, -3e-7, 3e-7);
  TH1D* hAbsDiffIm = new TH1D("hAbsDiffIm", "hAbsDiffIm;#Jgothic[Amp 1 - Amp 2];Count", 100000, -3e-7, 3e-7);
  TH1D* hRelDiffRe = new TH1D("hRelDiffRe", "hRelDiffRe;(#Rgothic[Ampl 1] - #Rgothic[Ampl 2]) / #Rgothic[Ampl 1];Count", 100000, -3e-7, 3e-7);
  TH1D* hRelDiffIm = new TH1D("hRelDiffIm", "hRelDiffIm;(#Jgothic[Ampl 1] - #Jgothic[Ampl 2]) / #Jgothic[Ampl 1];Count", 100000, -3e-7, 3e-7);
  TH2D* hCorrRe    = new TH2D("hCorrRe",    "hCorrRe;#Rgothic[Amp 1];#Rgothic[Amp 2]", 1000, -2, 2, 1000, -2, 2);
  TH2D* hCorrIm    = new TH2D("hCorrIm",    "hCorrIm;#Jgothic[Amp 1];#Jgothic[Amp 2]", 1000, -2, 2, 1000, -2, 2);

  // compare amplitudes
  double maxAbsDiff = 0;
  double maxRelDiff = 0;
  for (unsigned int i = 0; i < nmbAmps; ++i) {
    const complex<double> absDiff = amps[0][i] - amps[1][i];
    const complex<double> relDiff = complex<double>(absDiff.real() / amps[0][i].real(),
						    absDiff.imag() / amps[0][i].imag());
    // fill histograms
    if (outFile) {
      hAmps1Re->SetBinContent(i + 1, amps[0][i].real());
      hAmps1Im->SetBinContent(i + 1, amps[0][i].imag());
      hAmps2Re->SetBinContent(i + 1, amps[1][i].real());
      hAmps2Im->SetBinContent(i + 1, amps[1][i].imag());
      hAbsDiffRe->Fill(absDiff.real());
      hAbsDiffIm->Fill(absDiff.imag());
      hRelDiffRe->Fill(relDiff.real());
      hRelDiffIm->Fill(relDiff.imag());
      hCorrRe->Fill(amps[0][i].real(), amps[1][i].real());
      hCorrIm->Fill(amps[0][i].imag(), amps[1][i].imag());
    }
    // print amplitudes
    if (debug) {
      const unsigned int nmbDigits = numeric_limits<double>::digits10 + 1;
      ostringstream s;
      s.precision(nmbDigits);
      s.setf(ios_base::scientific, ios_base::floatfield);
      s << "    " << setw(49) << amps[0][i] << " - " << setw(49) << amps[1][i]
	<< " = " << setw(49) << absDiff	<< ", relative "
	<< "(" << setw(23) << relDiff.real() << ", " << setw(23) << relDiff.imag() << " )";
      cout << s.str() << endl;
    }
    // compute maximum deviations
    maxAbsDiff = max(fabs(absDiff.real()), fabs(absDiff.imag()));
    maxRelDiff = max(fabs(relDiff.real()), fabs(relDiff.imag()));
  }

  printInfo << "maximum observed deviation absolute = " << maxAbsDiff << ", "
	    << "relative = " << maxRelDiff << endl;
  if (outFile) {
    outFile->Write();
    outFile->Close();
    delete outFile;
    printInfo << "wrote difference plots to '" << outFileName << "'" << endl;
  }
  return 0;
}
