#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>

#include "TTree.h"
#include "TMatrixT.h"
#include "TComplex.h"

#include "../TFitBin.h"
#include "../utilities.h"


using namespace std;


void
testNewTFitBinInterface(TTree* tree)
{
  const bool verbose = true;
  //const bool verbose = false;
  TFitBin*   massBin = new TFitBin();
  tree->SetBranchAddress("fitbin", &massBin);
  tree->GetEntry(30);
  const unsigned int n = massBin->nmbWaves();

  if (0) {
    massBin->printWaveNames();
    cout << endl;
    massBin->printProdAmpNames();
    cout << endl;
    massBin->printProdAmps();
    cout << endl;
  }
  
  if (0) {
    complex<double> maxDelta = 0;
    for (unsigned int i = 0; i < n; ++i)
      for (unsigned int j = 0; j < n; ++j) {
	const TComplex        temp   = massBin->spinDens(i, j);
	const complex<double> oldVal = complex<double>(temp.Re(), temp.Im());
	const complex<double> newVal = massBin->spinDensityMatrixElem(i, j);
	const complex<double> delta  = oldVal - newVal;
	maxDelta.real() = (fabs(maxDelta.real()) < fabs(delta.real())) ? delta.real() : maxDelta.real();
	maxDelta.imag() = (fabs(maxDelta.imag()) < fabs(delta.imag())) ? delta.imag() : maxDelta.imag();
	if (verbose)
	  cout << "spinDensityMatrixElem(" << massBin->waveName(i) << ", "  << massBin->waveName(j) << "): "
	       << setprecision(12) << newVal << " vs. " << oldVal << ", delta = " << delta << endl;
      }
    cout << "spinDensityMatrixElem() max. deviation = " << maxDelta << endl << endl;
  }

  if (1) {
    double maxDelta = 0;
    for (unsigned int i = 0; i < n; ++i) {
      const double oldVal = massBin->intens(i);
      const double newVal = massBin->intensity(i);
      const double delta  = oldVal - newVal;
      maxDelta = (fabs(maxDelta) < fabs(delta)) ? delta : maxDelta;
      if (verbose)
	cout << "intensity(" << massBin->waveName(i) << "): "
	     << setprecision(12) << newVal << " vs. " << oldVal << ", delta = " << oldVal - newVal << endl;
    }
    cout << "intensity() max. deviation = " << maxDelta << endl << endl;
  }

  if (1) {
    double maxDelta = 0;
    for (unsigned int i = 0; i < n; ++i) {
      const double oldVal = massBin->err(i);
      const double newVal = massBin->intensityErr(i);
      const double delta  = oldVal - newVal;
      maxDelta = (fabs(maxDelta) < fabs(delta)) ? delta : maxDelta;
      if (verbose)
	cout << "intensityErr(" << massBin->waveName(i) << "): "
	     << setprecision(12) << newVal << " vs. " << oldVal << ", delta = " << oldVal - newVal << endl;
    }
    cout << "intensityErr() max. deviation = " << maxDelta << endl << endl;
  }

  if (1) {
    double maxDelta = 0;
    const string waveNamePatterns[] = {"",  // total intensity
				       "flat",
				       "0++0-",
				       "0-+0+",
				       "1++0+",
				       "2-+0+",
				       "2++0-",
				       "2++1+",
				       "1-+0-",
				       "1-+1+",
				       "3++0+",
				       "4++1+",
				       "3-+1+",
				       "3-+1-",
				       "3-+0-"};
    for (unsigned int i = 0; i < sizeof(waveNamePatterns) / sizeof(string); ++i) {
      const double oldVal = massBin->intens(waveNamePatterns[i].c_str());
      const double newVal = massBin->intensity(waveNamePatterns[i]);
      const double delta  = oldVal - newVal;
      maxDelta = (fabs(maxDelta) < fabs(delta)) ? delta : maxDelta;
      if (verbose)
	cout << "intensity(" << waveNamePatterns[i] << "): "
	     << setprecision(12) << newVal << " vs. " << oldVal << ", delta = " << oldVal - newVal << endl;
    }
    cout << "intensity() max. deviation = " << maxDelta << endl << endl;
  }

  if (1) {
    double maxDelta = 0;
    const string waveNamePatterns[] = {"",  // total intensity
				       "flat",
				       "0++0-",
				       "0-+0+",
				       "1++0+",
				       "2-+0+",
				       "2++0-",
				       "2++1+",
				       "1-+0-",
				       "1-+1+",
				       "3++0+",
				       "4++1+",
				       "3-+1+",
				       "3-+1-",
				       "3-+0-"};
    for (unsigned int i = 0; i < sizeof(waveNamePatterns) / sizeof(string); ++i) {
      const double oldVal = massBin->err(waveNamePatterns[i].c_str());
      const double newVal = massBin->intensityErr(waveNamePatterns[i]);
      const double delta  = oldVal - newVal;
      maxDelta = (fabs(maxDelta) < fabs(delta)) ? delta : maxDelta;
      if (verbose)
	cout << "intensityErr(" << waveNamePatterns[i] << "): "
	     << setprecision(12) << newVal << " vs. " << oldVal << ", delta = " << oldVal - newVal << endl;
    }
    cout << "intensityErr() max. deviation = " << maxDelta << endl << endl;
  }

  if (0) {
    double maxDelta = 0;
    for (unsigned int i = 0; i < n; ++i)
      for (unsigned int j = 0; j < n; ++j) {
	const double oldVal = massBin->phase(i, j);
	const double newVal = massBin->phaseNew(i, j);
	const double delta  = oldVal - newVal;
	maxDelta = (fabs(maxDelta) < fabs(delta)) ? delta : maxDelta;
	if (verbose)
	  cout << "phaseNew(" << massBin->waveName(i) << ", "  << massBin->waveName(j) << "): "
	       << setprecision(12) << newVal << " vs. " << oldVal << ", delta = " << delta << endl;
      }
    cout << "phaseNew() max. deviation = " << maxDelta << endl << endl;
  }

}
