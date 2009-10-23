#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>

#include "TTree.h"
#include "TMatrixT.h"
#include "TComplex.h"

#include "../TFitBin.h"
#include "../TFitResult.h"
#include "../utilities.h"


using namespace std;


void
testTFitResult(TTree* tree)
{
  const bool verbose = true;
  //const bool verbose = false;
  const bool copyBin = false;

  TFitBin* massBin = new TFitBin();
  tree->SetBranchAddress("fitbin", &massBin);
  TFitResult* newMassBin = 0;
  if (!copyBin) {
    newMassBin = new TFitResult();
    tree->SetBranchAddress("fitResult", &newMassBin);
  }
  //tree->GetEntry(30);
  tree->GetEntry(0);
  
  if (copyBin)
    newMassBin = new TFitResult(*massBin);
  const unsigned int n = newMassBin->nmbWaves();

  if (1)
    cout << *newMassBin << endl;
  
  if (0) {
    complex<double> maxDelta = 0;
    for (unsigned int i = 0; i < n; ++i)
      for (unsigned int j = 0; j < n; ++j) {
	const TComplex        temp   = massBin->spinDens(i, j);
	const complex<double> oldVal = complex<double>(temp.Re(), temp.Im());
	const complex<double> newVal = newMassBin->spinDensityMatrixElem(i, j);
	const complex<double> delta  = oldVal - newVal;
	maxDelta.real() = (fabs(maxDelta.real()) < fabs(delta.real())) ? delta.real() : maxDelta.real();
	maxDelta.imag() = (fabs(maxDelta.imag()) < fabs(delta.imag())) ? delta.imag() : maxDelta.imag();
	if (verbose)
	  cout << "spinDensityMatrixElem(" << newMassBin->waveName(i) << ", "  << newMassBin->waveName(j) << "): "
	       << setprecision(12) << newVal << " vs. " << oldVal << ", delta = " << delta << endl;
      }
    cout << "spinDensityMatrixElem() max. deviation = " << maxDelta << endl << endl;
  }

  if (1) {
    double maxDelta = 0;
    for (unsigned int i = 0; i < n; ++i) {
      const double oldVal = massBin->intens(i);
      const double newVal = newMassBin->intensity(i);
      const double delta  = oldVal - newVal;
      maxDelta = (fabs(maxDelta) < fabs(delta)) ? delta : maxDelta;
      if (verbose)
	cout << "intensity(" << newMassBin->waveName(i) << "): "
	     << setprecision(12) << newVal << " vs. " << oldVal << ", delta = " << oldVal - newVal << endl;
    }
    cout << "intensity() max. deviation = " << maxDelta << endl << endl;
  }

  if (1) {
    double maxDelta = 0;
    for (unsigned int i = 0; i < n; ++i) {
      const double oldVal = massBin->err(i);
      const double newVal = newMassBin->intensityErr(i);
      const double delta  = oldVal - newVal;
      maxDelta = (fabs(maxDelta) < fabs(delta)) ? delta : maxDelta;
      if (verbose)
	cout << "intensityErr(" << newMassBin->waveName(i) << "): "
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
      const double newVal = newMassBin->intensity(waveNamePatterns[i].c_str());
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
      const double newVal = newMassBin->intensityErr(waveNamePatterns[i].c_str());
      const double delta  = oldVal - newVal;
      maxDelta = (fabs(maxDelta) < fabs(delta)) ? delta : maxDelta;
      if (verbose)
	cout << "intensityErr(" << waveNamePatterns[i] << "): "
	     << setprecision(12) << newVal << " vs. " << oldVal << ", delta = " << oldVal - newVal << endl;
    }
    cout << "intensityErr() max. deviation = " << maxDelta << endl << endl;
  }

  if (1) {
    double maxDelta = 0;
    for (unsigned int i = 0; i < n; ++i)
      for (unsigned int j = 0; j < n; ++j) {
	const double oldVal = massBin->phase(i, j);
	const double newVal = newMassBin->phase(i, j);
	const double delta  = oldVal - newVal;
	maxDelta = (fabs(maxDelta) < fabs(delta)) ? delta : maxDelta;
	if (verbose)
	  cout << "phase(" << newMassBin->waveName(i) << ", "  << newMassBin->waveName(j) << "): "
	       << setprecision(12) << newVal << " vs. " << oldVal << ", delta = " << delta << endl;
      }
    cout << "phase() max. deviation = " << maxDelta << endl << endl;
  }

  if (0) {
    for (unsigned int i = 0; i < n; ++i)
      for (unsigned int j = 0; j < n; ++j) {
	const double coh    = newMassBin->coherence(i, j);
	const double cohErr = newMassBin->coherenceErr(i, j);
	cout << "coh(" << newMassBin->waveName(i) << ", "  << newMassBin->waveName(j) << "): "
	       << setprecision(12) << coh << " +- " << cohErr << endl;
      }
  }

}
