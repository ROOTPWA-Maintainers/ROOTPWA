#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>

#include "TTree.h"
#include "TMatrixT.h"
#include "TComplex.h"

#include "reportingUtils.hpp"
#include "TFitResult.h"
#include "fitResult.h"


using namespace std;
using namespace rpwa;


#if TFITRESULT_ENABLED


void
testFitResult(TTree*        oldTree,
	      TTree*        newTree,
	      const bool    verbose       = true,
	      const string& oldBranchName = "fitResult",
	      const string& newBranchName = "fitResult_v2")
{

  const bool copyBin = false;

  TFitResult* oldResult = 0;
  oldTree->SetBranchAddress(oldBranchName.c_str(), &oldResult);
  oldTree->GetEntry(0);
  fitResult* newResult = 0;
  if (!copyBin) {
    newTree->SetBranchAddress(newBranchName.c_str(), &newResult);
    newTree->GetEntry(0);
  } else
    newResult = new fitResult(*oldResult);
  const unsigned int n = newResult->nmbWaves();

  if (0)
    cout << *newResult << endl;
  
  if (1) {
    complex<double> maxDelta = 0;
    for (unsigned int i = 0; i < n; ++i)
      for (unsigned int j = 0; j < n; ++j) {
	const complex<double> oldVal = oldResult->spinDensityMatrixElem(i, j);
	const complex<double> newVal = newResult->spinDensityMatrixElem(i, j);
	const complex<double> delta  = oldVal - newVal;
	maxDelta.real() = (fabs(maxDelta.real()) < fabs(delta.real())) ? delta.real() : maxDelta.real();
	maxDelta.imag() = (fabs(maxDelta.imag()) < fabs(delta.imag())) ? delta.imag() : maxDelta.imag();
	if (verbose)
	  cout << "spinDensityMatrixElem(" << newResult->waveName(i) << ", "  << newResult->waveName(j) << "): "
	       << setprecision(12) << newVal << " vs. " << oldVal << ", delta = " << delta << endl;
      }
    cout << "spinDensityMatrixElem() max. deviation = " << maxDelta << endl << endl;
  }

  if (1) {
    double maxDelta = 0;
    for (unsigned int i = 0; i < n; ++i) {
      const double oldVal = oldResult->intensity(i);
      const double newVal = newResult->intensity(i);
      const double delta  = oldVal - newVal;
      maxDelta = (fabs(maxDelta) < fabs(delta)) ? delta : maxDelta;
      if (verbose)
	cout << "intensity(" << newResult->waveName(i) << "): "
	     << setprecision(12) << newVal << " vs. " << oldVal << ", delta = " << oldVal - newVal << endl;
    }
    cout << "intensity() max. deviation = " << maxDelta << endl << endl;
  }

  if (1) {
    double maxDelta = 0;
    for (unsigned int i = 0; i < n; ++i) {
      const double oldVal = oldResult->intensityErr(i);
      const double newVal = newResult->intensityErr(i);
      const double delta  = oldVal - newVal;
      maxDelta = (fabs(maxDelta) < fabs(delta)) ? delta : maxDelta;
      if (verbose)
	cout << "intensityErr(" << newResult->waveName(i) << "): "
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
      const double oldVal = oldResult->intensity(waveNamePatterns[i].c_str());
      const double newVal = newResult->intensity(waveNamePatterns[i].c_str());
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
      const double oldVal = oldResult->intensityErr(waveNamePatterns[i].c_str());
      const double newVal = newResult->intensityErr(waveNamePatterns[i].c_str());
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
	const double oldVal = oldResult->phase(i, j);
	const double newVal = newResult->phase(i, j);
	const double delta  = oldVal - newVal;
	maxDelta = (fabs(maxDelta) < fabs(delta)) ? delta : maxDelta;
	if (verbose)
	  cout << "phase(" << newResult->waveName(i) << ", "  << newResult->waveName(j) << "): "
	       << setprecision(12) << newVal << " vs. " << oldVal << ", delta = " << delta << endl;
      }
    cout << "phase() max. deviation = " << maxDelta << endl << endl;
  }

  if (1) {
    double maxDelta = 0;
    for (unsigned int i = 0; i < n; ++i)
      for (unsigned int j = 0; j < n; ++j) {
	const double oldVal = oldResult->phaseErr(i, j);
	const double newVal = newResult->phaseErr(i, j);
	const double delta  = oldVal - newVal;
	maxDelta = (fabs(maxDelta) < fabs(delta)) ? delta : maxDelta;
	if (verbose)
	  cout << "phaseErr(" << newResult->waveName(i) << ", "  << newResult->waveName(j) << "): "
	       << setprecision(12) << newVal << " vs. " << oldVal << ", delta = " << delta << endl;
      }
    cout << "phaseErr() max. deviation = " << maxDelta << endl << endl;
  }

  if (1) {
    double maxDelta = 0;
    for (unsigned int i = 0; i < n; ++i)
      for (unsigned int j = 0; j < n; ++j) {
	const double oldVal = oldResult->coherence(i, j);
	const double newVal = newResult->coherence(i, j);
	const double delta  = oldVal - newVal;
	maxDelta = (fabs(maxDelta) < fabs(delta)) ? delta : maxDelta;
	if (verbose)
	  cout << "coherence(" << newResult->waveName(i) << ", "  << newResult->waveName(j) << "): "
	       << setprecision(12) << newVal << " vs. " << oldVal << ", delta = " << delta << endl;
      }
    cout << "coherence() max. deviation = " << maxDelta << endl << endl;
  }

  if (1) {
    double maxDelta = 0;
    for (unsigned int i = 0; i < n; ++i)
      for (unsigned int j = 0; j < n; ++j) {
	const double oldVal = oldResult->coherenceErr(i, j);
	const double newVal = newResult->coherenceErr(i, j);
	const double delta  = oldVal - newVal;
	maxDelta = (fabs(maxDelta) < fabs(delta)) ? delta : maxDelta;
	if (verbose)
	  cout << "coherenceErr(" << newResult->waveName(i) << ", "  << newResult->waveName(j) << "): "
	       << setprecision(12) << newVal << " vs. " << oldVal << ", delta = " << delta << endl;
      }
    cout << "coherenceErr() max. deviation = " << maxDelta << endl << endl;
  }
}

  
#endif  // TFITRESULT_ENABLED
