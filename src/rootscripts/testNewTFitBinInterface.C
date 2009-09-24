#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>

#include "TTree.h"
#include "TMatrixT.h"
#include "TComplex.h"

#include "../TFitBin.h"


using namespace std;


//#define OLD


inline
ostream&
operator << (ostream&        o,
	     const TComplex& c)
{
  o << "(" << c.Re() << ", " << c.Im() << ")";
  return o;
}

template <typename T>
ostream&
operator << (ostream&           o,
	     const TMatrixT<T>& A)
{
  o << "(";
  for (int row = 0; row < A.GetNrows(); ++row) {
    o << "(";
    for (int col = 0; col < A.GetNcols(); ++col) {
      o << A[row][col];
      if (col < A.GetNcols() - 1)
	o << ", ";
    }
    if (row < A.GetNrows() - 1)
      o << "), ";
    else
      o << ")";
  }
  o << ")";
  return o;
}


void
testNewTFitBinInterface(TTree* tree)
{
  const bool verbose = true;
  TFitBin*   massBin = new TFitBin();
  tree->SetBranchAddress("fitbin", &massBin);
  tree->GetEntry(30);
  //massBin->printWaveNames();
  //massBin->printProdAmpNames();
  const unsigned int n = massBin->nmbWaves();

  if (1) {
    complex<double> maxDelta = 0;
    for (unsigned int i = 0; i < n; ++i)
      for (unsigned int j = 0; j < n; ++j) {
	const TComplex        temp   = massBin->spinDens(i, j);
	const complex<double> oldVal = complex<double>(temp.Re(), temp.Im());
	const complex<double> newVal = massBin->spinDensityMatrixElem(i, j);
	const complex<double> delta  = oldVal - newVal;
	maxDelta.real() = (maxDelta.real() < delta.real()) ? delta.real() : maxDelta.real();
	maxDelta.imag() = (maxDelta.imag() < delta.imag()) ? delta.imag() : maxDelta.imag();
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
      const double newVal = massBin->spinDensityMatrixElem(i, i).real();
      const double delta  = oldVal - newVal;
      maxDelta = (maxDelta < delta) ? delta : maxDelta;
      if (verbose)
	cout << "intensity(" << massBin->waveName(i) << ") via spinDensityMatrixElem() :"
	     << setprecision(12) << newVal << " vs. " << oldVal << ", delta = " << delta << endl;
    }
    cout << "intensity via spinDensityMatrixElem() max. deviation = " << maxDelta << endl << endl;
  }

  if (1) {
    double maxDelta = 0;
    for (unsigned int i = 0; i < n; ++i) {
      const double oldVal = massBin->err(i);
      const double newVal = sqrt(massBin->spinDensityMatrixElemCov(i, i)[0][0]);
      const double delta  = oldVal - newVal;
      maxDelta = (maxDelta < delta) ? delta : maxDelta;
      if (verbose)
	cout << "intensityErr(" << massBin->waveName(i) << ") via spinDensityMatrixElemCov() :"
	     << setprecision(12) << newVal << " vs. " << oldVal << ", delta = " << delta << endl;
    }
    cout << "intensityErr via spinDensityMatrixElemCov() max. deviation = " << maxDelta << endl << endl;
  }

  if (1) {
    double maxDelta = 0;
    for (unsigned int i = 0; i < n; ++i) {
      const double oldVal = massBin->intens(i);
      const double newVal = massBin->intensity(i);
      const double delta  = oldVal - newVal;
      maxDelta = (maxDelta < delta) ? delta : maxDelta;
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
      maxDelta = (maxDelta < delta) ? delta : maxDelta;
      if (verbose)
      cout << "intensityErr(" << massBin->waveName(i) << "): "
	   << setprecision(12) << newVal << " vs. " << oldVal << ", delta = " << oldVal - newVal << endl;
    }
    cout << "intensityErr() max. deviation = " << maxDelta << endl << endl;
  }
}
