///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
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
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Likelihood function Object to use with ROOT minimizers
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef TPWALIKELIHOOD_HH
#define TPWALIKELIHOOD_HH


#include <vector>
#include <string>
#include <complex>

#include "Math/IFunction.h"

// PWA2000 classes
#include "integral.h"
#include "matrix.h"


class TString;
class TCMatrix;


//typedef matrix<complex<double> > intmat;


class TPWALikelihood : public ROOT::Math::IGradientFunctionMultiDim {

public:

  // enum for function call counters
  enum functionCallEnum {FDF          = 0,
			 GRADIENT     = 1,
			 DOEVAL       = 2,
			 DODERIVATIVE = 3};

  // Constructors/Destructors ---------
  TPWALikelihood();
  ~TPWALikelihood();

  // overload public IGradientFunctionMultiDim member functions
  virtual TPWALikelihood* Clone() const { return new TPWALikelihood(*this); }  // using default copy constructor
  virtual void FdF(const double* par,  // evaluate function and gradient at the same time
		   double&       funcVal,
		   double*       gradient) const;
  virtual void Gradient(const double* par,  // evaluate the full gradient vector at the vector value x
 			double*       gradient) const;
  virtual unsigned int NDim() const;

  // Accessors ------------------------
  std::string parname(const unsigned int i) const { return _parNames[i]; }
  const std::vector<std::string>& wavetitles() const { return _waveNames; }
  double parthreshold(const unsigned int i) const { return _parThresholds[i]; }
  double dLcache(const unsigned int i) const { return _derivCache[i]; }
  unsigned int ncalls(const functionCallEnum callType = FDF) const { return _nmbCalls[FDF]; }
  double Ltime() const { return _Ltime; }
  double Ntime() const { return _Ntime; }
  unsigned int nevents() const { return _nmbEvents; }
  const integral& normInt() const { return _normInt; }

  // Modifiers -----------------------
  void UseNormalizedAmps(const bool useNorm = true) { _useNorm = useNorm; }
  void SetWavelist(const TString& wavelist);
  void SetRank(const unsigned int rank);
  void SetQuiet(const bool flag = false) { _debug = !flag; }
  void SetMaxSampDL(const unsigned int samp);

  // Operations ----------------------
  // load Amplitudes into memory
  void LoadIntegrals(const TString& norm,
		     const TString& acceptance);
  void LoadAmplitudes();
  void getIntCMatrix(TCMatrix& integr,
		     TCMatrix& acceptance);

  // note: amps which do not exist in higher ranks are NOT built!
  void buildCAmps(const double*                       x,
		  std::vector<std::complex<double> >& V,
		  std::vector<std::pair<int,int> >&   indices,
		  std::vector<std::string>&           names,
		  const bool                          withFlat = false);

private:

  // overload private IGradientFunctionMultiDim member functions
  virtual double DoEval(const double* par) const;
  virtual double DoDerivative(const double* par,
			      unsigned int  derivativeIndex) const;

  void clearCache();
  int getReflectivity(const TString& waveName) const;

  matrix<complex<double> > reorderedIntegralMatrix(integral& integral) const;
  vector<vector<complex<double> > > copyFromParArray(const double* inPar,              // input parameter array
						     double&       outFlatVal) const;  // output value corresponding to flat wave
  void copyToParArray(const vector<vector<complex<double> > >& inVal,          // values corresponding to production amplitudes
		      const double                             inFlatVal,      // value corresponding to flat wave
		      double*                                  outPar) const;  // output parameter array

  unsigned int         _rank;             // rank of the spin density matrix
  unsigned int         _dim;              // number of function parameters
  unsigned int         _nmbEvents;        // number of events
  unsigned int         _nmbWaves;         // number of waves
  unsigned int         _nmbWavesPosRefl;  // number of positive reflectivity waves 
  unsigned int         _nmbWavesNegRefl;  // number of negative reflectivity waves
  mutable unsigned int _nmbCalls[4];      // function call counters
  unsigned int         _nmbEventsGrad;    // number of events used to calculate derivative
  mutable double       _Ltime;            // total time spent calculating L
  mutable double       _Ntime;            // total time spent calculating normalization

  bool _debug;    // if set debug messages are suppressed
  bool _useNorm;  // use normalized amplitudes

  std::vector<std::string>    _waveNames;       // wave names
  std::vector<std::string>    _parNames;        // function parameter names
  std::vector<int>            _waveRefl;        // reflectivities of waves
  std::vector<double>         _waveThresholds;  // mass thresholds of waves
  std::vector<double>         _parThresholds;   // mass thresholds of parameters
  mutable std::vector<double> _parCache;        // parameter cache for derivative calc.
  mutable std::vector<double> _derivCache;      // cache for derivatives

  std::vector<std::vector<std::complex<double> > > _decayAmps;  // data cache

  // normalization integrals 
  integral _normInt;
  matrix<std::complex<double> > _normMatrix;
  integral _accInt;
  mutable matrix<std::complex<double> > _accMatrix;

};


#endif  // TPWALIKELIHOOD_HH


//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
