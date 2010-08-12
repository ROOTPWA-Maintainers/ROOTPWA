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
#include <utility>
#include <iostream>

#define BOOST_DISABLE_ASSERTS
#include "boost/multi_array.hpp"
#include "boost/tuple/tuple.hpp"

#include "Math/IFunction.h"

// PWA2000 classes
#include "integral.h"
#include "matrix.h"
#include "../cuda/cuda_helper.h"

class TString;
class TCMatrix;


template<typename T = double>  // template for type of internal variables used for intermediate results
class TPWALikelihood : public ROOT::Math::IGradientFunctionMultiDim {

private:

	// define array types
	typedef boost::multi_array<std::string,            2> waveNameArrayType;    // array for wave names
	typedef boost::multi_array<double,                 2> waveThrArrayType;     // array for wave thresholds
	typedef boost::multi_array<unsigned int,           2> waveToIntMapType;     // array for mapping of waves to integral indices
	typedef boost::multi_array<unsigned int,           2> waveToListMapType;    // array for mapping of waves to position in wave list
	typedef boost::multi_array<boost::tuple<int, int>, 3> ampToParMapType;      // array for mapping of amplitudes to parameters
	typedef boost::multi_array<std::complex<double>,   3> ampsArrayType;        // array for production and decay amplitudes
	typedef boost::multi_array<std::complex<double>,   4> normMatrixArrayType;  // array for normalization matrices


public:

	// enum for function call counters
	enum functionCallEnum {
		FDF                  = 0,
		GRADIENT             = 1,
		DOEVAL               = 2,
		DODERIVATIVE         = 3,
		NMB_FUNCTIONCALLENUM = 4
	};

	TPWALikelihood();
	~TPWALikelihood();

	// overload public IGradientFunctionMultiDim member functions:
	/// clones the function using the default copy constructor
	virtual TPWALikelihood* Clone() const { return new TPWALikelihood(*this); }
	/// returns total number of function parameters (= dimension of the function)
	virtual unsigned int NDim() const { return nmbPars(); }
	/// optimized method to evaluate function value and derivative at a point defined by par at the same time
	virtual void FdF(const double* par,
	                 double&       funcVal,
	                 double*       gradient) const;
	/// calculates gradient (vector of partial derivatives) of function at point defined by par
	virtual void Gradient(const double* par,
	                      double*       gradient) const;
	
	unsigned int             nmbEvents   ()                                    const { return _nmbEvents;               }  ///< returns number of events that enter in the likelihood
	unsigned int             rank        ()                                    const { return _rank;                    }  ///< returns rank of spin density matrix
	inline unsigned int      nmbWaves    (const int          reflectivity = 0) const;                                      ///< returns total number of waves (reflectivity == 0) or number or number of waves with positive/negative reflectivity; flat wave is not counted!
	unsigned int             nmbPars     ()                                    const { return _nmbPars;                 }  ///< returns total number of parameters
	std::string              waveName    (const unsigned int waveIndex)        const { return _waveNames[waveIndex];    }  ///< returns name of wave at waveIndex
	std::vector<std::string> waveNames   ()                                    const;  ///< returns vector with all wave names ordered like in input wave list
	std::string              parName     (const unsigned int parIndex)         const { return _parNames[parIndex];      }  ///< returns name of likelihood parameter at parIndex
	double                   parThreshold(const unsigned int parIndex)         const { return _parThresholds[parIndex]; }  ///< returns threshold in GeV/c^2 above which likelihood parameter at parIndex becomes free
  
	double dLcache(const unsigned int i) const { return _derivCache[i]; }
	unsigned int ncalls(const functionCallEnum callType = FDF) const { return _nmbCalls[FDF]; }
	double Ltime() const { return _Ltime; }
	double Ntime() const { return _Ntime; }
	//const integral& normInt() const { return _normInt; }

	// modifiers
	void        useNormalizedAmps(const bool useNorm = true) { _useNormalizedAmps = useNorm; }
	static void setQuiet         (const bool flag    = true) { _debug             = !flag;   }

	// operations
	void init(const unsigned int rank,
	          const std::string& waveListFileName,
	          const std::string& normIntFileName,
	          const std::string& accIntFileName,
	          const std::string& ampDirName    = ".",
	          const unsigned int numbAccEvents = 0);  ///< prepares all internal data structures

	void getIntCMatrix(TCMatrix& normMatrix,
	                   TCMatrix& accMatrix) const;
  
	// note: amplitudes which do not exist in higher ranks are NOT built!
	void buildCAmps(const double*                       inPar,
	                std::vector<std::complex<double> >& prodAmps,
	                std::vector<std::pair<int,int> >&   parIndices,
	                std::vector<std::string>&           prodAmpNames,
	                const bool                          withFlat = false) const;
	
	std::ostream& print(std::ostream& out = std::cout) const;
	friend std::ostream& operator << (std::ostream&         out,
	                                  const TPWALikelihood& func) { return func.print(out); }
	
	// overload private IGradientFunctionMultiDim member functions
	virtual double DoEval      (const double* par) const;
	virtual double DoEvalCuda  (const double* par) const;
	virtual double DoDerivative(const double* par,
	                            unsigned int  derivativeIndex) const;

	std::vector<unsigned int> orderedParIndices() const;  // helper function for backwards-compatibility

private:

	// helper functions
	void readWaveList       (const std::string& waveListFileName);  ///< reads wave names and thresholds from wave list file
	void buildParDataStruct (const unsigned int rank);              ///< builds parameter data structures
	void readIntegrals      (const std::string& normIntFileName,
	                         const std::string& accIntFileName);    ///< reads normalization and acceptance integrals from file
	void readDecayAmplitudes(const std::string& ampDirName = ".");  ///< reads decay amplitudes from files in specified directory
  
  
	void clear();
	static int getReflectivity(const TString& waveName);
  
	void reorderIntegralMatrix(integral&            integral,
	                           normMatrixArrayType& reorderedMatrix) const;
	void copyFromParArray(const double*  inPar,              // input parameter array
	                      ampsArrayType& outVal,             // output values organized as 3D array of complex numbers with [rank][reflectivity][wave index]
	                      T&             outFlatVal) const;  // output value corresponding to flat wave
	void copyToParArray(const ampsArrayType& inVal,          // values corresponding to production amplitudes [rank][reflectivity][wave index]
	                    const T              inFlatVal,      // value corresponding to flat wave
	                    double*              outPar) const;  // output parameter array

	unsigned int _nmbEvents;        // number of events
	unsigned int _rank;             // rank of spin density matrix
	unsigned int _nmbWaves;         // number of waves
	unsigned int _nmbWavesRefl[2];  // number of negative (= 0) and positive (= 1) reflectivity waves 
	unsigned int _nmbPars;          // number of function parameters

	mutable unsigned int _nmbCalls[4];  // function call counters
	mutable double       _Ltime;        // total time spent calculating L
	mutable double       _Ntime;        // total time spent calculating normalization

	bool        _useNormalizedAmps;  // if true normalized amplitudes are used
	static bool _debug;              // if true debug messages are printed

	unsigned int _numbAccEvents; // number of input events used for acceptance integrals (accepted + rejected!)

	waveNameArrayType        _waveNames;            // wave names [reflectivity][wave index]
	waveThrArrayType         _waveThresholds;       // mass thresholds of waves
	waveToListMapType        _waveToWaveList;       // maps wave to its index in wave list
	std::vector<std::string> _parNames;             // function parameter names
	std::vector<double>      _parThresholds;        // mass thresholds of parameters
	ampToParMapType          _prodAmpToFuncParMap;  // maps each production amplitude to the indices of
	// its real and imginary part in the parameter
	// array; negative indices mean that the parameter
	// is not existing due to rank restrictions
  
	ampsArrayType _decayAmps;  // precalculated decay amplitudes [event index][reflectivity][wave index]
  
	rpwa::complex<double>* _cudaDecayAmps;  // precalculated decay amplitudes [event index][reflectivity][wave index]
	unsigned int _num_threads,_num_blocks;

	mutable std::vector<double> _parCache;    // parameter cache for derivative calc.
	mutable std::vector<double> _derivCache;  // cache for derivatives
  
	// normalization integrals 
	normMatrixArrayType _normMatrix;  // normalization matrix w/o acceptance [reflectivity 1][wave index 1][reflectivity 2][wave index 2]
	normMatrixArrayType _accMatrix;   // normalization matrix with acceptance [reflectivity 1][wave index 1][reflectivity 2][wave index 2]
  
};


template<typename T>
unsigned int
TPWALikelihood<T>::nmbWaves(const int reflectivity) const
{
	if (reflectivity == 0)
		return _nmbWaves;
	else if (reflectivity > 0)
		return _nmbWavesRefl[1];  // positive reflectivity
	else
		return _nmbWavesRefl[0];  // negative reflectivity
}


#include "TPWALikelihood.cc"


#endif  // TPWALIKELIHOOD_HH


//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
