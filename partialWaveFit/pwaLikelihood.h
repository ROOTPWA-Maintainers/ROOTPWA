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


#ifndef PWALIKELIHOOD_H
#define PWALIKELIHOOD_H


#include <vector>
#include <string>
#include <iostream>

#define BOOST_DISABLE_ASSERTS
#include "boost/multi_array.hpp"
#include "boost/tuple/tuple.hpp"

#include "Math/IFunction.h"
#include "TFile.h"
#include "TH1.h"
#include "TMatrixT.h"

#include "sumAccumulators.hpp"
#include "ampIntegralMatrix.h"


class TString;


namespace rpwa {


	class complexMatrix;


	template<typename complexT>  // type of internal variables used for intermediate results
	class pwaLikelihood : public ROOT::Math::IGradientFunctionMultiDim {

	public:

		typedef typename complexT::value_type value_type;

		// define array types
		typedef boost::multi_array<std::string,                    2> waveNameArrayType;    // array for wave names
		typedef boost::multi_array<double,                         2> waveThrArrayType;     // array for wave thresholds
		typedef boost::multi_array<unsigned int,                   2> waveToIntMapType;     // array for mapping of waves to integral indices
		typedef boost::multi_array<boost::tuples::tuple<int, int>, 3> ampToParMapType;      // array for mapping of amplitudes to parameters
		typedef boost::multi_array<complexT,                       3> ampsArrayType;        // array for production and decay amplitudes
		typedef boost::multi_array<complexT,                       4> normMatrixArrayType;  // array for normalization matrices
		typedef boost::multi_array<value_type,                     2> phaseSpaceIntType;    // array for phase space integrals


	public:

		// enum for function call counters
		enum functionCallEnum {
			FDF                  = 0,
			GRADIENT             = 1,
			DOEVAL               = 2,
			DODERIVATIVE         = 3,
			HESSIAN              = 4,
			NMB_FUNCTIONCALLENUM = 5
		};

		struct functionCallInfo {
			typedef boost::accumulators::accumulator_set
			  <double, boost::accumulators::stats
				<boost::accumulators::tag::sum(boost::accumulators::compensated)> > timeAccType;
			unsigned int nmbCalls;   // number of times function was called
			timeAccType  funcTime;   // time needed to calculate function value(s) (w/o normalization)
			timeAccType  normTime;   // time needed to normalize function value(s)
			timeAccType  totalTime;  // total execution time of function
		};

		enum priorEnum {
			FLAT,
			HALF_CAUCHY
		};

		pwaLikelihood();
		~pwaLikelihood();

		// overload public IGradientFunctionMultiDim member functions:
		/// clones the function using the default copy constructor
		virtual pwaLikelihood* Clone() const { return new pwaLikelihood(*this); }
		/// returns total number of function parameters (= dimension of the function)
		virtual unsigned int NDim() const { return nmbPars(); }
		/// calculates gradient (vector of partial derivatives) of function at point defined by par
		virtual void Gradient(const double* par,
		                      double*       gradient) const;
		/// optimized method to evaluate function value and derivative at a point defined by par at the same time
		virtual void FdF(const double* par,
		                 double&       funcVal,
		                 double*       gradient) const;

		// overload private IGradientFunctionMultiDim member functions
		virtual double DoEval      (const double* par) const;
		virtual double DoDerivative(const double* par,
		                            unsigned int  derivativeIndex) const;

		/// calculates Hessian of function at point defined by par
		TMatrixT<double> Hessian(const double* par) const;

		/// calculates covariance matrix of function at point defined by par
		TMatrixT<double> CovarianceMatrix(const double* par) const;
		/// turns hessian into covariance matrix
		TMatrixT<double> CovarianceMatrix(const TMatrixT<double>& hessian) const;

		/// flips the signs of the paramaters according to conventions (amplitudes of each anchor wave and the flat wave are real and positive)
		std::vector<double> CorrectParamSigns(const double* par) const;

		unsigned int             nmbEvents   ()                                    const { return _nmbEvents;               }  ///< returns number of events that enter in the likelihood
		unsigned int             rank        ()                                    const { return _rank;                    }  ///< returns rank of spin density matrix
		unsigned int             nmbWaves    (const int          reflectivity = 0) const;                                      ///< returns total number of waves (reflectivity == 0) or number or number of waves with positive/negative reflectivity; flat wave is not counted!
		unsigned int             nmbPars     ()                                    const { return _nmbPars;                 }  ///< returns total number of parameters
		unsigned int             nmbParsFixed()                                    const { return _nmbParsFixed;            }  ///< returns number of fixed parameters
		std::string              parName     (const unsigned int parIndex)         const { return _parNames[parIndex];      }  ///< returns name of likelihood parameter at parIndex
		double                   parThreshold(const unsigned int parIndex)         const { return _parThresholds[parIndex]; }  ///< returns threshold in GeV/c^2 above which likelihood parameter at parIndex becomes free
		bool                     parFixed    (const unsigned int parIndex)         const { return _parFixed[parIndex];      }  ///< returns whether likelihood parameter at parIndex is fixed due to mass threshold

		double dLcache(const unsigned int i) const { return _derivCache[i]; }
		unsigned int ncalls(const functionCallEnum func = FDF) const
		{ return _funcCallInfo[func].nmbCalls; }
		double Ltime(const functionCallEnum func = FDF) const
		{ return boost::accumulators::sum(_funcCallInfo[func].funcTime); }
		double Ntime(const functionCallEnum func = FDF) const
		{ return boost::accumulators::sum(_funcCallInfo[func].normTime); }
		//const integral& normInt() const { return _normInt; }

		// modifiers
		void          enableCuda       (const bool      enableCuda = true);
		bool          cudaEnabled      () const;
		void          useNormalizedAmps(const bool      useNorm    = true) { _useNormalizedAmps = useNorm;   }
		void          setPriorType     (const priorEnum priorType  = FLAT) { _priorType         = priorType; }
		priorEnum     priorType        () const                            { return _priorType;              }
		void          setCauchyWidth   (const double&   cauchyWidth)       { _cauchyWidth = cauchyWidth;     }
		const double& cauchyWidth      ()                                  { return _cauchyWidth;            }
		static void   setQuiet         (const bool      flag       = true) { _debug             = !flag;     }

		// operations
		void init(const unsigned int rank,
		          const std::map<std::string, std::string>& ampFileList,
		          const double       massBinCenter,
		          const std::string& waveListFileName,
		          const std::string& normIntFileName,
		          const std::string& accIntFileName,
		          const unsigned int numbAccEvents = 0);  ///< prepares all internal data structures

		void getIntegralMatrices(rpwa::complexMatrix&       normMatrix,
		                         rpwa::complexMatrix&       accMatrix,
		                         std::vector<double>&       phaseSpaceIntegral) const;

		// note: amplitudes which do not exist in higher ranks are NOT built!
		void buildProdAmpArrays(const double*                       inPar,
		                        std::vector<std::complex<double> >& prodAmps,
		                        std::vector<std::pair<int,int> >&   parIndices,
		                        std::vector<std::string>&           prodAmpNames,
		                        const bool                          withFlat = false) const;

		std::ostream& print(std::ostream& out = std::cout) const;
		std::ostream& printFuncInfo(std::ostream& out = std::cout) const;
		friend std::ostream& operator << (std::ostream&         out,
		                                  const pwaLikelihood& func) { return func.print(out); }


	private:

		// helper functions
		void readWaveList       (const std::string& waveListFileName);  ///< reads wave names and thresholds from wave list file
		void buildParDataStruct (const unsigned int rank,
		                         const double       massBinCenter);     ///< builds parameter data structures
		void readIntegrals      (const std::string& normIntFileName,
		                         const std::string& accIntFileName,
		                         const std::string& integralTKeyName = "integral");  ///< reads normalization and acceptance integrals from file

		void readDecayAmplitudes(const std::map<std::string, std::string>& ampFileList,
		                         const std::string& ampLeafName = "amplitude");  ///< reads decay amplitudes from files in specified directory


		void clear();

		void reorderIntegralMatrix(const rpwa::ampIntegralMatrix& integral,
		                           normMatrixArrayType&           reorderedMatrix) const;

	public:

		void copyFromParArray(const double*  inPar,              // input parameter array
		                      ampsArrayType& outVal,             // output values organized as 3D array of complex numbers with [rank][reflectivity][wave index]
		                      value_type&    outFlatVal) const;  // output value corresponding to flat wave
		void copyToParArray(const ampsArrayType& inVal,          // values corresponding to production amplitudes [rank][reflectivity][wave index]
		                    const value_type     inFlatVal,      // value corresponding to flat wave
		                    double*              outPar) const;  // output parameter array

	private:

		void resetFuncCallInfo() const;

		unsigned int _nmbEvents;        // number of events
		unsigned int _rank;             // rank of spin density matrix
		unsigned int _nmbWaves;         // number of waves
		unsigned int _nmbWavesRefl[2];  // number of negative (= 0) and positive (= 1) reflectivity waves
		unsigned int _nmbWavesReflMax;  // maximum of number of negative and positive reflectivity waves
		unsigned int _nmbPars;          // number of function parameters
		unsigned int _nmbParsFixed;     // number of fixed function parameters

	#ifdef USE_CUDA
		bool                _cudaEnabled;        // if true CUDA kernels are used for some calculations
	#endif
		bool                _useNormalizedAmps;  // if true normalized amplitudes are used
		priorEnum           _priorType;          // which prior to apply to parameters
		double              _cauchyWidth;        // width for the half-Cauchy prior
		static bool         _debug;              // if true debug messages are printed

		unsigned int _numbAccEvents; // number of input events used for acceptance integrals (accepted + rejected!)
		double       _totAcc;        // total acceptance in this bin

		waveNameArrayType        _waveNames;            // wave names [reflectivity][wave index]
		waveThrArrayType         _waveThresholds;       // mass thresholds of waves
		std::vector<std::string> _parNames;             // function parameter names
		std::vector<double>      _parThresholds;        // mass thresholds of parameters
		std::vector<bool>        _parFixed;             // parameter fixed due to mass thresholds
		ampToParMapType          _prodAmpToFuncParMap;  // maps each production amplitude to the indices
		                                                // of its real and imginary part in the parameter
		                                                // array; negative indices mean that the parameter
		                                                // is not existing due to rank restrictions

		ampsArrayType _decayAmps;  // precalculated decay amplitudes [event index][reflectivity][wave index]

		mutable std::vector<double> _parCache;    // parameter cache for derivative calc.
		mutable std::vector<double> _derivCache;  // cache for derivatives

		// normalization integrals
		normMatrixArrayType _normMatrix;          // normalization matrix w/o acceptance [reflectivity 1][wave index 1][reflectivity 2][wave index 2]
		normMatrixArrayType _accMatrix;           // normalization matrix with acceptance [reflectivity 1][wave index 1][reflectivity 2][wave index 2]
		phaseSpaceIntType   _phaseSpaceIntegral;  // phase space integrals

		mutable functionCallInfo _funcCallInfo[NMB_FUNCTIONCALLENUM];  // collects function call statistics

	};


}

#endif  // PWALIKELIHOOD_H
