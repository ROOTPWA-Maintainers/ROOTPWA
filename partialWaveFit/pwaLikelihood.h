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


#include <iostream>
#include <map>
#include <string>
#include <vector>

#define BOOST_DISABLE_ASSERTS
#include "boost/multi_array.hpp"
#include "boost/tuple/tuple.hpp"

#include "Math/IFunction.h"
#include "TMatrixT.h"
#include "TVectorT.h"

#include "ampIntegralMatrix.h"
#include "sumAccumulators.hpp"


class TString;


namespace rpwa {


	class complexMatrix;
	class eventMetadata;


	template<typename complexT>  // type of internal variables used for intermediate results
	class pwaLikelihood : public ROOT::Math::IGradientFunctionMultiDim {

	public:

		typedef typename complexT::value_type value_type;

		// define array types
		typedef boost::multi_array<std::string,                               2> waveNameArrayType;     // array for wave names
		typedef boost::multi_array<double,                                    2> waveThrArrayType;      // array for wave thresholds
		typedef boost::multi_array<boost::tuples::tuple<int, int>,            3> ampToParMapType;       // array for mapping of amplitudes to parameters
		typedef boost::multi_array<complexT,                                  3> prodAmpsArrayType;     // array for production and decay amplitudes
		typedef boost::multi_array<complexT,                                  2> decayAmpsArrayType;    // with memory layout to save memory
		typedef boost::multi_array<complexT,                                  4> normMatrixArrayType;   // array for normalization matrices
		typedef boost::multi_array<value_type,                                2> phaseSpaceIntType;     // array for phase space integrals
		typedef boost::multi_array<bool,                                      2> waveAmpAddedArrayType; // array for wave amplitudes read
		typedef std::map<std::string, std::pair<unsigned int, unsigned int>    > waveParamsType;        // map wave names to reflectivity and index in reflectivity
		typedef boost::tuples::tuple<std::string, rpwa::waveDescription, double> waveDescThresType;     // tuple for wave name, wave description and threshold


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

		// class to store fit parameter information
		class fitParameter {

		public:

			fitParameter()
				: _waveName(""), _rank(0), _threshold(0.), _fixed(true), _realPart(true) {}
			fitParameter(const std::string& waveName,
			             const unsigned int rank,
			             const double       threshold,
			             const bool         fixed,
			             const bool         realPart)
				: _waveName(waveName), _rank(rank), _threshold(threshold), _fixed(fixed), _realPart(realPart) {}

			bool operator==(const fitParameter& rhs) const {
				return (_waveName == rhs.waveName() and _rank == rhs.rank() and _threshold == rhs.threshold()
				        and _fixed == rhs.fixed() and _realPart == rhs.realPart());
			}

			const std::string& waveName () const { return _waveName;  }
			unsigned int       rank     () const { return _rank;      }
			double             threshold() const { return _threshold; }
			bool               fixed    () const { return _fixed;     }
			bool               realPart () const { return _realPart;  }

			std::string        parName  () const;

                private:

			std::string  _waveName;
			unsigned int _rank;
			double       _threshold;
			bool         _fixed;
			bool         _realPart;

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
		virtual double DoDerivative(const double*      par,
		                            const unsigned int derivativeIndex) const;

		/// calculates Hessian of function at point defined by par
		TMatrixT<double> Hessian(const double* par) const;
		/// calculates eigenvectors/-values of Hessian
		std::vector<std::pair<TVectorT<double>, double> > HessianEigenVectors(const TMatrixT<double>& hessian) const;

		/// calculates covariance matrix of function at point defined by par
		TMatrixT<double> CovarianceMatrix(const double* par) const;
		/// turns hessian into covariance matrix
		TMatrixT<double> CovarianceMatrix(const TMatrixT<double>& hessian) const;

		/// flips the signs of the paramaters according to conventions (amplitudes of each anchor wave and the flat wave are real and positive)
		std::vector<double> CorrectParamSigns(const double* par) const;

		unsigned int                     nmbEvents   ()                                    const { return _nmbEvents;               }  ///< returns number of events that enter in the likelihood
		unsigned int                     rank        ()                                    const { return _rank;                    }  ///< returns rank of spin density matrix
		unsigned int                     nmbWaves    (const int          reflectivity = 0) const;                                      ///< returns total number of waves (reflectivity == 0) or number or number of waves with positive/negative reflectivity; flat wave is not counted!
		unsigned int                     nmbPars     ()                                    const { return _nmbPars;                 }  ///< returns total number of parameters
		unsigned int                     nmbParsFixed()                                    const { return _nmbParsFixed;            }  ///< returns number of fixed parameters
		const fitParameter&              parameter   (const unsigned int parIndex)         const { return _parameters[parIndex];    }  ///< returns information about likelihood parameter at parIndex
		const std::vector<fitParameter>& parameters  ()                                    const { return _parameters;              }  ///< returns information about likelihood parameters

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
		void          setCauchyWidth   (const double    cauchyWidth)       { _cauchyWidth = cauchyWidth;     }
		double        cauchyWidth      () const                            { return _cauchyWidth;            }
		static void   setQuiet         (const bool      flag       = true) { _debug             = !flag;     }

		// operations
		bool init(const std::vector<waveDescThresType>& waveDescThres,
		          const unsigned int                    rank = 1,
		          const double                          massBinCenter = 0.);  ///< prepares all internal data structures

		bool addNormIntegral(const rpwa::ampIntegralMatrix& normMatrix);

		bool addAccIntegral(rpwa::ampIntegralMatrix& accMatrix, const unsigned int accEventsOverride = 0);

		bool addAmplitude(const std::vector<const rpwa::amplitudeMetadata*>& meta);

		bool finishInit();

		bool setOnTheFlyBinning(const std::map<std::string, std::pair<double, double> >& binningMap,
		                        const std::vector<const eventMetadata*>&                 evtMeta);

		void getIntegralMatrices(rpwa::complexMatrix&       normMatrix,
		                         rpwa::complexMatrix&       accMatrix,
		                         std::vector<double>&       phaseSpaceIntegral,
		                         const bool                 withFlat = false) const;

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

		void clear();

		// helper functions
		bool readWaveList      (const std::vector<waveDescThresType>& waveDescThres);  ///< reads wave names and thresholds from wave list file
		bool buildParDataStruct(const unsigned int rank,
		                        const double       massBinCenter);                     ///< builds parameter data structures

		void reorderIntegralMatrix(const rpwa::ampIntegralMatrix& integral,
		                           normMatrixArrayType&           reorderedMatrix) const;

	public:

		void copyFromParArray(const double*      inPar,              // input parameter array
		                      prodAmpsArrayType& outVal,             // output values organized as 3D array of complex numbers with [rank][reflectivity][wave index]
		                      value_type&        outFlatVal) const;  // output value corresponding to flat wave
		void copyToParArray(const prodAmpsArrayType& inVal,          // values corresponding to production amplitudes [rank][reflectivity][wave index]
		                    const value_type         inFlatVal,      // value corresponding to flat wave
		                    double*                  outPar) const;  // output parameter array

	private:

		void resetFuncCallInfo() const;

		unsigned int _nmbEvents;        // number of events
		unsigned int _rank;             // rank of spin density matrix
		unsigned int _nmbWaves;         // number of waves
		unsigned int _nmbWavesRefl[2];  // number of negative (= 0) and positive (= 1) reflectivity waves
		unsigned int _nmbWavesReflMax;  // maximum of number of negative and positive reflectivity waves
		unsigned int _nmbPars;          // number of function parameters
		unsigned int _nmbParsFixed;     // number of fixed function parameters
		bool         _initialized;      // was init method called?
		bool         _normIntAdded;     // was normalization integral matrix added?
		bool         _accIntAdded;      // was acceptance integral matrix added?
		bool         _initFinished;     // was initialization finished?

	#ifdef USE_CUDA
		bool                _cudaEnabled;        // if true CUDA kernels are used for some calculations
	#endif
		bool                _useNormalizedAmps;  // if true normalized amplitudes are used
		priorEnum           _priorType;          // which prior to apply to parameters
		double              _cauchyWidth;        // width for the half-Cauchy prior
		static bool         _debug;              // if true debug messages are printed

		unsigned int _numbAccEvents; // number of input events used for acceptance integrals (accepted + rejected!)
		double       _totAcc;        // total acceptance in this bin

		waveNameArrayType         _waveNames;            // wave names [reflectivity][wave index]
		waveThrArrayType          _waveThresholds;       // mass thresholds of waves
		waveAmpAddedArrayType     _waveAmpAdded;         // amplitude read for waves
		std::vector<fitParameter> _parameters;           // information about function parameters
		waveParamsType            _waveParams;           // map wave name to reflectivity and index in
		                                                 // reflectivity
		ampToParMapType           _prodAmpToFuncParMap;  // maps each production amplitude to the indices
		                                                // of its real and imginary part in the parameter
		                                                // array; negative indices mean that the parameter
		                                                // is not existing due to rank restrictions

		decayAmpsArrayType _decayAmps[2];  // precalculated decay amplitudes [event index][reflectivity][wave index]

		mutable std::vector<double> _parCache;    // parameter cache for derivative calc.
		mutable std::vector<double> _derivCache;  // cache for derivatives

		// normalization integrals
		normMatrixArrayType _normMatrix;          // normalization matrix w/o acceptance [reflectivity 1][wave index 1][reflectivity 2][wave index 2]
		normMatrixArrayType _accMatrix;           // normalization matrix with acceptance [reflectivity 1][wave index 1][reflectivity 2][wave index 2]
		phaseSpaceIntType   _phaseSpaceIntegral;  // phase space integrals

		mutable functionCallInfo _funcCallInfo[NMB_FUNCTIONCALLENUM];  // collects function call statistics

		std::vector<std::string> _eventFileHashOrder;
		std::map<std::string, std::pair<size_t, std::vector<size_t> > > _eventFileProperties;
	};


}

#endif  // PWALIKELIHOOD_H
