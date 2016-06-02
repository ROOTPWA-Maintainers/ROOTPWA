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
//      Data storage class for PWA fit result of one kinematic bin
//
// Environment:
//      Software developed for the COMPASS experiment at CERN
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#ifndef FITRESULT_H
#define FITRESULT_H


#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <complex>

#ifndef __CINT__
#include <boost/shared_ptr.hpp>
#endif

#include "TObject.h"
#include "TComplex.h"
#include "TMatrixT.h"
#include "TString.h"
#include "TPRegexp.h"

#include "complexMatrix.h"
#include "conversionUtils.hpp"
#include "partialWaveFitHelper.h"
#include "reportingUtils.hpp"
#include "reportingUtilsRoot.hpp"


namespace rpwa {


	inline
	std::string
	escapeRegExpSpecialChar(const std::string& s)  ///< escapes all special characters used in regular expressions
	{
		TString escapedS(s);
		const char specialChars[] = {'^', '$', '.', '[', ']', '*', '+', '?'};
		for (unsigned int i = 0; i < sizeof(specialChars) / sizeof(specialChars[0]); ++i) {
			const TString escapedChar = TString("\\").Append(specialChars[i]);
			escapedS.ReplaceAll(specialChars[i], escapedChar);
			const TString escapedTwiceChar = TString("\\").Append(escapedChar);
			while(escapedS.Contains(escapedTwiceChar)) {
				escapedS.ReplaceAll(escapedTwiceChar, escapedChar);
			}
		}
		return std::string(escapedS.Data());
	}


	inline
	std::string
	unescapeRegExpSpecialChar(const std::string& s)  ///< unescapes all special characters used in regular expressions
	{
		TString escapedS(s);
		const char specialChars[] = {'^', '$', '.', '[', ']', '*', '+', '?'};
		for (unsigned int i = 0; i < sizeof(specialChars) / sizeof(specialChars[0]); ++i) {
			const TString escapedChar = TString("\\").Append(specialChars[i]);
			while(escapedS.Contains(escapedChar)) {
				escapedS.ReplaceAll(escapedChar, specialChars[i]);
			}
		}
		return std::string(escapedS.Data());
	}


#ifndef __CINT__
	class fitResult;
	typedef boost::shared_ptr<fitResult> fitResultPtr;
#endif


	/// \brief data storage class for PWA fit result of one kinematic bin
	class fitResult : public TObject {

	public:

		fitResult();
		fitResult(const fitResult& result);
		virtual ~fitResult();

		fitResult* variedProdAmps() const;  ///< create a copy with production amplitudes varied according to covariance matrix

		void reset();
		void fill(const unsigned int                        nmbEvents,               // number of events in bin
		          const unsigned int                        normNmbEvents,           // number of events to normalize to
		          const double                              massBinCenter,           // center value of mass bin
		          const double                              logLikelihood,           // log(likelihood) at maximum
		          const int                                 rank,                    // rank of fit
		          const std::vector<std::complex<double> >& prodAmps,                // production amplitudes
		          const std::vector<std::string>&           prodAmpNames,            // names of production amplitudes used in fit
		          const TMatrixT<double>&                   fitParCovMatrix,         // covariance matrix of fit parameters
		          const std::vector<std::pair<int, int> >&  fitParCovMatrixIndices,  // indices of fit parameters for real and imaginary part in covariance matrix
		          const rpwa::complexMatrix&                normIntegral,            // normalization integral matrix
		          const rpwa::complexMatrix&                acceptedNormIntegral,    // normalization integral matrix with acceptance
		          const std::vector<double>&                phaseSpaceIntegral,      // normalization integral over full phase space without acceptance
		          const bool                                converged,               // indicates whether fit has converged (according to minimizer)
		          const bool                                hasHessian);             // indicates whether Hessian matrix has been calculated successfully
		void fill(const fitResult& result);

		unsigned int        nmbEvents         () const { return _nmbEvents;         }  ///< returns number of events in bin
		unsigned int        normNmbEvents     () const { return _normNmbEvents;     }  ///< returns number of events to normalize to
		double              massBinCenter     () const { return _massBinCenter;     }  ///< returns center value of mass bin
		double              logLikelihood     () const { return _logLikelihood;     }  ///< returns log(likelihood) at maximum
		double              evidence          () const;                                ///< returns the model evidence (OccamFactorMethod)
		std::vector<double> evidenceComponents() const;                                ///< returns a vector { maxLogL, ln(sqrt((2\pi)^m*cov)), -ln(V_A^k), \sum(S_\alpha) }, i.e. evidence = sum(evidenceComponents[i]) for i in {1,2,3,4}
		unsigned int        rank              () const { return _rank;              }  ///< returns rank of fit
		bool                covMatrixValid    () const { return _covMatrixValid;    }
		bool                converged         () const { return _converged;         }  ///< returns whether fit has converged (according to minimizer)
		bool                hasHessian        () const { return _hasHessian;        }  ///< returns whether Hessian matrix has been calculated successfully
		unsigned int        nmbWaves          () const { return _waveNames.size();  }  ///< returns number of waves in fit
		unsigned int        nmbProdAmps       () const { return _prodAmps.size();   }  ///< returns number of production amplitudes

		const std::string& waveName          (const unsigned int waveIndex)    const { return _waveNames[waveIndex];                                }  ///< returns name of wave at index
		std::string        waveNameEsc       (const unsigned int waveIndex)    const { return escapeRegExpSpecialChar(waveName(waveIndex));         }  ///< returns name of wave at index with special regexp characters escaped
		const std::string& prodAmpName       (const unsigned int prodAmpIndex) const { return _prodAmpNames[prodAmpIndex];                          }  ///< returns name of production amplitude at index
		std::string        prodAmpNameEsc    (const unsigned int prodAmpIndex) const { return escapeRegExpSpecialChar(_prodAmpNames[prodAmpIndex]); }  ///< returns name of production amplitude at index with special regexp characters escaped
		inline std::string waveNameForProdAmp(const unsigned int prodAmpIndex) const;
		inline int         rankOfProdAmp     (const unsigned int prodAmpIndex) const;

		int waveIndex   (const std::string& waveName   ) const;  ///< returns wave index corresponding to wave name
		int prodAmpIndex(const std::string& prodAmpName) const;  ///< returns production amplitude index corresponding to production amplitude name

		double fitParameter   (const std::string& parName  ) const;  ///< returns value of fit parameter with name
		double fitParameterErr(const std::string& parName  ) const;  ///< returns error of fit parameter with name
		/// returns covariance of fit parameters at index A and B
		double fitParameterCov(const unsigned int parIndexA,
		                       const unsigned int parIndexB) const { return _fitParCovMatrix[parIndexA][parIndexB]; }

		/// returns production amplitude value at index
		std::complex<double>    prodAmp   (const unsigned int prodAmpIndex) const { return std::complex<double>(_prodAmps[prodAmpIndex].Re(), _prodAmps[prodAmpIndex].Im()); }
		/// returns covariance matrix of production amplitude value at index
		inline TMatrixT<double> prodAmpCov(const unsigned int prodAmpIndex) const;
		///< returns covariance matrix for a set of production amplitudes given by index list
		TMatrixT<double>        prodAmpCov(const std::vector<unsigned int>& prodAmpIndices) const;

		/// returns normalization integral for pair of waves at index A and B
		std::complex<double> normIntegral        (const unsigned int waveIndexA,
		                                          const unsigned int waveIndexB) const { return _normIntegral(waveIndexA, waveIndexB); }
		/// returns accepted normalization integral for pair of waves at index A and B
		std::complex<double> acceptedNormIntegral(const unsigned int waveIndexA,
		                                          const unsigned int waveIndexB) const { return _acceptedNormIntegral(waveIndexA, waveIndexB); }

		/// returns the sqrt(!) of the phase space integral for given wave
		double phaseSpaceIntegral(const unsigned int waveIndex) const { return _phaseSpaceIntegral[waveIndex];          }
		double phaseSpaceIntegral(const std::string& waveName ) const { return phaseSpaceIntegral(waveIndex(waveName)); }


		// accessors to information from spin-density matrix
		// * access by wave indices

		/// returns spin density matrix element for pair of waves at index A and B
		std::complex<double> spinDensityMatrixElem   (const unsigned int waveIndexA,
		                                              const unsigned int waveIndexB) const;
		/// returns covariance matrix of spin density matrix element for pair of waves at index A and B
		TMatrixT<double>     spinDensityMatrixElemCov(const unsigned int waveIndexA,
		                                              const unsigned int waveIndexB) const;

		/// returns phase difference between pair of waves at index A and B
		double phase   (const unsigned int waveIndexA,
		                const unsigned int waveIndexB) const;
		/// returns error of phase difference between pair of waves at index A and B
		double phaseErr(const unsigned int waveIndexA,
		                const unsigned int waveIndexB) const;

		/// returns coherence of pair of waves at index A and B
		double coherence   (const unsigned int waveIndexA,
		                    const unsigned int waveIndexB) const;
		/// returns error of coherence of pair of waves at index A and B
		double coherenceErr(const unsigned int waveIndexA,
		                    const unsigned int waveIndexB) const;

		/// returns overlap of pair of waves at index A and B
		double overlap   (const unsigned int waveIndexA,
		                  const unsigned int waveIndexB) const;
		/// returns error of overlap of pair of waves at index A and B
		double overlapErr(const unsigned int waveIndexA,
		                  const unsigned int waveIndexB) const;

		// accessors to information from spin-density matrix
		// * access by wave names

		/// returns spin density matrix element for pair of waves A and B
		std::complex<double> spinDensityMatrixElem   (const std::string& waveNameA,
		                                              const std::string& waveNameB) const { return spinDensityMatrixElem   (waveIndex(waveNameA), waveIndex(waveNameB)); }
		/// returns covariance matrix of spin density matrix element for pair of waves A and B
		TMatrixT<double>     spinDensityMatrixElemCov(const std::string& waveNameA,
		                                              const std::string& waveNameB) const { return spinDensityMatrixElemCov(waveIndex(waveNameA), waveIndex(waveNameB)); }

		/// returns phase difference between pair of waves A and B
		double phase   (const std::string& waveNameA,
		                const std::string& waveNameB) const { return phase   (waveIndex(waveNameA), waveIndex(waveNameB)); }
		/// returns error of phase difference between pair of waves A and B
		double phaseErr(const std::string& waveNameA,
		                const std::string& waveNameB) const { return phaseErr(waveIndex(waveNameA), waveIndex(waveNameB)); }

		/// returns coherence of pair of waves A and B
		double coherence   (const std::string& waveNameA,
		                    const std::string& waveNameB) const { return coherence   (waveIndex(waveNameA), waveIndex(waveNameB)); }
		/// returns error of coherence of pair of waves A and B
		double coherenceErr(const std::string& waveNameA,
		                    const std::string& waveNameB) const { return coherenceErr(waveIndex(waveNameA), waveIndex(waveNameB)); }

		/// returns overlap of pair of waves A and B
		double overlap   (const std::string& waveNameA,
		                  const std::string& waveNameB) const { return overlap   (waveIndex(waveNameA), waveIndex(waveNameB)); }
		/// returns error of overlap of pair of waves A and B
		double overlapErr(const std::string& waveNameA,
		                  const std::string& waveNameB) const { return overlapErr(waveIndex(waveNameA), waveIndex(waveNameB)); }


		// accessors to intensities
		// * access by wave indices

		/// returns intensity of single wave at index
		double intensity   (const unsigned int waveIndex) const { return spinDensityMatrixElem(waveIndex, waveIndex).real();         }
		/// returns error of intensity of single wave at index
		double intensityErr(const unsigned int waveIndex) const { return sqrt(spinDensityMatrixElemCov(waveIndex, waveIndex)[0][0]); }

		// accessors to intensities
		// * access by wave name patterns

		/// returns intensity of sum of waves matching name pattern
		double intensity   (const std::string& waveNamePattern) const;
		/// returns error of intensity of sum of waves matching name pattern
		double intensityErr(const std::string& waveNamePattern) const;

		// accessors to intensities
		// * total intensity

		/// returns total intensity
		double intensity   () const { return intensity   (".*"); }
		/// returns error of total intensity
		double intensityErr() const { return intensityErr(".*"); }


		// low level interface to make copying easier
		const std::vector<TComplex>&                 prodAmps                  () const { return _prodAmps;               }
		const std::vector<std::string>&              prodAmpNames              () const { return _prodAmpNames;           }
		const std::vector<std::string>&              waveNames                 () const { return _waveNames;              }
		const TMatrixT<Double_t>&                    fitParCovMatrix           () const { return _fitParCovMatrix;        }
		const std::vector<std::pair<Int_t, Int_t> >& fitParCovIndices          () const { return _fitParCovMatrixIndices; }
		const rpwa::complexMatrix&                   normIntegralMatrix        () const { return _normIntegral;           }
		const rpwa::complexMatrix&                   acceptedNormIntegralMatrix() const { return _acceptedNormIntegral;   }
		const std::vector<double>&                   phaseSpaceIntegralVector  () const { return _phaseSpaceIntegral;     }
		const std::map<Int_t, Int_t>&                normIntIndexMap           () const { return _normIntIndexMap;        }


		inline std::ostream& printProdAmpNames(std::ostream& out = std::cout) const;  ///< prints all production amplitude names
		inline std::ostream& printWaveNames   (std::ostream& out = std::cout) const;  ///< prints all wave names
		inline std::ostream& printProdAmps    (std::ostream& out = std::cout) const;  ///< prints all production amplitudes and their covariance matrix
		inline std::ostream& printWaves       (std::ostream& out = std::cout) const;  ///< prints all wave intensities and their errors

		virtual inline std::ostream& print(std::ostream& out = std::cout) const;
		friend std::ostream& operator << (std::ostream&    out,
		                                  const fitResult& result) { return result.print(out); }

	private:

		// helper functions

		std::complex<double> normIntegralForProdAmp(const unsigned int prodAmpIndexA,
		                                            const unsigned int prodAmpIndexB) const;

		/// returns indices of waves that match the regular expression
		inline std::vector<unsigned int> waveIndicesMatchingPattern   (const std::string& waveNamePattern   ) const;
		/// returns indices of production amplitudes that match the regular expression
		inline std::vector<unsigned int> prodAmpIndicesMatchingPattern(const std::string& prodAmpNamePattern) const;
		/// returns indices of production amplitudes that belong to wave at index
		std::vector<unsigned int>        prodAmpIndicesForWave        (const unsigned int waveIndex         ) const
		{ return prodAmpIndicesMatchingPattern(waveNameEsc(waveIndex)); }
		/// returns pair of indices of production amplitudes that belong to the waves and have the same rank
		inline std::vector<std::pair<unsigned int, unsigned int> > prodAmpIndexPairsForWaves(const unsigned int waveIndexA,
		                                                                                     const unsigned int waveIndexB) const;

		///< returns covariance matrix for a set of production amplitudes given by a list of index pairs
		inline TMatrixT<double> prodAmpCov(const std::vector<std::pair<unsigned int, unsigned int> >& prodAmpIndexPairs) const;
		///< returns covariance matrix for a set of production amplitudes given by index lists A and B
		inline TMatrixT<double> prodAmpCov(const std::vector<unsigned int>& prodAmpIndicesA,
		                                   const std::vector<unsigned int>& prodAmpIndicesB) const;

		inline double realValVariance(const unsigned int      waveIndexA,
		                              const unsigned int      waveIndexB,
		                              const TMatrixT<double>& jacobian) const;


		// stored data
		UInt_t                                _nmbEvents;                 ///< number of events in bin
		UInt_t                                _normNmbEvents;             ///< number of events to normalize to
		Double_t                              _massBinCenter;             ///< center value of mass bin
		Double_t                              _logLikelihood;             ///< log(likelihood) at maximum
		Int_t                                 _rank;                      ///< rank of fit
		std::vector<TComplex>                 _prodAmps;                  ///< production amplitudes
		std::vector<std::string>              _prodAmpNames;              ///< names of production amplitudes used in fit
		std::vector<std::string>              _waveNames;                 ///< names of waves used in fit
		Bool_t                                _covMatrixValid;            ///< indicates whether bin has a valid covariance matrix
		TMatrixT<Double_t>                    _fitParCovMatrix;           ///< covariance matrix of fit parameters
		std::vector<std::pair<Int_t, Int_t> > _fitParCovMatrixIndices;    ///< indices of fit parameters for real and imaginary part in covariance matrix matrix
		rpwa::complexMatrix                   _normIntegral;         //|| ///< normalization integral over full phase space without acceptance
		rpwa::complexMatrix                   _acceptedNormIntegral; //|| ///< normalization integral over accepted phase space
		std::map<Int_t, Int_t>                _normIntIndexMap;           ///< maps production amplitude indices to indices in normalization integral
		std::vector<double>                   _phaseSpaceIntegral;        ///< diagonals of phase space integrals (without acceptance)
		bool                                  _converged;                 ///< indicates whether fit has converged (according to minimizer)
		bool                                  _hasHessian;                ///< indicates whether Hessian matrix has been calculated successfully
		// add more info about fit: quality of fit information, ndf, list of fixed parameters, ...

	public:

		ClassDef(fitResult, 6)

	};  // class fitResult


	/// returns covariance matrix for a single production amplitude
	inline
	TMatrixT<double>
	fitResult::prodAmpCov(const unsigned int prodAmpIndex) const {
		TMatrixT<double> cov(2, 2);
		if(not covMatrixValid()) {
			printWarn << "no valid covariance matrix to return, return 0-matrix." << std::endl;
			cov = 0.;
			return cov;
		}
		// get parameter indices
		const int i = _fitParCovMatrixIndices[prodAmpIndex].first;
		const int j = _fitParCovMatrixIndices[prodAmpIndex].second;
		cov[0][0] = _fitParCovMatrix[i][i];
		if (j >= 0) {
			cov[0][1] = _fitParCovMatrix[i][j];
			cov[1][0] = _fitParCovMatrix[j][i];
			cov[1][1] = _fitParCovMatrix[j][j];
		}
		return cov;
	}


	/// \brief constructs covariance matrix for production amplitudes specified by index pair list
	///
	/// layout:
	///         cov(first,  first)        cov(first,  second)
	///         cov(second, first)        cov(second, second)
	inline
	TMatrixT<double>
	fitResult::prodAmpCov(const std::vector<std::pair<unsigned int, unsigned int> >& prodAmpIndexPairs) const
	{
		std::vector<unsigned int> prodAmpIndices;
		// copy first production amplitude indices
		for (unsigned int i = 0; i < prodAmpIndexPairs.size(); ++i)
			prodAmpIndices.push_back(prodAmpIndexPairs[i].first);
		// copy second production amplitude indices
		for (unsigned int i = 0; i < prodAmpIndexPairs.size(); ++i)
			prodAmpIndices.push_back(prodAmpIndexPairs[i].second);
		return prodAmpCov(prodAmpIndices);
	}


	/// \brief constructs covariance matrix of production amplitudes specified by two index lists
	///
	/// layout:
	///         cov(A, A)        cov(A, B)
	///         cov(B, A)        cov(B, B)
	inline
	TMatrixT<double>
	fitResult::prodAmpCov(const std::vector<unsigned int>& prodAmpIndicesA,
	                      const std::vector<unsigned int>& prodAmpIndicesB) const
	{
		std::vector<unsigned int> prodAmpIndices;
		// copy wave A production amplitude indices
		prodAmpIndices.assign(prodAmpIndicesA.begin(), prodAmpIndicesA.end());
		// copy wave B production amplitude indices
		prodAmpIndices.insert(prodAmpIndices.end(), prodAmpIndicesB.begin(), prodAmpIndicesB.end());
		return prodAmpCov(prodAmpIndices);
	}


	// prints all production amplitude names
	inline
	std::ostream&
	fitResult::printProdAmpNames(std::ostream& out) const
	{
		out << "    Production amplitude names:" << std::endl;
		for (unsigned int i = 0; i < nmbProdAmps(); ++i)
			out << "        " << std::setw(3) << i << " " << _prodAmpNames[i] << std::endl;
		return out;
	}


	// prints all wave names
	inline
	std::ostream&
	fitResult::printWaveNames(std::ostream& out) const
	{
		out << "    Wave names:" << std::endl;
		for (unsigned int i = 0; i < nmbWaves(); ++i)
			out << "        " << std::setw(3) << i << " " << waveName(i) << std::endl;
		return out;
	}


	// prints all production amplitudes and their covariance matrix
	inline
	std::ostream&
	fitResult::printProdAmps(std::ostream& out) const
	{
		out << "Production amplitudes:" << std::endl;
		for (unsigned int i = 0; i < nmbProdAmps(); ++i) {
			out << "    " << std::setw(3) << i << " " << _prodAmpNames[i] << " = "  << prodAmp(i)
			    << ", cov = " << prodAmpCov(i) << std::endl;
		}
		return out;
	}


	// prints all wave intensities and their errors
	inline
	std::ostream&
	fitResult::printWaves(std::ostream& out) const
	{
		out << "Waves:" << std::endl;
		for (unsigned int i = 0; i < nmbWaves(); ++i) {
			out << "    " << std::setw(3) << i << " " << waveName(i) << " = "  << intensity(i)
			    << " +- " << intensityErr(i) << std::endl;
		}
		return out;
	}


	/// dumps all raw data stored in object
	inline
	std::ostream&
	fitResult::print(std::ostream& out) const
	{
		out << "fitResult dump:" << std::endl
		    << "    number of events ..................... " << nmbEvents()                   << std::endl
		    << "    number of events to normalize to ..... " << normNmbEvents()               << std::endl
		    << "    center value of mass bin ............. " << massBinCenter()               << std::endl
		    << "    log(likelihood) at maximum ........... " << logLikelihood()               << std::endl
		    << "    rank of fit .......................... " << rank()                        << std::endl
		    << "    bin has a valid covariance matrix .... " << rpwa::yesNo(covMatrixValid()) << std::endl
		    << "    fit has converged .................... " << rpwa::yesNo(converged())      << std::endl
		    << "    Hessian matrix has been calculated ... " << rpwa::yesNo(hasHessian())     << std::endl;
		printProdAmps(out);
		printWaveNames(out);
		out << "    covariance matrix:" << std::endl << _fitParCovMatrix << std::endl;
		out << "    covariance matrix indices:" << std::endl;
		for (unsigned int i = 0; i < _fitParCovMatrixIndices.size(); ++i)
			out << "        index " << std::setw(3) << i << " = (" << std::setw(3) << _fitParCovMatrixIndices[i].first
			    << ", " << std::setw(3) << _fitParCovMatrixIndices[i].second << ")" << std::endl;
		out << "    normalization integral (w/o acceptance):" << std::endl << _normIntegral << std::endl;
		out << "    normalization integral (with acceptance):" << std::endl << _acceptedNormIntegral << std::endl;
		out << "    map of production amplitude indices to indices in normalization integral:" << std::endl;
		for (std::map<Int_t, Int_t>::const_iterator i = _normIntIndexMap.begin(); i != _normIntIndexMap.end(); ++i)
			out << "        prod. amp [" << std::setw(3) << i->first << "] "
			    << "-> norm. int. [" << std::setw(3) << i->second << "]" << std::endl;
		return out;
	}


	/// extracts rank of production amplitude from its name
	inline
	int
	fitResult::rankOfProdAmp(const unsigned int prodAmpIndex) const
	{
		const std::string& ampName = _prodAmpNames[prodAmpIndex];
		if (ampName == "V_flat")
			return -1;
		else
			// the rank is encoded between the V and the next underscore,
			// so starting from the second character to the first character
			// which cannot be converted
			return atoi(&ampName.c_str()[1]);
	}


	/// generates list of wave indices that match name pattern
	inline
	std::vector<unsigned int>
	fitResult::waveIndicesMatchingPattern(const std::string& waveNamePattern) const
	{
		TPRegexp Regexp(waveNamePattern);
		std::vector<unsigned int> waveIndices;
		for (unsigned int waveIndex = 0; waveIndex < nmbWaves(); ++waveIndex) {
			if (TString(waveName(waveIndex)).Contains(Regexp))
				waveIndices.push_back(waveIndex);
		}
		return waveIndices;
	}


	/// generates list of production amplitude indices that match name pattern
	inline
	std::vector<unsigned int>
	fitResult::prodAmpIndicesMatchingPattern(const std::string& ampNamePattern) const
	{
		TPRegexp Regexp(ampNamePattern);
		std::vector<unsigned int> prodAmpIndices;
		for (unsigned int prodAmpIndex = 0; prodAmpIndex < nmbProdAmps(); ++prodAmpIndex)
			if (TString(prodAmpName(prodAmpIndex)).Contains(Regexp))
				prodAmpIndices.push_back(prodAmpIndex);
		return prodAmpIndices;
	}


	/// generates list of pairs of production amplitude indices of same rank for a pair of waves
	inline
	std::vector<std::pair<unsigned int, unsigned int> >
	fitResult::prodAmpIndexPairsForWaves(const unsigned int waveIndexA,
	                                     const unsigned int waveIndexB) const
	{
		std::vector<std::pair<unsigned int, unsigned int> > prodAmpIndexPairs;
		// return empty vector if the two waves have different reflectivities
		if (rpwa::partialWaveFitHelper::getReflectivity(waveName(waveIndexA))
		    != rpwa::partialWaveFitHelper::getReflectivity(waveName(waveIndexB))) {
			return prodAmpIndexPairs;
		}

		const std::vector<unsigned int> prodAmpIndicesA = prodAmpIndicesForWave(waveIndexA);
		const std::vector<unsigned int> prodAmpIndicesB = prodAmpIndicesForWave(waveIndexB);
		for (unsigned int countAmpA = 0; countAmpA < prodAmpIndicesA.size(); ++countAmpA) {
			const unsigned int ampIndexA = prodAmpIndicesA[countAmpA];
			const int          ampRankA  = rankOfProdAmp(ampIndexA);
			// find production amplitude of wave B with same rank
			for (unsigned int countAmpB = 0; countAmpB < prodAmpIndicesB.size(); ++countAmpB) {
				const unsigned int ampIndexB = prodAmpIndicesB[countAmpB];
				const int          ampRankB  = rankOfProdAmp(ampIndexB);
				if (ampRankA == ampRankB) {
					prodAmpIndexPairs.push_back(std::make_pair(ampIndexA, ampIndexB));
					break;
				}
			}
		}
		return prodAmpIndexPairs;
	}


	/// calculates variance of a real-valued function of a spin density matrix element for wave A and wave B
	inline
	double
	fitResult::realValVariance(const unsigned int      waveIndexA,
	                           const unsigned int      waveIndexB,
	                           const TMatrixT<double>& jacobian) const  // Jacobian of real valued function (d f/ d Re[rho]   d f / d Im[rho])
	{
		if (!covMatrixValid())
			return 0;

		const TMatrixT<double> spinDensCov = spinDensityMatrixElemCov(waveIndexA, waveIndexB);  // 2 x 2 matrix
		const TMatrixT<double> jacobianT(TMatrixT<double>::kTransposed, jacobian);              // 2 x 1 matrix
		const TMatrixT<double> spinDensCovJT = spinDensCov * jacobianT;                         // 2 x 1 matrix
		const TMatrixT<double> cov           = jacobian * spinDensCovJT;                        // 1 x 1 matrix
		return cov[0][0];
	}


	/// \brief returns wave name for given production amplitude
	///
	/// naming scheme for production amplitudes:
	/// "V<rank>_<wave name>" with exception "V_flat" for flat wave
	/// <rank> is assumed to be a single-digit number
	inline
	std::string
	fitResult::waveNameForProdAmp(const unsigned int prodAmpIndex) const
	{
		const std::string& prodAmpName = _prodAmpNames[prodAmpIndex];
		if (prodAmpName.length() == 0 or prodAmpName[0] != 'V' or prodAmpName.find('_') == std::string::npos) {
			printErr << "production amplitude name '" << prodAmpName << "' does not follow the naming convention. "
			         << "cannot deduce corresponding wave name." << std::endl;
			return "";
		}
		return prodAmpName.substr(prodAmpName.find('_')+1);
	}


}  // namespace rpwa


#endif  // FITRESULT_H
