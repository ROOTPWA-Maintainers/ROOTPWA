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


//
//  !!! deprecated class --- use fitResult instead !!!
//


#ifndef TFITRESULT_HH
#define TFITRESULT_HH


#ifdef USE_TFITRESULT


#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <complex>

#include "TObject.h"
#include "TComplex.h"
#include "TMatrixT.h"
#include "TString.h"

#include "reportingUtils.hpp"
#include "reportingUtilsRoot.hpp"
#include "TCMatrix.h"
#include "TFitBin.h"


/// \brief data storage class for PWA fit result of one kinematic bin
class TFitResult : public TObject {

public:
    
	TFitResult();
	TFitResult(const TFitResult& result);
	TFitResult(const TFitBin&    fitBin);
	virtual ~TFitResult();
  
	void reset();
	void fill(const unsigned int                        nmbEvents,               // number of events in bin                      
	          const unsigned int                        normNmbEvents,	         // number of events to normalize to		 
	          const double                              massBinCenter,	         // center value of mass bin			 
	          const double                              logLikelihood,	         // log(likelihood) at maximum		 
	          const int                                 rank,		         // rank of fit				 
	          const std::vector<std::complex<double> >& prodAmps,	         // production amplitudes			 
	          const std::vector<std::string>&           prodAmpNames,	         // names of production amplitudes used in fit
	          const TMatrixT<double>&                   fitParCovMatrix,         // covariance matrix of fit parameters
	          const std::vector<std::pair<int, int> >&  fitParCovMatrixIndices,  // indices of fit parameters for real and imaginary part in covariance matrix matrix
	          const TCMatrix&                           normIntegral);           // normalization integral over full phase space without acceptance

	double       massBinCenter() const { return _massBinCenter;    }  ///< returns center value of mass bin
	double       logLikelihood() const { return _logLikelihood;    }  ///< returns log(likelihood) at maximum
	double       evidence     () const ;                              ///< return the model evidence (OccamFactorMethod)
	unsigned int rank         () const { return _rank;             }  ///< returns rank of fit
	unsigned int nmbEvents    () const { return _nmbEvents;        }  ///< returns number of events in bin
	unsigned int normNmbEvents() const { return _normNmbEvents;    }  ///< returns number of events to normalize to
	unsigned int nmbWaves     () const { return _waveNames.size(); }  ///< returns number of waves in fit
	unsigned int nmbProdAmps  () const { return _prodAmps.size();  }  ///< returns number of production amplitudes
  
	TString waveName    (const unsigned int waveIndex)    const { return _waveNames[waveIndex];       }  ///< returns name of wave at index
	TString prodAmpName (const unsigned int prodAmpIndex) const { return _prodAmpNames[prodAmpIndex]; }  ///< returns name of production amplitude at index

	double      fitParameter   (const std::string& parName)   const;  ///< returns value of fit parameter with name
	/// returns covariance of fit parameters at index A and B
	double      fitParameterCov(const unsigned int parIndexA,
	                            const unsigned int parIndexB) const { return _fitParCovMatrix[parIndexA][parIndexB]; }
	inline void fitParameters  (double*            parArray)  const;  ///< copies fit parameters into array

	/// returns production amplitude value at index
	std::complex<double>    prodAmp   (const unsigned int prodAmpIndex) const { return std::complex<double>(_prodAmps[prodAmpIndex].Re(), _prodAmps[prodAmpIndex].Im()); }
	inline TMatrixT<double> prodAmpCov(const unsigned int prodAmpIndex) const;   ///< returns covariance matrix of production amplitude value at index
	///< returns covariance matrix for a set of production amplitudes given by index list
	TMatrixT<double>        prodAmpCov(const std::vector<unsigned int>&                           prodAmpIndices)    const;
	///< returns covariance matrix for a set of production amplitudes given by a list of index pairs
	inline TMatrixT<double> prodAmpCov(const std::vector<std::pair<unsigned int, unsigned int> >& prodAmpIndexPairs) const;
	///< returns covariance matrix for a set of production amplitudes given by index lists A and B
	inline TMatrixT<double> prodAmpCov(const std::vector<unsigned int>& prodAmpIndicesA,
	                                   const std::vector<unsigned int>& prodAmpIndicesB) const;
	bool                    covMatrixValid() const { return _hasErrors; }

	/// returns normalization integral for pair of waves at index A and B
	inline std::complex<double> normIntegral(const unsigned int waveIndexA,
	                                         const unsigned int waveIndexB) const;

	/// returns spin density matrix element for pair of waves at index A and B
	std::complex<double> spinDensityMatrixElem   (const unsigned int waveIndexA,
	                                              const unsigned int waveIndexB) const;
	/// returns covariance matrix of spin density matrix element for pair of waves at index A and B
	TMatrixT<double>     spinDensityMatrixElemCov(const unsigned int waveIndexA,
	                                              const unsigned int waveIndexB) const;

	/// returns intensity of single wave at index
	double intensity   (const unsigned int waveIndex)       const { return spinDensityMatrixElem(waveIndex, waveIndex).real();         }
	/// returns error of intensity of single wave at index
	double intensityErr(const unsigned int waveIndex)       const { return sqrt(spinDensityMatrixElemCov(waveIndex, waveIndex)[0][0]); }
	double intensity   (const char*        waveNamePattern) const;                              ///< returns intensity of sum of waves matching name pattern
	double intensityErr(const char*        waveNamePattern) const;                              ///< returns error of intensity of sum of waves matching name pattern
	double intensity   ()                                   const { return intensity("");    }  ///< returns total intensity
	double intensityErr()                                   const { return intensityErr(""); }  ///< returns error of total intensity

	double phase       (const unsigned int waveIndexA,
	                    const unsigned int waveIndexB) const;  ///< returns phase difference between two waves at index A and B
	double phaseErr    (const unsigned int waveIndexA,
	                    const unsigned int waveIndexB) const;  ///< returns error of phase difference between two waves at index A and B
	double coherence   (const unsigned int waveIndexA,
	                    const unsigned int waveIndexB) const;  ///< returns coherence of two waves at index A and B
	double coherenceErr(const unsigned int waveIndexA,
	                    const unsigned int waveIndexB) const;  ///< returns error of coherence of two waves at index A and B
	double overlap     (const unsigned int waveIndexA,
	                    const unsigned int waveIndexB) const;  ///< returns overlap of two waves at index A and B
	double overlapErr  (const unsigned int waveIndexA,
	                    const unsigned int waveIndexB) const;  ///< returns error of overlap of two waves at index A and B

	// low level interface to make copying easier
	const std::vector<TComplex>&                 prodAmps          () const { return _prodAmps;               }
	const std::vector<std::string>&              prodAmpNames      () const { return _prodAmpNames;           }
	const std::vector<std::string>&              waveNames         () const { return _waveNames;              }
	const TMatrixT<Double_t>&                    fitParCovMatrix   () const { return _fitParCovMatrix;        }
	const std::vector<std::pair<Int_t, Int_t> >& fitParCovIndices  () const { return _fitParCovMatrixIndices; }
	const TCMatrix&                              normIntegralMatrix() const { return _normIntegral;           }
	const std::map<Int_t, Int_t>&                normIntIndexMap   () const { return _normIntIndexMap;        }

	inline std::ostream& printProdAmpNames(std::ostream& out = std::cout) const;  ///< prints all production amplitude names
	inline std::ostream& printWaveNames   (std::ostream& out = std::cout) const;  ///< prints all wave names
	inline std::ostream& printProdAmps    (std::ostream& out = std::cout) const;  ///< prints all production amplitudes and their covariance matrix
	inline std::ostream& printWaves       (std::ostream& out = std::cout) const;  ///< prints all wave intensities and their errors

	virtual inline std::ostream& print(std::ostream& out = std::cout) const;
	friend std::ostream& operator << (std::ostream&     out,
	                                  const TFitResult& fitResult) { return fitResult.print(out); }

private:

	UInt_t                                _nmbEvents;               ///< number of events in bin			 
	UInt_t                                _normNmbEvents;           ///< number of events to normalize to		 
	Double_t                              _massBinCenter;           ///< center value of mass bin			 
	Double_t                              _logLikelihood;           ///< log(likelihood) at maximum		 
	Int_t                                 _rank;                    ///< rank of fit				 
	std::vector<TComplex>                 _prodAmps;                ///< production amplitudes                   
	std::vector<std::string>              _prodAmpNames;            ///< names of production amplitudes used in fit
	std::vector<std::string>              _waveNames;               ///< names of waves used in fit
	Bool_t                                _hasErrors;               ///< indicates whether bin has a valid covariance matrix
	TMatrixT<Double_t>                    _fitParCovMatrix;         ///< covariance matrix of fit parameters
	std::vector<std::pair<Int_t, Int_t> > _fitParCovMatrixIndices;  ///< indices of fit parameters for real and imaginary part in covariance matrix matrix
	TCMatrix                              _normIntegral;            ///< normalization integral over full phase space without acceptance
	std::map<Int_t, Int_t>                _normIntIndexMap;         ///< maps production amplitude indices to indices in normalization integral

	// add more info about fit: quality of fit information, ndf, list of fixed parameters, ...
  
	// helper functions
	inline static TMatrixT<double> matrixRepr(const std::complex<double>& c);

	inline int rankOfProdAmp(const unsigned int prodAmpIndex) const;

	std::complex<double> normIntegralForProdAmp(const unsigned int prodAmpIndexA,
	                                            const unsigned int prodAmpIndexB) const;

	inline std::vector<unsigned int> waveIndicesMatchingPattern   (const std::string& waveNamePattern) const;
	inline std::vector<unsigned int> prodAmpIndicesMatchingPattern(const std::string& ampNamePattern ) const;
	std::vector<unsigned int>        prodAmpIndicesForWave        (const unsigned int waveIndex)       const
	{ return prodAmpIndicesMatchingPattern(_waveNames[waveIndex]); }  ///< returns indices of production amplitudes that belong to wave at index
	inline std::vector<std::pair<unsigned int, unsigned int> > prodAmpIndexPairsForWaves(const unsigned int waveIndexA,
		               const unsigned int waveIndexB) const;

	inline double realValVariance(const unsigned int      waveIndexA,
	                              const unsigned int      waveIndexB,
	                              const TMatrixT<double>& jacobian) const;


	TString wavetitle(int i) const
	{
		const TString ampName = _prodAmpNames[i];
		if (ampName.Contains("V_"))
			return ampName(2, 1000);
		else
			return ampName(3, 1000);
	}
	void buildWaveMap();

public:

	ClassDef(TFitResult,1)
  
};
  

/// copies fit parameters into array; assumes that array has sufficient size
inline
void
TFitResult::fitParameters(double* parArray) const
{
	unsigned int countPar = 0;
	for (unsigned int i = 0; i < nmbProdAmps(); ++i) {
		parArray[countPar++] = prodAmp(i).real();
		// check if amplitude has imaginary part
		if (_fitParCovMatrixIndices[i].second >= 0)
			parArray[countPar++] = prodAmp(i).imag();
	}
}


/// returns covariance matrix for a single production amplitude
inline
TMatrixT<double>
TFitResult::prodAmpCov(const unsigned int prodAmpIndex) const {
	// get parameter indices
	const int i = _fitParCovMatrixIndices[prodAmpIndex].first;
	const int j = _fitParCovMatrixIndices[prodAmpIndex].second;
	TMatrixT<double> cov(2, 2);
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
TFitResult::prodAmpCov(const std::vector<std::pair<unsigned int, unsigned int> >& prodAmpIndexPairs) const
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
TFitResult::prodAmpCov(const std::vector<unsigned int>& prodAmpIndicesA,
                       const std::vector<unsigned int>& prodAmpIndicesB) const
{
	std::vector<unsigned int> prodAmpIndices;
	// copy wave A production amplitude indices
	prodAmpIndices.assign(prodAmpIndicesA.begin(), prodAmpIndicesA.end());
	// copy wave B production amplitude indices
	prodAmpIndices.insert(prodAmpIndices.end(), prodAmpIndicesB.begin(), prodAmpIndicesB.end());
	return prodAmpCov(prodAmpIndices);
}


// returns normalization integral for pair of waves at index A and B
inline
std::complex<double>
TFitResult::normIntegral(const unsigned int waveIndexA,
                         const unsigned int waveIndexB) const
{
	const TComplex norm = _normIntegral(waveIndexA, waveIndexB);
	return std::complex<double>(norm.Re(), norm.Im());
}


// prints all production amplitude names
inline
std::ostream&
TFitResult::printProdAmpNames(std::ostream& out) const
{
	out << "    Production amplitude names:" << std::endl;
	for (unsigned int i = 0; i < nmbProdAmps(); ++i)
		out << "        " << std::setw(3) << i << " " << _prodAmpNames[i] << std::endl;
	return out;
}


// prints all wave names
inline
std::ostream&
TFitResult::printWaveNames(std::ostream& out) const
{
	out << "    Wave names:" << std::endl;
	for (unsigned int i = 0; i < nmbWaves(); ++i)
		out << "        " << std::setw(3) << i << " " << _waveNames[i] << std::endl;
	return out;
}


// prints all production amplitudes and their covariance matrix
inline
std::ostream&
TFitResult::printProdAmps(std::ostream& out) const
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
TFitResult::printWaves(std::ostream& out) const
{
	out << "Waves:" << std::endl;
	for (unsigned int i = 0; i < nmbWaves(); ++i) {
		out << "    " << std::setw(3) << i << " " << _waveNames[i] << " = "  << intensity(i)
		    << " +- " << intensityErr(i) << std::endl;
	}
	return out;
}


/// dumps all raw data stored in object
inline
std::ostream&
TFitResult::print(std::ostream& out) const
{
	out << "TFitResult dump:" << std::endl
	    << "    number of events .................... " << _nmbEvents     << std::endl
	    << "    number of events to normalize to .... " << _normNmbEvents << std::endl
	    << "    center value of mass bin ............ " << _massBinCenter << std::endl
	    << "    log(likelihood) at maximum .......... " << _logLikelihood << std::endl
	    << "    rank of fit ......................... " << _rank          << std::endl
	    << "    bin has a valid covariance matrix ... " << _hasErrors     << std::endl;
	printProdAmps(out);
	printWaveNames(out);
	out << "    covariance matrix:" << std::endl << _fitParCovMatrix << std::endl;
	out << "    covariance matrix indices:" << std::endl;
	for (unsigned int i = 0; i < _fitParCovMatrixIndices.size(); ++i)
		out << "        index " << std::setw(3) << i << " = (" << std::setw(3) << _fitParCovMatrixIndices[i].first
		    << ", " << std::setw(3) << _fitParCovMatrixIndices[i].second << ")" << std::endl;
	out << "    normalization integral (w/o acceptance):" << std::endl << _normIntegral << std::endl;
	out << "    map of production amplitude indices to indices in normalization integral:" << std::endl;
	for (std::map<Int_t, Int_t>::const_iterator i = _normIntIndexMap.begin(); i != _normIntIndexMap.end(); ++i)
		out << "        prod. amp [" << std::setw(3) << i->first << "] "
		    << "-> norm. int. [" << std::setw(3) << i->second << "]" << std::endl;
	return out;
}


/// \brief returns matrix representation of complex number
///
///   c.re  -c.im
///   c.im   c.re
// !!! possible optimization: return TMatrixTLazy
inline
TMatrixT<double>
TFitResult::matrixRepr(const std::complex<double>& c)
{
	TMatrixT<double> m(2, 2);
	m[0][0] = m[1][1] = c.real();
	m[0][1] = -c.imag();
	m[1][0] =  c.imag();
	return m;
}


/// extracts rank of production amplitude from its name
inline
int
TFitResult::rankOfProdAmp(const unsigned int prodAmpIndex) const
{
	const TString ampName = _prodAmpNames[prodAmpIndex];
	if (ampName.Contains("flat"))
		return -1;
	else
		// the rank is encoded in second character of parameter name
		return TString(ampName(1, 1)).Atoi();
}


/// generates list of wave indices that match name pattern
inline
std::vector<unsigned int>
TFitResult::waveIndicesMatchingPattern(const std::string& waveNamePattern) const
{
	std::vector<unsigned int> waveIndices;
	for (unsigned int waveIndex = 0; waveIndex < nmbWaves(); ++waveIndex)
		if (waveName(waveIndex).Contains(waveNamePattern))
			waveIndices.push_back(waveIndex);
	return waveIndices;
}


/// generates list of production amplitude indices that match name pattern
inline
std::vector<unsigned int>
TFitResult::prodAmpIndicesMatchingPattern(const std::string& ampNamePattern) const
{
	std::vector<unsigned int> prodAmpIndices;
	for (unsigned int prodAmpIndex = 0; prodAmpIndex < nmbProdAmps(); ++prodAmpIndex)
		if (prodAmpName(prodAmpIndex).Contains(ampNamePattern))
			prodAmpIndices.push_back(prodAmpIndex);
	return prodAmpIndices;
}


/// generates list of pairs of production amplitude indices of same rank for a pair of waves
inline
std::vector<std::pair<unsigned int, unsigned int> >
TFitResult::prodAmpIndexPairsForWaves(const unsigned int waveIndexA,
                                      const unsigned int waveIndexB) const
{
	std::vector<std::pair<unsigned int, unsigned int> > prodAmpIndexPairs;
	const std::vector<unsigned int> prodAmpIndicesA = prodAmpIndicesForWave(waveIndexA);
	const std::vector<unsigned int> prodAmpIndicesB = prodAmpIndicesForWave(waveIndexB);
	for (unsigned int countAmpA = 0; countAmpA < prodAmpIndicesA.size(); ++countAmpA) {
		const unsigned int ampIndexA = prodAmpIndicesA[countAmpA];
		const int          ampRankA  = rankOfProdAmp(ampIndexA);
		// find production amplitude of wave B with same rank
		int ampIndexB = -1;
		for (unsigned int countAmpB = 0; countAmpB < prodAmpIndicesB.size(); ++countAmpB)
			if (rankOfProdAmp(prodAmpIndicesB[countAmpB]) == ampRankA) {
				ampIndexB = prodAmpIndicesB[countAmpB];
				break;
			}
		if (ampIndexB >=0)
			prodAmpIndexPairs.push_back(std::make_pair(ampIndexA, (unsigned int)ampIndexB));
	}
	return prodAmpIndexPairs;
}


/// calculates variance of a real-valued function of a spin density matrix element for wave A and wave B
inline
double
TFitResult::realValVariance(const unsigned int      waveIndexA,
                            const unsigned int      waveIndexB,
                            const TMatrixT<double>& jacobian) const  // Jacobian of real valued function (d f/ d Re[rho]   d f / d Im[rho])
{
	if (!_hasErrors)
		return 0;
	const TMatrixT<double> spinDensCov = spinDensityMatrixElemCov(waveIndexA, waveIndexB);  // 2 x 2 matrix
	const TMatrixT<double> jacobianT(TMatrixT<double>::kTransposed, jacobian);              // 2 x 1 matrix
	const TMatrixT<double> spinDensCovJT = spinDensCov * jacobianT;                         // 2 x 1 matrix
	const TMatrixT<double> cov           = jacobian * spinDensCovJT;                        // 1 x 1 matrix
	return cov[0][0];
}


#endif // USE_TFITRESULT


#endif  // TFITRESULT_HH


//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
