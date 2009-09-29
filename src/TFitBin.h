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
//      Data storage class for PWA fit result of one mass bin
//
//
// Environment:
//      Software developed for the COMPASS experiment at CERN
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#ifndef TFITBIN_HH
#define TFITBIN_HH


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

#include "utilities.h"
#include "TCMatrix.h"


/// \brief data storage class for PWA fit result of one kinematic bin
class TFitBin : public TObject {
public:

  TFitBin();
  virtual ~TFitBin();

  void fill(const std::vector<TComplex>&             prodAmps,
	    const std::vector<std::pair<int, int> >& indices,
	    const std::vector<TString>&              prodAmpNames,
	    const int                                nevt,
	    const unsigned int                       nmbEvents,
	    const double                             massBinCenter,
	    const TCMatrix&                          normIntegral,
	    const TMatrixD&                          fitParCovMatrix,
	    const double                             logLikelihood,
	    const int                                rank);

  // proposed new interface
  double       massBinCenter() const { return _mass;              }  ///< returns center value of mass bin
  double       logLikelihood() const { return _logli;             }  ///< returns log(likelihood) at minimum
  unsigned int rank         () const { return _rank;              }  ///< returns rank of fit
  unsigned int nmbEvents    () const { return _rawevents;         }  ///< returns number of events in bin
  unsigned int nmbWaves     () const { return _wavetitles.size(); }  ///< returns number of waves in fit
  unsigned int nmbProdAmps  () const { return _amps.size();       }  ///< returns number of production amplitudes

  TString waveName    (const unsigned int waveIndex)    const { return _wavetitles[waveIndex];   }  ///< returns name of wave at index
  TString prodAmpName (const unsigned int prodAmpIndex) const { return _wavenames[prodAmpIndex]; }  ///< returns name of production amplitude at index

  double      fitParameter   (const std::string& parName)   const;
  /// returns covariance of fit parameters at index A and B
  double      fitParameterCov(const unsigned int parIndexA,
			      const unsigned int parIndexB) const { return _errm[parIndexA][parIndexB]; }
  inline void fitParameters  (double*            parArray)  const;

  /// returns production amplitude value at index
  std::complex<double>    prodAmp   (const unsigned int prodAmpIndex) const { return std::complex<double>(_amps[prodAmpIndex].Re(), _amps[prodAmpIndex].Im()); }
  inline TMatrixT<double> prodAmpCov(const unsigned int prodAmpIndex) const;
  TMatrixT<double>        prodAmpCov(const std::vector<unsigned int>&                           prodAmpIndices)    const;
  inline TMatrixT<double> prodAmpCov(const std::vector<std::pair<unsigned int, unsigned int> >& prodAmpIndexPairs) const;
  inline TMatrixT<double> prodAmpCov(const std::vector<unsigned int>& prodAmpIndicesA,
				     const std::vector<unsigned int>& prodAmpIndicesB) const;

  inline std::complex<double> normIntegral(const unsigned int waveIndexA,
					   const unsigned int waveIndexB) const;

  std::complex<double> spinDensityMatrixElem   (const unsigned int waveIndexA,
						const unsigned int waveIndexB) const;
  TMatrixT<double>     spinDensityMatrixElemCov(const unsigned int waveIndexA,
 						const unsigned int waveIndexB) const;

  /// returns intensity of single wave at index
  double intensity   (const unsigned int waveIndex)       const { return spinDensityMatrixElem(waveIndex, waveIndex).real();         }
  /// returns error of intensity of single wave at index
  double intensityErr(const unsigned int waveIndex)       const { return sqrt(spinDensityMatrixElemCov(waveIndex, waveIndex)[0][0]); }
  double intensity   (const char*        waveNamePattern) const;
  double intensityErr(const char*        waveNamePattern) const;
  double intensity   ()                                   const { return intens(""); }  ///< returns total intensity
  double intensityErr()                                   const { return err("");    }  ///< returns error of total intensity

  double phaseNew    (const unsigned int waveIndexA,
		      const unsigned int waveIndexB) const;
  double phaseErrNew (const unsigned int waveIndexA,
		      const unsigned int waveIndexB) const;
  double coherence   (const unsigned int waveIndexA,
		      const unsigned int waveIndexB) const;
  double coherenceErr(const unsigned int waveIndexA,
		      const unsigned int waveIndexB) const;
  double overlap     (const unsigned int waveIndexA,
		      const unsigned int waveIndexB) const;
  double overlapErr  (const unsigned int waveIndexA,
		      const unsigned int waveIndexB) const;

  inline std::ostream& printProdAmpNames(std::ostream& out = std::cout) const;
  inline std::ostream& printWaveNames   (std::ostream& out = std::cout) const;
  inline std::ostream& printProdAmps    (std::ostream& out = std::cout) const;


  // Operations ----------------------
  Double_t norm(const char* tag) const;
  Double_t normI(Int_t i) const {return norm(_wavenames[i].Data());}
  Double_t intens() const;  // total intensity
  Double_t intens(const char* tag) const; // added intensity of waves containing tag
  Double_t intens(Int_t i) const;
  Double_t phase(Int_t i, Int_t j) const; // phase difference between wave i and j
  Double_t phaseErr(Int_t i, Int_t j) const; 
  Double_t coh(Int_t i, Int_t j) const; // coherence between wave i and j
  Double_t mass() const {return _mass;}
  Double_t logli() const {return _logli;}
  Double_t getInt(Int_t i, Int_t j, bool re) const {if(re)return getInt(i,j).Re();else return getInt(i,j).Im();}
  UInt_t rawEvents() const {return _rawevents;}

  Int_t nwaves() const {return _wavetitles.size();}
  TString wavename(unsigned int i) const {return _wavenames[i];}
  TString waveDesignator(unsigned int i) const {return _wavetitles[i];}
  unsigned int namps() const { return _amps.size();}
  TComplex amp(unsigned int i)  const { return _amps.at(i);}
  void getParameters(double* par) const; // returns by filling the par array
  double getParameter(const char* name) const;

  Double_t err(const char* tag) const;
  Double_t err(Int_t i) const;

  void listwaves() const;
  void Reset();
  void PrintParameters() const;

private:

  // helper functions for proposed new interface
  inline static TMatrixT<double> matrixRepr(const TComplex& c);

  inline int rankOfProdAmp(const unsigned int prodAmpIndex) const;

  std::complex<double> normIntegralForProdAmp(const unsigned int prodAmpIndexA,
					      const unsigned int prodAmpIndexB) const;

  inline std::vector<unsigned int> waveIndicesMatchingPattern   (const std::string& waveNamePattern) const;
  inline std::vector<unsigned int> prodAmpIndicesMatchingPattern(const std::string& ampNamePattern ) const;
  std::vector<unsigned int>        prodAmpIndicesForWave        (const unsigned int waveIndex)       const
  { return prodAmpIndicesMatchingPattern(waveName(waveIndex).Data()); }  ///< returns indices of production amplitudes that belong to wave at index
  inline std::vector<std::pair<unsigned int, unsigned int> > prodAmpIndexPairsForWaves(const unsigned int waveIndexA,
										       const unsigned int waveIndexB) const;

  inline double realValVariance(const unsigned int      waveIndexA,
				const unsigned int      waveIndexB,
				const TMatrixT<double>& jacobian) const;

  std::vector<TComplex>             _amps;        ///< production amplitudes
  std::vector<std::pair<int, int> > _indices;     ///< indices of fit parameters in error matrix;
  std::vector<TString>              _wavenames;   ///< names of production amplitudes used in fit
  std::vector<TString>              _wavetitles;  ///< names of waves used in fit
  std::map<int,int>                 _wavemap;     ///< maps production amplitude indices to indices in normalization integral
  
  int          _nevt;       ///< number of events normalized (?)
  unsigned int _rawevents;  ///< number of events in bin
  Double_t     _mass;       ///< center value of mass bin
  TCMatrix     _int;        ///< normalization integral over full phase space without acceptance
  TMatrixD     _errm;       ///< covariance matrix of fit parameters
  Double_t     _logli;      ///< log(likelihood) at minimum
  Int_t        _rank;       ///< rank of fit
  Bool_t       _hasErrors;  ///< indicates whether bin has a valid covariance matrix
  

  // Private Methods -----------------
  TMatrixD getErr(unsigned int i) const {return getErr(_indices[i]);}
  TMatrixD getErr(std::pair<int,int>) const; // returns cov matrix for complex parameter i
  void getCov(const char* tag, TMatrixD& C, std::vector<int>& cpar) const;
  void getCov(int i, int j, TMatrixD& C) const;

  // Waves indices run over all ranks. In the integral each wave appears only
  // once (without notion of rank). So we have to map this:
  TComplex getInt(int i, int j) const;

  TComplex spinDens(int i, int j) const;
  TMatrixD spinDensErr(int i, int j) const;
  TMatrixD M(const TComplex& c) const;
  // the rank should be encoded into the second parameter of the wave
  int rankofwave(int i) const ;
  TString wavetitle(int i) const {if(_wavenames[i].Contains("V_"))return _wavenames[i](2,1000); else return _wavenames[i](3,1000);}

  void buildWaveMap();

public:
  ClassDef(TFitBin,8)

};


/// copies fit parameters into array; assumes that array has sufficient size
inline
void
TFitBin::fitParameters(double* parArray) const
{
  unsigned int countPar = 0;
  for (unsigned int i = 0; i < nmbProdAmps(); ++i) {
    parArray[countPar++] = _amps[i].Re();
    // check if amplitude has imaginary part
    if (_indices[i].second >= 0)
      parArray[countPar++] = _amps[i].Im();
  }
}


/// returns covariance matrix for a single production amplitude
inline
TMatrixT<double>
TFitBin::prodAmpCov(const unsigned int prodAmpIndex) const {
  // get parameter indices
  const int i = _indices[prodAmpIndex].first;
  const int j = _indices[prodAmpIndex].second;
  TMatrixT<double> cov(2, 2);
  cov[0][0] = _errm[i][i];
  if (j >= 0) {
    cov[0][1] = _errm[i][j];
    cov[1][0] = _errm[j][i];
    cov[1][1] = _errm[j][j];
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
TFitBin::prodAmpCov(const std::vector<std::pair<unsigned int, unsigned int> >& prodAmpIndexPairs) const
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
TFitBin::prodAmpCov(const std::vector<unsigned int>& prodAmpIndicesA,
		    const std::vector<unsigned int>& prodAmpIndicesB) const
{
  std::vector<unsigned int> prodAmpIndices;
  // copy wave A production amplitude indices
  prodAmpIndices.assign(prodAmpIndicesA.begin(), prodAmpIndicesA.end());
  // copy wave B production amplitude indices
  prodAmpIndices.insert(prodAmpIndices.end(), prodAmpIndicesB.begin(), prodAmpIndicesB.end());
  return prodAmpCov(prodAmpIndices);
}


/// returns normalization integral of wave A and wave B
inline
std::complex<double>
TFitBin::normIntegral(const unsigned int waveIndexA,
		      const unsigned int waveIndexB) const
{
  const TComplex norm = _int(waveIndexA, waveIndexB);
  return std::complex<double>(norm.Re(), norm.Im());
}


/// prints list of all production amplitudes
inline
std::ostream&
TFitBin::printProdAmpNames(std::ostream& out) const
{
  out << "Production amplitude names:" << std::endl;
  for (unsigned int i = 0; i < nmbProdAmps(); ++i)
    out << "    " << std::setw(3) << i << " " << prodAmpName(i) << std::endl;
  return out;
}


/// prints list of all production amplitudes
inline
std::ostream&
TFitBin::printWaveNames(std::ostream& out) const
{
  out << "Wave names:" << std::endl;
  for (unsigned int i = 0; i < nmbWaves(); ++i)
    out << "    " << std::setw(3) << i << " " << waveName(i) << std::endl;
  return out;
}


/// prints values and covariance matrices of all production amplitudes
inline
std::ostream&
TFitBin::printProdAmps(std::ostream& out) const
{
  out << "Production amplitudes:" << std::endl;
  for (unsigned int i = 0; i < nmbProdAmps(); ++i) {
    out << "    " << std::setw(3) << i << " " << prodAmpName(i) << " = "  << prodAmp(i)
	<< ", cov = " << prodAmpCov(i) << std::endl;
  }
  return out;
}


/// \brief returns matrix representation of complex number
///
///   c.re  -c.im
///   c.im   c.re
// !!! possible optimization: return TMatrixTLazy
inline
TMatrixT<double>
TFitBin::matrixRepr(const TComplex& c)
{
  TMatrixT<double> m(2, 2);
  m[0][0] = m[1][1] = c.Re();
  m[0][1] = -c.Im();
  m[1][0] =  c.Im();
  return m;
}


/// extracts rank of production amplitude from its name
inline
int
TFitBin::rankOfProdAmp(const unsigned int prodAmpIndex) const
{
  const TString ampName = prodAmpName(prodAmpIndex);
  if (ampName.Contains("flat"))
    return -1;
  else
    // the rank is encoded in second character of parameter name
    return TString(ampName(1, 1)).Atoi();
}


/// generates list of wave indices that match name pattern
inline
std::vector<unsigned int>
TFitBin::waveIndicesMatchingPattern(const std::string& waveNamePattern) const
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
TFitBin::prodAmpIndicesMatchingPattern(const std::string& ampNamePattern) const
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
TFitBin::prodAmpIndexPairsForWaves(const unsigned int waveIndexA,
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
TFitBin::realValVariance(const unsigned int      waveIndexA,
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


#endif  // TFITBIN_HH


//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
