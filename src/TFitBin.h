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
#include <vector>
#include <map>
#include <string>
#include <complex>

#include "TObject.h"
#include "TComplex.h"
#include "TMatrixT.h"
#include "TString.h"

#include "TCMatrix.h"


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
  double       massBinCenter() const { return _mass;              }
  double       logLikelihood() const { return _logli;             }
  unsigned int rank         () const { return _rank;              }
  unsigned int nmbEvents    () const { return _rawevents;         }
  unsigned int nmbWaves     () const { return _wavetitles.size(); }
  unsigned int nmbProdAmps  () const { return _amps.size();       }

  TString waveName    (const unsigned int waveIndex)    const { return _wavetitles[waveIndex];   }
  TString prodAmpName (const unsigned int prodAmpIndex) const { return _wavenames[prodAmpIndex]; }

  double      fitParameter   (const std::string& parName)   const;
  inline void fitParameters  (double*            parArray)  const;
  double      fitParameterCov(const unsigned int parIndexA,
			      const unsigned int parIndexB) const { return _errm[parIndexA][parIndexB]; }

  std::complex<double> spinDensityMatrixElem   (const unsigned int waveIndexA,
						const unsigned int waveIndexB) const;
  TMatrixT<double>     spinDensityMatrixElemCov(const unsigned int waveIndexA,
 						const unsigned int waveIndexB) const;

  double intensity   (const unsigned int waveIndex)       const { return intens(waveIndex);       }  // intensity of single wave
  double intensityErr(const unsigned int waveIndex)       const { return err(waveIndex);          }  // corresponding error
  double intensity   (const char*        waveNamePattern) const { return intens(waveNamePattern); }  // intensity sum of waves matching name pattern
  double intensityErr(const char*        waveNamePattern) const { return err(waveNamePattern);    }  // corresponding error
  double intensity   ()                                   const { return intens("");              }  // total intensity
  double intensityErr()                                   const { return err("");                 }  // corresponding error

  double coherence   (const unsigned int waveIndexA,  // coherence of wave A and wave B
		      const unsigned int waveIndexB) const { return coh(waveIndexA, waveIndexB); }
  double coherenceErr(const unsigned int waveIndexA,  // corresponding error
		      const unsigned int waveIndexB) const { return 0; }
  double overlap     (const unsigned int waveIndexA,  // overlap between wave A and wave B
		      const unsigned int waveIndexB) const { return 0; }
  double overlapErr  (const unsigned int waveIndexA,  // corresponding error
		      const unsigned int waveIndexB) const { return 0; }

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
  inline static TMatrixT<double> matrixRepr(const TComplex& c);  // matrix representation of complex number

  inline int rankOfProdAmp(const unsigned int prodAmpIndex) const;

  inline std::vector<unsigned int>
  prodAmpIndicesMatchingPattern(const std::string& ampNamePattern ) const;
  std::vector<unsigned int>
  prodAmpIndicesForWave(const unsigned int waveIndex) const
  { return prodAmpIndicesMatchingPattern(waveName(waveIndex).Data()); }
  inline std::vector<std::pair<unsigned int, unsigned int> >
  prodAmpIndexPairsForWaves(const unsigned int waveIndexA,
			    const unsigned int waveIndexB) const;

  TMatrixT<double>        prodAmpCovariance(const std::vector<unsigned int>& prodAmpIndices) const;
  inline TMatrixT<double> prodAmpCovariance(const std::vector<std::pair<unsigned int, unsigned int> >& prodAmpIndexPairs) const;


  // Private Data Members ------------
  std::vector<TComplex> _amps; // Fitted amplitudes
  std::vector<std::pair<int,int> > _indices; // indices of parameters in error matrix;
  std::vector<TString> _wavenames; // rank included!!!
  std::vector<TString> _wavetitles; // without rank
  std::map<int,int> _wavemap; // maps wave indices to indices in integral
  
  
  int _nevt; // number of events normalized;
  unsigned int _rawevents; // raw number of events in this bin
  Double_t _mass; // bin center;
  TCMatrix _int; // normalization
  TMatrixD _errm; // errormatrix;
  Double_t _logli; // log likelyhood
  Int_t _rank;
  
  Bool_t _hasErrors;

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


// returns parameter array; assumes that array has sufficient size
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


// returns matrix representation of complex number
//   c.re  -c.im
//   c.im   c.re
//!!! possible optimization: return TMatrixTLazy
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


// extracts rank of production amplitude from its name
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


// generates list of production amplitude indices that match name pattern
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


// generates list of pairs of production amplitude indices of same rank for a pair of waves
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


// constructs covariance matrix for production amplitudes specified by index pair list
inline
TMatrixT<double>
TFitBin::prodAmpCovariance(const std::vector<std::pair<unsigned int, unsigned int> >& prodAmpIndexPairs) const
{
  std::vector<unsigned int> prodAmpIndices;
  // copy wave A production amplitude indices
  for (unsigned int i = 0; i < prodAmpIndexPairs.size(); ++i)
    prodAmpIndices.push_back(prodAmpIndexPairs[i].first);
  // copy wave B production amplitude indices
  for (unsigned int i = 0; i < prodAmpIndexPairs.size(); ++i)
    prodAmpIndices.push_back(prodAmpIndexPairs[i].second);
  return prodAmpCovariance(prodAmpIndices);
}


#endif  // TFITBIN_HH


//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
