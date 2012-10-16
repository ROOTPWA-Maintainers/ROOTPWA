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
  void printAmpsGenPW(std::ostream& s) const;

  // accessors that allow copying of TFitBin
  const std::vector<TComplex>&             prodAmps()         const { return _amps;       }
  const std::vector<std::pair<int, int> >& fitParCovIndices() const { return _indices;    }
  const std::vector<TString>&              prodAmpNames()     const { return _wavenames;  }
  const std::vector<TString>&              waveNames()        const { return _wavetitles; }
  const std::map<int,int>&                 prodAmpIndexMap()  const { return _wavemap;    }

  int             normNmbEvents()        const { return _nevt;      }
  unsigned int    nmbEvents()            const { return _rawevents; }
  Double_t        massBinCenter()        const { return _mass;      }
  const TCMatrix& normIntegral()         const { return _int;       }
  const TMatrixD& fitParCovMatrix()      const { return _errm;      }
  Double_t        logLikelihood()        const { return _logli;     }
  Int_t           rank()                 const { return _rank;      }
  Bool_t          fitParCovMatrixValid() const { return _hasErrors; }

private:

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


#endif  // TFITBIN_HH


//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
