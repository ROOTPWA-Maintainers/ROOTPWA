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
//      Fit summary for one mass bin
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

// Base Class Headers ----------------
#include "TObject.h"

// Collaborating Class Headers -------
#include <ostream> // remove if you do not need streaming op
#include "TCMatrix.h"
#include <vector>
#include <map>
#include "TComplex.h"
#include "TString.h"

// Collaborating Class Declarations --

using std::vector;
using std::map;
using std::pair;

class TFitBin : public TObject{
public:

  // Constructors/Destructors ---------
  TFitBin();
  virtual ~TFitBin();

  // Accessors -----------------------


  // Modifiers -----------------------
  void fill(const vector<TComplex>& amps,
	    const vector<pair<int,int> >& indices,
	    const vector<TString>& wavenames,
	    int nevt,
	    unsigned int rawnevt,
	    Double_t mass,
	    const TCMatrix& integr,
	    const TMatrixD& errm,
	    Double_t logli,
	    Int_t rank);

  // Operations ----------------------

  Double_t norm(const char* tag);
  Double_t normI(Int_t i){return norm(_wavenames[i].Data());}
  Double_t intens();  // total intensity
  Double_t intens(const char* tag); // added intensity of waves containing tag
  Double_t intens(Int_t i);
  Double_t phase(Int_t i, Int_t j); // phase difference between wave i and j
  Double_t phaseErr(Int_t i, Int_t j); 
  Double_t coh(Int_t i, Int_t j); // coherence between wave i and j
  Double_t mass() const {return _mass;}
  Double_t logli() const {return _logli;}
  Double_t getInt(Int_t i, Int_t j, bool re){if(re)return getInt(i,j).Re();else return getInt(i,j).Im();}
  UInt_t rawEvents()const {return _rawevents;}

  Int_t nwaves(){return _wavetitles.size();}
  TString wavename(unsigned int i)const {return _wavenames[i];}
  TString waveDesignator(unsigned int i)const {return _wavetitles[i];}
  unsigned int namps() const { return _amps.size();}
  TComplex amp(unsigned int i) { return _amps.at(i);}
  void getParameters(double* par); // returns by filling the par array
  double getParameter(const char* name);

  Double_t err(const char* tag);
  Double_t err(Int_t i);

  void listwaves();
  void Reset();
  void PrintParameters();

private:

  // Private Data Members ------------
  vector<TComplex> _amps; // Fitted amplitudes
  vector<pair<int,int> > _indices; // indices of parameters in error matrix;
  vector<TString> _wavenames; // rank included!!!
  vector<TString> _wavetitles; // without rank
  map<int,int> _wavemap; // maps wave indices to indices in integral
  
  
  int _nevt; // number of events normalized;
  unsigned int _rawevents; // raw number of events in this bin
  Double_t _mass; // bin center;
  TCMatrix _int; // normalization
  TMatrixD _errm; // errormatrix;
  Double_t _logli; // log likelyhood
  Int_t _rank;
  
  Bool_t _hasErrors;

  // Private Methods -----------------
  TMatrixD getErr(unsigned int i){return getErr(_indices[i]);}
  TMatrixD getErr(pair<int,int>); // returns cov matrix for complex parameter i
  void getCov(const char* tag, TMatrixD& C, std::vector<int>& cpar);
  void getCov(int i, int j, TMatrixD& C);

  // Waves indices run over all ranks. In the integral each wave appears only
  // once (without notion of rank). So we have to map this:
  TComplex getInt(int i, int j); 

  TComplex spinDens(int i, int j);
  TMatrixD spinDensErr(int i, int j);
  TMatrixD M(const TComplex& c);
  // the rank should be encoded into the second parameter of the wave
  int rankofwave(int i);
  TString wavetitle(int i){if(_wavenames[i].Contains("V_"))return _wavenames[i](2,1000); else return _wavenames[i](3,1000);}

  void buildWaveMap();

public:
  ClassDef(TFitBin,8)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
