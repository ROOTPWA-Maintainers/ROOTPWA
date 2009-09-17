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
//      Supports Constraint amplitudes!
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef TPWALIKELIHOODC_HH
#define TPWALIKELIHOODC_HH

// Base Class Headers ----------------


// Collaborating Class Headers -------
#include <vector>
using std::vector;
#include <string>
using std::string;
#include <complex>
using std::complex;


#include "TString.h"
#include "TMatrixD.h"
#include "Math/IFunction.h"
using namespace ROOT;

#include "integral.h"
#include "matrix.h"
#include "TPWAAmp.h"


// Collaborating Class Declarations --
class TString;
class TCMatrix;

typedef matrix<complex<double> > intmat;


class TPWALikelihoodC : public Math::IGradientFunctionMultiDim {
public:

  // Constructors/Destructors ---------
  TPWALikelihoodC();
  ~TPWALikelihoodC();
  TPWALikelihoodC* Clone() const {return new TPWALikelihoodC(*this);}
  // Accessors -----------------------
  
  unsigned int NDim() const;
  string parname(unsigned int i)const {return _parnames.at(i);}
  const vector<string>& wavetitles() const {return _wavenames;}
  double parthreshold(unsigned int i) const {return _parthresholds[i];}
  double dLcache(unsigned int i) const {return _dLcache[i];}
  unsigned int ncalls() const {return _ncalls;}
  double Ltime() const {return _Ltime;}
  double Ntime() const {return _Ntime;}
  unsigned int nevents() const {return _nevents;}
  const integral& normInt()const {return _I;}
  

  // Modifiers -----------------------
  //void SetRank(unsigned int rank);
  void UseNormalizedAmps(bool flag=true){_useNorm=flag;}
  // Set wavelist filename:
  void Init(const TString& wavelist, unsigned int rank,
	    const TString& norm, const TString& acceptance,
	    int scaleevents, const TString constraints="");
  void SetQuiet(bool flag=true){_quiet=flag;}


  // Operations ----------------------
  // Load Amplitudes into memory
  void LoadAmplitudes();
  //void LoadIntegrals(const TString& norm, const TString& acceptance);
  void getIntCMatrix(TCMatrix& integr,
		     TCMatrix& acceptance);

  // convert parameter array into amplitudes:
  void partoamp(const double* x) const;
  void amptopar(double* x) const;
  // Note: amps which do not exist in higher ranks are NOT built!
// convert parameters error matrix into error matrix of amplitudes
  void buildCAmps(const double* x,
		  vector<complex<double> >& V,
		  vector<pair<int,int> >& indices,
		  vector<string>& names,
		  const TMatrixD& errpar,
		  TMatrixD& erramp,
		  bool withFlat=false);
  


  virtual double DoEval(const double*) const;
  virtual double DoDerivative(const double*, unsigned int) const;

  // Calculate Function f and Derivatives df in one go
  virtual void FdF(const double* x, double& f, double* df) const;
  void Print() const ;
private:

  // Private Data Members ------------
  unsigned int _rank;
  unsigned int _dim;
  unsigned int _nposrefl; // number of positive reflectivity waves 
  unsigned int _nnegrefl; // number of negative reflectivity waves
  unsigned int _nevents;
  unsigned int _nwaves;
  mutable unsigned int _ncalls;
  mutable double _Ltime; // total time spent calculating L
  mutable double _Ntime; // total time spent calculating Normalization

  bool _quiet;
  bool _useNorm; // use normalized amplitudes

  // production amplitude vectors per rank
  vector<vector<TPWAAmp>*> _V; // production amplitudes
  mutable TPWAAmp _Vflat; // amplitude for flat wave
  unsigned int _NumNullConstraints; // number of null constraint waves
  vector<string> _wavenames;
  vector<string> _parnames;
  vector<int> _reflect; // reflectivity of parameter
  vector<double> _thresholds; // mass threshold for a wave
  vector<double> _parthresholds; // mass threshold for a parameter
  vector<unsigned int> _wmap; // map: waveindex -> id in normalization integral
  vector<unsigned int> _accmap; // map: waveindex -> id in accpetance integral
  mutable vector<double> _parcache; // parameter cache for derivative calc.
  mutable vector<double> _dLcache;  // cache for derivatives

  // derivative with resp to real/imaginary part resp.
  // the dL are NOT well defined complex numbers!
  mutable vector<vector<complex<double> >* > _dL; // complex derivatives
  // derivative contribution of single event:
  mutable vector<vector<complex<double> >* > _dl; // contribution of each event
  // normalization contribution of derivative:
  mutable vector<vector<complex<double> >* > _dLN;


  // data cache:
  vector<vector<complex<double> >* > _data; 

  // normalization integrals 
  integral _I;
  matrix<complex<double> > _mat;

  integral _Acc;
  mutable matrix<complex<double> > _accmat; // has to be mutable because
                                            // of matrix accessor 
                                            // declaration not const

  // Private Methods -----------------
  void ClearCache();
  int geteps(TString wavename);
  // add constraints return number of constraints that have been added
  // used in Init
  unsigned int addConstraints(TString ConstraintsFile);

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
