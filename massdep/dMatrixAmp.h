///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2011 Sebastian Neubert
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
//      Resonance Mixing Amplitude in the D-Matrix Formalism
//      for References see 
//      arXiv:hep-ex/0706.1341v2
//      arXiv:hep-ph/9702339v1 
//
// Environment:
//      Software developed for the COMPASS experiment at CERN
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------
#include <complex>
#include <vector>
#include <map>
#include <fstream>
#include <string>

#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"

#include "Math/GSLIntegrator.h"
#include "Math/IFunction.h"
#include <complex>
#include "TCanvas.h"
#include "TF1.h"


 #include <boost/numeric/ublas/vector.hpp>
 #include <boost/numeric/ublas/vector_proxy.hpp>
 #include <boost/numeric/ublas/matrix.hpp>
 #include <boost/numeric/ublas/triangular.hpp>
 #include <boost/numeric/ublas/lu.hpp>
 #include <boost/numeric/ublas/io.hpp>

 namespace ublas = boost::numeric::ublas;
typedef std::complex<double> cnum;
 typedef ublas::matrix<cnum> cmatrix;
typedef ublas::matrix<double> rmatrix;

class dMatrixPole {
 public:
 dMatrixPole(): fBkg(false){}
 dMatrixPole(bool bkg): fBkg(bkg){}
  //dMatrixPole(double m, cnum prodAmp):fm(m), fgProd(prodAmp) {}
  ~dMatrixPole(){}

  // Modifiers
  void setChannels(std::vector<TF1*>* psp, const rmatrix& gammas)
  {fpsp=psp;fgamma=gammas;}
  //void setProdAmp(cnum a){fgProd=a;}
  void setMass(double m){fm=m;}
  void setGamma(double gamma, unsigned int i){fgamma(0,i)=gamma;}
  void setBkg(bool flag=true){fBkg=flag;}

  // Accessors
  cnum M2(double m) const ; // return complex pole position
  //cnum gProd()const {return fgProd;} // return production coupling
  cnum gDec(unsigned int i) const ;   // return decay coupling (real) for channel i
  cnum gamma(unsigned int i)const {return cnum(fgamma(0,i),0);} // return partial width for channel i
  double psp(double m, unsigned int i)const ; // return phase space element normalized to pole position
  double m() const {return fm;} 
  double gammaTot() const ; // total width (at resonance)
  double gammaTot(double m) const; // total mass dependent width
  
 private:
  double fm;
  rmatrix fgamma; // partial widths
  const std::vector<TF1*>* fpsp; // phase space functions for different channels
  //cnum fgProd; // production coupling of this state
  bool fBkg;   // is background? -> purely real M2
};


class dMatrixAmp {
 public:
  dMatrixAmp();
  ~dMatrixAmp();

  
  void setNPoles(unsigned int n);
  void setNBkg(unsigned int n);
  void addChannel(TF1*);

  void Setup( const rmatrix& mbare, // (nPoles+nBkg) row  vector of bare masses
	      const rmatrix& gamma,  // (nPoles+nBkg x nChannels) matrix of couplings
	      const cmatrix& production, // (nPoles+nBkg)column vector of production couplings
	      const rmatrix& mixing); // (nPoles+nBkg)^2 matrix of mixing
                                      // only off-diagonals of upper triangle count!
  

  // processor:
  cnum amp(double m, unsigned int channel); /// amplitude in a channel
  
  // parameter mapping:
  unsigned int getNPar() const;
  void setPar(const double* par);
  void getPar(double* par) const;
  

  // helpers:
  unsigned int nPoles()const {return fPoles.size();}
  const dMatrixPole& getPole(unsigned int i)const{return fPoles[i];}
  unsigned int nBkg()const {return fBkg.size();}
  const dMatrixPole& getBkg(unsigned int i) const {return fBkg[i];}
  unsigned int nChannels() const {return fChannels.size();}
  TF1* getPS(unsigned int i)const {return fChannels[i];}

 private:
  std::vector<dMatrixPole> fPoles;
  std::vector<dMatrixPole> fBkg; // background terms
  std::vector<TF1*> fChannels;
  cmatrix fprod; // production amplitudes
  rmatrix fmixing;

};
