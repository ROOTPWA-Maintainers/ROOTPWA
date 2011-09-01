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
// File and Version Information:
// $Id$
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

using namespace std;
 #include <boost/numeric/ublas/vector.hpp>
 #include <boost/numeric/ublas/vector_proxy.hpp>
 #include <boost/numeric/ublas/matrix.hpp>
 #include <boost/numeric/ublas/triangular.hpp>
 #include <boost/numeric/ublas/lu.hpp>
 #include <boost/numeric/ublas/io.hpp>

 namespace ublas = boost::numeric::ublas;
 typedef complex<double> cnum;
 typedef ublas::matrix<cnum> cmatrix;

class dMatrixPole {
 public:
 dMatrixPole(double m):fm(m){};
  ~dMatrixPole(){};

  // Modifiers
  void addChannel(TF1* psp, double gamma);
  void setProdAmp(cnum a){fgProd=a;}
  void setMass(double m){fm=m;}

  // Accessors
  cnum M2(double m); // return complex pole position
  cnum gProd(){return fgProd;} // return production coupling
  cnum gDec(unsigned int i);   // return decay coupling (real) for channel i
  cnum gamma(unsigned int i){return cnum(fgamma[i],0);} // return partial width for channel i
  double psp(double m, unsigned int i); // return phase space element normalized to pole position
  double gammaTot(); // total width (at resonance)
  double gammaTot(double m); // total mass dependent width
  
 private:
  double fm;
  vector<double> fgamma; // partial widths
  vector<TF1*> fpsp; // phase space functions for different channels
  cnum fgProd; // production coupling of this state

};




class dMatrixAmp {
 public:
  dMatrixAmp();
  ~dMatrixAmp();
  
  // first add all channels
  void addChannels(map<string, TF1*> channels);

  // then add poles
  void addPole(const dMatrixPole& pole);
  
  // processor:
  cnum amp(double m, unsigned int channel); /// amplitude in a channel
  
  // helpers:
  unsigned int nPoles(){return fPoles.size();}
  dMatrixPole& getPole(unsigned int i){return fPoles[i];}

 private:
  vector<dMatrixPole> fPoles;


};
