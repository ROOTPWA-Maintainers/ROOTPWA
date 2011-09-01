
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

#include "dMatrixAmp.h"
#include "cnumTools.h"


cnum
dMatrixPole::gDec(unsigned int i) {
  return cnum(sqrt(fm*fgamma[i]),0);
}

 double
 dMatrixPole::psp(double m, unsigned int i){
   return fpsp[i]->Eval(m)/fpsp[i]->Eval(fm);
 }

 
void 
dMatrixPole::addChannel(TF1* psp, double gamma){
  fpsp.push_back(psp);
  fgamma.push_back(gamma);
}

double 
dMatrixPole::gammaTot(){
  unsigned int n=fgamma.size();
  double result=0;
  for(unsigned int i=0;i<n;++i){
    result+=fgamma[i];
  }
  return result;
}

double 
dMatrixPole::gammaTot(double m){
  unsigned int n=fgamma.size();
  double result=0;
  for(unsigned int i=0;i<n;++i){
    result+=fgamma[i]*psp(m,i);
  }
  return result;
}


cnum 
dMatrixPole::M2(double m){
   double m2=fm*fm;
   return cnum(m2,gammaTot(m));
}






dMatrixAmp::dMatrixAmp(){};

dMatrixAmp::~dMatrixAmp(){};


// first add all channels
void 
dMatrixAmp::addChannels(map<string, TF1*> channels){}

// then add poles
void 
dMatrixAmp::addPole(const dMatrixPole& pole){}


cnum 
dMatrixAmp::amp(double m, unsigned int channel){
  cnum s(m*m,0);
  // build propagator matrix:
  unsigned int n=nPoles();
  cmatrix Dinv(n,n);
  cmatrix AProd(n,1);
  cmatrix ADec(1,n);
  for(unsigned int i=0;i<n;++i){
    for(unsigned int j=0;j<n;++j){
      if(i==j){
	Dinv(i,j)=getPole(i).M2(m)-s;
      } // end setting 
      else Dinv(i,j)=0; // no mixing at the moment
    }
    // build collumn vector of production amplitudes 
    // (complex, independent of channel)
    AProd(i,0)=getPole(i).gProd();
    // build row-vector of decay amplitudes
    // (real, depend on channel)
    ADec(0,i)=getPole(i).gDec(channel);
  } // end loop over poles
  
  cmatrix D;
  InvertMatrix(Dinv,D);

 

  return 0;
}; /// full amplitude
