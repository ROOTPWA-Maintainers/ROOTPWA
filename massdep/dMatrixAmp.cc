
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

dMatrixAmp::dMatrixAmp(){};

dMatrixAmp::~dMatrixAmp(){};


// first add all channels
void 
dMatrixAmp::addChannels(map<string, TF1*> channels){}

// then add poles
void 
dMatrixAmp::addPole(double m, double width, map<string, double>& branchings){}

unsigned int 
dMatrixAmp::nPoles(){return 0;}

cnum
dMatrixAmp::PoleM2(unsigned int i){
  return 0;
}

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
	Dinv(i,j)=GetPole(i).M2()-s;
      } // end setting 
      else Dinv(i,j)=0; // no mixing at the moment
    }
    // build collumn vector of production amplitudes 
    // (complex, independent of channel)
    AProd(i,0)=GetPole(i).gProd();
    // build row-vector of decay amplitudes
    // (real, depend on channel)
    ADec(0,i)=GetPole(i).gDec(channel);
  } // end loop over poles
  
  cmatrix D;
  InvertMatrix(Dinv,D);

 

  return 0;
}; /// full amplitude
