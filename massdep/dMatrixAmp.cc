
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

#include "dMatrixAmp.h"
#include "cnumTools.h"
#include <iostream>
using namespace ublas;
using namespace std;

cnum
dMatrixPole::gDec(unsigned int i) const {
  return cnum(sqrt(fm*fgamma(0,i)),0);
}

 double
 dMatrixPole::psp(double m, unsigned int i) const {
   return (*fpsp)[i]->Eval(m)/(*fpsp)[i]->Eval(fm);
 }

 
double 
dMatrixPole::gammaTot() const {
  unsigned int n=fgamma.size2();
  double result=0;
  for(unsigned int i=0;i<n;++i){
    result+=fgamma(0,i);
  }
  return result;
}

double 
dMatrixPole::gammaTot(double m) const {
  unsigned int n=fgamma.size2();
  double result=0;
  for(unsigned int i=0;i<n;++i){
    result+=fgamma(0,i)*psp(m,i);
  }
  return result;
}


cnum 
dMatrixPole::M2(double m) const {
   double m2=fm*fm;
   if(fBkg)return cnum(m2,0);
   else return cnum(m2,-gammaTot(m));
}



dMatrixAmp::dMatrixAmp(){};

dMatrixAmp::~dMatrixAmp(){};

void 
dMatrixAmp::addChannel(TF1* ch){
  fChannels.push_back(ch);
}

void 
dMatrixAmp::setNPoles(unsigned int n){
  fPoles.clear();
  fPoles.resize(n);
}
void 
dMatrixAmp::setNBkg(unsigned int n){
  fBkg.clear();
  fBkg.resize(n,dMatrixPole(true));
}



void 
dMatrixAmp::Setup( const rmatrix& mbare, 
	     const rmatrix& gamma, 
	     const cmatrix& production,
	     const rmatrix& mixing){
  unsigned int np=nPoles();
  unsigned int nb=nBkg();
  unsigned int nc=nChannels();
  fprod=production;
   fmixing=mixing;
   for(unsigned int ip=0;ip<np+nb;++ip){
     cerr << "ip="<< ip << endl;
     if(ip<np)fPoles[ip].setMass(mbare(0,ip));
     else fBkg[ip-np].setMass(mbare(0,ip));
       matrix_range<const rmatrix> gam(gamma,range(ip,ip+1),range(0,nc));
       cerr << gam << endl;
     if(ip<np)fPoles[ip].setChannels(&fChannels,gam);
     else fBkg[ip-np].setChannels(&fChannels,gam);
   }// end loop over poles
 }


cnum 
dMatrixAmp::amp(double m, unsigned int channel){
  cnum s(m*m,0);
  // build propagator matrix:
  unsigned int n=nPoles()+nBkg();
  cmatrix Dinv(n,n);
  cmatrix ADec(1,n);
  for(unsigned int i=0;i<n;++i){
    for(unsigned int j=0;j<n;++j){
      if(i==j){
	if(i<nPoles())Dinv(i,j)=getPole(i).M2(m)-s;
	else Dinv(i,j)=-getBkg(i-nPoles()).M2(m);
      } // end setting 
      else Dinv(i,j)=-fmixing(i,j); 
    }
  
    // build row-vector of decay amplitudes
    // (real, depend on channel)
    	if(i<nPoles())ADec(0,i)=getPole(i).gDec(channel);
	else ADec(0,i)=getBkg(i-nPoles()).gDec(channel);
  } // end loop over poles
  
  cmatrix D(n,n);
  ///cerr << "Inverting Matrix now" << endl;
  InvertMatrix(Dinv,D);
  //cerr << "Matrix inverted" << endl;

  cmatrix DAProd=prod(D,fprod);
  cmatrix ADecDAProd=prod(ADec,DAProd);

  return ADecDAProd(0,0);
}; /// full amplitude


// parameter mapping:
unsigned int 
dMatrixAmp::getNPar() const {
  unsigned int nbare=fPoles.size()+fBkg.size(); // bare masses
  // widths/decay couplings
  unsigned int ngamma=nbare*fChannels.size();
  // production constants (complex!) 
  unsigned int nprod = nbare*2;
  // mixing terms
  unsigned int nmix=0.5*nbare*(nbare - 1); // diagonal terms are unphysical!
                                           // symmetrie!!!
  return nbare+ngamma+nprod+nmix;

}

void
dMatrixAmp::getPar(double* par) const {
  unsigned int counter=0;
  // bare masses
  // loop over poles
  for(unsigned int ip=0;ip<nPoles();++ip){
    par[counter++]=fPoles[ip].m();
  }
  // loop over bkg
  for(unsigned int ip=0;ip<nBkg();++ip){
    par[counter++]=fBkg[ip].m();
  }
  // get gammas
  // loop over channels
  for(unsigned int ic=0;ic<nChannels();++ic){
    // loop over poles
    for(unsigned int ip=0;ip<nPoles();++ip){
      par[counter++]=fPoles[ip].gamma(ic).real();
    }
    // loop over bkg
    for(unsigned int ip=0;ip<nBkg();++ip){
      par[counter++]=fBkg[ip].gamma(ic).real();
    }
  }// end loop over channels
  // get production vector
  unsigned int nbare=nBkg()+nPoles();
  for(unsigned int ip=0;ip<nbare;++ip){
    par[counter++]=fprod(ip,0).real();
    par[counter++]=fprod(ip,0).imag();
  }
  // get mixing off diagonal elements
  for(unsigned int i=0;i<nbare;++i){
    for(unsigned int j=i+1;j<nbare;++j){
      if(i!=j)par[counter++]=fmixing(i,j);
    }
  }
  assert(counter==getNPar());

}

void
dMatrixAmp::setPar(const double* par) {
  unsigned int counter=0;
  // bare masses
  unsigned int nbare=nBkg()+nPoles();
  rmatrix mbare(1,nbare); 
  // loop over poles&bkg
  for(unsigned int ip=0;ip<nbare;++ip){
    mbare(0,ip)=par[counter++];
  }
  // get gammas
  rmatrix gamma(nbare,nChannels());
  // loop over channels
  for(unsigned int ic=0;ic<nChannels();++ic){
    // loop over poles&bkg
    for(unsigned int ip=0;ip<nbare;++ip){
      gamma(ip,ic)=par[counter++];
    }
  }// end loop over channels
  // get production vector
  cmatrix prod(nbare,1);
  for(unsigned int ip=0;ip<nbare;++ip){
    double real=par[counter++];
    double imag=par[counter++];
    prod(ip,0)=cnum(real,imag);
  }
  // get mixing  -- off diagonal elements
  rmatrix mix(nbare,nbare);
  for(unsigned int i=0;i<nbare;++i){
    for(unsigned int j=i;j<nbare;++j){
      if(i!=j){
	mix(i,j)=par[counter++];
	mix(j,i)=mix(i,j);
      }
      else mix(i,j)=0;
    }
  }
  assert(counter==getNPar());

  Setup(mbare,gamma,prod,mix);
}
