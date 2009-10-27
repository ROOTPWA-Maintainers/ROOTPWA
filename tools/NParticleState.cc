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
/////////////////////////////////////////////////////////////////////////////-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class NParticleState
//      see NParticleState.hh for details
//
// Environment:
//      Software developed for the COMPASS Experiment at CERN
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

// This Class' Header ------------------
#include "NParticleState.h"

// C/C++ Headers ----------------------
#include <iostream>

// Collaborating Class Headers --------
#include "FSParticle.h"


// Class Member definitions -----------


NParticleState::NParticleState()
  : _n(0), _q(0), _p(0,0,0,0), _beam(0,0,0,0)
{}

NParticleState::~NParticleState()
{}

bool
NParticleState::addParticle(FSParticle* pi)
{
   _fspart.push_back(pi);
  _n=_fspart.size();
  _q+=pi->q();
  _p+=pi->p();
  return true;
}


void NParticleState::setBeam(const TLorentzVector& beam){
  _beam=beam;
}


TLorentzVector 
NParticleState::p() const 
{
  TLorentzVector result;
  //std::cout<<_fspart.GetEntriesFast()<<"   "<<_n<<std::endl;
  for(unsigned int i=0;i<_n;++i){
    FSParticle* pa=_fspart.at(i);
    if(pa!=NULL)result+=pa->p();
    else std::cout<<"pion not found!"<<std::endl;
  }
  return result;
}

double
NParticleState::Q2()
{
  TLorentzVector beam=_beam;
  // recalibrate beam -- assumes exclusivity!
  TVector3 dir=beam.Vect();
  double const mpi=0.13957;
  double k=sqrt(_p.E()*_p.E()-mpi*mpi)/dir.Mag();
  dir*=k;
  beam.SetVectM(dir,mpi);
  return (beam-_p).M2();
}


double
NParticleState::rapidity()
{
  return 0.5*TMath::Log((_p.E()+_p.Vect().Z())/(_p.E()-_p.Vect().Z()));
}

TLorentzVector 
NParticleState::pfs(unsigned int i) const 
{
  if(i<_n)return _fspart.at(i)->p();
  else return TLorentzVector();
}

int
NParticleState::qabs() const 
{
  int result=0;
  for(unsigned int i=0;i<_n;++i){
    FSParticle* pa=_fspart.at(i);
    if(pa!=NULL)result+=abs(pa->q());
    else std::cout<<"pion not found!"<<std::endl;
  }
  return result;
}

TVector3 
NParticleState::vertex() const
{
  TVector3 result;
  for(unsigned int i=0;i<_n;++i){
    FSParticle* pa=_fspart.at(i);
    if(pa!=NULL)result+=pa->v();
  }
  if(_n!=0)result*=1./_n;
  return result;
}

bool
NParticleState::Exclusive(double d){
  if(fabs(_p.E()-_beam.E())<d)return true;
  else return false;
}



bool 
NParticleState::isSubstate(NParticleState* motherstate){
  // check if all fspart in this state are also part of the mother
  for(unsigned int i=0;i<_n;++i){//loop over fspart
    FSParticle* myp=getParticle(i);
    //loop over fspart in mother
    unsigned int nm=motherstate->n();
    bool found=false;
    for(unsigned int j=0;j<nm;++j){ // loop mother fspart
      if(myp->IsEqual(motherstate->getParticle(j))){
	found=true;
	//std::cout<<"found pion "<<i<<" as pion "<<j<<" in mother"<<std::endl;
	break;
      }
    } // end loop over mother-fspart
    if(!found){
      //std::cout<<"did not find pion number "<<i<<std::endl;
      return false;
    }
    
  } // end loop over my fspart
  // if we survive here we found all fspart in mother
  return true;
}

bool 
NParticleState::isDisjunctFrom(NParticleState* isobar){
  // check if all fspart in this state are not equal from partner
  for(unsigned int i=0;i<_n;++i){//loop over fspart
    FSParticle* myp=getParticle(i);
    //loop over fspart in mother
    unsigned int nm=isobar->n();
    for(unsigned int j=0;j<nm;++j){ // loop mother fspart
      if(myp->IsEqual(isobar->getParticle(j))){
	return false;
      }
    }// end loop over isobar fspart
  }// end loop over my fspart
  return true;
}
  
