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
//      Implementation of class TPWAAmp
//      see TPWAAmp.h for details
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


// This Class' Header ------------------
#include "TPWAAmp.h"

// C/C++ Headers ----------------------
#include<iostream>
#include<sstream>

// Collaborating Class Headers --------


// Class Member definitions -----------


TPWAAmp::TPWAAmp(string name, int rank, double threshold, 
		 unsigned int normindex,
		 unsigned int accindex)
  : _name(name), _rank(rank), _threshold(threshold),
    _integralindex(normindex), _acceptanceindex(accindex),
    _constr(NULL)
{
  // set reflectivity 
  // Relevant format of wave-names: IGJPCM\eps... \eps=reflectivity
  if(name[6] == '-') _reflectivity= -1;
  else if(name[6] == '+') _reflectivity= 1;
  else {
    std::cout<<"Wrong wavename format. "
	     <<"Could not determine reflectivity! Aborting!"<<std::endl;
    throw;
  }
}


TPWAAmp::TPWAAmp(const TPWAAmp& amp)
  : _name(amp._name),_reflectivity(amp._reflectivity), _rank(amp._rank),_threshold(amp._threshold), _integralindex(amp._integralindex), _acceptanceindex(amp._acceptanceindex),_constr(NULL){
  if(amp._constr!=NULL)setConstraint(amp._constr->clone());
}

TPWAAmp::~TPWAAmp(){
  if(_constr!=0)delete _constr;
  _constr=0;
}

void
TPWAAmp::operator=(const TPWAAmp& amp){
  _name=amp._name;
  _reflectivity=amp._reflectivity;
  _rank=amp._rank;
  _threshold=amp._threshold;
  _integralindex=amp._integralindex;
  _acceptanceindex=amp._acceptanceindex;
  setConstraint(amp._constr->clone());
}


complex<double> operator*(TPWAAmp& lhs, TPWAAmp& rhs){
  return lhs.amp() * rhs.amp();
}

complex<double> operator*(TPWAAmp& lhs, const complex<double>& rhs){
  return lhs.amp() * rhs;
}

complex<double> operator*(const complex<double>& lhs,TPWAAmp& rhs){
  return lhs * rhs.amp();
}

std::ostream& operator<< (std::ostream& s, const TPWAAmp& me){
  s << me.name() << " " << me.amp();
  return s;
}

const complex<double>&
TPWAAmp::amp() const {
  return _cached;
}

complex<double> 
TPWAAmp::updateAmp() {
  if(_constr!=NULL)_cached=_constr->cAmp(_amp);
  else _cached=_amp;
  return _cached;
}


complex<double> 
TPWAAmp::setPar(double* par){
  _amp.real()=par[0];
  _amp.imag()=par[1];
  updateAmp();
  return amp();
}

int
TPWAAmp::npar() const {
  if(_constr==NULL)return 2;
  else return _constr->npar();
}

string 
TPWAAmp::type() const {
  if(_constr==NULL)return "UnConstraint";
  else return _constr->type();
}

string
TPWAAmp::parname(unsigned int i) const {
  std::stringstream name;
  name << "V" << _rank << "_" << _name;
  if(_constr==NULL){
    if(i==0)name << "_RE";
    else if(i==1)name << "_IM";
  }
  else name << _constr->parname(i);
  return name.str();
}


double 
TPWAAmp::par(unsigned int i) const {
  return i==0 ? _amp.real() : _amp.imag();
}




void 
TPWAAmp::setConstraint(TPWAConstraint* c){
  if(_constr!=NULL){
    std::cerr << "*** " << this->name() << std::endl;
    std::cerr << "*** OVERRIDING PREVIOUSLY INSTALLED CONSTRAINT! " 
	      << std::endl;
    delete _constr;
  }
  _constr=c;
  updateAmp();
} 

complex<double> 
TPWAAmp::dampdpar(unsigned int i) const {
   if(_constr==NULL){
     if(i==0)return complex<double>(1,0);
     else if(i==1)return complex<double>(0,1);
     else return complex<double>(0,0);
   }
   else {
     return _constr->dampdpar(i);
   }
}
