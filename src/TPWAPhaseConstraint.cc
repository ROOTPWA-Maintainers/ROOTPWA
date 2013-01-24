//-----------------------------------------------------------
//
// Description:
//      Implementation of class TPWAPhaseConstraint
//      see TPWAPhaseConstraint.hh for details
//
// Environment:
//      rootpwa
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------



// This Class' Header ------------------
#include "TPWAPhaseConstraint.h"


TPWAPhaseConstraint::TPWAPhaseConstraint(double phase, TPWAAmp* master)
  : _phi(phase), _master(master)
{}


std::complex<double> 
TPWAPhaseConstraint::cAmp(const std::complex<double>& amp){
  // parameter 1 is the length of the vector in C-plane
  double phi1=std::arg(_master->amp());
  phi1+=_phi;
  std::complex<double> result(cos(phi1),sin(phi1));
  return amp.real()*result;
}

std::complex<double> 
TPWAPhaseConstraint::dampdpar(unsigned int i){
  double phi1=std::arg(_master->amp());
  phi1+=_phi;
  std::complex<double> result(cos(phi1),sin(phi1));
  return result;  
}
