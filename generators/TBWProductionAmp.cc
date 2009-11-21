//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class TBWProductionAmp
//      see TBWProductionAmp.hh for details
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

// This Class' Header ------------------
#include "TBWProductionAmp.h"

// C/C++ Headers ----------------------


// Collaborating Class Headers --------


// Class Member definitions -----------


TBWProductionAmp::TBWProductionAmp(double mass, double width)
  : _mass(mass),_m2(mass*mass), _width(width),_mw(mass*width)
{}


std::complex<double>
TBWProductionAmp::amp(double mass){

  std::complex<double> denom(_m2-mass*mass,-_mw);
  std::complex<double> nom(_mw);
  
  return nom/denom;
}
