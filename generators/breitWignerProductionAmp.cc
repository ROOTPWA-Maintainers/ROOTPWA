//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class breitWignerProductionAmp
//      see breitWignerProductionAmp.hh for details
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#include "breitWignerProductionAmp.h"


using namespace rpwa;

breitWignerProductionAmp::breitWignerProductionAmp(double mass,
                                                   double width,
                                                   std::complex<double> coupling)
	: _mass(mass),
	  _m2(mass*mass),
	  _width(width),
	  _mw(mass*width),
	  _coupling(coupling) { }


std::complex<double>
breitWignerProductionAmp::amp(double mass) {

	std::complex<double> denom(_m2-mass*mass,-_mw);
	std::complex<double> nom(_mw);

	return _coupling*nom/denom;

}

