//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Parameterization of Production amplitude
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef TBWPRODUCTIONAMP_HH
#define TBWPRODUCTIONAMP_HH

// Base Class Headers ----------------
#include "TProductionAmp.h"

// Collaborating Class Headers -------
#include <complex>

// Collaborating Class Declarations --



class TBWProductionAmp : public TProductionAmp {
public:

  // Constructors/Destructors ---------
  TBWProductionAmp(double mass, double width,
		   std::complex<double> coupling=std::complex<double>(1,0));
  virtual ~TBWProductionAmp(){;}

  // Accessors -----------------------
  virtual std::complex<double> amp(double mass); //> simple one parameter case

private:

  // Private Data Members ------------
  double _mass;
  double _m2;
  double _width;
  double _mw;
  std::complex<double> _coupling;

  // Private Methods -----------------

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
