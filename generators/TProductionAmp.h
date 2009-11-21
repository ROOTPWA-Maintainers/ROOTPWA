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

#ifndef TPRODUCTIONAMP_HH
#define TPRODUCTIONAMP_HH

// Base Class Headers ----------------


// Collaborating Class Headers -------
#include <complex>

// Collaborating Class Declarations --



class TProductionAmp {
public:

  // Constructors/Destructors ---------
  TProductionAmp(){};
  TProductionAmp(const std::complex<double>& a);
  virtual ~TProductionAmp(){;}

  // Accessors -----------------------
  virtual std::complex<double> amp(double mass) {return _amp;} //> simple one parameter case

private:

  // Private Data Members ------------
  std::complex<double> _amp;

  // Private Methods -----------------

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
