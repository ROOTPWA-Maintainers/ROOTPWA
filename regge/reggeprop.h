//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Reggeized propagator
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef REGGEPROP_HH
#define REGGEPROP_HH

// Base Class Headers ----------------


// Collaborating Class Headers -------
#include <ostream> // remove if you do not need streaming op

// Collaborating Class Declarations --

#include <complex>

class reggeprop {
public:

  // Constructors/Destructors ---------
  reggeprop(){};


  // Operators
 

  // Accessors -----------------------
  std::complex<double> amp(double t, double s);

  // Modifiers -----------------------


  // Operations ----------------------

private:

  // Private Data Members ------------


  // Private Methods -----------------

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
