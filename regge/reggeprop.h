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
//#include <ostream> // remove if you do not need streaming op

// Collaborating Class Declarations --

#include <complex>

class reggeprop {
public:

  // Constructors/Destructors ---------
 reggeprop(): spi(0.019479835) {};


  // Operators
 

  // Accessors -----------------------
  /// t : mass^2 of reggeon
  /// s : mass^2 of dipion system with reggeon between
  /// tout : momentum transfer below reggeon
  /// sin : mass^2 of system above
  /// s1,s2 : mass^2 of forward out going particles
  /* Sketch:

     sin --
           --
             --*-------s1   +
               |            |
	       t            s
	       |            |
	     --*-------s2   +
	   --
    tout --


   */


  std::complex<double> amp(double t, double s, double sin, double tout, double s1, double s2);

  double alphapi(double t); // pion trajectory
  double alphapiprime(); // slope of pion trajectory at 1GeV
  double S(double tin,double t, double tout, double s1, double s, double s2); // multiparticle kinematic function

  // Modifiers -----------------------


  // Operations ----------------------

private:

  // Private Data Members ------------
  const double spi;// squared charged pion mass

  // Private Methods -----------------

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
