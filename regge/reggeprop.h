//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Reggeized propagator
//      Literatur:
//      Berger Phys Rev 166 (1968) 1525
//      Ascoli Phys Rev D9 (1974) 1963
//      Harris Z. Phys C9 (1981) 275
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
 reggeprop(): spi(0.019479835), _a(2.) {
    _alphapiprime = 1.+1./(1.-spi);
    
}; // m_pi^2


  // Operators
 

  // Accessors -----------------------
  /// t : mass^2 of reggeon t<0 for scattering
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


 
 std::complex<double> ampHarris(double t, double s, double sin, double tout, double s1, double s2);
 // Nomenclature following Ascoli
  std::complex<double> ampPDK(double t, double s, double sin, double tout, double s1, double s2);
   std::complex<double> ampFULL(double t, double s, double sin, double tout, double s1, double s2);
std::complex<double> ampSMU(double t, double s, double sin, double tout, double s1, double s2);
   

  double alphapi(double t) const; // pion trajectory
  double alphapiprime() const; // slope of pion trajectory at 1GeV
  double S(double tin,double t, double tout, double s1, double s, double s2); // multiparticle kinematic function

  // Modifiers -----------------------
  void setA(double a){_a=a;}

  // Operations ----------------------

private:

  // Private Data Members ------------
  const double spi;// squared charged pion mass
  double _alphapiprime; // slope of piontrajectory at a scale of 1 GeV
  double _a;            // Ascoli style damping factor
  // Private Methods -----------------

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
