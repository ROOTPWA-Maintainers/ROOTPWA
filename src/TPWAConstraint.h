//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Abstract Base class for constraints
//
//
// Environment:
//      Software developed for the PANDA Detector at FAIR.
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef TPWACONSTRAINT_HH
#define TPWACONSTRAINT_HH

// Base Class Headers ----------------


// Collaborating Class Headers -------
#include <complex>
#include <string>

// Collaborating Class Declarations --

class TPWAConstraint {
public:

  // Constructors/Destructors ---------
  TPWAConstraint(){}
  virtual ~TPWAConstraint(){}

  virtual TPWAConstraint* clone()=0;

   // Accessors -----------------------
  virtual int npar() const =0; // returns number of free parameters 0,1 or even 2
  virtual std::string type() const =0;
  virtual std::string parname(unsigned int i) const =0;
  // Modifiers -----------------------


  // Operations ----------------------
  virtual std::complex<double> cAmp(const std::complex<double>& amp)=0;
  virtual std::complex<double> dampdpar(unsigned int i)=0;


private:

  // Private Data Members ------------


  // Private Methods -----------------

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
