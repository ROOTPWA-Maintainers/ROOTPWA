//-----------------------------------------------------------
//
// Description:
//      Contraints an amplitude to real values
//
//
// Environment:
//      rootpwa
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef TPWAREALCONSTRAINT_HH
#define TPWAREALCONSTRAINT_HH

// Base Class Headers ----------------
#include "TPWAConstraint.h"

// Collaborating Class Headers -------
#include <complex>

// Collaborating Class Declarations --

class TPWARealConstraint : public TPWAConstraint {
public:

  // Constructors/Destructors ---------
  TPWARealConstraint(){}
  virtual ~TPWARealConstraint(){}

  virtual TPWAConstraint* clone(){return new TPWARealConstraint();}

   // Accessors -----------------------
  virtual int npar()const {return 1;} // returns number of free parameters 0,1 or even 2
 virtual std::string type()const {return "RealConstraint";}
  virtual std::string parname(unsigned int i) const {return "_RE";}

   // Operations ----------------------
  virtual std::complex<double> cAmp(const std::complex<double>& amp)
  {return amp.real();}
  virtual std::complex<double> dampdpar(unsigned int i)
  {return std::complex<double>(1,0);}
};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
