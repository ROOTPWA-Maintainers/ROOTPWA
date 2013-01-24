//-----------------------------------------------------------
//
// Description:
//      Contraints an amplitude to a defined phase difference
//      with respect to another "master" amplitude
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

#ifndef TPWAPHASECONSTRAINT_HH
#define TPWAPHASECONSTRAINT_HH

// Base Class Headers ----------------
#include "TPWAConstraint.h"

// Collaborating Class Headers -------
#include <complex>
#include "TPWAAmp.h"

// Collaborating Class Declarations --
class TPWAAmp;


class TPWAPhaseConstraint : public TPWAConstraint {
public:

  // Constructors/Destructors ---------
  // Make sure master is on the heap!!!
  TPWAPhaseConstraint(double phase, TPWAAmp* master); 
  virtual ~TPWAPhaseConstraint(){}

  virtual TPWAConstraint* clone(){return new TPWAPhaseConstraint(*this);}

   // Accessors -----------------------
  virtual int npar() const {return 1;} // returns number of free parameters 0,1 or even 2
 virtual std::string type() const {return "PhaseConstraint";}
  virtual std::string parname(unsigned int i) const {return "_ABS";}


   // Operations ----------------------
  virtual std::complex<double> cAmp(const std::complex<double>& amp);
  virtual std::complex<double> dampdpar(unsigned int i);
private:
  double _phi;
  TPWAAmp* _master;
  

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
