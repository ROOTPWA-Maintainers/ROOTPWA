///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      A class to encapsulate PWA production amplitudes
//
//
// Environment:
//      ROOTPWA
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef TPWAAMP_HH
#define TPWAAMP_HH

// Base Class Headers ----------------


// Collaborating Class Headers -------
#include <ostream> // remove if you do not need streaming op
#include <string>
#include <complex>
#include "TPWAConstraint.h"
#include "TPWARealConstraint.h"

using std::string;
using std::complex;

// Collaborating Class Declarations --

class TPWAAmp {
public:

  // Constructors/Destructors ---------
  TPWAAmp(string name, int rank, 
	  double threshold=0, 
	  unsigned int normindex=0,
	  unsigned int accindex=0);
  TPWAAmp(const TPWAAmp&); // nontrivial CPCTOR! deep copy constraint!
  ~TPWAAmp();
  
  void operator =(const TPWAAmp&);

  // Operators
   friend std::ostream& operator<< (std::ostream& s, const TPWAAmp& me);

  // Arithmetic operators
  friend complex<double> operator*(TPWAAmp& lhs, TPWAAmp& rhs);
  friend complex<double> operator*(TPWAAmp& lhs, const complex<double>& rhs);
  friend complex<double> operator*(const complex<double>& lhs, TPWAAmp& rhs);

  // Accessors -----------------------
  string name() const {return _name;}
  string parname(unsigned int i) const; // return paramter name 
                                        // i=0,1 Real/imaginary
  double threshold() const {return _threshold;}
  string type() const;
  double par(unsigned int i) const;     // returns parameter (w/o constraint)
  const complex<double>& amp() const; // returns cached amplitude 
                               // (takes into account evtl. constraints)
  // derivative of amp wrt parameter i:
  complex<double> dampdpar(unsigned int i) const;

  complex<double> updateAmp();

  int rank() const {return _rank;}
  int reflectivity() const {return _reflectivity;}

  int npar() const ;  // returns number of parameters 2,1,0 depends on constraint
  unsigned int normindex() const {return _integralindex;} 
  unsigned int accindex() const {return _acceptanceindex;} 
  // index of normalization integral


  // Modifiers -----------------------
  // returns amplitude with contraints!
  complex<double> setPar(double* par); // recalculates and caches amp
  void setIntegralIndices(unsigned int inorm, unsigned int iacc)
  {_integralindex=inorm;_acceptanceindex=iacc;}
 

  // Operations ----------------------
  // recalculates and caches amp
  void setConstraint(TPWAConstraint* c);

  private:

  // Private Data Members ------------
  string _name;
  int _reflectivity;
  int _rank;
  double _threshold;
  complex<double> _amp;
  complex<double> _cached;
  // indices in normalization/acceptance integral
  unsigned int _integralindex;
  unsigned int _acceptanceindex;

  // a constraint is called by the object and
  // alters the amplitude. This usually reduces the
  // number of free parameters
  TPWAConstraint* _constr;

  // Private Methods -----------------

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
