//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Pion-Pion scattering amplitude
//
//      Literature:
//      Pelaez, Yndurain              Phys.Rev.D71 (2005) 074016
//      Kaminski, Pelaez, Yndurain    Phys.Rev.D74 (2006) 014001
//      Kaminski, Pelaez, Yndurain    Phys.Rev.D77 (2008) 054015
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

class pipiamp {
public:
 // Constructors/Destructors ---------
  pipiamp();
  ~pipiamp(){;}

 // Operators

 // Accessors -----------------------
 // Modifiers -----------------------
 // Operations ----------------------

  
  /// PionPion scattering amplitude
  /// amp = F^i(s,t)= 2 \sum_l (2l+2)P_l(cosTheta)f^I_l(s)
  std::complex<double> amp(std::complex<double> s, 
			   std::complex<double> t, 
			   unsigned int I) const;

  std::complex<double> amp(double s, 
			   double t, 
			   double cosTheta,
			   unsigned int I) const 
    {return amp(std::complex<double>(s,0), 
		std::complex<double>(t,0), 
		I);}

  /// Partial wave amplitudes
  /// f^I_l(s)=2sqrt(s)/k[(eta_l exp(2i\delta_l)-1)/(2i)]
  /// k=sqrt(s/4-mpi^2)
  std::complex<double> f(std::complex<double> s,
			 unsigned  int l=0,
			 unsigned int I=0)const;

  std::complex<double> f(double s,
			 unsigned  int l=0,
			 unsigned int I=0) const {return f(std::complex<double>(s,0),l,I);}
  
  /// Phase shifts
  double cotDelta(std::complex<double> s,
		unsigned int l=0,
		unsigned int I=0) const;


  double cotDelta(double s,
		unsigned int l=0,
		  unsigned int I=0) const  {return cotDelta(s,l,I);}


  double eta(std::complex<double> s,
		  unsigned int l=0,
		  unsigned int I=0) const;


  double eta(double s,
	     unsigned int l=0,
	     unsigned int I=0) const {return eta(s,l,I);}


private:

  // Private Data Members ------------
  unsigned int _lmax;   /// parameterizations are available up to this
  double _mpic2;         /// charged pion mass^2
  double _mpi02;         /// neutral pion mass^2

  // Private Methods -----------------
  std::complex<double> kpipi(const std::complex<double>& s) const {return sqrt(s*0.25-_mpic2);} /// breakup momentum

};

#endif
