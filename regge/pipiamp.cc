//-----------------------------------------------------------
//
// Description:
//      Implementation of class reggeprop
//      see pipiamp.h for details
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


// This Class' Header ------------------
#include "pipiamp.h"


// C/C++ Headers ----------------------
#include "TMath.h"
#include <iostream>

using namespace std;



// Collaborating Class Headers --------
#include <boost/math/special_functions/legendre.hpp>

// Class Member definitions -----------

pipiamp::pipiamp() : _lmax(4),
		     _mpic2(0.01947983515),
		     _mpi02(0.01821868255)
{;}


 std::complex<double>
 pipiamp::amp(std::complex<double> s,
	      std::complex<double> t,
	      unsigned int I) const {
   unsigned int lmin=1;
   if(I%2==0)lmin=0;

   // charged pion-pion scattering incoming 3-momentum magnitude
   // for kinematics see Martin/Spearman 4.3.1
   complex<double> q2=0.25*(s-4.0*_mpic2);

   std::complex<double> costheta=1.0 + 0.5*t/q2;

   std::complex<double> result(0,0);
   // loop over partial waves
   for(unsigned int l=lmin;lmin<_lmax;lmin+=2){
     result += (2.*l + 1.) * boost::math::legendre_p(l,costheta.real())*f(s,l,I);
   }

   complex<double> sqrts=sqrt(s);
   complex<double> k=kpipi(s);
   return 4.*sqrts/k*0.318398862 * result;
 }


// partial wave amplitudes in terms of phase shifts and elasitcity

 std::complex<double>
 pipiamp::f(std::complex<double> s,
	    unsigned  int l,
	    unsigned int I)const
 {

   double e=eta(s,l,I);
   std::complex<double> result=e/complex<double>(cotDelta(s,l,I),-1);
   result+= 0.5*complex<double>(0,1-e);
   return result;
}

// phase-shifts
double
pipiamp::cotDelta(std::complex<double> s,
		  unsigned int l,
		  unsigned int I) const
 {
   return 0;
 }


// elasticities
double
pipiamp::eta(std::complex<double> s,
		  unsigned int l,
	       unsigned int I) const
{
  return 1;
}
