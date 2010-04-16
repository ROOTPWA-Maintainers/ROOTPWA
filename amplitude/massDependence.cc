///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
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
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      functor class hierarchy for mass-dependent part of the amplitude
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "utilities.h"
#include "isobarDecayVertex.h"
#include "isobarHelicityAmplitude.h"
#include "massDependence.h"


using namespace std;
using namespace rpwa;


bool massDependence::_debug = false;


ostream&
massDependence::print(ostream& out) const
{
  out << "mass dependence";
  return out;
}


bool flatMassDependence::_debug = false;


complex<double>
flatMassDependence::amp(const isobarDecayVertex&)
{
  return 1;
}


ostream&
flatMassDependence::print(ostream& out) const
{
  out << "flat mass dependence";
  return out;
}


bool relativisticBreitWigner::_debug = false;


complex<double>
relativisticBreitWigner::amp(const isobarDecayVertex& v)
{
  const particlePtr& mother = v.mother();

  const double M      = mother->lzVec().M();          // mother mass
  const double m1     = v.daughter1()->lzVec().M();   // daughter 1 mass
  const double m2     = v.daughter2()->lzVec().M();   // daughter 2 mass
  const double M0     = mother->mass();               // resonance peak position
  const double Gamma0 = mother->width();              // resonance peak width
  const double q      = breakupMomentum(M,  m1, m2);  // breakup momentum
  const double q0     = breakupMomentum(M0, m1, m2);  // breakup momentum at peak position
  const double L      = v.L();

  return breitWigner(M, M0, Gamma0, L, q, q0);
}


ostream&
relativisticBreitWigner::print(ostream& out) const
{
  out << "relativistic Breit-Wigner";
  return out;
}


// AMP_M::AMP_M() {
//   ves_sheet = 0;
//   _Pmax     = 2;
//   _Nmax     = 5;

//   _rho = matrix<complex<double> >(2, 2);
//   _M   = matrix<complex<double> >(2, 2);
//   _T   = matrix<complex<double> >(2, 2);

//   _f = matrix<complex<double> >(1, 2);
//   _f.el(0, 0) = 0.1968;
//   _f.el(0, 1) = -0.0154;

//   _a = vector<matrix<complex<double> > >(2, matrix<complex<double> >(2, 2));
//   _a[0].el(0, 0) = 0.1131;
//   _a[0].el(0, 1) = 0.0150;
//   _a[0].el(1, 0) = 0.0150;
//   _a[0].el(1, 1) = -0.3216;
//   _a[1].el(0, 0) = _f.el(0, 0) * _f.el(0, 0);
//   _a[1].el(0, 1) = _f.el(0, 0) * _f.el(0, 1);
//   _a[1].el(1, 0) = _f.el(0, 1) * _f.el(0, 0);
//   _a[1].el(1, 1) = _f.el(0, 1) * _f.el(0, 1);

	
//   _c = vector<matrix<complex<double> > >(5, matrix<complex<double> >(2, 2));
//   _c[0].el(0, 0) = 0.0337;
//   _c[1].el(0, 0) = -0.3185;
//   _c[2].el(0, 0) = -0.0942;
//   _c[3].el(0, 0) = -0.5927;
//   _c[4].el(0, 0) = 0.1957; 
//   _c[0].el(0, 1) = _c[0].el(1, 0) = -0.2826;
//   _c[1].el(0, 1) = _c[1].el(1, 0) = 0.0918;
//   _c[2].el(0, 1) = _c[2].el(1, 0) = 0.1669;
//   _c[3].el(0, 1) = _c[3].el(1, 0) = -0.2082;
//   _c[4].el(0, 1) = _c[4].el(1, 0) = -0.1386;
//   _c[0].el(1, 1) = 0.3010;
//   _c[1].el(1, 1) = -0.5140;
//   _c[2].el(1, 1) = 0.1176;
//   _c[3].el(1, 1) = 0.5204;
//   _c[4].el(1, 1) = -0.3977; 

//   _sP = matrix<double>(1, 2);
//   _sP.el(0, 0) = -0.0074;
//   _sP.el(0, 1) = 0.9828;
// }


// complex<double> AMP_M::val(const particle& p) 
// {
//   extern particleDataTable PDGtable;

//   _M.el(0, 0) = 0;
//   _M.el(0, 1) = 0;
//   _M.el(1, 0) = 0;
//   _M.el(0, 1) = 0;

//   complex <double> ret;
//   complex <double> i(0, 1);

//   double pi_mass    = PDGtable.get("pi").Mass();
//   double pi0_mass   = PDGtable.get("pi0").Mass();
//   double K_mass     = PDGtable.get("K").Mass();
//   double K0_mass    = PDGtable.get("K0").Mass();
//   double Kmean_mass = 0.5 * (K_mass + K0_mass);
//   double My_mass    = ~(p.get4P());
//   double s          = My_mass * My_mass;
//   if (fabs(s - _sP.el(0, 1)) < 1e-6) {
//     My_mass += 1e-6;
//     s = My_mass * My_mass;
//   }

//   complex<double> q_pipi   = q(My_mass, pi_mass,    pi_mass);
//   complex<double> q_pi0pi0 = q(My_mass, pi0_mass,   pi0_mass);
//   complex<double> q_KK     = q(My_mass, K_mass,     K_mass);
//   complex<double> q_K0K0   = q(My_mass, K0_mass,    K0_mass);
//   complex<double> q_KmKm   = q(My_mass, Kmean_mass, Kmean_mass);
//   _rho.el(0, 0) = 0.5 * ((2.0 * q_pipi) / My_mass + (2.0 * q_pi0pi0) / My_mass);
//   _rho.el(1, 1) = 0.5 * ((2.0 * q_KK)   / My_mass + (2.0 * q_K0K0)   / My_mass);

//   if (ves_sheet) {
//     if (q_KmKm.imag() > 0.0)
//       q_KmKm *= -1;
//     _rho.el(0, 0) = (2.0 * q_pipi) / My_mass;
//     _rho.el(1, 1) = (2.0 * q_KmKm) / My_mass;
//   }
//   //_rho.print();

//   double scale = (s / (4.0 * Kmean_mass * Kmean_mass)) - 1;
//   //cout << "scale=" << scale << endl;

//   for (int p = 0; p < _Pmax; ++p) {
//     complex<double> fa = 1. / (s - _sP.el(0, p));
//     //cout << "fa" << p <<"="<< fa << endl;
//     _M += fa * _a[p];
//   }
//   for (int n = 0; n < _Nmax; ++n) {
//     complex<double> sc = pow(scale, n);
//     //cout << "sc" << n <<"="<< sc << endl;
//     _M += sc *_c[n];
//   }
	
//   // Modification: Here we set off-diagonal terms to 0
//   _M.el(0, 1) = 0;
//   _M.el(1, 0) = 0;
//   //_M.print();

//   _T = (_M - i * _rho).inv();
//   //_T.print();

//   ret = _T.el(0, 0);
//   return ret;
// }


// complex<double> AMP_ves::val(const particle& p) {
//   complex<double>        bw, amp_m;
//   static complex<double> coupling(-0.3743, 0.3197);
//   static double          massf0 = 0.9837, widthf0 = 0.0376;

//   double pi_mass = PDGtable.get("pi").Mass();
//   double My_mass = ~(p.get4P());

//   ves_sheet = 1;
//   amp_m     = AMP_M::val(p);

//   if (My_mass > 2.0 * pi_mass) {
//     double p, p0, gam, denom, A, B, C;
//     p     = q(My_mass, pi_mass, pi_mass).real();
//     p0    = q(massf0,  pi_mass, pi_mass).real();
//     gam   = widthf0 * (p / p0);
//     A     = massf0 * massf0 - My_mass * My_mass;
//     B     = massf0 * gam;
//     C     = B * (My_mass / p);
//     denom = C / (A * A + B * B);
//     bw    = denom * complex<double>(A, B);
//   }

//   return amp_m - coupling * bw;
// }


// AMP_kach::AMP_kach(): AMP_M()
// {
//   // Change parameters according to Kachaev's prescription
//   _c[4].el(0, 0) = 0; // was 0.1957;
//   _c[4].el(1, 1) = 0; // was -0.3977;

//   _a[0].el(0, 1) = 0; // setting off-diagonal values to 0
//   _a[0].el(1, 0) = 0; // 
 
//   // _a[1] are the f's from the AMP paper! Setting to 0!
//   _a[1].el(0, 0) = 0;
//   _a[1].el(0, 1) = 0;
//   _a[1].el(1, 0) = 0;
//   _a[1].el(1, 1) = 0;
// }
