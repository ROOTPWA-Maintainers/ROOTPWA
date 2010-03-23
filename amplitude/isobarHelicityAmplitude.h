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
//      prototype implementation for general isobar decay amplitude in
//      the helicity framework
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef ISOBARHELICITYAMPLITUDE_H
#define ISOBARHELICITYAMPLITUDE_H


#include <complex>

#include "pputil.h"

#include "isobarDecayTopology.h"


namespace rpwa {	

  class isobarHelicityAmplitude {
	
  public:
			
    isobarHelicityAmplitude();
    isobarHelicityAmplitude(isobarDecayTopology& decay);
    virtual ~isobarHelicityAmplitude();

    void setDecayTopology(isobarDecayTopology& decay);

    std::complex<double> amplitude();                            ///< computes amplitude
    std::complex<double> operator () () { return amplitude(); }  ///< computes amplitude

    std::ostream& print(std::ostream& out) const { return out; }  ///< prints amplitude parameters in human-readable form

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


    static TLorentzRotation hfTransform(const TLorentzVector& daughterLv);  ///< constructs Lorentz-transformation to helicity RF of daughter particle

    static TLorentzRotation gjTransform(const TLorentzVector& beamLv,
					const TLorentzVector& XLv);  ///< constructs Lorentz-transformation to X Gottfried-Jackson frame

    void transformDaughters();  ///< boosts Lorentz-vectors of decay daughters into frames where angular distributions are defined

    std::complex<double> twoBodyDecayAmplitude(const isobarDecayVertex& vertex,
					       const bool               topVertex) const;  ///< calculates amplitude for two-body decay a -> b + c; where b and c are stable
    
    std::complex<double> twoBodyDecayAmplitudeSum(isobarDecayVertex& vertex,
						  const bool         topVertex = false);  ///< recursively sums up decay amplitudes for all allowed helicitities for all vertices below given vertex
    

  private:

    isobarDecayTopology* _decay;  ///< isobar decay topology with all external information
    
    static bool _debug;  ///< if set to true, debug messages are printed
    
  };
  

  inline
  std::ostream&
  operator << (std::ostream&                  out,
	       const isobarHelicityAmplitude& amp) { return amp.print(out); }


  // some wrappers for libpp functions
  // !NOTE! all angular momenta and spin projections are in units of hbar/2
  inline
  double
  normFactor(const int J) { return sqrt(J + 1); }


  inline
  std::complex<double>
  DFuncConj(const int    J,
	    const int    M,
	    const int    lambda,
	    const double phi,
	    const double theta) { return conj(D(phi, theta, 0, J, M, lambda)); }


  inline
  double
  cgCoeff(const int J1,
	  const int M1,
	  const int J2,
	  const int M2,
	  const int J,
	  const int M) { return clebsch(J1, J2, J, M1, M2, M); }

  
  inline
  double barrierFactor(const int    L,
		       const double breakupMom) { return F(L, breakupMom); }


  inline
  std::complex<double> breitWigner(const double m,
				   const double m0,
				   const double Gamma0,
				   const int    L,
				   const double q,
				   const double q0)
  {
    const double Gamma  = Gamma0 * (m0 / m) * (q / q0)
                          * (pow(F(L, q), 2) / pow(F(L, q0), 2));
    return (m0 * Gamma0) / (m0 * m0 - m * m - imag * m0 * Gamma);
  }


} // namespace rpwa


#endif  // ISOBARHELICITYAMPLITUDE_H
