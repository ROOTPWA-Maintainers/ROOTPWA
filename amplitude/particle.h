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
//      container class for all particle related external information
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef PARTICLE_H
#define PARTICLE_H


#include <string>

#include "TVector3.h"
#include "TLorentzVector.h"

#include "particleProperties.h"


namespace rpwa {	

  class particle : public particleProperties {
	
  public:
			
    particle();
    particle(const particle&           part);
    particle(const int                 charge,
	     const TVector3&           momentum,
	     const particleProperties& partProp);
    particle(const std::string&        partName,
	     const TVector3&           momentum);
    virtual ~particle();

    particle& operator = (const particle& part);

    int            charge() const { return _charge; }  ///< returns particle's charge
    TLorentzVector lzVec()  const { return _lzVec;  }  ///< returns Lorentz vector of particle

    void setCharge  (const int       charge)   { _charge = charge;                                                            }  ///< sets particle's charge
    void setMomentum(const TVector3& momentum) { _lzVec  = TLorentzVector(momentum, sqrt(momentum.Mag2() + mass() * mass())); }  ///< sets particle's Lorentz vector

    void print(std::ostream& out) const;  ///< prints particle data in human-readable form
    friend std::ostream& operator << (std::ostream&   out,
				      const particle& part);

    static bool debug() const { return _debug; }
    static void setDebug(const bool debug = true) { _debug = debug; }


  private:
			
    int            _charge;  ///< charge
    TLorentzVector _lzVec;   ///< Lorentz vector

    static bool _debug;  ///< if set to true, debug messages are printed

  };


} // namespace rpwa


#endif  // PARTICLE_H
