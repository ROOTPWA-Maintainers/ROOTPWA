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

#include <boost/shared_ptr.hpp>

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "particleProperties.h"


namespace rpwa {	

  class particle : public particleProperties {
	
  public:
			
    particle();
    particle(const particle&           part);
    particle(const particleProperties& partProp,
	     const int                 charge,
	     const TVector3&           momentum = TVector3(),
	     const int                 spinProj = 0);
    particle(const std::string&        partName,
	     const TVector3&           momentum = TVector3(),
	     const int                 spinProj = 0);
    particle(const std::string&        partName,
	     const int                 isospin,
	     const int                 G,
	     const int                 J,
	     const int                 P,
	     const int                 C,
	     const int                 spinProj);
    virtual ~particle();

    virtual particle& operator =(const particle& part);
    virtual particle& clone() const;

    std::string           name()     const;                       ///< returns particle name including charge
    std::string           summary()  const;                       ///< returns particle summary in short form
    int                   charge()   const { return _charge;   }  ///< returns particle's charge
    int                   spinProj() const { return _spinProj; }  ///< returns particle's spin projection quantum number
    const TLorentzVector& lzVec()    const { return _lzVec;    }  ///< returns Lorentz vector of particle

    void setCharge  (const int             charge)   { _charge   = charge;                                                            }  ///< sets particle's charge
    void setSpinProj(const int             spinProj) { _spinProj = spinProj;                                                          }  ///< sets particle's spin projection quantum number
    void setMomentum(const TVector3&       momentum) { _lzVec    = TLorentzVector(momentum, sqrt(momentum.Mag2() + mass() * mass())); }  ///< sets particle's Lorentz vector
    void setLzVec   (const TLorentzVector& lzVec)    { _lzVec    = lzVec;                                                             }  ///< sets particle's Lorentz vector; if this is used to inject external data the mass values likely become inconsistent

    const TLorentzVector& transform(const TLorentzRotation& L) { return (_lzVec *= L); }  ///< applies Lorentz-transformation to particle

    virtual std::ostream& print(std::ostream& out) const;  ///< prints particle parameters in human-readable form

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag
			

  private:

    int            _charge;    ///< charge
    int            _spinProj;  ///< spin projection quantum number; can be either M or helicity
    TLorentzVector _lzVec;     ///< Lorentz vector [GeV]

    static bool _debug;  ///< if set to true, debug messages are printed

  };


  typedef boost::shared_ptr<particle> particlePtr;


  inline
  particlePtr
  createParticle(const std::string& partName,
		 const TVector3&    momentum = TVector3(),
		 const int          spinProj = 0)
  {
    particlePtr p(new particle(partName, momentum, spinProj));
    return p;
  }


  inline
  particlePtr
  createParticle(const std::string& partName,
		 const int          isospin,
		 const int          G,
		 const int          J,
		 const int          P,
		 const int          C,
		 const int          spinProj)
  {
    particlePtr p(new particle(partName, isospin, G, J, P, C, spinProj));
    return p;
  }


  inline
  std::ostream&
  operator <<(std::ostream&   out,
	      const particle& part)
  { return part.print(out); }


}  // namespace rpwa


#endif  // PARTICLE_H
