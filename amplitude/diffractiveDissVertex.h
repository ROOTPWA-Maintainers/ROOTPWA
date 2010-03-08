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
//      class that describes production vertex in diffractive dissociation
//      beam-Reggeon-(X-system) vertex has exactly one incoming beam and
//      one outgoing X particle, which unambiguously defines the Reggeon
//      kinematics
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef DIFFRACTIVEDISSVERTEX_H
#define DIFFRACTIVEDISSVERTEX_H


#include "vertex.h"


namespace rpwa {

  class diffractiveDissVertex : public vertex {

  public:
  
    diffractiveDissVertex(particle& beam,
			  particle& XSystem);  ///< vertex makes sense only if both beam and X-system are specified
    diffractiveDissVertex(const diffractiveDissVertex& vert);
    virtual ~diffractiveDissVertex();
		
    virtual bool addInParticle (particle&) { return false; }  ///< disabled; only 1 incoming particle (beam) is allowed
    virtual bool addOutParticle(particle&) { return false; }  ///< disabled; only 1 outgoing particle (X-system) is allowed

    // diffractive dissociation specific accessors
    particle& beam()    { return *(inParticles() [0]); }  ///< returns beam particle
    particle& XSystem() { return *(outParticles()[0]); }  ///< returns X particle

    virtual std::ostream& print(std::ostream& out) const;  ///< prints vertex parameters in human-readable form

  };


  inline
  std::ostream&
  operator << (std::ostream&                out,
	       const diffractiveDissVertex& vert) { return vert.print(out); }


} // namespace rpwa


#endif  // DIFFRACTIVEDISSVERTEX_H
