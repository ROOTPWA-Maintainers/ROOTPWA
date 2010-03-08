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
//      virtual base class that desbribes general interaction vertex between particles
//      only pointers to partciles are stored
//      vertex does not "own" its particles; creation and destruction
//      of particles has to be done in calling code
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef VERTEX_H
#define VERTEX_H


#include <iostream>
#include <vector>

#include "particle.h"


namespace rpwa {

  class vertex {

  public:
  
    vertex();
    vertex(const vertex& vert);
    virtual ~vertex();
		
    virtual vertex& operator = (const vertex& vert);
    //vertex& operator *= (const lorentzTransform& L);

    virtual void addInParticle (particle* part) { _inParticles.push_back (part); }
    virtual void addOutParticle(particle* part) { _outParticles.push_back(part); }

    virtual const std::vector<particle*>&       inParticles ()       { return _inParticles;  }
    virtual const std::vector<const particle*>& inParticles () const { return _inParticles;  }
    virtual const std::vector<particle*>&       outParticles()       { return _outParticles; }
    virtual const std::vector<const particle*>& outParticles() const { return _outParticles; }

    virtual std::complex<double> amplitude() = 0;
		
    virtual void print(std::ostream& out) const;
    friend std::ostream& operator << (std::ostream&   out,
				      const particle& part);

    static bool debug() const { return _debug; }
    static void setDebug(const bool debug = true) { _debug = debug; }


  private:

    std::vector<particle*> _inParticles;    ///< array of pointers to incoming particles
    std::vector<particle*> _outParticles;   ///< array of pointers to outgoing particles

    static bool _debug;  ///< if set to true, debug messages are printed
	
  };


} // namespace rpwa


#endif  // VERTEX_H
