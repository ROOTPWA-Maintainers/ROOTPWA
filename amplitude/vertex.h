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
#include <complex>

#include "particle.h"


namespace rpwa {

  class vertex {

  public:
  
    vertex();
    vertex(const vertex& vert);
    virtual ~vertex();
		
    virtual vertex& operator = (const vertex& vert);
    //vertex& operator *= (const lorentzTransform& L);

    virtual void addInParticle (particle* part) { _inParticles.push_back (part); }  ///< adds an incoming particle to vertex
    virtual void addOutParticle(particle* part) { _outParticles.push_back(part); }  ///< adds an outgoing particle to vertex

    virtual const std::vector<particle*>& inParticles () { return _inParticles;  }  ///< returns array of incoming particles
    virtual const std::vector<particle*>& outParticles() { return _outParticles; }  ///< returns array of outgoing particles

    virtual std::complex<double> amplitude() = 0;  ///< returns vertex amplitude
		
    virtual void print(std::ostream& out) const;  ///< prints vertex parameters in human-readable form
    friend std::ostream& operator << (std::ostream&   out,
				      const particle& part);

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


  private:

    std::vector<particle*> _inParticles;    ///< array of pointers to incoming particles
    std::vector<particle*> _outParticles;   ///< array of pointers to outgoing particles

    static bool _debug;  ///< if set to true, debug messages are printed
	
  };


} // namespace rpwa


#endif  // VERTEX_H
