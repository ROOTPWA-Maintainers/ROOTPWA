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
//      base class that desbribes general interaction vertex between particles
//      !NOTE! class stores pointers to particles, which are inserted via
//             references in order to ensure existence of objects; the calling
//             code has to ensure that lifetime of the particle instances is
//             longer than life time of the vertex instances they are assigned to
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef INTERACTIONVERTEX_H
#define INTERACTIONVERTEX_H


#include <iostream>
#include <vector>
#include <complex>

#include "particle.h"


namespace rpwa {

  class interactionVertex {

  public:
  
    interactionVertex();
    interactionVertex(const interactionVertex& vert);
    virtual ~interactionVertex();
		
    virtual interactionVertex& operator = (const interactionVertex& vert);
    //interactionVertex& operator *= (const lorentzTransform& L);
    virtual interactionVertex& clone() const;

    virtual bool addInParticle (particle& part);  ///< adds an incoming particle to vertex
    virtual bool addOutParticle(particle& part);  ///< adds an outgoing particle to vertex

    virtual std::vector<particle*>& inParticles () { return _inParticles;  }  ///< returns array of incoming particles
    virtual std::vector<particle*>& outParticles() { return _outParticles; }  ///< returns array of outgoing particles
    virtual const std::vector<particle*>& inParticles () const { return _inParticles;  }  ///< returns array of incoming particles
    virtual const std::vector<particle*>& outParticles() const { return _outParticles; }  ///< returns array of outgoing particles

    virtual unsigned int nmbInParticles () const { return _inParticles.size();  }  ///< returns number of incoming particles
    virtual unsigned int nmbOutParticles() const { return _outParticles.size(); }  ///< returns number of outgoing particles

    virtual bool dataAreValid() const { return _dataValid; }  ///< indicates whether vertex data are complete and valid
    virtual void setDataValid(const bool dataValid = true) { _dataValid = dataValid; }  ///< set internal flag that indicates that external data are complete and valid

    virtual std::ostream& print(std::ostream& out) const;  ///< prints vertex parameters in human-readable form

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


  protected:

    std::vector<particle*> _inParticles;   ///< array of pointers to incoming particles
    std::vector<particle*> _outParticles;  ///< array of pointers to outgoing particles
    bool                   _dataValid;     ///< indicates whether vertex data are complete and valid


  private:

    static bool _debug;  ///< if set to true, debug messages are printed
	
  };


  inline
  std::ostream&
  operator << (std::ostream&            out,
	       const interactionVertex& vert) { return vert.print(out); }


} // namespace rpwa


#endif  // INTERACTIONVERTEX_H
