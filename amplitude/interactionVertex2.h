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
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef INTERACTIONVERTEX2_H
#define INTERACTIONVERTEX2_H


#include <iostream>
#include <vector>
#include <complex>

#include <boost/shared_ptr.hpp>

#include "TLorentzRotation.h"

#include "particle.h"


namespace rpwa {

  class interactionVertex2 {

  public:
  
    interactionVertex2();
    interactionVertex2(const interactionVertex2& vert);
    virtual ~interactionVertex2();
		
    virtual interactionVertex2& operator =(const interactionVertex2& vert);
    virtual interactionVertex2* clone(const bool cloneInParticles  = false,
				      const bool cloneOutParticles = false) const;
    virtual void clear();

    virtual bool addInParticle (const particlePtr& part);  ///< adds an incoming particle to vertex
    virtual bool addOutParticle(const particlePtr& part);  ///< adds an outgoing particle to vertex

    void transformOutParticles(const TLorentzRotation& L);  ///< applies Lorentz-transformation to outgoing particles

    inline unsigned int nmbInParticles () const { return _inParticles.size();  }  ///< returns number of incoming particles
    inline unsigned int nmbOutParticles() const { return _outParticles.size(); }  ///< returns number of outgoing particles

    inline std::vector<particlePtr>&       inParticles ()       { return _inParticles;  }  ///< returns array of incoming particles
    inline std::vector<particlePtr>&       outParticles()       { return _outParticles; }  ///< returns array of outgoing particles
    inline const std::vector<particlePtr>& inParticles () const { return _inParticles;  }  ///< returns array of incoming particles
    inline const std::vector<particlePtr>& outParticles() const { return _outParticles; }  ///< returns array of outgoing particles

    virtual std::ostream& print(std::ostream& out) const;  ///< prints vertex parameters in human-readable form
    virtual std::ostream& dump (std::ostream& out) const;  ///< prints all vertex data in human-readable form

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


  protected:

    void cloneInParticles ();  ///< clones all incoming particles
    void cloneOutParticles();  ///< clones all outgoing particles

    std::vector<particlePtr> _inParticles;   ///< array of pointers to incoming particles
    std::vector<particlePtr> _outParticles;  ///< array of pointers to outgoing particles


  private:

    static bool _debug;  ///< if set to true, debug messages are printed
	
  };


  typedef boost::shared_ptr<interactionVertex2> interactionVertexPtr;


  inline
  interactionVertexPtr
  createInteractionVertex()
  {
    interactionVertexPtr v(new interactionVertex2());
    return v;
  }


  inline
  std::ostream&
  operator <<(std::ostream&             out,
	      const interactionVertex2& vert)
  {
    return vert.print(out);
  }


}  // namespace rpwa


#endif  // INTERACTIONVERTEX_H
