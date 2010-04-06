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
//      class that describes final state vertex decay topology
//      class is just used for internal book keeping
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef FSVERTEX_H
#define FSVERTEX_H


#include <boost/shared_ptr.hpp>

#include "interactionVertex2.h"


namespace rpwa {


  class fsVertex : public interactionVertex2 {

  public:
  
    fsVertex(const particlePtr& fsParticle);  ///< force vertex to have exactly one incoming and no outgoing particles
    fsVertex(const fsVertex& vert);
    virtual ~fsVertex();
		
    inline virtual bool addInParticle (const particlePtr&) { return false; }  ///< disabled; only 1 incoming particle (final state particle) is allowed
    inline virtual bool addOutParticle(const particlePtr&) { return false; }  ///< disabled; no outgoing particles are allowed

    // final-state specific accessors
    inline particlePtr&       fsParticle()       { return inParticles()[0]; }  ///< returns final state particle
    inline const particlePtr& fsParticle() const { return inParticles()[0]; }  ///< returns final state particle

    virtual std::ostream& print(std::ostream& out) const;  ///< prints vertex parameters in human-readable form
    virtual std::ostream& dump (std::ostream& out) const;  ///< prints all vertex data in human-readable form

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


  private:

    static bool _debug;  ///< if set to true, debug messages are printed

  };


  typedef boost::shared_ptr<fsVertex> fsVertexPtr;


  inline
  fsVertexPtr
  createFsVertex(const particlePtr& fsParticle)
  {
    fsVertexPtr v(new fsVertex(fsParticle));
    return v;
  }


  inline
  std::ostream&
  operator <<(std::ostream&   out,
	      const fsVertex& vert)
  {
    return vert.print(out);
  }


}  // namespace rpwa


#endif  // DIFFRACTIVEDISSVERTEX_H
