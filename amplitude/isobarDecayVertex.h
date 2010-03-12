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
//      class that describes decay vertex of isobar into two particles
//      the isobar -> particle1 + particle 2 vertex has exactly one
//      incoming mother and two outgoing daughter particle
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef ISOBARDECAYVERTEX_H
#define ISOBARDECAYVERTEX_H


#include "interactionVertex.h"


namespace rpwa {

  class isobarDecayVertex : public interactionVertex {

  public:
  
    isobarDecayVertex(particle&          mother,
		      particle&          daughter1,
		      particle&          daughter2,
		      const unsigned int L = 0,
		      const unsigned int S = 0);  ///< vertex makes sense only mother and the two daughter particles are specified
    isobarDecayVertex(const isobarDecayVertex& vert);
    virtual ~isobarDecayVertex();
		
    virtual isobarDecayVertex& operator = (const isobarDecayVertex& vert);

    virtual bool addInParticle (particle&) { return false; }  ///< disabled; only 1 incoming particle (mother) is allowed
    virtual bool addOutParticle(particle&) { return false; }  ///< disabled; only 2 outgoing particle (daughters) are allowed

    // isobar decay specific accessors
    particle& mother()    { return *(inParticles() [0]); }  ///< returns mother particle
    particle& daughter1() { return *(outParticles()[0]); }  ///< returns first daughter particle
    particle& daughter2() { return *(outParticles()[1]); }  ///< returns second daughter particle

    const TLorentzVector& updateMotherLzVec();  ///< sets mother Lorentz-vector to sum of daughter Lorentz-vetcors

    
    unsigned int L() const { return _L; }  ///< returns the relative orbital angular momentum between the two daughters
    unsigned int S() const { return _S; }  ///< returns the total spin of the two daughters * 2 (!!!)

    void setL(const unsigned int L) { _L = L; }  ///< sets the relative orbital angular momentum between the two daughters
    void setS(const unsigned int S) { _S = S; }  ///< sets the total spin of the two daughters * 2 (!!!)

    bool checkConsistency(); ///< checks quantum decomposition of in-particle to outparticles

    virtual std::ostream& print(std::ostream& out) const;  ///< prints vertex parameters in human-readable form


  private:
    
    unsigned int _L;  ///< relative orbital angular momentum between the two daughters
    unsigned int _S;  ///< total spin of the two daughters * 2 (!!!)

  };


  inline
  std::ostream&
  operator << (std::ostream&            out,
	       const isobarDecayVertex& vert) { return vert.print(out); }


} // namespace rpwa


#endif  // ISOBARDECAYVERTEX_H
