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


#ifndef ISOBARDECAYVERTEX2_H
#define ISOBARDECAYVERTEX2_H


#include <boost/shared_ptr.hpp>

#include "interactionVertex2.h"


namespace rpwa {

  class isobarDecayVertex2 : public interactionVertex2 {

  public:
  
    isobarDecayVertex2(const particlePtr& mother,
		       const particlePtr& daughter1,
		       const particlePtr& daughter2,
		       const unsigned int L = 0,
		       const unsigned int S = 0);  ///< force vertex to have exactly one incoming (mother) and two outgoing particles (daughters)
    isobarDecayVertex2(const isobarDecayVertex2& vert);
    virtual ~isobarDecayVertex2();
		
    virtual isobarDecayVertex2& operator =(const isobarDecayVertex2& vert);
    virtual isobarDecayVertex2* clone(const bool cloneInParticles  = false,
				      const bool cloneOutParticles = false) const;

    inline virtual bool addInParticle (const particlePtr&) { return false; }  ///< disabled; only 1 incoming particle (mother) is allowed
    inline virtual bool addOutParticle(const particlePtr&) { return false; }  ///< disabled; only 2 outgoing particle (daughters) are allowed

    // isobar decay specific accessors
    inline particlePtr&       mother   ()       { return inParticles ()[0]; }  ///< returns mother particle
    inline particlePtr&       daughter1()       { return outParticles()[0]; }  ///< returns first daughter particle
    inline particlePtr&       daughter2()       { return outParticles()[1]; }  ///< returns second daughter particle
    inline const particlePtr& mother   () const { return inParticles ()[0]; }  ///< returns mother particle
    inline const particlePtr& daughter1() const { return outParticles()[0]; }  ///< returns first daughter particle
    inline const particlePtr& daughter2() const { return outParticles()[1]; }  ///< returns second daughter particle

    const TLorentzVector& calcMotherLzVec();  ///< (re)calculates mother Lorentz-vector from daughter Lorentz-vectors
    
    inline unsigned int L() const { return _L; }  ///< returns the relative orbital angular momentum between the two daughters * 2 (!!!)
    inline unsigned int S() const { return _S; }  ///< returns the total spin of the two daughters * 2 (!!!)

    inline void setL(const unsigned int L) { _L = L; }  ///< sets the relative orbital angular momentum between the two daughters * 2 (!!!)
    inline void setS(const unsigned int S) { _S = S; }  ///< sets the total spin of the two daughters * 2 (!!!)

    bool checkConsistency();  ///< checks quantum decomposition of in-particle to outparticles

    virtual std::ostream& print        (std::ostream& out) const;  ///< prints vertex parameters in human-readable form
    virtual std::ostream& dump         (std::ostream& out) const;  ///< prints all vertex data in human-readable form
    virtual std::ostream& printPointers(std::ostream& out) const;  ///< prints particle pointers strored in vertex

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


  private:

    bool checkMultiplicativeQn(const int          motherQn,
			       const int          daughter1Qn,
			       const int          daughter2Qn,
			       const std::string& qnName = "");  ///< checks consistency of a multiplicative quantum number
    bool checkAdditiveQn      (const int          motherQn,
			       const int          daughter1Qn,
			       const int          daughter2Qn,
			       const std::string& qnName = "");  ///< checks consistency of an additive quantum number
    
    unsigned int _L;  ///< relative orbital angular momentum between the two daughters * 2 (!!!)
    unsigned int _S;  ///< total spin of the two daughters * 2 (!!!)

    static bool _debug;  ///< if set to true, debug messages are printed

  };


  typedef boost::shared_ptr<isobarDecayVertex2> isobarDecayVertexPtr;


  inline
  isobarDecayVertexPtr
  createIsobarDecayVertex(const particlePtr& mother,
			  const particlePtr& daughter1,
			  const particlePtr& daughter2,
			  const unsigned int L = 0,
			  const unsigned int S = 0)
  {
    isobarDecayVertexPtr v(new isobarDecayVertex2(mother, daughter1, daughter2, L, S));
    return v;
  }


  inline
  std::ostream&
  operator <<(std::ostream&             out,
	      const isobarDecayVertex2& vert)
  {
    return vert.print(out);
  }


}  // namespace rpwa


#endif  // ISOBARDECAYVERTEX2_H
