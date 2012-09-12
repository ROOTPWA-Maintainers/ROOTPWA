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
//      incoming parent and two outgoing daughter particle
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
#include "massDependence.h"


namespace rpwa {


	class isobarDecayVertex;
	typedef boost::shared_ptr<isobarDecayVertex> isobarDecayVertexPtr;


	class isobarDecayVertex : public interactionVertex {

	public:
  
		isobarDecayVertex(const particlePtr&       parent,
		                  const particlePtr&       daughter1,
		                  const particlePtr&       daughter2,
		                  const unsigned int       L       = 0,
		                  const unsigned int       S       = 0,
		                  const massDependencePtr& massDep = massDependencePtr());  ///< force vertex to have exactly one incoming (parent) and two outgoing particles (daughters)
		isobarDecayVertex(const isobarDecayVertex& vert);
		virtual ~isobarDecayVertex();
		
		isobarDecayVertex& operator =(const isobarDecayVertex& vert);
		isobarDecayVertexPtr clone(const bool cloneInParticles  = false,
		                           const bool cloneOutParticles = false) const  ///< creates copy of isobar decay vertex; must not be virtual
		{ return isobarDecayVertexPtr(doClone(cloneInParticles, cloneOutParticles)); }

		// comparison operators that check equality of all fields
		friend bool operator ==(const isobarDecayVertex& lhsVert,
		                        const isobarDecayVertex& rhsVert) { return lhsVert.isEqualTo(rhsVert); }
		friend bool operator !=(const isobarDecayVertex& lhsVert,
		                        const isobarDecayVertex& rhsVert) { return not(lhsVert == rhsVert); }

		virtual bool addInParticle (const particlePtr&);  ///< disabled; parent particle has to be specified at construction
		virtual bool addOutParticle(const particlePtr&);  ///< disabled; daughter particles have to be specified at construction

		// isobar decay specific accessors
		inline particlePtr&       parent   ()       { return inParticles ()[0]; }  ///< returns parent particle
		inline particlePtr&       daughter1()       { return outParticles()[0]; }  ///< returns first daughter particle
		inline particlePtr&       daughter2()       { return outParticles()[1]; }  ///< returns second daughter particle
		inline const particlePtr& parent   () const { return inParticles ()[0]; }  ///< returns parent particle
		inline const particlePtr& daughter1() const { return outParticles()[0]; }  ///< returns first daughter particle
		inline const particlePtr& daughter2() const { return outParticles()[1]; }  ///< returns second daughter particle

		const TLorentzVector& calcParentLzVec();  ///< (re)calculates parent Lorentz-vector from daughter Lorentz-vectors

		int calcParentCharge   ();  ///< sets parent charge to sum of daughter charges
		int calcParentBaryonNmb();  ///< sets parent baryon number to sum of daughter baryon numbers
    
		inline unsigned int L() const { return _L; }  ///< returns the relative orbital angular momentum between the two daughters * 2 (!!!)
		inline unsigned int S() const { return _S; }  ///< returns the total spin of the two daughters * 2 (!!!)

		inline void setL(const unsigned int L) { _L = L; }  ///< sets the relative orbital angular momentum between the two daughters * 2 (!!!)
		inline void setS(const unsigned int S) { _S = S; }  ///< sets the total spin of the two daughters * 2 (!!!)

		inline std::complex<double>     massDepAmplitude() const { return _massDep->amp(*this); }  ///< returns mass-dependent amplitude
		inline const massDependencePtr& massDependence  () const { return _massDep;             }  ///< returns mass-dependence
		inline void setMassDependence(const massDependencePtr& massDep) { _massDep = massDep; }    ///< sets mass dependence

		bool checkConsistency();  ///< checks quantum decomposition of in-particle to outparticles

		virtual std::ostream& print        (std::ostream& out) const;  ///< prints vertex parameters in human-readable form
		virtual std::ostream& dump         (std::ostream& out) const;  ///< prints all vertex data in human-readable form
		virtual std::ostream& printPointers(std::ostream& out) const;  ///< prints particle pointers strored in vertex

		virtual std::string name() const { return "isobarDecayVertex"; }  ///< returns label used in graph visualization, reporting, and key file

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	protected:

		virtual isobarDecayVertex* doClone(const bool cloneInParticles,
		                                   const bool cloneOutParticles) const;  ///< helper function to use covariant return types with smart pointers; needed for public clone()

		bool isEqualTo(const isobarDecayVertex& vert) const;  ///< returns whether vertex is equal to this by checking equality of all member variables

	private:

		bool checkMultiplicativeQn(const int          parentQn,
		                           const int          daughter1Qn,
		                           const int          daughter2Qn,
		                           const std::string& qnName = "");  ///< checks consistency of a multiplicative quantum number
		bool checkAdditiveQn      (const int          parentQn,
		                           const int          daughter1Qn,
		                           const int          daughter2Qn,
		                           const std::string& qnName = "");  ///< checks consistency of an additive quantum number
    
		unsigned int _L;  ///< relative orbital angular momentum between the two daughters * 2 (!!!)
		unsigned int _S;  ///< total spin of the two daughters * 2 (!!!)

		massDependencePtr _massDep;  ///< functor with mass dependence

		static bool _debug;  ///< if set to true, debug messages are printed

	};


	inline
	isobarDecayVertexPtr
	createIsobarDecayVertex(const particlePtr&       parent,
	                        const particlePtr&       daughter1,
	                        const particlePtr&       daughter2,
	                        const unsigned int       L       = 0,
	                        const unsigned int       S       = 0,
	                        const massDependencePtr& massDep = massDependencePtr())
	{
		isobarDecayVertexPtr vert(new isobarDecayVertex(parent, daughter1, daughter2, L, S, massDep));
		return vert;
	}


	inline
	std::ostream&
	operator <<(std::ostream&            out,
	            const isobarDecayVertex& vert)
	{
		return vert.print(out);
	}


}  // namespace rpwa


#endif  // ISOBARDECAYVERTEX_H
