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


#ifndef INTERACTIONVERTEX_H
#define INTERACTIONVERTEX_H


#include <iostream>
#include <vector>
#include <complex>

#include <boost/weak_ptr.hpp>

#include "TLorentzRotation.h"

#include "particle.h"


class TClonesArray;


namespace rpwa {

	class decayTopology;
	typedef boost::shared_ptr<decayTopology> decayTopologyPtr;
  
	class interactionVertex;
	typedef boost::shared_ptr<interactionVertex> interactionVertexPtr;


	class interactionVertex {

	public:
  
		interactionVertex();
		interactionVertex(const interactionVertex& vert);
		virtual ~interactionVertex();
		
		interactionVertex& operator =(const interactionVertex& vert);
		interactionVertexPtr clone(const bool cloneInParticles  = false,
		                           const bool cloneOutParticles = false) const  ///< creates deep copy of interaction vertex; must not be virtual
		{ return interactionVertexPtr(doClone(cloneInParticles, cloneOutParticles)); }
		virtual void clear();

		// comparison operators that check equality of all member variables
		bool operator ==(const interactionVertex& rhsVert) { return this->isEqualTo(rhsVert); }
		bool operator !=(const interactionVertex& rhsVert) { return not (*this == rhsVert);   }

		virtual bool addInParticle (const particlePtr& part);  ///< adds an incoming particle to vertex
		virtual bool addOutParticle(const particlePtr& part);  ///< adds an outgoing particle to vertex

		void transformOutParticles(const std::vector<TLorentzRotation>& L);  ///< applies Lorentz-transformation to outgoing particles

		inline unsigned int nmbInParticles () const { return _inParticles.size();  }  ///< returns number of incoming particles
		inline unsigned int nmbOutParticles() const { return _outParticles.size(); }  ///< returns number of outgoing particles

		virtual inline std::vector<particlePtr>&       inParticles ()       { return _inParticles;  }  ///< returns array of incoming particles
		virtual inline std::vector<particlePtr>&       outParticles()       { return _outParticles; }  ///< returns array of outgoing particles
		virtual inline const std::vector<particlePtr>& inParticles () const { return _inParticles;  }  ///< returns array of incoming particles
		virtual inline const std::vector<particlePtr>& outParticles() const { return _outParticles; }  ///< returns array of outgoing particles

		virtual const decayTopologyPtr decay() const { return _decay.lock(); } ///< Caution: this might return an empty shared_ptr, as _decay is only set when a topology is added to an isobarAmplitude
		virtual void setDecay(const decayTopologyPtr& decay) { _decay = boost::weak_ptr<decayTopology>(decay); }

		virtual std::ostream& print        (std::ostream& out) const;  ///< prints vertex parameters in human-readable form
		virtual std::ostream& dump         (std::ostream& out) const;  ///< prints all vertex data in human-readable form
		virtual std::ostream& printPointers(std::ostream& out) const;  ///< prints particle pointers strored in vertex

		virtual std::string name() const { return "interactionVertex"; }  ///< returns label used in graph visualization, reporting, and key file

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	protected:

		virtual interactionVertex* doClone(const bool cloneInParticles,
		                                   const bool cloneOutParticles) const;  ///< helper function to use covariant return types with smart pointers; needed for public clone()

		virtual bool isEqualTo(const interactionVertex& vert) const;  ///< returns whether vert is equal to this by checking equality of all member variables

		void cloneInParticles ();  ///< clones all incoming particles
		void cloneOutParticles();  ///< clones all outgoing particles

		std::vector<particlePtr> _inParticles;   ///< array of pointers to incoming particles
		std::vector<particlePtr> _outParticles;  ///< array of pointers to outgoing particles

		boost::weak_ptr<decayTopology> _decay; ///< pointer to the decay the vertex is a "member" of, has to be a weak_ptr to avoid cycles (decayTopology -> interactionVertexPtr -> decayTopologyPtr).


	private:

		static bool _debug;  ///< if set to true, debug messages are printed
	
	};


	inline
	interactionVertexPtr
	createInteractionVertex()
	{
		interactionVertexPtr vert(new interactionVertex());
		return vert;
	}


	inline
	std::ostream&
	operator <<(std::ostream&            out,
	            const interactionVertex& vert)
	{
		return vert.print(out);
	}


}  // namespace rpwa


#endif  // INTERACTIONVERTEX_H
