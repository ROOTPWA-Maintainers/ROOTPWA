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

#include "interactionVertex.h"


namespace rpwa {


	class fsVertex;
	typedef boost::shared_ptr<fsVertex> fsVertexPtr;


	class fsVertex : public interactionVertex {

	public:

		fsVertex(const particlePtr& fsParticle);  ///< force vertex to have exactly one incoming and no outgoing particles
		fsVertex(const fsVertex& vert);
		virtual ~fsVertex();

		fsVertexPtr clone(const bool cloneInParticles  = false,
		                  const bool                   = false) const  ///< creates deep copy of final state vertex; must not be virtual
		{ return fsVertexPtr(doClone(cloneInParticles, false)); }

		virtual bool addInParticle (const particlePtr&);  ///< disabled; final state particle has to be specified at construction
		virtual bool addOutParticle(const particlePtr&);  ///< disabled; no outgoing particles are allowed

		// final-state specific accessors
		inline particlePtr&       fsParticle()       { return inParticles()[0]; }  ///< returns final state particle
		inline const particlePtr& fsParticle() const { return inParticles()[0]; }  ///< returns final state particle

		virtual std::ostream& print        (std::ostream& out) const;  ///< prints vertex parameters in human-readable form
		virtual std::ostream& dump         (std::ostream& out) const;  ///< prints all vertex data in human-readable form
		virtual std::ostream& printPointers(std::ostream& out) const;  ///< prints particle pointers strored in vertex

		virtual std::string name() const { return "fsVertex"; }  ///< returns label used in graph visualization, reporting, and key file

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	protected:

		virtual fsVertex* doClone(const bool cloneInParticles,
		                          const bool cloneOutParticles) const;  ///< helper function to use covariant return types with smart pointers; needed for public clone()


	private:

		static bool _debug;  ///< if set to true, debug messages are printed

	};


	inline
	fsVertexPtr
	createFsVertex(const particlePtr& fsParticle)
	{
		fsVertexPtr vert(new fsVertex(fsParticle));
		return vert;
	}


	inline
	std::ostream&
	operator <<(std::ostream&   out,
	            const fsVertex& vert)
	{
		return vert.print(out);
	}


}  // namespace rpwa


#endif  // FSVERTEX_H
