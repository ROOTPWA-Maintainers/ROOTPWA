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
//      class that describes production vertex in diffractive dissociation
//      beam-Reggeon-(X-system) vertex has exactly one incoming beam and
//      one outgoing X particle, which unambiguously defines the Reggeon
//      kinematics
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef DIFFRACTIVEDISSVERTEX_H
#define DIFFRACTIVEDISSVERTEX_H


#include <boost/shared_ptr.hpp>

#include "productionVertex.h"


class TClonesArray;


namespace rpwa {


	class diffractiveDissVertex;
	typedef boost::shared_ptr<diffractiveDissVertex> diffractiveDissVertexPtr;


	class diffractiveDissVertex : public productionVertex {

	public:
  
		diffractiveDissVertex(const particlePtr& beam,
		                      const particlePtr& target,
		                      const particlePtr& XParticle,
		                      const particlePtr& recoil = particlePtr());  ///< force vertex to have two incoming (beam + target) and two outgoing particles (X + recoil); recoil is optional
		diffractiveDissVertex(const diffractiveDissVertex& vert);
		virtual ~diffractiveDissVertex();
		
		virtual diffractiveDissVertex& operator =(const diffractiveDissVertex& vert);
		diffractiveDissVertexPtr clone(const bool cloneInParticles  = false,
		                               const bool cloneOutParticles = false) const  ///< creates deep copy of diffractive dissociation vertex; must not be virtual
		{ return diffractiveDissVertexPtr(doClone(cloneInParticles, cloneOutParticles)); }

		virtual bool addInParticle (const particlePtr&);  ///< disabled; all incoming particles have to be specified at construction
		virtual bool addOutParticle(const particlePtr&);  ///< disabled; all outgoing particles have to be specified at construction

		// production specific accessors
		virtual TVector3             zAxis    ()     const { return beam()->lzVec().Vect().Unit(); }  ///< returns z-axis defined by production process
		virtual const particlePtr&   XParticle()     const { return outParticles()[0];             }  ///< returns X particle
		virtual std::complex<double> productionAmp() const;                                           ///< returns production amplitude

		// diffractive dissociation specific accessors
		inline const particlePtr& beam  () const { return inParticles ()[0]; }  ///< returns beam particle
		inline const particlePtr& target() const { return inParticles ()[1]; }  ///< returns target particle
		inline const particlePtr& recoil() const { return outParticles()[1]; }  ///< returns recoil particle
    
		virtual bool readData(const TClonesArray& prodKinParticles,
		                      const TClonesArray& prodKinMomenta);  ///< reads data from TClonesArrays

		virtual bool revertMomenta();  ///< resets momenta to the values of last event read

		virtual std::ostream& print        (std::ostream& out) const;  ///< prints vertex parameters in human-readable form
		virtual std::ostream& dump         (std::ostream& out) const;  ///< prints all vertex data in human-readable form
		virtual std::ostream& printPointers(std::ostream& out) const;  ///< prints particle pointers strored in vertex

		virtual std::string label() const { return "diffractive dissociation vertex"; }  ///< returns label used in graph visualization and reporting

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	protected:

		virtual diffractiveDissVertex* doClone(const bool cloneInParticles,
		                                       const bool cloneOutParticles) const;  ///< helper function to use covariant return types with smart pointers; needed for public clone()


	private:

		TVector3 _beamMomCache;    ///< caches beam momentum of last event read from input data; allows to "reset" kinematics for multiple passes over the same data
		TVector3 _targetMomCache;  ///< caches target momentum of last event read from input data; allows to "reset" kinematics for multiple passes over the same data
		TVector3 _recoilMomCache;  ///< caches recoil momentum of last event read from input data; allows to "reset" kinematics for multiple passes over the same data

		static bool _debug;  ///< if set to true, debug messages are printed

	};


	inline
	diffractiveDissVertexPtr
	createDiffractiveDissVertex(const particlePtr& beam,
	                            const particlePtr& target,
	                            const particlePtr& XParticle,
	                            const particlePtr& recoil = particlePtr())
	{
		diffractiveDissVertexPtr vert(new diffractiveDissVertex(beam, target, XParticle, recoil));
		return vert;
	}


	inline
	std::ostream&
	operator <<(std::ostream&                out,
	            const diffractiveDissVertex& vert)
	{
		return vert.print(out);
	}


}  // namespace rpwa


#endif  // DIFFRACTIVEDISSVERTEX_H
