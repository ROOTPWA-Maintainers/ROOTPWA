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
//      class that describes production vertex in diffractive
//      dissociation the beam-Reggeon-X vertex has exactly one
//      incoming beam and one outgoing particle (X), which
//      unambiguously defines the Reggeon kinematics; if the target
//      momentum is not specified, a fixed target is assumed; if the
//      recoil particle is not specified, elastic scattering is
//      assumed
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

		diffractiveDissVertex& operator =(const diffractiveDissVertex& vert);
		diffractiveDissVertexPtr clone(const bool cloneInParticles  = false,
		                               const bool cloneOutParticles = false) const  ///< creates deep copy of diffractive dissociation vertex; must not be virtual
		{ return diffractiveDissVertexPtr(doClone(cloneInParticles, cloneOutParticles)); }

		virtual bool addInParticle (const particlePtr&);  ///< disabled; all incoming particles have to be specified at construction
		virtual bool addOutParticle(const particlePtr&);  ///< disabled; all outgoing particles have to be specified at construction

		// production specific accessors
		virtual const ParVector<LorentzVector>& referenceLzVecs() const { return beam()->lzVecs();  }  ///< returns Lorentz-vector that defines z-axis for angular distributions
		virtual const particlePtr&              XParticle      () const { return outParticles()[0]; }  ///< returns X particle

		virtual ParVector<Complex> productionAmps() const;  ///< returns production amplitudes all events stored in particles

		virtual void setXFlavorQN();  ///< sets flavor quantum numbers of X (baryon nmb., S, C, B) to that of incoming beam particle (assumes Pomeron exchange)

		// diffractive dissociation specific accessors
		inline const particlePtr& beam  () const { return inParticles ()[0]; }  ///< returns beam particle
		inline const particlePtr& target() const { return inParticles ()[1]; }  ///< returns target particle
		inline const particlePtr& recoil() const { return outParticles()[1]; }  ///< returns recoil particle

		virtual bool initKinematicsData(const TClonesArray& prodKinPartNames                    );  ///< initializes input data
		virtual bool readKinematicsData(const std::vector<ParVector<Vector3> >& prodKinMomenta);  ///< reads multiple input data event

		virtual bool revertMomenta();  ///< resets momenta to the values of last read event block

		virtual std::ostream& print        (std::ostream& out) const;  ///< prints vertex parameters in human-readable form
		virtual std::ostream& dump         (std::ostream& out) const;  ///< prints all vertex data in human-readable form
		virtual std::ostream& printPointers(std::ostream& out) const;  ///< prints particle pointers stored in vertex

		virtual std::string name() const { return "diffractiveDissVertex"; }  ///< returns label used in graph visualization, reporting, and key file

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	protected:

		virtual diffractiveDissVertex* doClone(const bool cloneInParticles,
		                                       const bool cloneOutParticles) const;  ///< helper function to use covariant return types with smart pointers; needed for public clone()


	private:

		size_t             _nmbProdKinPart;  ///< number of production kinematics particles in input data arrays
		ParVector<Vector3> _beamMomCache;    ///< caches beam momenta of last block of events read from input data; allows to "reset" kinematics for multiple passes over the same data
		ParVector<Vector3> _recoilMomCache;  ///< caches recoil momenta of last block of events read from input data; allows to "reset" kinematics for multiple passes over the same data
		ParVector<Vector3> _targetMomCache;  ///< caches target momenta of last block of events read from input data; allows to "reset" kinematics for multiple passes over the same data

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
