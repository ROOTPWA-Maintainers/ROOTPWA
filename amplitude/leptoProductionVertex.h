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
//      class that describes leptoproduction vertex
//      the kinematics is defined by the incoming beam lepton, the
//      scattered lepton and the target; if the recoil particle is not
//      specified, elastic scattering is assumed
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef LEPTOPRODUCTIONVERTEX_H
#define LEPTOPRODUCTIONVERTEX_H


#include <cmath>

#include <boost/shared_ptr.hpp>

#include "productionVertex.h"


class TClonesArray;


namespace rpwa {


	class leptoProductionVertex;
	typedef boost::shared_ptr<leptoProductionVertex> leptoProductionVertexPtr;


	class leptoProductionVertex : public productionVertex {

	public:
  
		leptoProductionVertex(const particlePtr& beamLepton,
		                      const particlePtr& target,
		                      const particlePtr& XParticle,
		                      const particlePtr& recoil = particlePtr());  ///< force vertex to have two incoming (beam lepton + target) and three outgoing particles (X + recoil + scattered lepton); recoil is optional
		leptoProductionVertex(const leptoProductionVertex& vert);
		virtual ~leptoProductionVertex();
		
		leptoProductionVertex& operator =(const leptoProductionVertex& vert);
		leptoProductionVertexPtr clone(const bool cloneInParticles  = false,
		                               const bool cloneOutParticles = false) const  ///< creates deep copy of leptoproduction vertex; must not be virtual
		{ return leptoProductionVertexPtr(doClone(cloneInParticles, cloneOutParticles)); }

		virtual bool addInParticle (const particlePtr&);  ///< disabled; all incoming particles have to be specified at construction
		virtual bool addOutParticle(const particlePtr&);  ///< disabled; all outgoing particles have to be specified at construction

		// production specific accessors
		virtual const TLorentzVector& referenceLzVec() const { return virtPhoton()->lzVec(); }  ///< returns Lorentz-vector that defines z-axis for angular distributions
		virtual const particlePtr&    XParticle     () const { return outParticles()[0];     }  ///< returns X particle

		virtual std::complex<double> productionAmp() const;  ///< returns production amplitude

		virtual void setXFlavorQN();  ///< sets flavor quantum numbers of X (baryon nmb., S, C, B) to that of incoming beam particle (assumes Pomeron exchange)

		// leptoproduction specific accessors
		inline const particlePtr& beamLepton     () const { return inParticles ()[0]; }  ///< returns incoming beam particle
		inline const particlePtr& target         () const { return inParticles ()[1]; }  ///< returns target particle
		inline const particlePtr& virtPhoton     () const { return inParticles ()[2]; }  ///< returns virtual photon
		inline const particlePtr& scatteredLepton() const { return outParticles()[2]; }  ///< returns outgoing beam particle
		inline const particlePtr& recoil         () const { return outParticles()[1]; }  ///< returns recoil particle

		inline double beamPol() const { return _longPol; }  ///< returns (longitudinal) beam polarization
		inline void   setBeamPol(const double longPol = 0) { _longPol = longPol; }  ///< sets (longitudinal) beam polarization

		// leptoproduction kinematic variables
		inline double Q2     () const { return -virtPhoton()->lzVec().Mag2();                                }  ///< returns Q^2 of virtual photon
		inline double nu     () const { return target()->lzVec() * virtPhoton()->lzVec() / target()->mass(); }  ///< returns energy of virtual photon
		inline double y      () const;                                                                          ///< returns relative energy loss of virtual photon
		double        epsilon() const;                                                                          ///< returns photon's polarization parameter
		inline double delta  () const { return (2 * beamLepton()->mass2() / Q2()) * (1 - epsilon());         }  ///< returns photon's mass correction parameter
		inline double xBj    () const { return Q2() / (2 * (target()->lzVec() * virtPhoton()->lzVec()));     }  ///< returns Bjorken x
		inline double s      () const { return (target()->lzVec() + virtPhoton()->lzVec()).Mag2();           }  ///< returns total energy squared in (virtual photon, target) CM system
		inline double W      () const { return sqrt(s());                                                    }  ///< returns total energy in (virtual photon, target) CM system
		inline double delta(const double epsilon) const { return (2 * beamLepton()->mass2() / Q2()) * (1 - epsilon); }  ///< returns photon's mass correction parameter for known epsilon
    
		virtual bool initKinematicsData(const TClonesArray& prodKinPartNames);  ///< initializes input data
		virtual bool readKinematicsData(const TClonesArray& prodKinMomenta);    ///< reads input data

		virtual bool revertMomenta();  ///< resets momenta to the values of last event read

		virtual std::ostream& print        (std::ostream& out) const;  ///< prints vertex parameters in human-readable form
		virtual std::ostream& dump         (std::ostream& out) const;  ///< prints all vertex data in human-readable form
		virtual std::ostream& printPointers(std::ostream& out) const;  ///< prints particle pointers strored in vertex

		virtual std::string name() const { return "leptoProductionVertex"; }  ///< returns label used in graph visualization, reporting, and key file

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	protected:

		virtual leptoProductionVertex* doClone(const bool cloneInParticles,
		                                       const bool cloneOutParticles) const;  ///< helper function to use covariant return types with smart pointers; needed for public clone()


	private:

		double _longPol;  ///< longitudinal beam polarization

		int      _nmbProdKinPart;           ///< number of production kinematics particles in input data arrays
		TVector3 _beamLeptonMomCache;       ///< caches beam momentum of last event read from input data; allows to "reset" kinematics for multiple passes over the same data
		TVector3 _scatteredLeptonMomCache;  ///< caches beam momentum of last event read from input data; allows to "reset" kinematics for multiple passes over the same data
		TVector3 _recoilMomCache;           ///< caches recoil momentum of last event read from input data; allows to "reset" kinematics for multiple passes over the same data
		TVector3 _targetMomCache;           ///< caches target momentum of last event read from input data; allows to "reset" kinematics for multiple passes over the same data

		static bool _debug;  ///< if set to true, debug messages are printed

	};


	inline
	double
	leptoProductionVertex::y() const
	{
		return (target()->lzVec() * virtPhoton()->lzVec()) / (target()->lzVec() * beamLepton()->lzVec());
	}


	inline
	leptoProductionVertexPtr
	createLeptoProductionVertex(const particlePtr& beam,
	                            const particlePtr& target,
	                            const particlePtr& XParticle,
	                            const particlePtr& recoil = particlePtr())
	{
		leptoProductionVertexPtr vert(new leptoProductionVertex(beam, target, XParticle, recoil));
		return vert;
	}


	inline
	std::ostream&
	operator <<(std::ostream&                out,
	            const leptoProductionVertex& vert)
	{
		return vert.print(out);
	}


}  // namespace rpwa


#endif  // LEPTOPRODUCTIONVERTEX_H
