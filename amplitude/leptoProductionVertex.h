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
		virtual const std::vector<TLorentzVector>& referenceLzVecs() const { return virtPhoton()->lzVecs(); }  ///< returns Lorentz-vectors for a number of events that defines z-axis for angular distributions
		virtual const particlePtr&                 XParticle      () const { return outParticles()[0];      }  ///< returns X particle

		virtual std::vector<std::complex<double> > productionAmps() const;  ///< returns production amplitudes all events stored in particles

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
		void Q2     (std::vector<double>& result) const; ///< returns Q^2 of virtual photon
		void nu     (std::vector<double>& result) const; ///< returns energy of virtual photon
		void y      (std::vector<double>& result) const; ///< returns relative energy loss of virtual photon
		void epsilon(std::vector<double>& result) const; ///< returns photon's polarization parameter
		void delta  (std::vector<double>& result) const; ///< returns photon's mass correction parameter
		void xBj    (std::vector<double>& result) const; ///< returns Bjorken x
		void s      (std::vector<double>& result) const; ///< returns total energy squared in (virtual photon, target) CM system
		void W      (std::vector<double>& result) const; ///< returns total energy in (virtual photon, target) CM system
		void delta  (const std::vector<double>& epsilon, std::vector<double>& result) const; ///< returns photon's mass correction parameter for known epsilon

		virtual bool initKinematicsData(const TClonesArray& prodKinPartNames                     );  ///< initializes input data
		virtual bool readKinematicsData(const std::vector<std::vector<TVector3> >& prodKinMomenta);  ///< reads multiple input data event

		virtual bool revertMomenta();  ///< resets momenta to the values of last read event block

		virtual std::ostream& print        (std::ostream& out) const;  ///< prints vertex parameters in human-readable form
		virtual std::ostream& dump         (std::ostream& out) const;  ///< prints all vertex data in human-readable form
		virtual std::ostream& printPointers(std::ostream& out) const;  ///< prints particle pointers stored in vertex

		virtual std::string name() const { return "leptoProductionVertex"; }  ///< returns label used in graph visualization, reporting, and key file

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	protected:

		virtual leptoProductionVertex* doClone(const bool cloneInParticles,
		                                       const bool cloneOutParticles) const;  ///< helper function to use covariant return types with smart pointers; needed for public clone()


	private:

		double _longPol;  ///< longitudinal beam polarization

		int                   _nmbProdKinPart;           ///< number of production kinematics particles in input data arrays
		std::vector<TVector3> _beamLeptonMomCache;       ///< caches beam momenta of last block of events read from input data; allows to "reset" kinematics for multiple passes over the same data
		std::vector<TVector3> _scatteredLeptonMomCache;  ///< caches beam momenta of last block of events read from input data; allows to "reset" kinematics for multiple passes over the same data
		std::vector<TVector3> _recoilMomCache;           ///< caches recoil momenta of last block of events read from input data; allows to "reset" kinematics for multiple passes over the same data
		std::vector<TVector3> _targetMomCache;           ///< caches target momenta of last block of events read from input data; allows to "reset" kinematics for multiple passes over the same data

		static bool _debug;  ///< if set to true, debug messages are printed

	};

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
