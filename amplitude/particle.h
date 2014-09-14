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
//      container class for all particle related external information
//
//      particle momenta are stored in a vector that can hold data
//      from multiple events
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef PARTICLE_H
#define PARTICLE_H


#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "conversionUtils.hpp"
#include "parallelUtils.hpp"
#include "particleProperties.h"


namespace rpwa {

	class particle;
	typedef boost::shared_ptr<particle> particlePtr;


	class particle : public particleProperties {

	public:

		particle();
		particle(const particle&              part);
		particle(const particleProperties&    partProp,
		         const int                    index    = -1,
		         const int                    spinProj = 0,
		         const int                    refl     = 0,
		         const std::vector<TVector3>& momenta  = std::vector<TVector3>());
		particle(const std::string&           partName,
		         const bool                   requirePartInTable = true,
		         const int                    index              = -1,
		         const int                    spinProj           = 0,
		         const int                    refl               = 0,
		         const std::vector<TVector3>& momenta            = std::vector<TVector3>());
		particle(const std::string&           partName,
		         const int                    isospin,
		         const int                    G,
		         const int                    J,
		         const int                    P,
		         const int                    C,
		         const int                    spinProj,
		         const int                    refl  = 0,
		         const int                    index = -1);
		virtual ~particle();

		particle& operator =(const particle& part);
		particlePtr clone() const { return particlePtr(doClone()); }  ///< creates deep copy of particle; must not be virtual

		int                                spinProj     () const { return _spinProj; }  ///< returns particle's spin projection quantum number
		//const std::vector<TVector3>&       momenta      () const;                       ///< returns particle's three-momenta
		const std::vector<TLorentzVector>& lzVecs       () const { return _lzVecs;   }  ///< returns particle's Lorentz vectors
		std::vector<TLorentzVector>&       mutableLzVecs()       { return _lzVecs;   }  ///< returns particle's Lorentz vectors as modifiable refrence; use with care
		int                                index        () const { return _index;    }  ///< returns index label assigned to particle; -1 means undefined
		int                                reflectivity () const { return _refl;     }  ///< returns particle's reflectivity; 0 means undefined

		void setSpinProj    (const int                          spinProj) { _spinProj = spinProj;     }  ///< sets particle's spin projection quantum number
		void setMomenta     (const std::vector<TVector3>&       momenta );                               ///< sets particle's Lorentz vectors
		void setLzVecs      (const std::vector<TLorentzVector>& lzVecs  ) { _lzVecs   = lzVecs;       }  ///< sets particle's Lorentz vectors; if used to inject external data the mass values likely will become inconsistent
		void setIndex       (const int                          index   ) { _index    = index;        }  ///< sets particle's index label
		void setReflectivity(const int                          refl    ) { _refl     = signum(refl); }  ///< sets particle's reflectivity

		void setProperties(const particleProperties& prop);  ///< sets particle's poperties to those given by argument

		const std::vector<TLorentzVector>& transform  (const std::vector<TLorentzRotation>& lorentzTransforms);      ///< applies different Lorentz transformation to each Lorentz vector
		const std::vector<TLorentzVector>& transform  (const std::vector<TVector3>&         boosts           );      ///< applies different Lorentz boost to each Lorentz vector
		void                               scaleLzVecs(double scaleX, double scaleY, double scaleZ, double scaleE);  ///< multiplies the components of all Lorentz vectors by scaling values

		virtual std::string qnSummary() const;  ///< returns particle's quantum number summary in form name[IG(JPC)M]

		virtual std::ostream& print(std::ostream& out) const;  ///< prints particle parameters in human-readable form

		virtual std::string label() const;  ///< returns particle label

		size_t numParallelEvents() const { return _lzVecs.size(); }  ///< returns number of events handled in parallel (= size of event data arrays)

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	protected:

		virtual particle* doClone() const;  ///< helper function to use covariant return types with smart pointers; needed for public clone()

		virtual bool isEqualTo(const particleProperties& partProp) const;  ///< returns whether partProp is equal to this by checking equality of all member variables (except Lorentz vector)


	private:

		int                         _spinProj;  ///< spin projection quantum number; can be either M or helicity
		std::vector<TLorentzVector> _lzVecs;    ///< Lorentz vector [GeV]
		int                         _index;     ///< index that can be used to label indistinguishable particles
		int                         _refl;      ///< reflectivity

		static bool _debug;  ///< if set to true, debug messages are printed

	};


	inline
	particlePtr
	createParticle(const particle& part)
	{
		particlePtr partCopy(new particle(part));
		return partCopy;
	}


	inline
	particlePtr
	createParticle(const particleProperties&    partProp,
	               const int                    index    = -1,
	               const int                    spinProj = 0,
	               const int                    refl     = 0,
	               const std::vector<TVector3>& momenta  = std::vector<TVector3>())
	{
		particlePtr part(new particle(partProp, index, spinProj, refl, momenta));
		return part;
	}


	inline
	particlePtr
	createParticle(const std::string&           partName,
	               const bool                   requirePartInTable = true,
	               const int                    index              = -1,
	               const int                    spinProj           = 0,
	               const int                    refl               = 0,
	               const std::vector<TVector3>& momenta            = std::vector<TVector3>())
	{
		particlePtr part(new particle(partName, requirePartInTable, index, spinProj, refl, momenta));
		return part;
	}


	inline
	particlePtr
	createParticle(const std::string& partName,
	               const int          isospin,
	               const int          G,
	               const int          J,
	               const int          P,
	               const int          C,
	               const int          spinProj,
	               const int          refl  = 0,
	               const int          index = -1)
	{
		particlePtr part(new particle(partName, isospin, G, J, P, C, spinProj, refl, index));
		return part;
	}


	// predicate for ascending sort
	inline
	bool
	compareIndicesAsc(const particlePtr& a,
	                  const particlePtr& b)
	{
		if (not a or not b) {
			printWarn << "null pointer to particle. result undefined." << std::endl;
			return false;
		}
		return a->index() < b->index();
	}


	// predicate for descending sort
	inline
	bool
	compareIndicesDesc(const particlePtr& a,
	                   const particlePtr& b)
	{
		if (not a or not b) {
			printWarn << "null pointer to particle. result undefined." << std::endl;
			return false;
		}
		return a->index() > b->index();
	}


	inline
	std::ostream&
	operator <<(std::ostream&   out,
	            const particle& part)
	{
		return part.print(out);
	}


}  // namespace rpwa


#endif  // PARTICLE_H
