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
//      container class for all particle related external information
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

#include <boost/shared_ptr.hpp>

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "particleProperties.h"


namespace rpwa {

	class particle;
	typedef boost::shared_ptr<particle> particlePtr;


	class particle : public particleProperties {
	
	public:
			
		particle();
		particle(const particle&           part);
		particle(const particleProperties& partProp,
		         const int                 charge,
		         const int                 index    = -1,
		         const int                 spinProj = 0,
		         const int                 refl     = 0,
		         const TVector3&           momentum = TVector3());
		particle(const std::string&        partName,
		         const bool                requirePartInTable = true,
		         const int                 index              = -1,
		         const int                 spinProj           = 0,
		         const int                 refl               = 0,
		         const TVector3&           momentum           = TVector3());
		particle(const std::string&        partName,
		         const int                 isospin,
		         const int                 G,
		         const int                 J,
		         const int                 P,
		         const int                 C,
		         const int                 spinProj,
		         const int                 refl  = 0,
		         const int                 index = -1);
		virtual ~particle();

		particle& operator =(const particle& part);
		particlePtr clone() const { return particlePtr(doClone()); }  ///< creates deep copy of particle; must not be virtual

		std::string           name        () const;                                         ///< returns particle name including charge
		std::string           bareName    () const { return stripChargeFromName(name()); }  ///< returns particle name w/o charge
		int                   charge      () const { return _charge;                     }  ///< returns particle's charge
		int                   spinProj    () const { return _spinProj;                   }  ///< returns particle's spin projection quantum number
		TVector3              momentum    () const { return _lzVec.Vect();               }  ///< returns three-momentum of particle
		const TLorentzVector& lzVec       () const { return _lzVec;                      }  ///< returns Lorentz vector of particle
		int                   index       () const { return _index;                      }  ///< returns index label assigned to particle; -1 means undefined
		int                   reflectivity() const { return _refl;                       }  ///< returns particle's reflectivity; 0 means undefined
		int                   isospinProj () const
		{ return 2 * _charge - (baryonNmb() + strangeness() + charm() + beauty()); }  ///< returns z-component of isospin using Gell-Mann-Nishijima formula (see PDG 2008 eq. 14.1)

		void setCharge      (const int             charge  ) { _charge   = charge;                                                            }  ///< sets particle's charge
		void setSpinProj    (const int             spinProj) { _spinProj = spinProj;                                                          }  ///< sets particle's spin projection quantum number
		void setMomentum    (const TVector3&       momentum) { _lzVec    = TLorentzVector(momentum, sqrt(momentum.Mag2() + mass() * mass())); }  ///< sets particle's Lorentz vector
		void setLzVec       (const TLorentzVector& lzVec   ) { _lzVec    = lzVec;                                                             }  ///< sets particle's Lorentz vector; if this is used to inject external data the mass values likely become inconsistent
		void setIndex       (const int             index   ) { _index    = index;                                                             }  ///< sets particle's index label
		void setReflectivity(const int             refl    ) { _refl     = refl;                                                              }  ///< sets particle's reflectivity

		void setProperties(const particleProperties& prop);  ///< sets particle's poperties to those given by argument

		const TLorentzVector& transform(const TLorentzRotation& L)     { return _lzVec.Transform(L); }  ///< applies Lorentz-transformation to particle
		const TLorentzVector& transform(const TVector3&         boost)  ///< applies Lorentz-boost to particle
		{
			_lzVec.Boost(boost);
			return _lzVec;
		}

		virtual std::string qnSummary() const;  ///< returns particle's quantum number summary in form name[IG(JPC)M]

		virtual std::ostream& print(std::ostream& out) const;  ///< prints particle parameters in human-readable form

		virtual std::string label() const;  ///< returns graph label

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag

    
	protected:
    
		virtual particle* doClone() const;  ///< helper function to use covariant return types with smart pointers; needed for public clone()
			

	private:

		int            _charge;    ///< charge
		int            _spinProj;  ///< spin projection quantum number; can be either M or helicity
		TLorentzVector _lzVec;     ///< Lorentz vector [GeV]
		int            _index;     ///< index that can be used to label indistinguishable particles
		int            _refl;      ///< reflectivity

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
	createParticle(const particleProperties& partProp,
	               const int                 charge,
	               const int                 index    = -1,
	               const int                 spinProj = 0,
	               const int                 refl     = 0,
	               const TVector3&           momentum = TVector3())
	{
		particlePtr part(new particle(partProp, charge, index, spinProj, refl, momentum));
		return part;
	}


	inline
	particlePtr
	createParticle(const std::string& partName,
	               const bool         requirePartInTable = true,
	               const int          index              = -1,
	               const int          spinProj           = 0,
	               const int          refl               = 0,
	               const TVector3&    momentum           = TVector3())
	{
		particlePtr part(new particle(partName, requirePartInTable, index, spinProj, refl, momentum));
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
