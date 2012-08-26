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
//      virtual base class for isobar decay amplitude independent of formalism
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef ISOBARAMPLITUDE_H
#define ISOBARAMPLITUDE_H


#include <complex>
#include <map>
#include <vector>

#include "isobarDecayTopology.h"


namespace rpwa {  


	class isobarAmplitude;
	typedef boost::shared_ptr<rpwa::isobarAmplitude> isobarAmplitudePtr;


	class isobarAmplitude {
  
	public:
      
		isobarAmplitude();
		isobarAmplitude(const isobarDecayTopologyPtr& decay);
		virtual ~isobarAmplitude();

		const isobarDecayTopologyPtr& decayTopology   () const { return _decay; }              ///< returns pointer to decay topology
		void                          setDecayTopology(const isobarDecayTopologyPtr& decay);   ///< sets decay topology

		virtual void init();  ///< initializes amplitude; needs to be called when decay topology was changed

		bool reflectivityBasis    () const { return _useReflectivityBasis; }  ///< returns whether reflectivity basis is used
		bool boseSymmetrization   () const { return _boseSymmetrize;       }  ///< returns whether Bose symmetrization is used
		bool isospinSymmetrization() const { return _isospinSymmetrize;    }  ///< returns whether isospin symmetrization is used
		void enableReflectivityBasis    (const bool flag = true) { _useReflectivityBasis = flag; }  ///< en/disables use of reflectivity basis
		void enableBoseSymmetrization   (const bool flag = true) { _boseSymmetrize       = flag; }  ///< en/disables use of Bose symmetrization
		void enableIsospinSymmetrization(const bool flag = true) { _isospinSymmetrize    = flag; }  ///< en/disables use of isospin symmetrization

		bool doSpaceInversion() const { return _doSpaceInversion; }  ///< returns whether parity transformation is performed on decay
		bool doReflection    () const { return _doReflection;     }  ///< returns whether decay is reflected through production plane
		void enableSpaceInversion(const bool flag = true) { _doSpaceInversion = flag; }  ///< en/disables parity transformation of decay
		void enableReflection    (const bool flag = true) { _doReflection     = flag; }  ///< en/disables reflection of decay through production plane

		static TLorentzRotation gjTransform(const TLorentzVector& beamLv,
		                                    const TLorentzVector& XLv);  ///< constructs Lorentz-transformation to X Gottfried-Jackson frame
    
		std::complex<double> amplitude()   const;                         ///< computes amplitude
		std::complex<double> operator ()() const { return amplitude(); }  ///< computes amplitude

		virtual std::string   name           ()                  const { return "isobarAmplitude"; }
		virtual std::ostream& printParameters(std::ostream& out) const;  ///< prints amplitude parameters in human-readable form
		virtual std::ostream& print          (std::ostream& out) const;  ///< prints amplitude in human-readable form
    
		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag
    

	protected:

		void spaceInvertDecay() const;  ///< performs parity transformation on all decay three-momenta
		void reflectDecay    () const;  ///< performs reflection through production plane on all decay three-momenta

		virtual void transformDaughters() const = 0;  ///< boosts Lorentz-vectors of decay daughters into frames where angular distributions are defined

		virtual std::complex<double> twoBodyDecayAmplitude
		(const isobarDecayVertexPtr& vertex,
		 const bool                  topVertex) const = 0;  ///< calculates amplitude for two-body decay a -> b + c; where b and c are stable
    
		virtual std::complex<double> twoBodyDecayAmplitudeSum
		(const isobarDecayVertexPtr& vertex,
		 const bool                  topVertex = false) const;  ///< recursive function that sums up decay amplitudes for all allowed helicitities for all vertices below the given vertex

		virtual std::complex<double> symTermAmp(const std::vector<unsigned int>& fsPartPermMap) const;  ///< returns decay amplitude for a certain permutation of final-state particles

		virtual void genBoseSymTermMaps
		(const std::map<std::string, std::vector<unsigned int> >&     origFsPartIndices,
		 const std::map<std::string, std::vector<unsigned int> >&     newFsPartIndices,
		 std::map<std::string, std::vector<unsigned int> >::iterator& newFsPartIndicesEntry,
		 const std::vector<unsigned int>&                             baseFsPartPermMap,
		 std::vector<symTermMap>&                                     symTermMaps) const;  ///< recursive function that generates all permutation maps of indistinguishable final state particles
		virtual void initBoseSymTermMaps();  ///< generates final-state permutation maps for Bose symmetrization

		virtual void initIsospinSymTermMaps();  ///< generates final-state permutation maps for isospin symmetrization


		isobarDecayTopologyPtr  _decay;                 ///< isobar decay topology with all external information
		bool                    _useReflectivityBasis;  ///< if set, reflectivity basis is used to calculate the X decay node
		bool                    _boseSymmetrize;        ///< if set, amplitudes are Bose symmetrized
		bool                    _isospinSymmetrize;     ///< if set, amplitudes are isospin symmetrized
		bool                    _doSpaceInversion;      ///< is set, all three-momenta of the decay particles are parity transformed (for test purposes)
		bool                    _doReflection;          ///< is set, all three-momenta of the decay particles are reflected through production plane (for test purposes)
		std::vector<symTermMap> _symTermMaps;           ///< array of factors and permutation maps for symmetrization terms
    
		static bool _debug;  ///< if set to true, debug messages are printed
    
	};
  

	inline
	std::ostream&
	operator <<(std::ostream&          out,
	            const isobarAmplitude& amp)
	{
		return amp.print(out);
	}
  
  
} // namespace rpwa


#endif  // ISOBARAMPLITUDE_H
