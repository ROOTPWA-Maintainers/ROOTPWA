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
//      general isobar decay amplitude in helicity formalism
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef ISOBARHELICITYAMPLITUDE_H
#define ISOBARHELICITYAMPLITUDE_H

#include "isobarAmplitude.h"


namespace rpwa {  


	class isobarHelicityAmplitude;
	typedef boost::shared_ptr<isobarHelicityAmplitude> isobarHelicityAmplitudePtr;


	class isobarHelicityAmplitude : public isobarAmplitude {
  
	public:
      
		isobarHelicityAmplitude();
		isobarHelicityAmplitude(const isobarDecayTopologyPtr& decay);
		virtual ~isobarHelicityAmplitude();

		static ParVector<LorentzRotation> hfTransform(const ParVector<LorentzVector>& daughterLv);  ///< constructs Lorentz-transformation to helicity RF of daughter particle

		std::string name() const { return "isobarHelicityAmplitude"; }
    
		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag
    

	private:

		void transformDaughters() const;  ///< boosts Lorentz-vectors of decay daughters into frames where angular distributions are defined

		ParVector<Complex> twoBodyDecayAmplitude
		(const isobarDecayVertexPtr& vertex,
		 const bool                  topVertex) const;  ///< calculates amplitude for two-body decay a -> b + c; where b and c are stable
    
		static bool _debug;  ///< if set to true, debug messages are printed
    
	};
  
  
	inline
	isobarHelicityAmplitudePtr
	createIsobarHelicityAmplitude(const isobarDecayTopologyPtr& decay)
	{
		isobarHelicityAmplitudePtr amp(new isobarHelicityAmplitude(decay));
		return amp;
	}
  
  
} // namespace rpwa


#endif  // ISOBARHELICITYAMPLITUDE_H
