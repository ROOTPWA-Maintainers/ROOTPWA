///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
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
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////


/** @addtogroup generators
 * @{
 */

#ifndef TDIFFRACTIVEPHASESPACE_HH
#define TDIFFRACTIVEPHASESPACE_HH


#include <iostream>
#include <vector>


#include "generator.h"
#include "generatorParameters.hpp"
#include "nBodyPhaseSpaceGenerator.h"
#include "particle.h"
#include "beamAndVertexGenerator.h"


class TH1;
class TLorentzVector;

namespace rpwa {

	/** @brief Phase Space generator for diffractive pion dissociation
	 *  @author Sebastian Neubert TUM (original author)
	 *
	 */
	class diffractivePhaseSpace : public rpwa::generator {

	  public:

		diffractivePhaseSpace();
		~diffractivePhaseSpace() { };

		void setDecayProducts(const std::vector<rpwa::particle>& particles);
		void addDecayProduct(const rpwa::particle& particle);

		/** @brief calculate kinematics of X system in lab frame (= target RF)
		 *
		 * calculate the Lorentz vector of the X system in the laboratory frame
		 * (target rest frame) using the Lorentz vector of the beam particle,
		 * the masses of target, recoil and X, and the Mandelstam variables
		 */
		static TLorentzVector calculateXSystemLab(const TLorentzVector& beamLorentzVector,
		                                          const double          targetMass,
		                                          const double          xMass,
		                                          const double          recoilMass,
		                                          const double          s,
		                                          const double          tPrime);

		/** @brief generates one event
		 *
		 * returns number of attempts to generate this event and beam
		 * the decay products can be fetched with GetDecay(i)
		 */
		unsigned int event();

	  private:

		void buildDaughterList();

		rpwa::nBodyPhaseSpaceGenerator _phaseSpace;

		std::vector<double> _maxXMassSlices;
		std::vector<double> _maxWeightsForXMasses;
		const static unsigned int _numberOfMassSlices = 10;

	};

}  // namespace rpwa

#endif  // TDIFFRACTIVEPHASESPACE_HH
/* @} **/
