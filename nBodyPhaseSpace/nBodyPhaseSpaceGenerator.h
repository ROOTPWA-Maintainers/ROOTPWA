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
//
// calculates n-body phase space (constant matrix element) using various algorithms
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef NBODYPHASESPACEGENERATOR_HH
#define NBODYPHASESPACEGENERATOR_HH


#include <iostream>
#include <vector>

#include <TLorentzVector.h>

#ifndef __CINT__
#include "reportingUtils.hpp"
#include "mathUtils.hpp"
#endif

#include "nBodyPhaseSpaceKinematics.h"
#include "randomNumberGenerator.h"


namespace rpwa {


	class nBodyPhaseSpaceGenerator : public nBodyPhaseSpaceKinematics {

	public:

		nBodyPhaseSpaceGenerator();
		virtual ~nBodyPhaseSpaceGenerator();

		//----------------------------------------------------------------------------
		// high-level generator interface
		/// generates full event with certain n-body mass and momentum and returns event weight
		double generateDecay        (const TLorentzVector& nBody);          // Lorentz vector of n-body system in lab frame
		/// \brief generates full event with certain n-body mass and momentum only when event is accepted (return value = true)
		/// this function is more efficient, if only weighted events are needed
		bool   generateDecayAccepted(const TLorentzVector& nBody,           // Lorentz vector of n-body system in lab frame
		                             const double          maxWeight = 0);  // if positive, given value is used as maximum weight, otherwise _maxWeight


		//----------------------------------------------------------------------------
		// low-level generator interface
		/// randomly choses the (n - 2) effective masses of the respective (i + 1)-body systems
		virtual void pickMasses(const double nBodyMass);  // total energy of n-body system in its RF

		/// randomly choses the (n - 1) polar and (n - 1) azimuthal angles in the respective (i + 1)-body RFs
		inline void pickAngles();


		//----------------------------------------------------------------------------
		// weight routines
		void   setMaxWeight          (const double maxWeight) { _maxWeight = maxWeight;    }  ///< sets maximum weight used for hit-miss MC
		double maxWeight             () const                 { return _maxWeight;         }  ///< returns maximum weight used for hit-miss MC

		/// estimates maximum weight for given n-body mass
		double estimateMaxWeight(const double       nBodyMass,                 // sic!
		                         const unsigned int nmbOfIterations = 10000);  // number of generated events

		/// \brief applies event weight in form of hit-miss MC
		/// assumes that event weight has been already calculated by calcWeight()
		/// if maxWeight > 0 value is used as maximum weight, otherwise _maxWeight value is used
		inline bool eventAccepted(const double maxWeight = 0);


		virtual std::ostream& print(std::ostream& out = std::cout) const;  ///< prints generator status
		friend std::ostream& operator << (std::ostream&                   out,
		                                  const nBodyPhaseSpaceGenerator& gen) { return gen.print(out); }

	private:

		// internal variables
		double _maxWeight;  ///< maximum weight used to weight events in hit-miss MC

		ClassDef(nBodyPhaseSpaceGenerator, 1)

	};

}  // namespace rpwa


inline
void
rpwa::nBodyPhaseSpaceGenerator::pickAngles()
{
	randomNumberGenerator* random = randomNumberGenerator::instance();
	for (unsigned int i = 1; i < nmbOfDaughters(); ++i) {  // loop over 2- to n-bodies
		_cosTheta[i] = 2 * random->rndm() - 1;  // range [-1,    1]
		_phi[i]      = rpwa::twoPi * random->rndm();  // range [ 0, 2 pi]
	}
}


inline
bool
rpwa::nBodyPhaseSpaceGenerator::eventAccepted(const double maxWeight)  // if maxWeight > 0, given value is used as maximum weight, otherwise _maxWeight
{
	if (weightType() == FLAT)
		return true;  // no weighting
	const double max = (maxWeight <= 0) ? _maxWeight : maxWeight;
	if (max <= 0) {
		printErr << "maximum weight = " << max << " does not make sense. rejecting event." << std::endl;
		return false;
	}
	if ((eventWeight() / max) > randomNumberGenerator::instance()->rndm())
		return true;
	return false;
}


#endif  // NBODYPHASESPACEGENERATOR_H
