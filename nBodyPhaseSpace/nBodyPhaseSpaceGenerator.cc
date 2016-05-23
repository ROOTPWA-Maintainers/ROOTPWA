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


#include <algorithm>

#include "nBodyPhaseSpaceGenerator.h"


using namespace std;
using namespace rpwa;


nBodyPhaseSpaceGenerator::nBodyPhaseSpaceGenerator()
	: _maxWeight(0)
{ }


nBodyPhaseSpaceGenerator::~nBodyPhaseSpaceGenerator()
{ }


// generates event with certain n-body mass and momentum and returns event weigth
// general purpose function
double
nBodyPhaseSpaceGenerator::generateDecay(const TLorentzVector& nBody)  // Lorentz vector of n-body system in lab frame
{
	const double nBodyMass = nBody.M();
	if (nmbOfDaughters() < 2) {
		printWarn << "number of daughter particles = " << nmbOfDaughters() << " is smaller than 2. weight is set to 0." << endl;
		return 0;
	} else if (nBodyMass < sumOfDaughterMasses(nmbOfDaughters() - 1)) {
		printWarn << "n-body mass = " << nBodyMass << " is smaller than sum of daughter masses = "
			<< sumOfDaughterMasses(nmbOfDaughters() - 1) << ". weight is set to 0." << endl;
		return 0;
	} else {
		pickMasses(nBodyMass);
		calcWeight();
		pickAngles();
		calcEventKinematics(nBody);
	}
	return eventWeight();
}


// generates full event with certain n-body mass and momentum only, when event is accepted (return value = true)
// this function is more efficient, if only weighted events are needed
bool
nBodyPhaseSpaceGenerator::generateDecayAccepted(const TLorentzVector& nBody,      // Lorentz vector of n-body system in lab frame
                                                const double          maxWeight)  // if positive, given value is used as maximum weight, otherwise _maxWeight
{
	const double nBodyMass = nBody.M();
	if (nmbOfDaughters() < 2) {
		printWarn << "number of daughter particles = " << nmbOfDaughters() << " is smaller than 2. no event generated." << endl;
		return false;
	} else if (nBodyMass < sumOfDaughterMasses(nmbOfDaughters() - 1)) {
		printWarn << "n-body mass = " << nBodyMass << " is smaller than sum of daughter masses = "
			<< sumOfDaughterMasses(nmbOfDaughters() - 1) << ". no event generated." << endl;
		return false;
	}
	pickMasses(nBodyMass);
	calcWeight();
	if (!eventAccepted(maxWeight)) {
		return false;
	}
	pickAngles();
	calcEventKinematics(nBody);
	return true;
}


// randomly choses the (n - 2) effective masses of the respective (i + 1)-body systems
void
nBodyPhaseSpaceGenerator::pickMasses(const double nBodyMass)  // total energy of the system in its RF
{
	_M[nmbOfDaughters() - 1] = nBodyMass;
	switch (weightType()) {
		case NUPHAZ:
			{
				// for better readability notation was slightly changed w.r.t to Block's paper
				// \mathcal{F} -> F
				// \xi         -> x
				// U           -> u
				// p_i         -> prob
				for (unsigned int i = nmbOfDaughters() - 1; i >= 2; --i) {  // loop over 3- to n-bodies
					// generate variable for (i + 1)-body decay that follows importance sampling distribution as defined in eq. (12)
					//
					// 1) calculate probability
					const double sqrtU  = daughterMass(i) / _M[i];                                                    // cf. eq. (39) and (51)
					const double u      = sqrtU * sqrtU;
					const double xMin   = sumOfDaughterMasses(i - 1) * sumOfDaughterMasses(i - 1) / (_M[i] * _M[i]);  // cf. eq. (8)
					const double xMax   = (1 - sqrtU) * (1 - sqrtU);
					const double deltaX = xMax - xMin;
					const double term   = 1 + u - xMin;
					const double prob   = 1 / (i - (i - 1) * deltaX / term);                                          // cf. eq. (20)
					// 2) calculate generator for distribution
					double x;
					randomNumberGenerator* random = randomNumberGenerator::instance();
					if (random->rndm() < prob) {
						x = xMin + deltaX * pow(random->rndm(), 1 / (double)i) * pow(random->rndm(), 1 / (double)(i - 1));  // cf. eq. (21)
					} else {
						x = xMin + deltaX * pow(random->rndm(), 1 / (double)i);                                             // cf. eq. (22)
					}
					// 3) set effective isobar mass of i-body using x_i = _M(i - 1)^2 / _Mi^2
					_M[i - 1] = _M[i] * sqrt(x);
				}
			}
			break;
		default:
			{
				// create vector of sorted random values
				vector<double> r(nmbOfDaughters() - 2, 0);  // (n - 2) values needed for 2- through (n - 1)-body systems
				for (unsigned int i = 0; i < (nmbOfDaughters() - 2); ++i) {
					r[i] = randomNumberGenerator::instance()->rndm();
				}
				sort(r.begin(), r.end());
				// set effective masses of (intermediate) two-body decays
				const double massInterval = nBodyMass - sumOfDaughterMasses(nmbOfDaughters() - 1);  // kinematically allowed mass interval
				for (unsigned int i = 1; i < (nmbOfDaughters() - 1); ++i) {                         // loop over intermediate 2- to (n - 1)-bodies
					_M[i] = sumOfDaughterMasses(i) + r[i - 1] * massInterval;                   // sumOfDaughterMasses(i) is minimum effective mass
				}
			} // end default mass picking
			break;
	}
}


// calculates maximum weight for given n-body mass
double
nBodyPhaseSpaceGenerator::estimateMaxWeight(const double       nBodyMass,        // sic!
                                            const unsigned int nmbOfIterations)  // number of generated events
{
	double maxWeight = 0;
	for (unsigned int i = 0; i < nmbOfIterations; ++i) {
		pickMasses(nBodyMass);
		calcWeight();
		maxWeight = max(eventWeight(), maxWeight);
	}
	return maxWeight;
}


ostream&
nBodyPhaseSpaceGenerator::print(ostream& out) const
{
	nBodyPhaseSpaceKinematics::print(out);
	out << "nBodyPhaseSpaceGenerator parameters:" << endl
	    << "    maximum weight used in hit-miss MC ......... " << _maxWeight << endl;
	return out;
}
