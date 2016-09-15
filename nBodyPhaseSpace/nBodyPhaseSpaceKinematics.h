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
// the n-body decay is split up into (n - 2) successive 2-body decays
// each 2-body decay is considered in its own center-of-mass frame thereby
// separating the mass from the (trivial) angular dependence
//
// the event is boosted into the same frame in which the n-body system is
// given
//
// !Note! some of the implemented weights can be applied only in
//        certain cases; if you are not sure, use the "S.U. Chung"
//        weight
//
// !Note! in case this class is used to calculate the absolute value
//        of the phase space integral, be sure to use the "S.U. Chung"
//        weight; if in the integral the phase space is weighted with
//        some amplitude with angular depdendence, the integral has to
//        be divided by the correct power of (2) pi (the "S.U. Chung"
//        already contains a factor (4 pi)^(n - 1) from (trivial)
//        integration over all angles)
//
// based on:
// GENBOD (CERNLIB W515), see F. James, "Monte Carlo Phase Space", CERN 68-15 (1968)
// NUPHAZ, see M. M. Block, "Monte Carlo phase space evaluation", Comp. Phys. Commun. 69, 459 (1992)
// S. U. Chung, "Spin Formalism", CERN Yellow Report
// S. U. Chung et. al., "Diffractive Dissociation for COMPASS"
//
// index convention:
// - all vectors have the same size (= number of decay daughters)
// - index i corresponds to the respective value in the (i + 1)-body system: effective mass M, break-up momentum, angles
// - thus some vector elements are not used like breakupMom[0], theta[0], phi[0], ...
//   this overhead is negligible compared to the ease of notation
//
// the following graph illustrates how the n-body decay is decomposed into a sequence of two-body decays
//
// n-body       ...                   3-body                 2-body                 single daughter
//
// m[n - 1]                           m[2]                   m[1]
//  ^                                  ^                      ^
//  |                                  |                      |
//  |                                  |                      |
// M[n - 1] --> ... -->               M[2] -->               M[1] -->               M    [0] = m[0]
// theta[n - 1] ...                   theta[2]               theta[1]               theta[0] = 0 (not used)
// phi  [n - 1] ...                   phi  [2]               phi  [1]               phi  [0] = 0 (not used)
// mSum [n - 1] ...                   mSum [2]               mSum [1]               mSum [0] = m[0]
// = sum_0^(n - 1) m[i]               = m[2] + m[1] + m[0]   = m[1] + m[0]
// breakUpMom[n - 1] ...              breakUpMom[2]          breakUpMom[1]          breakUpMom[0] = 0 (not used)
// = q(M[n - 1], m[n - 1], M[n - 2])  = q(M[2], m[2], M[1])  = q(M[1], m[1], m[0])
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef NBODYPHASESPACEKINEMATICS_HH
#define NBODYPHASESPACEKINEMATICS_HH


#include <iostream>
#include <vector>

#include <TLorentzVector.h>


namespace rpwa {


	class nBodyPhaseSpaceKinematics {

	public:

		nBodyPhaseSpaceKinematics();
		virtual ~nBodyPhaseSpaceKinematics();

		//----------------------------------------------------------------------------
		// generator setup
		/// sets decay constants and prepares internal variables
		bool setDecay(const std::vector<double>& daughterMasses);  // daughter particle masses
		bool setDecay(const unsigned int   nmbOfDaughters,   // number of daughter particles
		              const double*        daughterMasses);  // array of daughter particle masses


		//----------------------------------------------------------------------------
		// low-level generator interface

		// copy effective masses of (i + 1)-body systems
		void setMasses(const std::vector<double>& M);

		// copy angles of the 2-body decay of the (i + 1)-body systems
		void setAngles(const std::vector<double>& cosTheta,
		               const std::vector<double>& phi);

		enum weightTypeEnum {S_U_CHUNG = 1,   // gives physically correct mass dependence
		                     NUPHAZ    = 2,   // weight used in nuphaz
		                     GENBOD    = 3,   // default weight used in genbod (gives wrong mass dependence)
		                     FLAT      = 4};   // uniform mass distribution; warning: produces distorted angular distribution
		void           setWeightType(const weightTypeEnum weightType) { _weightType = weightType; }  ///< selects formula used for weight calculation
		weightTypeEnum weightType   () const                          { return _weightType;       }  ///< returns formula used for weight calculation

		/// \brief computes breakup momenta
		void calcBreakupMomenta();

		/// \brief computes event weight and breakup momenta
		double calcWeight();

		enum kinematicsTypeEnum {BLOCK         = 1,   // method for calculation of event kinematics used in nuphaz (faster)
		                         RAUBOLD_LYNCH = 2};  // method for calculation of event kinematics used in genbod
		void               setKinematicsType(const kinematicsTypeEnum kinematicsType) { _kinematicsType = kinematicsType; }  ///< selects algorithm used to calculate event kinematics
		kinematicsTypeEnum kinematicsType   () const                                  { return _kinematicsType;           }  ///< returns algorithm used to calculate event kinematics

		/// \brief calculates full event kinematics from the effective masses of the (i + 1)-body systems and the Lorentz vector of the decaying system
		/// uses the break-up momenta calculated by calcBreakupMomenta() (either called
		/// directly or implicitly via calcWeight()) and angles from setAngles()
		void calcEventKinematics(const TLorentzVector& nBody);  // Lorentz vector of n-body system in lab frame


		double normalization         () const                 { return _norm;              }  ///< returns normalization used in weight calculation
		double eventWeight           () const                 { return _weight;            }  ///< returns weight of generated event
		double maxWeightObserved     () const                 { return _maxWeightObserved; }  ///< returns maximum observed weight since instantiation
		void   resetMaxWeightObserved()                       { _maxWeightObserved = 0;    }  ///< sets maximum observed weight back to zero

		//----------------------------------------------------------------------------
		// trivial accessors
		const TLorentzVector&              daughter           (const int index) const { return _daughters[index];  }  ///< returns Lorentz vector of daughter at index
		const std::vector<TLorentzVector>& daughters          ()                const { return _daughters;         }  ///< returns Lorentz vectors of all daughters
		unsigned int                       nmbOfDaughters     ()                const { return _n;                 }  ///< returns number of daughters
		double                             daughterMass       (const int index) const { return _m[index];          }  ///< returns invariant mass of daughter at index
		double                             intermediateMass   (const int index) const { return _M[index];          }  ///< returns intermediate mass of (index + 1)-body system
		double                             sumOfDaughterMasses(const int index) const { return _mSum[index];       }  ///< returns sum of invariant masses of daughters up to index
		double                             breakupMom         (const int index) const { return _breakupMom[index]; }  ///< returns breakup momentum in (index + 1)-body RF
		double                             cosTheta           (const int index) const { return _cosTheta[index];   }  ///< returns polar angle in (index + 1)-body RF
		double                             phi                (const int index) const { return _phi[index];        }  ///< returns azimuth in (index + 1)-body RF


		virtual std::ostream& print(std::ostream& out = std::cout) const;  ///< prints generator status
		friend std::ostream& operator << (std::ostream&                    out,
		                                  const nBodyPhaseSpaceKinematics& gen) { return gen.print(out); }

	private:

		/// \brief reduced 2-body phase-space factor = 2 * q / M
		/// needed for NUPHAS weight calculation
		static inline double F(const double x,
		                       const double y);

		// parameters and constants
		weightTypeEnum              _weightType;         ///< switches between different weight formulas
		kinematicsTypeEnum          _kinematicsType;     ///< switches between different ways of calculating event kinematics

		unsigned int                _n;                  ///< number of daughter particles
		std::vector<double>         _m;                  ///< masses of daughter particles
		std::vector<double>         _mSum;               ///< sums of daughter particle masses
		double                      _norm;               ///< normalization value

	protected:

		// variables that can be changed by inheriting classes
		std::vector<double>         _M;                  ///< effective masses of (i + 1)-body systems
		std::vector<double>         _cosTheta;           ///< cosine of polar angle of the 2-body decay of the (i + 1)-body system
		std::vector<double>         _phi;                ///< azimuthal angle of the 2-body decay of the (i + 1)-body system

	private:

		// variables
		std::vector<double>         _breakupMom;         ///< breakup momenta for the two-body decays: (i + 1)-body --> daughter_(i + 1) + i-body
		std::vector<TLorentzVector> _daughters;          ///< Lorentz vectors of the daughter particles
		double                      _weight;             ///< phase space weight of generated event
		double                      _maxWeightObserved;  ///< maximum event weight calculated processing the input data

		ClassDef(nBodyPhaseSpaceKinematics, 1)

	};

}  // namespace rpwa


inline
double
rpwa::nBodyPhaseSpaceKinematics::F(const double x,
                                   const double y)
{
	const double val = 1 + (x - y) * (x - y) - 2 * (x + y);
	if (val < 0)
		return 0;
	else
		return sqrt(val);
}


#endif  // NBODYPHASESPACEKINEMATICS_H
