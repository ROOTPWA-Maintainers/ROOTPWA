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


#include <algorithm>

#include "TMath.h"

#include "reportingUtilsRoot.hpp"
#include "physUtils.hpp"
#include "factorial.hpp"
#include "nBodyPhaseSpaceGen.h"


using namespace std;
using namespace rpwa;


ClassImp(nBodyPhaseSpaceGen);


nBodyPhaseSpaceGen::nBodyPhaseSpaceGen()
  : _n                (0),
    _weightType       (S_U_CHUNG),
    _norm             (0),
    _weight           (0),
 _impweight(1),
    _maxWeightObserved(0),
    _maxWeight        (0),

    _kinematicsType   (BLOCK),

    _verbose(false),

    _isoBWMass(1.0),_isoBWWidth(0.01) // dummy values!!!
{ }


nBodyPhaseSpaceGen::~nBodyPhaseSpaceGen()
{ }


// sets decay constants and prepares internal variables
bool
nBodyPhaseSpaceGen::setDecay(const vector<double>& daughterMasses)  // array of daughter particle masses
{
  _n = daughterMasses.size();
  if (_n < 2) {
    printErr << "number of daughters = " << _n << " does not make sense." << endl;
    return false;
  }
  // copy daughter masses
  _m.clear();
  _m = daughterMasses;
  // prepare effective mass vector
  _M.clear();
  _M.resize(_n, 0);
  _M[0] = _m[0];
  // prepare angle vectors
  _cosTheta.clear();
  _cosTheta.resize(_n, 0);
  _phi.clear();
  _phi.resize(_n, 0);
  // calculate daughter mass sums
  _mSum.clear();
  _mSum.resize(_n, 0);
  _mSum[0] = _m[0];
  for (unsigned int i = 1; i < _n; ++i)
    _mSum[i] = _mSum[i - 1] + _m[i];
  // prepare breakup momentum vector
  _breakupMom.clear();
  _breakupMom.resize(_n, 0);
  // prepare vector for daughter Lorentz vectors
  _daughters.clear();
  _daughters.resize(_n, TLorentzVector(0, 0, 0, 0));
  // calculate normalization
  switch (_weightType) {
  case S_U_CHUNG: case IMPORTANCE:
    // S. U. Chung's normalization
	  _norm = 1 / (2 * pow(twoPi, 2 * (int)_n - 3) * rpwa::factorial<double>(_n - 2));
    break;
  case NUPHAZ:
    {  // NUPHAZ's normalization: calculates phase space for massless daughters
	    const double fact = rpwa::factorial<double>(_n - 2);
      _norm = 8 / (pow(fourPi, 2 * (int)_n + 1) * fact * fact * (_n - 1));
    }
    break;
  case GENBOD: case FLAT:
    _norm = 1;
    break;
  default:
    printWarn << "unknown weight type. setting normalization to 1." << endl;
    _norm = 1;
    break;
  }
  resetMaxWeightObserved();
  return true;
}


// set decay constants and prepare internal variables
bool
nBodyPhaseSpaceGen::setDecay(const unsigned int nmbOfDaughters,  // number of daughter particles
			     const double*      daughterMasses)  // array of daughter particle masses
{
  vector <double> m;
  m.resize(nmbOfDaughters, 0);
  for (unsigned int i = 0; i < nmbOfDaughters; ++i)
    m[i] = daughterMasses[i];
  return setDecay(m);
}


// generates event with certain n-body mass and momentum and returns event weigth
// general purpose function
double
nBodyPhaseSpaceGen::generateDecay(const TLorentzVector& nBody)  // Lorentz vector of n-body system in lab frame
{
  _weight = 1;
  const double nBodyMass = nBody.M();
  if (_n < 2) {
    printWarn << "number of daughter particles = " << _n << " is smaller than 2. weight is set to 0." << endl;
    _weight = 0;
  } else if (nBodyMass < _mSum[_n - 1]) {
    printWarn << "n-body mass = " << nBodyMass << " is smaller than sum of daughter masses = "
	      << _mSum[_n - 1] << ". weight is set to 0." << endl;
    _weight = 0;
  } else {
    pickMasses(nBodyMass);
    calcWeight();
    pickAngles();
    calcEventKinematics(nBody);
  }
#if DEBUG
  print();
#endif
  return _weight;
}


// generates full event with certain n-body mass and momentum only, when event is accepted (return value = true)
// this function is more efficient, if only weighted evens are needed
bool
nBodyPhaseSpaceGen::generateDecayAccepted(const TLorentzVector& nBody,      // Lorentz vector of n-body system in lab frame
					  const double          maxWeight)  // if positive, given value is used as maximum weight, otherwise _maxWeight
{
  const double nBodyMass = nBody.M();
  if (_n < 2) {
    printWarn << "number of daughter particles = " << _n << " is smaller than 2. no event generated." << endl;
    return false;
  } else if (nBodyMass < _mSum[_n - 1]) {
    printWarn << "n-body mass = " << nBodyMass << " is smaller than sum of daughter masses = "
	      << _mSum[_n - 1] << ". no event generated." << endl;
    return false;
  }
  pickMasses(nBodyMass);
  calcWeight();
  if (!eventAccepted(maxWeight))
    return false;
  pickAngles();
  calcEventKinematics(nBody);
  return true;
}


// randomly choses the (n - 2) effective masses of the respective (i + 1)-body systems
void
nBodyPhaseSpaceGen::pickMasses(const double nBodyMass)  // total energy of the system in its RF
{

  _M[_n - 1] = nBodyMass;
  switch (_weightType) {
  case NUPHAZ:
    {  // the NUPHAZ algorithm calculates part of the weight along with the effective isobar masses
      // for better readability notation was slightly changed w.r.t to Block's paper
      // \mathcal{F} -> F
      // \xi         -> x
      // U           -> u
      // p_i         -> prob
      _weight = 1;
      for (unsigned int i = _n - 1; i >= 2; --i) {  // loop over 3- to n-bodies
	// generate variable for (i + 1)-body decay that follows importance sampling distribution as defined in eq. (12)
	//
	// 1) calculate probability
	const double sqrtU  = _m[i] / _M[i];                                  // cf. eq. (39) and (51)
	const double u      = sqrtU * sqrtU;
	const double xMin   = _mSum[i - 1] * _mSum[i - 1] / (_M[i] * _M[i]);  // cf. eq. (8)
	const double xMax   = (1 - sqrtU) * (1 - sqrtU);
	const double deltaX = xMax - xMin;
	const double term   = 1 + u - xMin;
	const double prob   = 1 / (i - (i - 1) * deltaX / term);              // cf. eq. (20)
	// 2) calculate generator for distribution
	double x;
	randomNumberGenerator* random = randomNumberGenerator::instance();
	if (random->rndm() < prob)
	  x = xMin + deltaX * pow(random->rndm(), 1 / (double)i) * pow(random->rndm(), 1 / (double)(i - 1));  // cf. eq. (21)
	else
	  x = xMin + deltaX * pow(random->rndm(), 1 / (double)i);  // cf. eq. (22)
	// 3) calculate weight factor
	const double deltaZ = (i * term - (i - 1) * deltaX) * pow(deltaX, (int)i - 1);  // cf. eq. (17) and (18)
	_weight *= deltaZ * pow(x / (x - xMin), (int)i - 2) * F(x, u) / (1 + u - x);    // cf. eq. (24)
	// 4) set effective isobar mass of i-body using x_i = _M(i - 1)^2 / _Mi^2
	_M[i - 1] = _M[i] * sqrt(x);
      }
    }
    break;
    /*  case IMPORTANCE:
    {

      // set effective masses of (intermediate) two-body decays
      const double massInterval = nBodyMass - _mSum[_n - 1];  // kinematically allowed mass interval

      //  do first isobar breitwigner importance sampling:
      double m=0;//unsigned int count=0;
      //while(m<_mSum[_n-2] || m>_mSum[_n-2]+massInterval){
      //m=gRandom->BreitWigner(_isoBWMass,_isoBWWidth);
	m=gRandom->Uniform(_mSum[_n-2],_mSum[_n-2]+massInterval);
	//printErr <<  _M[_n-3] << " < " << m << " < " << _mSum[_n-2]+massInterval << endl;
	//}
      _M[_n-2]=m;
      _impweight=TMath::BreitWigner(m,_isoBWMass,_isoBWWidth);  // for de-weighting (has to be read out explicitely!!!)

      // now that the first isobar mass is fixed, generate the rest:
      // create vector of sorted random values
      vector<double> r(_n - 2, 0);  // (n - 2) values needed for 2- through (n - 1)-body systems
      r[_n-3]=(_M[_n-2]-_mSum[_n-2])/massInterval;


      for (unsigned int i = _n-3; i > 0; --i)
	r[i-1] = gRandom->Uniform(0,r[i]);

      for (unsigned int i = 1; i < (_n - 2); ++i)             // loop over intermediate 2- to (n - 1)-bodies
	_M[i] = _mSum[i] + r[i - 1] * massInterval;           // _mSum[i] is minimum effective mass



      if(_verbose){
	for(unsigned int i =0; i < (_n - 1) ; ++i){
	  cerr << "M["<<i<<"]="<<_M[i] << endl;
	}
      }// end ifverbose
    }// end if importance sampling


    break;*/
  default:
    {

      // create vector of sorted random values
      vector<double> r(_n - 2, 0);  // (n - 2) values needed for 2- through (n - 1)-body systems
      for (unsigned int i = 0; i < (_n - 2); ++i)
	r[i] = randomNumberGenerator::instance()->rndm();
      sort(r.begin(), r.end());
      // set effective masses of (intermediate) two-body decays
      const double massInterval = nBodyMass - _mSum[_n - 1];  // kinematically allowed mass interval
      for (unsigned int i = 1; i < (_n - 1); ++i)             // loop over intermediate 2- to (n - 1)-bodies
	  _M[i] = _mSum[i] + r[i - 1] * massInterval;           // _mSum[i] is minimum effective mass

      //cerr << _M[1] << endl;
      if(_weightType==IMPORTANCE){
	_impweight=TMath::BreitWigner(_M[_n-2],_isoBWMass,_isoBWWidth);
	// BE CAREFULL:::: hard coded a1 BW in 3pi mass
	_impweight*=TMath::BreitWigner(_M[_n-3],1.23,0.425);
	// BE CAREFULL:::: hard coded rho BW in 2pi mass
	_impweight*=TMath::BreitWigner(_M[_n-4],0.77,0.15);
      }
      //cerr << _impweight << endl;
    } // end default mass picking
    break;
  }
}


// computes event weight (= integrand value) and breakup momenta
// uses vector of intermediate two-body masses prepared by pickMasses()
double
nBodyPhaseSpaceGen::calcWeight()
{
  for (unsigned int i = 1; i < _n; ++i)  // loop over 2- to n-bodies
    _breakupMom[i] = breakupMomentum(_M[i], _M[i - 1], _m[i]);
  switch (_weightType) {
  case S_U_CHUNG: case IMPORTANCE:
    {  // S. U. Chung's weight
      double momProd = 1;                    // product of breakup momenta
      for (unsigned int i = 1; i < _n; ++i)  // loop over 2- to n-bodies
	momProd *= _breakupMom[i];
      const double massInterval = _M[_n - 1] - _mSum[_n - 1];  // kinematically allowed mass interval
      _weight = _norm * pow(massInterval, (int)_n - 2) * momProd / _M[_n - 1] * _impweight;
    }
    break;
  case NUPHAZ:
    {  // NUPHAZ's weight
      const double M2 = _M[1] * _M[1];
      _weight *= _norm * pow(_M[_n - 1], 2 * (int)_n - 4) * F(_m[1] * _m[1] / M2, _m[0] * _m[0] / M2);
      //_weight *= F(_m[1] * _m[1] / M2, _m[0] * _m[0] / M2);
    }
    break;
  case GENBOD:
    {  // GENBOD's weight; does not reproduce dependence on n-body mass correctly
      double momProd = 1;                    // product of breakup momenta
      for (unsigned int i = 1; i < _n; ++i)  // loop over 2- to n-bodies
	momProd *= _breakupMom[i];
      double motherMassMax = _M[_n - 1] - _mSum[_n - 1] + _m[0];  // maximum possible value of decaying effective mass
      double momProdMax    = 1;                                   // product of maximum breakup momenta
      for (unsigned int i = 1; i < _n; ++i) {  // loop over 2- to n-bodies
	motherMassMax += _m[i];
	momProdMax    *= breakupMomentum(motherMassMax, _mSum[i - 1], _m[i]);
      }
      _weight = momProd / momProdMax;
    }
    break;
  case FLAT:
    // no weighting
    //!!! warning: produces distorted angular distribution
    _weight = 1;
    break;
  default:
    printWarn << "unknown weight type. setting weight to 1." << endl;
    _weight = 1;
    break;
  }
  if (_weight > _maxWeightObserved)
    _maxWeightObserved = _weight;
  if (isnan(_weight))
    printWarn << "weight = " << _weight << endl;
  return _weight;
}


// calculates complete event from the effective masses of the (i + 1)-body
// systems, the Lorentz vector of the decaying system, and the decay angles
// uses the break-up momenta calculated by calcWeight()
void
nBodyPhaseSpaceGen::calcEventKinematics(const TLorentzVector& nBody)  // Lorentz vector of n-body system in lab frame
{
  switch (_kinematicsType) {
  case RAUBOLD_LYNCH:
    {
      // builds event starting in 2-body RF going up to n-body
      // is less efficicient, since to calculate kinematics in i-body RF
      // all the daughters of the (i - 1)-body have to bu boosted into the i-body RF
      // and rotated according to the chosen angles
      //
      // construct Lorentz vector of first daughter in 2-body RF
      _daughters[0].SetPxPyPzE(0, 0, _breakupMom[1], sqrt(_m[0] * _m[0] + _breakupMom[1] * _breakupMom[1]));
      for (unsigned int i = 1; i < _n; ++i) {  // loop over remaining (n - 1) daughters
	// construct daughter Lorentz vector in (i + 1)-body RF; i-body is along z-axis
	_daughters[i].SetPxPyPzE(0, 0, -_breakupMom[i], sqrt(_m[i] * _m[i] + _breakupMom[i] * _breakupMom[i]));
	// rotate all daughters of the (i + 1)-body system
	const double cT = _cosTheta[i];
	const double sT = sqrt(1 - cT * cT);
	const double cP = cos(_phi[i]);
	const double sP = sin(_phi[i]);
	for (unsigned int j = 0; j <= i; ++j) {
	  TLorentzVector* daughter = &(_daughters[j]);
	  {  // rotate by theta around y-axis
	    const double x = daughter->Px();
	    const double z = daughter->Pz();
	    daughter->SetPz(cT * z - sT * x);
	    daughter->SetPx(sT * z + cT * x);
	  }
	  {  // rotate by phi around z-axis
	    const double x = daughter->Px();
	    const double y = daughter->Py();
	    daughter->SetPx(cP * x - sP * y);
	    daughter->SetPy(sP * x + cP * y);
	  }
	}
	if (i < (_n - 1)) {
	  // boost all daughters of the (i + 1)-body system from (i + 1)-body RF to (i + 2)-body RF; (i + 2)-body is along z-axis
	  const double betaGammaInv = _M[i] / _breakupMom[i + 1];  // 1 / (beta * gamma) of (i + 1)-body system in the (i + 2) RF
	  const double beta         = 1 / sqrt(1 + betaGammaInv * betaGammaInv);
	  for (unsigned int j = 0; j <= i; ++j)
	    _daughters[j].Boost(0, 0, beta);
	} else {
	  // final boost of all daughters from n-body RF to lab system
	  const TVector3 boost = nBody.BoostVector();
	  for (unsigned int j = 0; j < _n; ++j)
	    _daughters[j].Boost(boost);
	}
      }
    }
    break;
  case BLOCK:
    {
      // builds event starting in n-body RF going down to 2-body RF
      // is more efficicient, since it requitres only one rotation and boost per daughter
      TLorentzVector P = nBody;  // Lorentz of (i + 1)-body system in lab frame
      for (unsigned int i = _n - 1; i >= 1; --i) {  // loop from n-body down to 2-body
	// construct Lorentz vector of daughter _m[i] in (i + 1)-body RF
	const double    sinTheta = sqrt(1 - _cosTheta[i] * _cosTheta[i]);
	const double    pT       = _breakupMom[i] * sinTheta;
	TLorentzVector* daughter = &(_daughters[i]);
	daughter->SetPxPyPzE(pT * cos(_phi[i]),
			     pT * sin(_phi[i]),
			     _breakupMom[i] * _cosTheta[i],
			     sqrt(_m[i] * _m[i] + _breakupMom[i] * _breakupMom[i]));
	// boost daughter into lab frame
	daughter->Boost(P.BoostVector());
	// calculate Lorentz vector of i-body system in lab frame
	P -= *daughter;
      }
      // set last daughter
      _daughters[0] = P;
    }
    break;
  default:
    printErr << "unknown kinematics type. cannot calculate event kinematics." << endl;
    break;
  }
}



// calculates maximum weight for given n-body mass
double
nBodyPhaseSpaceGen::estimateMaxWeight(const double       nBodyMass,        // sic!
                                      const unsigned int nmbOfIterations)  // number of generated events

{
  double maxWeight = 0;
  for (unsigned int i = 0; i < nmbOfIterations; ++i) {
    _weight=1;
    pickMasses(nBodyMass);
    calcWeight();
    maxWeight = max(_weight, maxWeight);
  }
  return maxWeight;
}


ostream&
nBodyPhaseSpaceGen::print(ostream& out) const
{
  out << "nBodyPhaseSpaceGen parameters:" << endl
      << "    number of daughter particles ............... " << _n                 << endl;
  //  << "    masses of the daughter particles ........... " << _m                 << endl;
  //  << "    sums of daughter particle masses ........... " << _mSum              << endl
  //  << "    effective masses of (i + 1)-body systems ... " << _M                 << endl
  //  << "    cos(polar angle) in (i + 1)-body systems ... " << _cosTheta          << endl
  //  << "    azimuth in (i + 1)-body systems ............ " << _phi               << endl
  //   << "    breakup momenta in (i + 1)-body systems .... " << _breakupMom        << endl
  out  << "    weight formula ............................. " << _weightType        << endl
       << "    normalization value ........................ " << _norm              << endl
       << "    weight of generated event .................. " << _weight            << endl
       << "    maximum weight used in hit-miss MC ......... " << _maxWeight         << endl
       << "    maximum weight since instantiation ......... " << _maxWeightObserved << endl
       << "    algorithm for kinematics calculation ....... " << _kinematicsType    << endl
       << "    daughter four-momenta:" << endl;

  for (unsigned int i = 0; i < _n; ++i)
    out << "        daughter " << i << ": " << _daughters[i] << endl;
  return out;
}
