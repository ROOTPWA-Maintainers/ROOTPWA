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

// Base Class Headers ----------------


// Collaborating Class Headers -------
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include <iostream>
#include <vector>

// Collaborating Class Declarations --
class TH1;
class TGenPhaseSpace;

/** @brief Small helper class for bookkeeping
 */
class particleinfo {
public:
  particleinfo(int id,int q,double m):gid(id),charge(q),mass(m){}
  int gid;
  int charge;
  double mass;
};



/** @brief Phase Space generator for diffractive pion dissociation
 *  @author Sebastian Neubert TUM (original author)
 *  Based on TGenPhaseSpace
 *
 */
class TDiffractivePhaseSpace {
public:

  // Constructors/Destructors ---------
  TDiffractivePhaseSpace();
  ~TDiffractivePhaseSpace(){}

 // Accessors -----------------------
  TLorentzVector* GetDecay(unsigned int i){return phaseSpace.GetDecay(i);}
  TLorentzVector* GetBeam(){return &gbeam;}
  // Modifiers -----------------------
  /** @brief Set beam parameters
   * 
   * 2004 COMPASS beam:
   * gBeamDxDz      = 0.00026; // tilt from Quirin was in mrad
   * gBeamDxDzSigma = 0.00010;
   * gBeamDyDz      = 0.00001; // tilt from Quirin was in mrad
   * gBeamDyDzSigma = 0.00018;
   * gBeamMom       = 189 [GeV/c]
   * gBeamMomsigma  = 1.2
   */
  void SetBeam(double Mom=190, double MomSigma=1.2,
	       double DxDz=0, double DxDzSigma=0,
	       double DyDz=0, double DyDzSigma=0);

  /** @brief Set beam parameters
   *  The target is assumed to be a single cylinder
   *  centered at zPos and r=0;
   *
   *  Example: SetTarget(-20,40,2)
   *           Defines a target cell with 4cm diameter extending 
   *           from z=-40cm to z=0cm
   */
  void SetTarget(double zPos,double length, double r, double mass)
  {gTargetZPos=zPos;gTargetZLength=length;gTargetR=r;gRecoilMass=mass;}

  //void SetThetaDistribution(TH1* distr){thetaDistribution=distr;}


  /** @brief Set the slope b of the t-prime distribution
   *
   *  \f[ \frac{d\sigma}{dt'} \propto e^{-bt'} \f]
   */
  void SetTPrimeSlope(double b){gBT=1./b;} // inverse for simple usage with TRandom
  

  /** @brief Set mass range of produced system X
   * 
   *  Events will be generated uniformly in mass
   */
  void SetMassRange(double min, double max){xMassMin=min;xMassMax=max;}
  void SetDecayProducts(const std::vector<particleinfo>& info);
  void AddDecayProduct(const particleinfo& info);
  void SetSeed(int seed);
  

  // Operations ----------------------

  /** @brief generates on event
   * 
   * returns number of attempts to generate this event and beam
   * the decay products can be fetched with GetDecay(i)
   */
  unsigned int event(TLorentzVector& beam);

/** @brief generates on event
   * 
   * returns number of attempts to generate this event;
   * writes event to stream
   *
   */
  unsigned int event(ostream&);



private:

  // Private Data Members ------------
  TGenPhaseSpace phaseSpace;

  // target position
  double gTargetZPos;
  double gTargetZLength;
  double gTargetR;
  double gRecoilMass;

  // beam parameters:
  double gBeamMomSigma;  // [GeV/c]
  double gBeamMom;  // [GeV/c]

  double gBeamDxDz;
  double gBeamDxDzSigma;
  double gBeamDyDz;
  double gBeamDyDzSigma;

  TLorentzVector gbeam; // cache for last generated beam

  //TH1* thetaDistribution;

  double gBT;

  // cut on t-distribution
  double tMin; // [(GeV/c)^2]

  double xMassMin;
  double xMassMax;
  
  double* daughterMasses;
  std::vector<particleinfo> decayProducts;
  
  // Private Methods -----------------

  TLorentzVector makeBeam();
  bool writePwa2000Ascii(std::ostream&              out,
			 const TLorentzVector& beam,
			 TGenPhaseSpace&       event);
  
  std::ostream& progressIndicator(const long currentPos,
				  const long nmbTotal,
				  const int  nmbSteps = 10,
				  std::ostream& out   = std::cout);
  void  BuildDaughterList();
 // particle masses
  double gProtonMass;
  double gPionMass;
  double gPionMass2;


};

#endif
/* @} **/
