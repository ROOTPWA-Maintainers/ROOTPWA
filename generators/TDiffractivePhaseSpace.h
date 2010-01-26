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
#include <iostream>
#include <vector>

#include "TLorentzVector.h"

#include "nBodyPhaseSpaceGen.h"


// Collaborating Class Declarations --
class TH1;


namespace rpwa {


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
   *
   */
  class TDiffractivePhaseSpace {
  public:

    // Constructors/Destructors ---------
    TDiffractivePhaseSpace();
    ~TDiffractivePhaseSpace(){}

    // Accessors -----------------------
    const TLorentzVector* const GetDecay(unsigned int i){return &_phaseSpace.daughter(i);}
    TLorentzVector* GetBeam(){return &_beam;}
    // Modifiers -----------------------
    /** @brief Set beam parameters
     * 
     * 2004 COMPASS beam:
     * _beamDxDz      = 0.00026; // tilt from Quirin was in mrad
     * _beamDxDzSigma = 0.00010;
     * _beamDyDz      = 0.00001; // tilt from Quirin was in mrad
     * _beamDyDzSigma = 0.00018;
     * _beamMom       = 189 [GeV/c]
     * _beamMomSigma  = 1.2
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
    void SetTarget(double zPos,
		   double length,
		   double r,
		   double mass)
    {
      _targetZPos=zPos;
      _targetZLength=length;
      _targetR=r;
      _recoilMass=mass;
    }

    //void SetThetaDistribution(TH1* distr){thetaDistribution=distr;}


    /** @brief Set the slope b of the t-prime distribution
     *
     *  \f[ \frac{d\sigma}{dt'} \propto e^{-bt'} \f]
     */
    void SetTPrimeSlope(double b){_invSlopePar=1./b;} // inverse for simple usage with TRandom
  

    /** @brief Set mass range of produced system X
     * 
     *  Events will be generated uniformly in mass
     */
    void SetMassRange(double min,
		      double max)
    {
      _xMassMin=min;
      _xMassMax=max;
    }
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
    rpwa::nBodyPhaseSpaceGen _phaseSpace;

    // target position
    double _targetZPos;
    double _targetZLength;
    double _targetR;

    double _recoilMass;  // [GeV/c^2]

    // beam parameters:
    double _beamMomSigma;  // [GeV/c]
    double _beamMom;       // [GeV/c]

    double _beamDxDz;
    double _beamDxDzSigma;
    double _beamDyDz;
    double _beamDyDzSigma;

    TLorentzVector _beam; // cache for last generated beam

    //TH1* thetaDistribution;

    double _invSlopePar;  // inverse slope parameter 1 / b of t' distribution

    // cut on t-distribution
    double _tMin; // [(GeV/c)^2]

    double _xMassMin;
    double _xMassMax;
  
    double* _daughterMasses;
    std::vector<particleinfo> _decayProducts;
  
    // Private Methods -----------------

    TLorentzVector makeBeam();
    bool writePwa2000Ascii(std::ostream&             out,
			   const TLorentzVector&     beam,
			   rpwa::nBodyPhaseSpaceGen& event);
  
    std::ostream& progressIndicator(const long    currentPos,
				    const long    nmbTotal,
				    const int     nmbSteps = 10,
				    std::ostream& out   = std::cout);
    void  BuildDaughterList();
    // particle masses
    double _protonMass;
    double _pionMass;
    double _pionMass2;

  };

}  // namespace rpwa


#endif  // TDIFFRACTIVEPHASESPACE_HH
/* @} **/
