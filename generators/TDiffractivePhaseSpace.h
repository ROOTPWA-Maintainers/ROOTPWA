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

#include "TLorentzVector.h"

#include "nBodyPhaseSpaceGen.h"

#include "TPrimaryVertexGen.h"


class TH1;


namespace rpwa {


  /** @brief Small helper class for bookkeeping
   */
  class particleInfo {
  public:
    particleInfo(int    gId,
		 int    charge,
		 double mass)
      : _gId   (gId),
	_charge(charge),
	_mass  (mass)
    { }
    int    _gId;     // GEANT ID
    int    _charge;  // charge
    double _mass;    // mass [GeV/c^2]
  };


  /** @brief Phase Space generator for diffractive pion dissociation
   *  @author Sebastian Neubert TUM (original author)
   *
   */
  class TDiffractivePhaseSpace {
  public:

    // Constructors/Destructors ---------
    TDiffractivePhaseSpace();
    ~TDiffractivePhaseSpace();

    // Accessors -----------------------
    const TLorentzVector* const GetDecay(unsigned int i){return &_phaseSpace.daughter(i);}
    TLorentzVector* GetBeam(){return &_beamLab;}
    TVector3* GetVertex(){return &_vertex;}
    double Gettprime(){return _tprime;}

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
      _targetMass=mass;
      _recoilMass=mass;
    }

    //void SetThetaDistribution(TH1* distr){thetaDistribution=distr;}


    /** @brief Set the slope b of the t-prime distribution
     *
     *  \f[ \frac{d\sigma}{dt'} \propto e^{-bt'} \f]
     */
    void SetTPrimeSlope(double slopePar) {
    	_invSlopePar = new double[1];
    	_invSlopePar[0] = 1. / slopePar;
    }  // inverse for simple usage with TRandom
  
    /** @brief Set the slopes b of the t-prime distribution depending on the invariant mass
     *
     *  \f[ \frac{d\sigma}{dt'} \propto e^{-bt'} \f]
     *
     *  in case of more than one value assuming sorted ascending input for interpolation
     */
    void SetTPrimeSlope(double* slopePar, double* inv_m = NULL, int nvalues = 1) {
    	// delete previous arrays if existing
    	if (_invSlopePar) delete [] _invSlopePar;
    	if (_invM)	      delete [] _invM;
    	_invSlopePar = new double [nvalues];
    	_invM		 = new double [nvalues];
    	_ninvSlopePar = nvalues;
    	for (int i = 0; i < _ninvSlopePar; i++){
    		if (inv_m) _invM[i] = inv_m[i]; else _invM[i] = 0;
    		_invSlopePar[i] = 1. / slopePar[i];
    		//cout << _invM[i] << " " << _invSlopePar[i] << endl;
    	}
    }  // inverse for simple usage with TRandom

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
    void SetDecayProducts(const std::vector<particleInfo>& info);
    void AddDecayProduct(const particleInfo& info);
    void SetSeed(int seed);

  void setVerbose(bool flag){_phaseSpace.setVerbose(flag);}
  void SetImportanceBW(double mass, double width){
    _phaseSpace.setProposalBW(mass,width);
    _phaseSpace.setWeightType(nBodyPhaseSpaceGen::IMPORTANCE);
  }

  void SettMin(double tMin){_tMin = tMin;};

  /*
   * If you set the Primary Vertex Generator (create it first)
   * Vertex position, Beam Energy and Direction will be
   * created by the primary Vertex Generator
   */
  void SetPrimaryVertexGen(TPrimaryVertexGen* primaryVertexGen){_primaryVertexGen = primaryVertexGen;};



    /** @brief generates on event
     * 
     * returns number of attempts to generate this event and beam
     * the decay products can be fetched with GetDecay(i)
     */
    unsigned int event();

    /** @brief generates on event
     * 
     * returns number of attempts to generate this event;
     * writes event to stream
     *
     */
    unsigned int event(ostream&);

    /** @brief generates on event
     *
     * returns number of attempts to generate this event;
     * writes event to VES formatted stream and
     * for ComGeant fort.26 input files
     *
     */
    unsigned int event(ostream&, ostream&);

  double impWeight() const {return _phaseSpace.impWeight();}

  private:

    // Private Data Members ------------
    rpwa::nBodyPhaseSpaceGen _phaseSpace;

    TPrimaryVertexGen* _primaryVertexGen;


    // target position
    double _targetZPos;     // [cm]
    double _targetZLength;  // [cm]
    double _targetR;        // [cm]

    double _targetMass;  // [GeV/c^2]
    double _recoilMass;  // [GeV/c^2]

    // beam parameters:
    double _beamMomSigma;  // [GeV/c]
    double _beamMom;       // [GeV/c]

    double _beamDxDz;
    double _beamDxDzSigma;
    double _beamDyDz;
    double _beamDyDzSigma;

    TLorentzVector _beamLab; // cache for last generated beam (in lab frame)
    TLorentzVector _recoilprotonLab; // cache for last generated recoil proton (in lab frame)
    TVector3 _vertex; 		 // cache for last generated vertex
    double _tprime;			 // cache for last generated t' (recalculated)

    //TH1* thetaDistribution;

    double* _invSlopePar;  // inverse slope parameter(s) 1 / b of t' distribution [(GeV/c)^{-2}]
    double* _invM;		   // invariant masses corresponding to _invSlopePar
    int		_ninvSlopePar;

    // cut on t-distribution
    double _tMin;  // [(GeV/c)^2]

    double _xMassMin;  // [GeV/c^2]
    double _xMassMax;  // [GeV/c^2]
  
    std::vector<particleInfo> _decayProducts;
  
    // Private Methods -----------------

    TLorentzVector makeBeam();
    bool writePwa2000Ascii(std::ostream& out,
			   const int     beamGeantId,
			   const int     beamCharge);

    // writes event to ascii file read by ComGeant fort.26 interface
    bool writeComGeantAscii(
    		ostream& out,
    		bool  formated = true); // true: text file ; false: binary file (not implemented yet)
  
    void  BuildDaughterList();
    // particle masses
    double _protonMass;
    double _pionMass;
    double _pionMass2;

    // calculate the t' by using the information of the incoming and outgoing particle in the vertex
    float Calc_t_prime(const TLorentzVector& particle_In, const TLorentzVector& particle_Out);


    // case invariant_M < 0. : first entry in _invSlopePar is taken (if any)
    // case invariant_M >= 0.: extrapolation or interpolation of given points
    // for t' over invariant Mass
    double Get_inv_SlopePar(double invariant_M = -1.);
  };

}  // namespace rpwa


#endif  // TDIFFRACTIVEPHASESPACE_HH
/* @} **/
