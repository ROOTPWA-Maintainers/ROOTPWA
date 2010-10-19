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
//      helper functions that convert between standard ASCII PWA2000
//      .evt files and the new ROOT tree format
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <string>
#include <vector>
#include <complex>

// cint has problems parsing Boost stuff
//#ifndef __CINT__
#if !defined (__CINT__) || !defined (__MAKECINT__)
#include "isobarAmplitude.h"
#endif


class TTree;
class TChain;


namespace rpwa {


	class isobarDecayTopology;
	class isobarAmplitude;


	std::string particleNameFromGeantId(const int id,
	                                    const int charge);


	void idAndChargeFromParticleName(const std::string& name,
	                                 int&               id,
	                                 int&               charge);


	double getParticleMass(const std::string& name);


	bool fillTreeFromEvt(std::istream&      inEvt,
	                     TTree&             outTree,
	                     const long int     maxNmbEvents              = -1,
	                     const std::string& prodKinParticlesLeafName  = "prodKinParticles",
	                     const std::string& prodKinMomentaLeafName    = "prodKinMomenta",
	                     const std::string& decayKinParticlesLeafName = "decayKinParticles",
	                     const std::string& decayKinMomentaLeafName   = "decayKinMomenta",
	                     const std::string& targetParticleName        = "p+",
	                     const bool         debug                     = false);
  

	bool writeEvtFromTree(TChain&            inTree,
	                      std::ostream&      outEvt,
	                      const long int     maxNmbEvents              = -1,
	                      const std::string& inTreeName                = "rootPwaEvtTree",
	                      const std::string& prodKinParticlesLeafName  = "prodKinParticles",
	                      const std::string& prodKinMomentaLeafName    = "prodKinMomenta",
	                      const std::string& decayKinParticlesLeafName = "decayKinParticles",
	                      const std::string& decayKinMomentaLeafName   = "decayKinMomenta",
	                      const bool         debug                     = false);


// cint has problems parsing Boost stuff
#if !defined (__CINT__) || !defined (__MAKECINT__)
	long int processTree(TTree&                              tree,
	                     const isobarAmplitudePtr&           amplitude,
	                     std::vector<std::complex<double> >& ampValues,
	                     const long int                      maxNmbEvents              = -1,
	                     const std::string&                  prodKinParticlesLeafName  = "prodKinParticles",
	                     const std::string&                  prodKinMomentaLeafName    = "prodKinMomenta",
	                     const std::string&                  decayKinParticlesLeafName = "decayKinParticles",
	                     const std::string&                  decayKinMomentaLeafName   = "decayKinMomenta",
	                     const bool                          printProgress             = true);
#endif


}  // namespace rpwa
