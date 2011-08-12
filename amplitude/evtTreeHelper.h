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


#if !defined (__CINT__) && !defined (ROOT_CINT)
// CINT and ACLiC have problems parsing Boost stuff; ROOT_CINT is set in rootlogon.C
#include "isobarAmplitude.h"
#endif


class TTree;
class TChain;
class TClonesArray;


namespace rpwa {


	class isobarDecayTopology;
	class isobarAmplitude;


	std::string particleNameFromGeantId(const int id,
	                                    const int charge);


	void idAndChargeFromParticleName(const std::string& name,
	                                 int&               id,
	                                 int&               charge);


	double getParticleMass(const std::string& name);


	bool getParticleNamesFromRootFile(TFile&             inFile,
	                                  TClonesArray*&     prodKinPartNames,   // array of particle names to be filled
	                                  TClonesArray*&     decayKinPartNames,  // array of particle names to be filled
	                                  const std::string& inTreeName                = "rootPwaEvtTree",
	                                  const std::string& prodKinParticlesLeafName  = "prodKinParticles",
	                                  const std::string& decayKinParticlesLeafName = "decayKinParticles");

	bool openRootEvtFiles(std::vector<TTree*>&            inTrees,            // array of trees from .root and .evt files
	                      TClonesArray*&                  prodKinPartNames,   // array of particle names to be filled
	                      TClonesArray*&                  decayKinPartNames,  // array of particle names to be filled
	                      const std::vector<std::string>& rootFileNames,      // .root files to be opened
	                      const std::vector<std::string>& evtFileNames,       // .evt files to be converted to trees
	                      const std::string&              inTreeName                = "rootPwaEvtTree",
	                      const std::string&              prodKinParticlesLeafName  = "prodKinParticles",
	                      const std::string&              prodKinMomentaLeafName    = "prodKinMomenta",
	                      const std::string&              decayKinParticlesLeafName = "decayKinParticles",
	                      const std::string&              decayKinMomentaLeafName   = "decayKinMomenta",
	                      const std::string&              targetParticleName        = "p+",
	                      const bool                      debug                     = false);


	bool fillTreeFromEvt(std::istream&      inEvt,
	                     TTree&             outTree,            // tree to be filled
	                     TClonesArray&      prodKinPartNames,   // array of particle names to be filled
	                     TClonesArray&      decayKinPartNames,  // array of particle names to be filled
	                     const long int     maxNmbEvents            = -1,
	                     const std::string& prodKinMomentaLeafName  = "prodKinMomenta",
	                     const std::string& decayKinMomentaLeafName = "decayKinMomenta",
	                     const std::string& targetParticleName      = "p+",
	                     const bool         debug                   = false,
	                     const long int     treeCacheSize           = 25000000);  // 25 MByte ROOT tree read cache


	bool writeEvtFromTree(TChain&            inTree,
	                      std::ostream&      outEvt,
	                      const long int     maxNmbEvents              = -1,
	                      const std::string& inTreeName                = "rootPwaEvtTree",
	                      const std::string& prodKinParticlesLeafName  = "prodKinParticles",
	                      const std::string& prodKinMomentaLeafName    = "prodKinMomenta",
	                      const std::string& decayKinParticlesLeafName = "decayKinParticles",
	                      const std::string& decayKinMomentaLeafName   = "decayKinMomenta",
	                      const bool         debug                     = false);


#if !defined (__CINT__) && !defined (ROOT_CINT)
	// CINT and ACLiC have problems parsing Boost stuff; ROOT_CINT is set in rootlogon.C
	bool processTree(TTree&                              tree,               // tree to be read
	                 const TClonesArray&                 prodKinPartNames,   // array of particle names
	                 const TClonesArray&                 decayKinPartNames,  // array of particle names
	                 const isobarAmplitudePtr&           amplitude,
	                 std::vector<std::complex<double> >& ampValues,
	                 const long int                      maxNmbEvents            = -1,
	                 const std::string&                  prodKinMomentaLeafName  = "prodKinMomenta",
	                 const std::string&                  decayKinMomentaLeafName = "decayKinMomenta",
	                 const bool                          printProgress           = true,
	                 const std::string&                  treePerfStatOutFileName = "",         // root file name for tree performance result
	                 const long int                      treeCacheSize           = 25000000);  // 25 MByte ROOT tree read cache
#endif


}  // namespace rpwa
