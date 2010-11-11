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
//      class that reads/writes wave description from/to keyfiles,
//      constructs decay topologies and amplitudes
//      can be saved into .root files
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef WAVEDESCRIPTION_H
#define WAVEDESCRIPTION_H


#include <string>
#include <vector>

#include "TObject.h"

#if !defined (__CINT__) || !defined (__MAKECINT__)
#include "libconfig.h++"

#include "isobarDecayTopology.h"
#include "isobarAmplitude.h"
#endif


namespace rpwa {


	class waveDescription : public TObject {

  public:

	  waveDescription ();
    virtual ~waveDescription();

#if !defined (__CINT__) || !defined (__MAKECINT__)

		// construction of decay topology and amplitude objects
		bool parseKeyFile(const std::string& keyFileName);  ///< parses key file
		std::string keyFileContents() const { return _keyFileContents; }  ///< returns contents of key file
		bool constructDecayTopology(isobarDecayTopologyPtr& topo,
		                            const bool              requireXQnKey = true) const;  ///< construct isobar decay topology from keyfile
	  bool constructAmplitude(isobarAmplitudePtr& amplitude) const;   ///< construct isobar decay amplitude from keyfile
	  bool constructAmplitude(isobarAmplitudePtr&           amplitude,
	                          const isobarDecayTopologyPtr& topo) const;  ///< construct isobar amplitude using existing decay topology

	  static bool writeKeyFile(const std::string&         keyFileName,
	                           const isobarDecayTopology& topo,
	                           const bool                 writeProdVert = true);  ///< creates key file from decay topology
	  static bool writeKeyFile(const std::string&         keyFileName,
	                           const isobarAmplitude&     amplitude,
	                           const bool                 writeProdVert = true);  ///< creates key file from amplitude

	  static std::string waveNameFromTopology
	  (const isobarDecayTopology&  topo,
	   const isobarDecayVertexPtr& currentVertex = isobarDecayVertexPtr());  ///< recursive function that generates unique wave name from decay topology
	  static std::string waveNameFromTopologyOld
	  (isobarDecayTopology&        topo,
	   const isobarDecayVertexPtr& currentVertex = isobarDecayVertexPtr());  ///< recursive function that generates unique wave name from decay topology (old convention)

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


  private:

		void parseKeyString();  ///< parses _keyFileContents string

		// helper functions for construction of decay topology and ampltiude
		static bool constructXParticle(const libconfig::Setting& XQnKey,
		                               particlePtr&              X);  ///< creates X particle with quantum numbers defined in X key
    static productionVertexPtr mapProductionVertexType(const libconfig::Setting& prodVertKey,
                                                       const std::string&        vertType,
                                                       const particlePtr&        X);  ///< creates production vertex according to given type
		static bool constructProductionVertex(const libconfig::Setting& rootKey,
		                                      const particlePtr&        X,
		                                      productionVertexPtr&      prodVert);  ///< creates production vertex
    static bool constructParticle(const libconfig::Setting& particleKey,
                                  particlePtr&              particle);  ///< creates particle using name in particle key
    static massDependencePtr mapMassDependenceType(const std::string& massDepType);  ///< creates mass dependence functor of specified type
		static bool constructDecayVertex(const libconfig::Setting&          parentKey,
		                                 const particlePtr&                 parentParticle,
		                                 std::vector<isobarDecayVertexPtr>& decayVertices,
		                                 std::vector<particlePtr>&          fsParticles);  ///< recursively traverses decay chain and creates decay vertices and final state particles
	  static isobarAmplitudePtr mapAmplitudeType(const std::string&            formalismType,
	                                             const isobarDecayTopologyPtr& topo);  ///< creates amplitude for specified formalism

		// helper functions for writing key files from decay topology and ampltiude
    static bool setProductionVertexKeys(libconfig::Setting&        prodVertKey,
                                        const productionVertexPtr& prodVert);  ///< puts production vertex info into keys
	  static bool setXQuantumNumbersKeys(libconfig::Setting& XQnKey,
	                                     const particle&     X);  ///< puts X quantum numbers into keys
		static bool setMassDependence(libconfig::Setting&   isobarMassDepKey,
		                              const massDependence& massDep);  ///< puts mass dependence into key
	  static bool setXDecayKeys(libconfig::Setting&        parentDecayKey,
	                            const isobarDecayTopology& topo,
	                            const isobarDecayVertex&   vert);  ///< recursive function that puts X decay chain into keys
	  static bool writeKeyFile(libconfig::Setting&        rootKey,
	                           const isobarDecayTopology& topo,
	                           const bool                 writeProdVert = true);  ///< writes decay topology to key file
	  static bool setAmplitude(libconfig::Setting&    amplitudeKey,
	                           const isobarAmplitude& amplitude);  ///< puts amplitude specification into key
	  static bool writeKeyFile(libconfig::Setting&        rootKey,
	                           const isobarAmplitude&     amplitude);  ///< writes amplitude parameters to key file

		libconfig::Config*        _key;           ///< libConfig date structure constructed from key file
		const libconfig::Setting* _amplitudeKey;  ///< pointer to amplitude options

	  static bool _debug;  ///< if set to true, debug messages are printed

#endif  // __CINT__

		std::string _keyFileContents;  ///< copy of keyfile contents; can be written to .root file
		
		ClassDef(waveDescription,1)

  };


}  // namespace rpwa


#endif  // WAVEDESCRIPTION_H
