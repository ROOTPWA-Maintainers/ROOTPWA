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
#include <map>

#include "TObject.h"

#ifndef __CINT__
#include "libconfig.h++"

#include "isobarDecayTopology.h"
#include "isobarAmplitude.h"
#else
namespace libconfig {
	class Config;
	class Setting;
}
#endif


namespace rpwa {

	class amplitudeMetadata;

	class waveDescription : public TObject {

	public:

		waveDescription();
		waveDescription(const amplitudeMetadata* amplitudeMeta);
		virtual ~waveDescription();

		void clear();

		waveDescription& operator =(const waveDescription& waveDesc);

#ifndef __CINT__

		// construction of decay topology and amplitude objects
		bool parseKeyFile       (const std::string& keyFileName   );    ///< parses key file
		bool parseKeyFileContent(const std::string& keyFileContent);    ///< parses key file
		bool keyFileParsed() const { return _keyFileParsed; }  ///< returns whether key file was successfully parsed

		std::string   keyFileContent() const { return _keyFileLocalCopy; }  ///< returns content of key file
		std::ostream& printKeyFileContent(std::ostream&      out,
		                                  const std::string& keyFileContent = "") const;  ///< prints key file content string with line numbers
		bool constructDecayTopology(isobarDecayTopologyPtr& topo,
		                            const bool              fromTemplate = false) const;  ///< construct isobar decay topology from keyfile
		bool constructAmplitude(isobarAmplitudePtr& amplitude) const;  ///< construct isobar decay amplitude from keyfile
		bool constructAmplitude(isobarAmplitudePtr&           amplitude,
		                        const isobarDecayTopologyPtr& topo) const;  ///< construct isobar amplitude using existing decay topology

		template<class T>
		static bool writeKeyFile(std::ostream& out,
		                         const T&      topoOrAmp,
		                         const bool    writeProdVert = true);  ///< writes keys from decay topology or amplitude to stream
		template<class T>
		static bool writeKeyFile(const std::string& keyFileName,
		                         const T&           topoOrAmp,
		                         const bool         writeProdVert = true);  ///< creates key file from decay topology or amplitude

		static std::string waveNameFromTopology
		(isobarDecayTopology         topo,
		 const bool                  newConvention = false,
		 const isobarDecayVertexPtr& currentVertex = isobarDecayVertexPtr());  ///< recursive function that generates unique wave name from decay topology

		static std::string waveLaTeXFromTopology
		(isobarDecayTopology         topo,
		 const isobarDecayVertexPtr& currentVertex = isobarDecayVertexPtr());  ///< recursive function that generates unique wave name from decay topology

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	private:

		bool readKeyFileIntoLocalCopy(const std::string& keyFileName);  ///< reads key file content into _keyFileLocalCopy string

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
		                              particlePtr&              particle,
		                              const bool                requirePartInTable = true);  ///< creates particle using name in particle key
		static massDependencePtr mapMassDependenceType(const std::string& massDepType);  ///< creates mass dependence functor of specified type
		static bool constructDecayVertex(const libconfig::Setting&          parentKey,
		                                 const particlePtr&                 parentParticle,
		                                 std::vector<isobarDecayVertexPtr>& decayVertices,
		                                 std::vector<particlePtr>&          fsParticles,
		                                 const bool                         fromTemplate = false);  ///< recursively traverses decay chain and creates decay vertices and final state particles
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
		static bool setKeysFromTopology(libconfig::Setting&        rootKey,
		                                const isobarDecayTopology& topo,
		                                const bool                 setProdVert = true);  ///< fills keys from decay topology
		static bool setAmplitude(libconfig::Setting&    amplitudeKey,
		                         const isobarAmplitude& amplitude);  ///< puts amplitude specification into key
		static bool setKeysFromAmplitude(libconfig::Setting&    rootKey,
		                                 const isobarAmplitude& amplitude,
		                                 const bool             setProdVert = true);  ///< fills keys from amplitude
		static bool writeKeyFile(FILE&                      outStream,
		                         const isobarDecayTopology& topo,
		                         const bool                 writeProdVert = true);  ///< creates key file from decay topology and writes it to output stream
		static bool writeKeyFile(FILE&                      outStream,
		                         const isobarAmplitude&     amplitude,
		                         const bool                 writeProdVert = true);  ///< creates key file from amplitude and writes it to output stream

#endif  // __CINT__

		libconfig::Config* _key;            //! ///< libConfig date structure constructed from key file
		bool               _keyFileParsed;  //! ///< indicates whether key file was successfully parsed

		std::string _keyFileLocalCopy;  ///< copy of keyfile content; is written to .root file

		static bool _debug;  ///< if set to true, debug messages are printed
		static std::map<std::string, std::string> isobars; ///< LaTeX names of isobars


		ClassDef(waveDescription,2)

	};


	template<class T>
	inline
	bool
	waveDescription::writeKeyFile(std::ostream& out,
	                              const T&      topoOrAmp,
	                              const bool    writeProdVert)
	{
		// create pipe
		int pipeFileDescriptors[2];
		if (pipe(pipeFileDescriptors) == -1) {
			printErr << "failed to create pipe. cannot write keys." << std::endl;
			return false;
		}
		// open write end of pipe
		FILE* pipeWriteEnd = fdopen(pipeFileDescriptors[1], "wt");
		if (!pipeWriteEnd) {
			printErr << "could not open write end of pipe. cannot write keys." << std::endl;
			return false;
		}
		// write keys to pipe and close write end
		if (not writeKeyFile(*pipeWriteEnd, topoOrAmp, writeProdVert)) {
			printWarn << "problems writing keys for decay topology. cannot write keys." << std::endl;
			fclose(pipeWriteEnd);
			return false;
		}
		fclose(pipeWriteEnd);
		// read keys from pipe
		char         buf;
		unsigned int countChar = 0;
		while (read(pipeFileDescriptors[0], &buf, 1) > 0) {
			out << buf;
			++countChar;
		}
		close(pipeFileDescriptors[0]);
		if (countChar > 0)
			return true;
		else {
			printWarn << "nothing was written" << std::endl;
			return false;
		}
	}


	template<class T>
	inline
	bool
	waveDescription::writeKeyFile(const std::string& keyFileName,
	                              const T&           topoOrAmp,
	                              const bool         writeProdVert)
	{
		if (_debug)
			printDebug << "writing key file '" << keyFileName << "'" << std::endl;
		std::ofstream outFile(keyFileName.c_str());
		if (not outFile) {
			printErr << "cannot create key file '" << keyFileName << "'" << std::endl;
			return false;
		}
		if (writeKeyFile(outFile, topoOrAmp, writeProdVert)) {
			outFile.close();
			return true;
		} else {
			printWarn << "problems writing keys for decay topology. cannot write key file." << std::endl;
			outFile.close();
			return false;
		}
	}


}  // namespace rpwa


#endif  // WAVEDESCRIPTION_H
