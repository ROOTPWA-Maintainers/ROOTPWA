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


#include <fstream>
#include <map>

#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

#include "TClass.h"

#include "libConfigUtils.hpp"
#include "conversionUtils.hpp"
#include "diffractiveDissVertex.h"
#include "leptoProductionVertex.h"
#include "isobarHelicityAmplitude.h"
#include "isobarCanonicalAmplitude.h"
#include "waveDescription.h"

  
using namespace std;
using namespace boost;
using namespace libconfig;
using namespace rpwa;


ClassImp(waveDescription);


bool waveDescription::_debug = false;


waveDescription::waveDescription()
	: TObject         (),
	  _key            (0),
	  _amplitudeKey   (0),
	  _keyFileContents("")
{
	waveDescription::Class()->IgnoreTObjectStreamer();  // don't store TObject's fBits and fUniqueID
}


waveDescription::~waveDescription()
{
	if (_key)
		delete _key;
}


bool
waveDescription::parseKeyFile(const string& keyFileName)
{
	if (not _key)
		_key = new Config();
	_amplitudeKey = 0;
	if (not parseLibConfigFile(keyFileName, *_key, _debug)) {
		printWarn << "problems reading key file '" << keyFileName << "'. "
		          << "cannot construct decay topology." << endl;
		delete _key;
		_key = 0;
		return false;
	}		
	// read key file contents into string
	{
		ifstream keyFile(keyFileName.c_str());
    if (not keyFile or not keyFile.good()) {
	    printWarn << "cannot read from file '"  << keyFileName << "'" << endl;
	    return false;
    }
		_keyFileContents = "";
    string line;
    while(getline(keyFile, line))
	    _keyFileContents += line + "\n";
	}
	// find amplitude group
	const Setting& rootKey = _key->getRoot();
	_amplitudeKey = findLibConfigGroup(rootKey, "amplitude", false);
	return true;
}


ostream&
waveDescription::printKeyFileContents(ostream& out) const
{
	if (_keyFileContents != "") {
		typedef tokenizer<char_separator<char> > tokenizer;
		char_separator<char> separator("\n");
		tokenizer            keyFileLines(_keyFileContents, separator);
		unsigned int         lineNumber = 0;
		for (tokenizer::iterator i = keyFileLines.begin(); i != keyFileLines.end(); ++i)
			out << setw(5) << ++lineNumber << "  " << *i << endl;
	} else
		out << "key file contents string is empty" << endl;
	return out;
}


bool
waveDescription::constructDecayTopology(isobarDecayTopologyPtr& topo,
                                        const bool              fromTemplate) const
{
	if (not _key) {
		printWarn << "null pointer to libconfig data structure. "
		          << "cannot construct decay topology. was parsing successful?" << endl;
	}

	if (topo)
		topo.reset();
	topo = isobarDecayTopologyPtr();  // null pointer

	const Setting& rootKey = _key->getRoot();
	printInfo << "constructing decay topology from key file" << endl;

	// find wave group
	const Setting* waveKey = findLibConfigGroup(rootKey, "wave");
	if (not waveKey) {
		printWarn << "cannot find 'wave' group. cannot construct decay topology." << endl;
		return false;
	}

	// find  X quantum numbers group
	const Setting* XQnKey = findLibConfigGroup(*waveKey, "XQuantumNumbers", not fromTemplate);
	particlePtr    X;
	if (not XQnKey)
		if (not fromTemplate) {
			printWarn << "cannot find 'XQuantumNumbers' group. cannot construct decay topology." << endl;
			return false;
		} else {
			// create default X particle
			if (_debug)
				printDebug << "cannot find 'XQuantumNumbers' group. "
				           << "creating X particle with default properties." << endl;
			X = createParticle("X", not fromTemplate);
		}
	else
		// create X particle from key
		if (not constructXParticle(*XQnKey, X)) {
			printWarn << "problems constructing X particle. cannot construct decay topology." << endl;
			return false;
		}
  
	// create production vertex
	productionVertexPtr prodVert = productionVertexPtr();
	if (not constructProductionVertex(rootKey, X, prodVert)) {
		printWarn << "problems constructing production vertex. cannot construct decay topology." << endl;
		return false;
	}
  
	// find X decay group
	const Setting* XDecayKey = findLibConfigGroup(*waveKey, "XDecay");
	if (not XDecayKey) {
		printWarn << "cannot find 'XDecay' group. cannot construct decay topology." << endl;
		return false;
	}
	// traverse decay chain and create final state particles and isobar decay vertices
	vector<isobarDecayVertexPtr> decayVertices;
	vector<particlePtr>          fsParticles;
	if (not constructDecayVertex(*XDecayKey, X, decayVertices, fsParticles, fromTemplate)) {
		printWarn << "problems constructing decay chain. cannot construct decay topology." << endl;
		return false;
	}

	// create decay isobar decay topology
	topo = createIsobarDecayTopology(prodVert, decayVertices, fsParticles);
	// backward compatibility: allow sloppy key files, where charges of
	// isobars are not explicitely defined
	topo->calcIsobarCharges();
	//!!! user should correctly define quantum numbers
	//topo->calcIsobarBaryonNmbs();
	//topo->productionVertex()->setXFlavorQN();  // sets baryon nmb, S, C, and B of X
  
	printSucc << "constructed decay topology from key file" << endl;
	return true;
}


bool
waveDescription::constructAmplitude(isobarAmplitudePtr& amplitude) const
{
	isobarDecayTopologyPtr topo;
	if (not constructDecayTopology(topo)) {
		printWarn << "problems constructing decay topology. cannot construct decay amplitude." << endl;
		return false;
	}
	return constructAmplitude(amplitude, topo);
}


bool
waveDescription::constructAmplitude(isobarAmplitudePtr&           amplitude,
                                    const isobarDecayTopologyPtr& topo) const
{
	if (not topo) {
		printWarn << "null pointer to decay topology. cannot construct decay amplitude." << endl;
		return false;
	}
	if (not topo->checkTopology() or not topo->checkConsistency()) {
		printWarn << "decay topology has issues. cannot construct decay amplitude." << endl;
		return false;
	}
	if (amplitude)
		amplitude.reset();
	// get amplitude parameters from key file
	// default values
	string formalism            = "helicity";
	bool   boseSymmetrize       = true;
	bool   useReflectivityBasis = true;
	if (_amplitudeKey) {
		if (_amplitudeKey->lookupValue("formalism", formalism) and _debug)
			printDebug << "setting amplitude formalism to '" << formalism << "'" << endl;
		if (_amplitudeKey->lookupValue("boseSymmetrize", boseSymmetrize) and _debug)
			printDebug << "setting amplitude option 'boseSymmetrize' to "
			           << ((boseSymmetrize) ? "true" : "false") << endl;
		if (_amplitudeKey->lookupValue("useReflectivityBasis", useReflectivityBasis) and _debug)
			printDebug << "setting amplitude option 'useReflectivityBasis' to "
			           << ((useReflectivityBasis) ? "true" : "false") << endl;
	}
	// construct amplitude
	amplitude = mapAmplitudeType(formalism, topo);
	if (_debug)
		printDebug << "constructed amplitude '"<< amplitude->name() << "': "
		           << ((boseSymmetrize      ) ? "en" : "dis") << "abled Bose symmetrization, "
		           << ((useReflectivityBasis) ? "en" : "dis") << "abled reflectivity basis" << endl;
	if (not amplitude) {
		printWarn << "problems constructing decay amplitude." << endl;
		return false;
	}
	amplitude->enableBoseSymmetrization(boseSymmetrize      );
	amplitude->enableReflectivityBasis (useReflectivityBasis);
	return true;
}


bool
waveDescription::writeKeyFile(const string&              keyFileName,
                              const isobarDecayTopology& topo,
                              const bool                 writeProdVert)
{
	if (_debug)
		printDebug << "writing key file '" << keyFileName << "'" << endl;
	Config   key;
	Setting& rootKey = key.getRoot();
	if (not writeKeyFile(rootKey, topo, writeProdVert)) {
		printWarn << "problems writing keys for decay topology. cannot write key file." << endl;
		return false;
	}
	try {
		key.writeFile(keyFileName.c_str());
	} catch (const FileIOException& ioEx) {
		printWarn << "I/O error while writing key file '" << keyFileName << "'" << endl;
		return false;
	}
	printSucc << "written key file '" << keyFileName << "'" << endl;
	return true;
}


bool
waveDescription::writeKeyFile(const string&          keyFileName,
                              const isobarAmplitude& amplitude,
                              const bool             writeProdVert)
{
	if (_debug)
		printDebug << "writing key file '" << keyFileName << "'" << endl;
	Config   key;
	Setting& rootKey = key.getRoot();
	if (not writeKeyFile(rootKey, *amplitude.decayTopology(), writeProdVert)) {
		printWarn << "problems writing keys for decay topology. cannot write key file." << endl;
		return false;
	}
	if (not writeKeyFile(rootKey, amplitude)) {
		printWarn << "problems writing keys for amplitude. cannot write key file." << endl;
		return false;
	}
	try {
		key.writeFile(keyFileName.c_str());
	} catch (const FileIOException& ioEx) {
		printWarn << "I/O error while writing key file '" << keyFileName << "'" << endl;
		return false;
	}
	printSucc << "written key file '" << keyFileName << "'" << endl;
	return true;
}


string
waveDescription::waveNameFromTopology(const isobarDecayTopology&  topo,
                                      const isobarDecayVertexPtr& currentVertex)
{
	ostringstream fileName;
	if (currentVertex == interactionVertexPtr()) {
		if (not topo.checkTopology() or not topo.checkConsistency()) {
			printWarn << "decay topology has issues. cannot construct key file name." << endl;
			return "";
		}
		// X quantum numbers
		const particle& X = *(topo.XParticle());
		fileName << spinQn(X.isospin()) << parityQn(X.G())
		         << spinQn(X.J()) << parityQn(X.P()) << parityQn(X.C())
		         << spinQn(X.spinProj()) << parityQn(X.reflectivity());
		// start traversing down decay chain
		fileName << waveNameFromTopology(topo, topo.XIsobarDecayVertex());
	} else {
		// recurse down decay chain
		// first daughter
		fileName << "=[" << currentVertex->daughter1()->name();
		if (not topo.isFsParticle(currentVertex->daughter1()))
			fileName << waveNameFromTopology
				(topo, static_pointer_cast<isobarDecayVertex>(topo.toVertex(currentVertex->daughter1())));
		// L, S
		fileName << "#" << spinQn(currentVertex->L()) << "," << spinQn(currentVertex->S()) << "#";
		// second daughter
		fileName << currentVertex->daughter2()->name();
		if (not topo.isFsParticle(currentVertex->daughter2()))
			fileName << waveNameFromTopology
				(topo, static_pointer_cast<isobarDecayVertex>(topo.toVertex(currentVertex->daughter2())));
		fileName << "]";
	}
	string keyFileName = fileName.str();
	replace_all(keyFileName, "(", "_");
	replace_all(keyFileName, ")", "_");
	return keyFileName;
}


string
waveDescription::waveNameFromTopologyOld(isobarDecayTopology&        topo,
                                         const isobarDecayVertexPtr& currentVertex)
{
	ostringstream fileName;
	if (currentVertex == interactionVertexPtr()) {  // null pointer
		if (not topo.checkTopology() or not topo.checkConsistency()) {
			printWarn << "decay topology has issues. cannot construct key file name." << endl;
			return "";
		}
		// X quantum numbers and first decay
		const particle& X = *(topo.XParticle());
		fileName << spinQn(X.isospin()) << parityQn(X.G())
		         << spinQn(X.J()) << parityQn(X.P()) << parityQn(X.C())
		         << spinQn(X.spinProj()) << parityQn(X.reflectivity())
		         << waveNameFromTopologyOld(topo, static_pointer_cast<isobarDecayVertex>
		                                    (topo.toVertex(topo.XIsobarDecayVertex()->daughter1())))
		         << "_" << spinQn(topo.XIsobarDecayVertex()->L())
		         << spinQn(topo.XIsobarDecayVertex()->S()) << "_"
		         << waveNameFromTopologyOld(topo, static_pointer_cast<isobarDecayVertex>
		                                    (topo.toVertex(topo.XIsobarDecayVertex()->daughter2())));
	} else {
		// recurse down decay chain
		fileName << ((currentVertex->parent()->charge() != 0) ? currentVertex->parent()->name()
		             : currentVertex->parent()->bareName());
		isobarDecayTopology subGraph = topo.subDecay(topo.node(currentVertex));
		if (not topo.isFsVertex(currentVertex) and subGraph.nmbFsParticles() > 2)
			fileName << "="
			         << waveNameFromTopologyOld(topo, static_pointer_cast<isobarDecayVertex>
			                                    (topo.toVertex(currentVertex->daughter1())))
			         << "_" << spinQn(currentVertex->L()) << spinQn(currentVertex->S()) << "_"
			         << waveNameFromTopologyOld(topo, static_pointer_cast<isobarDecayVertex>
			                                    (topo.toVertex(currentVertex->daughter2())));
	}
	string keyFileName = fileName.str();
	replace_all(keyFileName, "(", "");
	replace_all(keyFileName, ")", "");
	return keyFileName;
}


void
waveDescription::parseKeyString()
{
	if (not _key)
		_key = new Config();
	_amplitudeKey = 0;
	if (not parseLibConfigString(_keyFileContents, *_key, _debug)) {
		printWarn << "problems parsing key file string:" << endl;
		printKeyFileContents(cout);
		cout  << "    cannot construct decay topology." << endl;
		delete _key;
		_key = 0;
		return;
	}
	// find amplitude group
	const Setting& rootKey = _key->getRoot();
	_amplitudeKey = findLibConfigGroup(rootKey, "amplitude", false);
}


bool
waveDescription::constructXParticle(const Setting& XQnKey,
                                    particlePtr&   X)
{
	if (_debug)
		printDebug << "reading X quantum numbers from '" << XQnKey.getPath() << "':" << endl;
	// get X quantum numbers
	map<string, int> mandatoryXQn, optionalXQn;
	mandatoryXQn["isospin"    ];
	optionalXQn ["G"          ];
	mandatoryXQn["J"          ];
	mandatoryXQn["P"          ];
	optionalXQn ["C"          ];
	mandatoryXQn["M"          ];
	optionalXQn ["refl"       ];
	optionalXQn ["baryonNmb"  ];
	optionalXQn ["strangeness"];
	optionalXQn ["charm"      ];
	optionalXQn ["beauty"     ];
	bool success = true;
	for (map<string, int>::iterator i = mandatoryXQn.begin(); i != mandatoryXQn.end(); ++i)
		if (not XQnKey.lookupValue(i->first, i->second)) {
			printWarn << "cannot find integer field '" << i->first << "in "
			          << "'" << XQnKey.getPath() << "'" << endl;
			success = false;
		}
	for (map<string, int>::iterator i = optionalXQn.begin(); i != optionalXQn.end(); ++i)
		if (not XQnKey.lookupValue(i->first, i->second))
			i->second = 0;
	if (not success) {
		printWarn << "cannot find X quantum numbers. cannot construct decay topology." << endl;
		return false;
	}
	// create X particle
	X = createParticle("X",
	                   mandatoryXQn["isospin"], optionalXQn["G"],
	                   mandatoryXQn["J"], mandatoryXQn["P"], optionalXQn["C"],
	                   mandatoryXQn["M"], optionalXQn["refl"]);
	X->setBaryonNmb(mandatoryXQn["baryonNmb"]);
	X->setSCB(optionalXQn["strangeness"], optionalXQn["charm"], optionalXQn["beauty"]);
	if (_debug)
		printDebug << "constructed X particle: " << X->qnSummary() << endl;
	return true;
}


productionVertexPtr
waveDescription::mapProductionVertexType(const Setting&       prodVertKey,
                                         const string&        vertType,
                                         const particlePtr&   X)
{
	productionVertexPtr prodVert;
	bool                success = true;
	if ((vertType == "diffractiveDissVertex") or (vertType == "leptoProductionVertex")) {
		const string             prodKinParticleNames[] = {"beam", "target", "recoil"};
		map<string, particlePtr> prodKinParticles;
		for (unsigned int i = 0; i < sizeof(prodKinParticleNames) / sizeof(string); ++i)
			prodKinParticles[prodKinParticleNames[i]] = particlePtr();
		// construct particles
		const Setting* beamParticleKey = 0;  // needed for leptoproduction vertex
		for (map<string, particlePtr>::iterator i = prodKinParticles.begin();
		     i != prodKinParticles.end(); ++i) {
			const Setting* particleKey = findLibConfigGroup(prodVertKey, i->first,
			                                                (i->first != "recoil") ? true : false);
			if (particleKey) {
				if (i->first == "beam")
					beamParticleKey = particleKey;
				if (not constructParticle(*particleKey, i->second))
					success = false;
			} else if (i->first != "recoil")
				success = false;
		}
		if (success) {
			if (vertType == "diffractiveDissVertex")
				prodVert = createDiffractiveDissVertex(prodKinParticles["beam"], prodKinParticles["target"],
				                                       X, prodKinParticles["recoil"]);
			else if (vertType == "leptoProductionVertex" ) {
				prodVert = createLeptoProductionVertex(prodKinParticles["beam"], prodKinParticles["target"],
				                                       X, prodKinParticles["recoil"]);
				double beamLongPol = 0;
				if ((beamParticleKey->lookupValue("longPol", beamLongPol))) {
					if (_debug)
						printDebug << "setting polarization of beam " << prodKinParticles["beam"]->qnSummary()
						           << " to " << beamLongPol << endl;
				} else
					printWarn << "no polarization is given for beam " << prodKinParticles["beam"]->qnSummary()
					          << ". assuming unpolarized beam." << endl;
				static_pointer_cast<leptoProductionVertex>(prodVert)->setBeamPol(beamLongPol);
			}
		}
	} else
		printWarn << "unknown production vertex type '" << vertType << "'" << endl;
	if (_debug) {
		if (success and prodVert)
			printDebug << "constructed " << *prodVert << endl;
		else
			printWarn << "problems constructing production vertex of type '" << vertType << "'" << endl;
	}
	return prodVert;
}


bool
waveDescription::constructProductionVertex(const Setting&       rootKey,
                                           const particlePtr&   X,
                                           productionVertexPtr& prodVert)
{
	// find production vertex group
	const Setting* prodVertKey = findLibConfigGroup(rootKey, "productionVertex");
	if (not prodVertKey) {
		printWarn << "cannot find 'productionVertex' group" << endl;
		return false;
	}
	if (_debug)
		printDebug << "reading production vertex from '" << prodVertKey->getPath() << "':" << endl;
	bool success = true;
	// get vertex type
	string vertType;
	if (not prodVertKey->lookupValue("type", vertType)) {
		printWarn << "cannot find 'type' entry in '" << prodVertKey->getPath() << "'" << endl;
		success = false;
	}
	// create production vertex
	prodVert = mapProductionVertexType(*prodVertKey, vertType, X);
	return (prodVert) ? true : false;
}


bool
waveDescription::constructParticle(const Setting& particleKey,
                                   particlePtr&   particle,
                                   const bool     requirePartInTable)
{
	if (_debug)
		printDebug << "reading particle information from '" << particleKey.getPath() << "':" << endl;
	string name;
	if (particleKey.lookupValue("name", name)) {
		particle = createParticle(name, requirePartInTable);
		int spinProj;
		if (particleKey.lookupValue("spinProj", spinProj))
			particle->setSpinProj(spinProj);
		if (_debug)
			printDebug << "created particle " << particle->qnSummary() << flush;
		int index;
		if (particleKey.lookupValue("index", index)) {
			particle->setIndex(index);
			if (_debug)
				cout << " with index [" << index << "]" << flush;
		}
		if (_debug)
			cout << endl;
		return true;
	} else {
		printWarn << "cannot construct particle from entry '" << particleKey.getPath() << "'. "
		          << "name field is missing." << endl;
		particle = particlePtr();
		return false;
	}
}


massDependencePtr
waveDescription::mapMassDependenceType(const string& massDepType)
{
	massDependencePtr massDep;
	if (   (massDepType == "BreitWigner")
	    or (massDepType == ""))  // default mass dependence
		massDep = createRelativisticBreitWigner();
	else if (massDepType == "flat")
		massDep = createFlatMassDependence();
	else if (massDepType == "piPiSWaveAuMorganPenningtonM")
		massDep = createPiPiSWaveAuMorganPenningtonM();
	else if (massDepType == "piPiSWaveAuMorganPenningtonVes")
		massDep = createPiPiSWaveAuMorganPenningtonVes();
	else if (massDepType == "piPiSWaveAuMorganPenningtonKachaev")
		massDep = createPiPiSWaveAuMorganPenningtonKachaev();
	else {
		printWarn << "unknown mass dependence '" << massDepType << "'. using Breit-Wigner." << endl;
		massDep = createRelativisticBreitWigner();
	}
	return massDep;
}


bool
waveDescription::constructDecayVertex(const Setting&                parentKey,
                                      const particlePtr&            parentParticle,
                                      vector<isobarDecayVertexPtr>& decayVertices,
                                      vector<particlePtr>&          fsParticles,
                                      const bool                    fromTemplate)
{
	if (_debug)
		printDebug << "reading decay vertex information from '" << parentKey.getPath() << "':" << endl;
	bool success = true;

	const Setting*      fsPartKeys = findLibConfigList(parentKey, "fsParticles", false);
	vector<particlePtr> fsDaughters;
	if (fsPartKeys)
		for (int i = 0; i < fsPartKeys->getLength(); ++i) {
			particlePtr p;
			if (constructParticle((*fsPartKeys)[i], p)) {
				fsDaughters.push_back(p);
				fsParticles.push_back(p);
			} else
				success = false;
		}
  
	const Setting*      isobarKeys = findLibConfigList(parentKey, "isobars", false);
	vector<particlePtr> isobarDaughters;
	if (isobarKeys)
		for (int i = 0; i < isobarKeys->getLength(); ++i) {
			particlePtr p;
			if (constructParticle((*isobarKeys)[i], p, not fromTemplate))
				isobarDaughters.push_back(p);
			else
				success = false;
			success &= constructDecayVertex((*isobarKeys)[i], isobarDaughters.back(),
			                                decayVertices, fsParticles, fromTemplate);
		}
  
	const unsigned int nmbDaughters = fsDaughters.size() + isobarDaughters.size();
	if (nmbDaughters != 2) {
		printWarn << "cannot construct isobar vertex, because number of daughters "
		          << "(" << nmbDaughters << ") in '" << parentKey.getPath() << "' "
		          << "is not 2" << endl;
		success = false;
	}
	if (not success)
		return false;
  
	// get isobar vertex parameters
	int L = 0, S = 0;
	if ((   not parentKey.lookupValue("L", L)
	     or not parentKey.lookupValue("S", S)) and not fromTemplate)
		printWarn << "Either L or S are not specified in '" << parentKey.getPath() << "'. "
		          << "using zero." << endl;
  
	// get mass dependence
	massDependencePtr massDep;
	if (parentParticle->bareName() != "X") {
		string         massDepType = "";
		const Setting* massDepKey  = findLibConfigGroup(parentKey, "massDep", false);
		if (massDepKey)
			massDepKey->lookupValue("name", massDepType);
		massDep = mapMassDependenceType(massDepType);
	}
  
	// if there is 1 final state particle and 1 isobar put them in the
	// same order as in the key file
	vector<particlePtr> daughters(2, particlePtr());
	if ((fsDaughters.size() == 1) and (isobarDaughters.size() == 1)) {
		if (fsPartKeys->getIndex() < isobarKeys->getIndex()) {
			daughters[0] = fsDaughters    [0];
			daughters[1] = isobarDaughters[0];
		} else {
			daughters[0] = isobarDaughters[0];
			daughters[1] = fsDaughters    [0];
		}
	} else if (fsDaughters.size() == 2)
		daughters = fsDaughters;
	else if (isobarDaughters.size() == 2)
		daughters = isobarDaughters;
  
	// construct isobar decay vertex
	decayVertices.push_back(createIsobarDecayVertex(parentParticle, daughters[0],
	                                                daughters[1], L, S, massDep));
	if (_debug)
		printDebug << "constructed " << *decayVertices.back() << endl;
	return true;
}


isobarAmplitudePtr
waveDescription::mapAmplitudeType(const string&                 formalismType,
                                  const isobarDecayTopologyPtr& topo)
{
	isobarAmplitudePtr amplitude;
	if (formalismType == "helicity")
		amplitude = createIsobarHelicityAmplitude(topo);
	else if (formalismType == "canonical")
		amplitude = createIsobarCanonicalAmplitude(topo);
	else {
		printWarn << "unknown amplitude formalism '" << formalismType << "'. "
		          << "using helicity formalism." << endl;
		amplitude = createIsobarHelicityAmplitude(topo);
	}
	return amplitude;
}


bool
waveDescription::setProductionVertexKeys(Setting&                   prodVertKey,
                                         const productionVertexPtr& prodVert)
{
	const string prodVertName = prodVert->name();
	if ((prodVertName == "diffractiveDissVertex") or (prodVertName == "leptoProductionVertex")) {
		if (_debug)
			printDebug << "setting keys for " << *prodVert << endl;
		prodVertKey.add("type", Setting::TypeString) = prodVertName;
		map<string, particlePtr> prodKinParticles;
		if (prodVertName == "diffractiveDissVertex") {
			const diffractiveDissVertexPtr& vert = static_pointer_cast<diffractiveDissVertex>(prodVert);
			prodKinParticles["beam"  ] = vert->beam();
			prodKinParticles["target"] = vert->target();
			if (vert->recoil()->name() != vert->target()->name())
				prodKinParticles["recoil"] = vert->recoil();
		} else if (prodVertName == "leptoProductionVertex") {
			const leptoProductionVertexPtr& vert = static_pointer_cast<leptoProductionVertex>(prodVert);
			prodKinParticles["beam"  ] = vert->beamLepton();
			prodKinParticles["target"] = vert->target();
			if (vert->recoil()->name() != vert->target()->name())
				prodKinParticles["recoil"] = vert->recoil();
		}
		// write particles
		Setting* beamParticleKey = 0;  // needed for leptoproduction vertex
		for (map<string, particlePtr>::iterator i = prodKinParticles.begin();
		     i != prodKinParticles.end(); ++i) {
			Setting& particleKey = prodVertKey.add(i->first, Setting::TypeGroup);
			particleKey.add("name", Setting::TypeString) = i->second->name();
			if (i->first == "beam")
				beamParticleKey = &particleKey;
		}
		if ((prodVertName == "leptoProductionVertex") and beamParticleKey)
			beamParticleKey->add("longPol", Setting::TypeFloat)
				= static_pointer_cast<leptoProductionVertex>(prodVert)->beamPol();
		return true;
	} else {
		printWarn << "writing of keys for production vertex of this type is not yet implemented:"
		          << *prodVert << endl;
		return false;
	}
}


bool
waveDescription::setXQuantumNumbersKeys(Setting&        XQnKey,
                                        const particle& X)
{
	if (_debug)
		printDebug << "setting quantum number keys for " << X.qnSummary() << endl;
	XQnKey.add("isospin", Setting::TypeInt) = X.isospin();
	if (X.G() != 0)
		XQnKey.add("G",     Setting::TypeInt) = X.G();
	XQnKey.add("J",       Setting::TypeInt) = X.J();
	if (X.P() != 0)
		XQnKey.add("P",     Setting::TypeInt) = X.P();
	if (X.C() != 0)
		XQnKey.add("C",     Setting::TypeInt) = X.C();
	XQnKey.add("M",       Setting::TypeInt) = X.spinProj();
	if (X.reflectivity() != 0)
		XQnKey.add("refl",  Setting::TypeInt) = X.reflectivity();
	return true;
}


bool
waveDescription::setMassDependence(Setting&              isobarDecayKey,
                                   const massDependence& massDep)
{
	const string massDepName = massDep.name();
	if (massDepName == "flatMassDependence")
		// default for X
		return true;
	else if (massDepName == "relativisticBreitWigner")
		// default mass dependence for isobars
		return true;
	else {
		if (_debug)
			printDebug << "setting key for '" << massDep << "' mass dependence to "
			           << "'" << massDepName << "'" << endl;
		Setting& massDepKey = isobarDecayKey.add("massDep", Setting::TypeGroup);
		massDepKey.add("name", Setting::TypeString) = massDepName;
	}
	return true;
}


bool
waveDescription::setXDecayKeys(Setting&                   parentDecayKey,
                               const isobarDecayTopology& topo,
                               const isobarDecayVertex&   vert)
{
	if (_debug)
		printDebug << "setting keys for decay " << vert << endl;
	setMassDependence(parentDecayKey, *vert.massDependence());
	vector<particlePtr> fsParticles;
	vector<particlePtr> isobars;
	for (unsigned int i = 0; i < vert.nmbOutParticles(); ++i) {
		const particlePtr& part = vert.outParticles()[i];
		if (topo.isFsParticle(part))
			fsParticles.push_back(part);
		else
			isobars.push_back(part);
	}
	bool success = true;	
	if (isobars.size() > 0) {
		Setting& isobarsKey = parentDecayKey.add("isobars", Setting::TypeList);
		for (unsigned int i = 0; i < isobars.size(); ++i) {
			Setting& isobarKey = isobarsKey.add(Setting::TypeGroup);
			isobarKey.add("name", Setting::TypeString) = isobars[i]->name();
			if (not setXDecayKeys(isobarKey, topo,
			                      *static_pointer_cast<isobarDecayVertex>(topo.toVertex(isobars[i]))))
				success = false;
		}
	}
	parentDecayKey.add("L", Setting::TypeInt) = (int)vert.L();
	parentDecayKey.add("S", Setting::TypeInt) = (int)vert.S();
	if (fsParticles.size() > 0) {
		Setting& fsParticlesKey = parentDecayKey.add("fsParticles", Setting::TypeList);
		for (unsigned int i = 0; i < fsParticles.size(); ++i) {
			Setting& fsParticleKey = fsParticlesKey.add(Setting::TypeGroup);
			fsParticleKey.add("name", Setting::TypeString) = fsParticles[i]->name();
		}
	}
	return success;
}


bool
waveDescription::writeKeyFile(Setting&                   rootKey,
                              const isobarDecayTopology& topo,
                              const bool                 writeProdVert)
{
	if (not topo.checkTopology() or not topo.checkConsistency()) {
		printWarn << "decay topology has issues. cannot write key file." << endl;
		return false;
	}
	if (_debug)
		printDebug << "writing " << topo;
	if (writeProdVert) {
		Setting& prodVertKey = rootKey.add("productionVertex", Setting::TypeGroup);
		setProductionVertexKeys(prodVertKey, topo.productionVertex());
	}
	Setting& waveKey   = rootKey.add("wave",            Setting::TypeGroup);
	Setting& XQnKey    = waveKey.add("XQuantumNumbers", Setting::TypeGroup);
	Setting& XDecayKey = waveKey.add("XDecay",          Setting::TypeGroup);
	if (not setXQuantumNumbersKeys(XQnKey, *(topo.XParticle()))) {
		printWarn << "problems writing X quantum numbers" << endl;
		return false;
	}
	if (not setXDecayKeys(XDecayKey, topo, *(topo.XIsobarDecayVertex()))) {
		printWarn << "problems writing X decay" << endl;
		return false;
	}
	return true;
}


bool
waveDescription::setAmplitude(Setting&               amplitudeKey,
                              const isobarAmplitude& amplitude)
{
	string       formalism     = "";
	const string amplitudeName = amplitude.name();
	if (amplitudeName == "isobarCanonicalAmplitude") {
		formalism = "canonical";
	} else if (amplitudeName != "isobarHelicityAmplitude") {
		printWarn << "unknown amplitude type '" << amplitudeName << "'" << endl;
		return false;
	}
	if (formalism != "") {
		if (_debug)
			printDebug << "setting 'formalism' key for '" << amplitudeName << "' to "
			           << "'" << formalism << "'" << endl;
		amplitudeKey.add("formalism", Setting::TypeString) = formalism;
	}
	if (not amplitude.boseSymmetrization())
		amplitudeKey.add("boseSymmetrize", Setting::TypeBoolean) = amplitude.boseSymmetrization();
	if (not amplitude.reflectivityBasis())
		amplitudeKey.add("useReflectivityBasis", Setting::TypeBoolean) = amplitude.reflectivityBasis();
	return true;
}


bool
waveDescription::writeKeyFile(Setting&               rootKey,
                              const isobarAmplitude& amplitude)
{
	printInfo << "writing ";
	amplitude.printParameters(cout);
	Setting& amplitudeKey = rootKey.add("amplitude", Setting::TypeGroup);
	if (not setAmplitude(amplitudeKey, amplitude)) {
		printWarn << "problems writing amplitude paramaters" << endl;
		return false;
	}
	return true;
}
