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
//      singleton class that constructs decay topologies according to key files
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <fstream>
#include <map>

#include "utilities.h"
#include "diffractiveDissVertex.h"
#include "leptoProductionVertex.h"
#include "keyFileParser.h"

  
using namespace std;
using namespace boost;
using namespace libconfig;
using namespace rpwa;


keyFileParser keyFileParser::_instance;
Config        keyFileParser::_key;
bool          keyFileParser::_boseSymmetrize       = true;
bool          keyFileParser::_useReflectivityBasis = true;
bool          keyFileParser::_debug                = false;


bool
keyFileParser::parse(const string& keyFileName)
{
	printInfo << "parsing key file '" << keyFileName << "'" << endl;
	try {
		_key.readFile(keyFileName.c_str());
	} catch(const FileIOException& ioEx) {
		printWarn << "I/O error while reading key file '" << keyFileName << "'. "
		          << "cannot construct decay topology." << endl;
		return false;
	} catch(const ParseException&  parseEx) {
		printWarn << "parse error in '" << parseEx.getFile() << "' line " << parseEx.getLine()
		          << ": " << parseEx.getError() << ". cannot construct decay topology." << endl;
		return false;
	}
	return true;
}


bool
keyFileParser::constructDecayTopology(isobarDecayTopologyPtr& topo)
{
	_boseSymmetrize       = true;
	_useReflectivityBasis = true;
	if (topo)
		topo.reset();
	topo = isobarDecayTopologyPtr();  // null pointer

	const Setting& rootKey = _key.getRoot();
	printInfo << "constructing decay topology from key file" << endl;

	// find wave group
	const Setting* waveKey = findGroup(rootKey, "wave");
	if (not waveKey) {
		printWarn << "cannot find 'wave' group. cannot construct decay topology." << endl;
		return false;
	}

	// find wave options group
	const Setting* waveOptionsKey = findGroup(*waveKey, "options", false);
	if (waveOptionsKey) {
		if (waveOptionsKey->lookupValue("boseSymmetrize", _boseSymmetrize) and _debug)
			printInfo << "setting wave option 'boseSymmetrize' to "
			          << ((_boseSymmetrize) ? "true" : "false") << endl;
		if (waveOptionsKey->lookupValue("useReflectivityBasis", _useReflectivityBasis) and _debug)
			printInfo << "setting wave option 'useReflectivityBasis' to "
			          << ((_useReflectivityBasis) ? "true" : "false") << endl;
	}

	// find  X quantum numbers group
	const Setting* XQnKey = findGroup(*waveKey, "XQuantumNumbers");
	if (not XQnKey) {
		printWarn << "cannot find 'XQuantumNumbers' group. cannot construct decay topology." << endl;
		return false;
	}
	// create X particle
	particlePtr X;
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
	const Setting* XDecayKey = findGroup(*waveKey, "XDecay");
	if (not XDecayKey) {
		printWarn << "cannot find 'XDecay' group. cannot construct decay topology." << endl;
		return false;
	}
	// traverse decay chain and create final state particles and isobar decay vertices
	vector<isobarDecayVertexPtr> decayVertices;
	vector<particlePtr>          fsParticles;
	if (not constructDecayVertex(*XDecayKey, X, decayVertices, fsParticles)) {
		printWarn << "problems constructing decay chain. cannot construct decay topology." << endl;
		return false;
	}

	// create decay isobar decay topology
	topo = createIsobarDecayTopology(prodVert, decayVertices, fsParticles);
	topo->calcIsobarCharges();
  
	printInfo << "successfully constructed decay topology from key file" << endl;
	return true;
}


void
keyFileParser::setAmplitudeOptions(isobarHelicityAmplitude& amp)
{
	if (_debug)
		printInfo << "setting amplitude options: "
		          << ((_boseSymmetrize) ? "en" : "dis") << "abling Bose symmetrization, "
		          << ((_useReflectivityBasis) ? "en" : "dis") << "abling reflectivity basis" << endl;
	amp.enableBoseSymmetrization(_boseSymmetrize      );
	amp.enableReflectivityBasis (_useReflectivityBasis);
}


bool
keyFileParser::writeKeyFile(const string&                 keyFileName,
                            const isobarDecayTopologyPtr& topo,
                            const bool                    writeProdVert,
                            const bool)
{
	printInfo << "writing key file '" << keyFileName << "': " << *topo;
	Config key;
	Setting& rootKey = key.getRoot();

	if (writeProdVert) {
		Setting& prodVertKey = rootKey.add("productionVertex", Setting::TypeGroup);
		setProductionVertexKeys(prodVertKey, topo->productionVertex());
	}
	
	Setting& waveKey = rootKey.add("wave", Setting::TypeGroup);

	// if (writeWaveOpt) {
	// 	Setting& waveOptionsKey = waveKey.add("options", Setting::TypeGroup);
	// 	waveOptionsKey.add("boseSymmetrize",       Setting::TypeBoolean) = _boseSymmetrize;
	// 	waveOptionsKey.add("useReflectivityBasis", Setting::TypeBoolean) = _useReflectivityBasis;
	// }
	
	Setting& XQnKey = waveKey.add("XQuantumNumbers", Setting::TypeGroup);
	setXQuantumNumbersKeys(XQnKey, *(topo->XParticle()));

	Setting& XDecayKey = waveKey.add("XDecay", Setting::TypeGroup);
	setXDecayKeys(XDecayKey, *topo, *(topo->XIsobarDecayVertex()));

	try {
		key.writeFile(keyFileName.c_str());
	} catch(const FileIOException& ioEx) {
		printWarn << "I/O error while writing key file '" << keyFileName << "'" << endl;
		return false;
	}
	printInfo << "successfully written key file '" << keyFileName << "'" << endl;
	return true;
}


const Setting*
keyFileParser::findGroup(const Setting&     parent,
                         const std::string& groupName,
                         const bool         mustExist)
{
	// find field
	if (not parent.exists(groupName)) {
		if (mustExist)
			printWarn << "cannot find '" << groupName << "' field in '" << parent.getPath() << "' "
			          << "of key file '" << parent.getSourceFile() << "'" << endl;
		return 0;
	}
	const Setting* groupKey = &parent[groupName];
	// check that it is a group
	if (not groupKey->isGroup()) {
		printWarn << "'" << groupName << "' field in '" << parent.getPath() << "' "
		          << "of key file '" << parent.getSourceFile() << "' is not a group" << endl;
		return 0;
	}
	return groupKey;
}


const Setting*
keyFileParser::findList(const Setting&     parent,
                        const std::string& listName,
                        const bool         mustExist)
{
	// find field
	if (not parent.exists(listName)) {
		if (mustExist)
			printWarn << "cannot find '" << listName << "' field in '" << parent.getPath() << "' "
			          << "of key file '" << parent.getSourceFile() << "'" << endl;
		return 0;
	}
	const Setting* listKey = &parent[listName];
	// check that it is a list
	if (not listKey->isList()) {
		printWarn << "'" << listName << "' field in '" << parent.getPath() << "' "
		          << "of key file '" << parent.getSourceFile() << "' is not a list. "
		          << "check that braces are correct." << endl;
		return 0;
	}
	// check that it is not empty
	if (listKey->getLength() < 1) {
		printWarn << "list '" << listName << "' in '" << parent.getPath() << "' "
		          << "of key file '" << parent.getSourceFile() << "' is empty" << endl;
		return 0;
	}
	return listKey;
}


bool
keyFileParser::constructXParticle(const Setting& XQnKey,
                                  particlePtr&   X)
{
	if (_debug)
		printInfo << "reading X quantum numbers from '" << XQnKey.getPath() << "':" << endl;
	// get X quantum numbers
	map<string, int> XQn;
	XQn["isospin"];
	XQn["G"      ];
	XQn["J"      ];
	XQn["P"      ];
	XQn["C"      ];
	XQn["M"      ];
	if (_useReflectivityBasis)
		XQn["refl"];
	bool success = true;
	for (map<string, int>::iterator i = XQn.begin(); i != XQn.end(); ++i)
		if (not XQnKey.lookupValue(i->first, i->second)) {
			printWarn << "cannot find integer field '" << i->first << "in "
			          << "'" << XQnKey.getPath() << "'" << endl;
			success = false;
		}
	if (not success) {
		printWarn << "cannot find X quantum numbers. cannot construct decay topology." << endl;
		return false;
	}
	// create X particle
	X = createParticle("X",
	                   XQn["isospin"], XQn["G"],
	                   XQn["J"], XQn["P"], XQn["C"],
	                   XQn["M"], (_useReflectivityBasis) ? XQn["refl"]: 0);
	if (_debug)
		printInfo << "constructed X particle: " << X->qnSummary() << endl;
	return true;
}


bool
keyFileParser::setXQuantumNumbersKeys(Setting&        XQnKey,
                                      const particle& X)
{
	if (_debug)
		printInfo << "setting quantum number keys for " << X.qnSummary() << endl;
	XQnKey.add("isospin", Setting::TypeInt) = X.isospin();
	XQnKey.add("G",       Setting::TypeInt) = X.G();
	XQnKey.add("J",       Setting::TypeInt) = X.J();
	XQnKey.add("P",       Setting::TypeInt) = X.P();
	XQnKey.add("C",       Setting::TypeInt) = X.C();
	XQnKey.add("M",       Setting::TypeInt) = X.spinProj();
	if (_useReflectivityBasis)
		XQnKey.add("refl", Setting::TypeInt) = X.reflectivity();
	return true;
}


bool
keyFileParser::constructParticle(const Setting& particleKey,
                                 particlePtr&   particle)
{
	if (_debug)
		printInfo << "reading particle information from '" << particleKey.getPath() << "':" << endl;
	string name;
	if (particleKey.lookupValue("name", name)) {
		particle = createParticle(name);
		int spinProj;
		if (particleKey.lookupValue("spinProj", spinProj))
			particle->setSpinProj(spinProj);
		if (_debug)
			printInfo << "created particle " << particle->qnSummary() << flush;
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


bool
keyFileParser::constructDecayVertex(const Setting&                parentKey,
                                    const particlePtr&            parentParticle,
                                    vector<isobarDecayVertexPtr>& decayVertices,
                                    vector<particlePtr>&          fsParticles)
{
	if (_debug)
		printInfo << "reading decay vertex information from '" << parentKey.getPath() << "':" << endl;
	bool success = true;

	const Setting*      fsPartKeys = findList(parentKey, "fsParticles", false);
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
  
	const Setting*      isobarKeys = findList(parentKey, "isobars", false);
	vector<particlePtr> isobarDaughters;
	if (isobarKeys)
		for (int i = 0; i < isobarKeys->getLength(); ++i) {
			particlePtr p;
			if (constructParticle((*isobarKeys)[i], p))
				isobarDaughters.push_back(p);
			else
				success = false;
			success &= constructDecayVertex((*isobarKeys)[i], isobarDaughters.back(),
			                                decayVertices, fsParticles);
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
	if (   not parentKey.lookupValue("L", L)
	       or not parentKey.lookupValue("S", S))
		printWarn << "Either L or S are not specified in '" << parentKey.getPath() << "'. "
		          << "using zero." << endl;
  
	// get mass dependence
	massDependencePtr massDep;
	if (parentParticle->bareName() != "X") {
		string         massDepType = "";
		const Setting* massDepKey  = findGroup(parentKey, "massDep", false);
		if (massDepKey)
			massDepKey->lookupValue("name", massDepType);
		massDep = mapMassDependence(massDepType);
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
		printInfo << "constructed " << *decayVertices.back() << endl;
	return true;
}


massDependencePtr
keyFileParser::mapMassDependence(const string& massDepType)
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
keyFileParser::setMassDependence(Setting&              isobarDecayKey,
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
			printInfo << "setting key for " << massDep << " mass dependence to "
			          << "'" << massDepName << "'" << endl;
		Setting& massDepKey = isobarDecayKey.add("massDep", Setting::TypeGroup);
		massDepKey.add("name", Setting::TypeString) = massDepName;
	}
	return true;
}


bool
keyFileParser::constructProductionVertex(const Setting&       rootKey,
                                         const particlePtr&   X,
                                         productionVertexPtr& prodVert)
{
	// find production vertex group
	const Setting* prodVertKey = findGroup(rootKey, "productionVertex");
	if (not prodVertKey) {
		printWarn << "cannot find 'productionVertex' group" << endl;
		return false;
	}
	if (_debug)
		printInfo << "reading production vertex from '" << prodVertKey->getPath() << "':" << endl;
	bool success = true;
	// get vertex type
	string vertType;
	if (not prodVertKey->lookupValue("type", vertType)) {
		printWarn << "cannot find 'type' entry in '" << prodVertKey->getPath() << "'" << endl;
		success = false;
	}
	// create production vertex
	return mapProductionVertexType(*prodVertKey, vertType, X, prodVert);
}


bool
keyFileParser::mapProductionVertexType(const Setting&       prodVertKey,
                                       const string&        vertType,
                                       const particlePtr&   X,
                                       productionVertexPtr& prodVert)
{
	prodVert = productionVertexPtr();
	bool success = true;
	if ((vertType == "diffractiveDissVertex") or (vertType == "leptoProductionVertex")) {
		const string             prodKinParticleNames[] = {"beam", "target", "recoil"};
		map<string, particlePtr> prodKinParticles;
		for (unsigned int i = 0; i < sizeof(prodKinParticleNames) / sizeof(string); ++i)
			prodKinParticles[prodKinParticleNames[i]] = particlePtr();
		// construct particles
		const Setting* beamParticleKey = 0;  // needed for leptoproduction vertex
		for (map<string, particlePtr>::iterator i = prodKinParticles.begin();
		     i != prodKinParticles.end(); ++i) {
			const Setting* particleKey = findGroup(prodVertKey, i->first,
			                                       (i->first != "recoil") ? true : false);
			if (particleKey) {
				if (i->first != "beam")
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
			if (vertType == "leptoProductionVertex" ) {
				prodVert = createLeptoProductionVertex(prodKinParticles["beam"], prodKinParticles["target"],
				                                       X, prodKinParticles["recoil"]);
				double beamLongPol;
				if (beamParticleKey->lookupValue("longPol", beamLongPol)) {
					if (_debug)
						printInfo << "setting polarization of beam " << prodKinParticles["beam"]->qnSummary()
						          << " to " << beamLongPol << endl;
					static_pointer_cast<leptoProductionVertex>(prodVert)->setBeamPol(beamLongPol);
				} else
					printWarn << "no polarization is given for beam " << prodKinParticles["beam"]->qnSummary()
					          << ". assuming unpolarized beam." << endl;
			}
		}
	} else {
		printWarn << "unknown production vertex type '" << vertType << "'" << endl;
		return false;
	}
	if (_debug) {
		if (success and prodVert)
			printInfo << "constructed " << *prodVert << endl;
		else
			printWarn << "problems constructing production vertex of type '" << vertType << "'" << endl;
	}
	return success;
}


bool
keyFileParser::setProductionVertexKeys(Setting&                   prodVertKey,
                                       const productionVertexPtr& prodVert)
{
	const string prodVertName = prodVert->name();
	if ((prodVertName == "diffractiveDissVertex") or (prodVertName == "leptoProductionVertex")) {
		if (_debug)
			printInfo << "setting keys for " << *prodVert << endl;
		prodVertKey.add("type", Setting::TypeString) = prodVertName;
		map<string, particlePtr> prodKinParticles;
		if (prodVertName == "diffractiveDissVertex") {
			const diffractiveDissVertexPtr& vert = static_pointer_cast<diffractiveDissVertex>(prodVert);
			prodKinParticles["beam"  ] = vert->beam  ();
			prodKinParticles["target"] = vert->target();
			prodKinParticles["recoil"] = vert->recoil();
		} else if (prodVertName == "leptoProductionVertex") {
			const leptoProductionVertexPtr& vert = static_pointer_cast<leptoProductionVertex>(prodVert);
			prodKinParticles["beam"  ] = vert->beamLepton();
			prodKinParticles["target"] = vert->target    ();
			prodKinParticles["recoil"] = vert->recoil    ();
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
keyFileParser::setXDecayKeys(Setting&                   parentDecayKey,
                             const isobarDecayTopology& topo,
                             const isobarDecayVertex&   vert)
{
	if (_debug)
		printInfo << "setting keys for decay " << vert << endl;
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
