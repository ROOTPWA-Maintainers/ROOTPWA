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
#include "keyFileParser.h"

	
using namespace std;
using namespace libconfig;
using namespace rpwa;


keyFileParser                keyFileParser::_instance;
interactionVertexPtr         keyFileParser::_prodVert;
vector<isobarDecayVertexPtr> keyFileParser::_decayVertices;
vector<particlePtr>          keyFileParser::_fsParticles;
bool                         keyFileParser::_useReflectivityBasis = true;
bool                         keyFileParser::_boseSymmetrize       = true;
bool                         keyFileParser::_debug                = false;


bool
keyFileParser::parse(const string&           keyFileName,
		     isobarDecayTopologyPtr& topo)
{
  printInfo << "parsing key file '" << keyFileName << "'" << endl;

  _prodVert = interactionVertexPtr();
  _decayVertices.clear();
  _fsParticles.clear();
  _useReflectivityBasis = true;
  _boseSymmetrize       = true;
  if (topo)
    topo.reset();
  topo = isobarDecayTopologyPtr();  // null pointer

  // parse key file
  Config key;
  try {
    key.readFile(keyFileName.c_str());
  } catch(const FileIOException& ioEx) {
    printWarn << "I/O error while reading key file '" << keyFileName << "'. "
	      << "cannot construct decay topology." << endl;
    return false;
  } catch(const ParseException&  parseEx) {
    printWarn << "parse error in '" << parseEx.getFile() << "' line " << parseEx.getLine()
	      << ": " << parseEx.getError() << ". cannot construct decay topology." << endl;
    return false;
  }

  const Setting& rootKey = key.getRoot();

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

  // find X decay group
  const Setting* XDecayKey = findGroup(*waveKey, "XDecay");
  if (not XDecayKey) {
    printWarn << "cannot find 'XDecay' group. cannot construct decay topology." << endl;
    return false;
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
  if (not constructProductionVertex(rootKey, X)) {
    printWarn << "problems constructing production vertex. cannot construct decay topology." << endl;
    return false;
  }

  // traverse decay chain and create final state particles and isobar decay vertices
  if (not constructDecayVertex(*XDecayKey, X)) {
    printWarn << "problems constructing decay chain. cannot construct decay topology." << endl;
    return false;
  }

  // create decay isobar decay topology
  topo = createIsobarDecayTopology(_prodVert, _decayVertices, _fsParticles);
  topo->calcIsobarCharges();
  
  printInfo << "successfully constructed decay topology from key file "
	    << "'" << keyFileName << "'." << endl;
  return true;
}


void
keyFileParser::setAmplitudeOptions(isobarHelicityAmplitude& amp)
{
  if (_debug)
    printInfo << "setting amplitude options: "
	      << ((_boseSymmetrize) ? "en" : "dis") << "abling Bose symmetrization, "
	      << ((_useReflectivityBasis) ? "en" : "dis") << "abling reflectivity basis" << endl;
  amp.enableReflectivityBasis (_useReflectivityBasis);
  amp.enableBoseSymmetrization(_boseSymmetrize      );
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
keyFileParser::constructParticle(const Setting& particleKey,
				 particlePtr&   particle)
{
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
keyFileParser::constructDecayVertex(const Setting&     parentKey,
				    const particlePtr& parentParticle)
{
  bool                success = true;
  vector<particlePtr> daughters;

  if (_debug)
    printInfo << "reading final state particles from '" << parentKey.getPath() << "':" << endl;
  const Setting* fsPartKeys = findList(parentKey, "fsParticles", false);
  if (fsPartKeys)
    for (int i = 0; i < fsPartKeys->getLength(); ++i) {
      particlePtr p;
      if (constructParticle((*fsPartKeys)[i], p)) {
	daughters.push_back   (p);
	_fsParticles.push_back(p);
      }	else
	success = false;
    }

  if (_debug)
    printInfo << "reading isobars from '" << parentKey.getPath() << "':" << endl;
  const Setting* isobarKeys = findList(parentKey, "isobars", false);
  if (isobarKeys)
    for (int i = 0; i < isobarKeys->getLength(); ++i) {
      particlePtr p;
      if (constructParticle((*isobarKeys)[i], p))
	daughters.push_back(p);
      else
	success = false;
      success &= constructDecayVertex((*isobarKeys)[i], daughters.back());
    }

  if (daughters.size() < 2) {
    printWarn << "cannot construct isobar vertex, because number of daughters "
	      << "(" << daughters.size() << ") in '" << parentKey.getPath() << "' "
	      << "is too low" << endl;
    success = false;
  }
  if (daughters.size() > 2)
    printWarn << "number of daughters " << "(" << daughters.size() << ") "
	      << "in '" << parentKey.getPath() << "' is larger than 2. "
	      << "using daughters " << daughters[0]->qnSummary() << " and "
	      << daughters[1]->qnSummary() << endl;
  
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
  
  // construct isobar decay vertex
  _decayVertices.push_back(createIsobarDecayVertex(parentParticle, daughters[0],
						   daughters[1], L, S, massDep));
  if (_debug)
    printInfo << "constructed " << *_decayVertices.back() << endl;
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
keyFileParser::constructProductionVertex(const Setting&     rootKey,
					 const particlePtr& X)
{
  // find production vertex group
  const Setting* prodVertKey = findGroup(rootKey, "productionVertex");
  if (not prodVertKey) {
    printWarn << "cannot find 'productionVertex' group" << endl;
    return false;
  }
  bool success = true;
  // get vertex type
  string vertType;
  if (not prodVertKey->lookupValue("type", vertType)) {
    printWarn << "cannot find 'type' entry in '" << prodVertKey->getPath() << "'" << endl;
    success = false;
  }
  // get production kinematics particles
  const Setting* isPartKeys = findList(*prodVertKey, "particles");
  if (not isPartKeys)
    success = false;
  if (not success)
    return false;
  // create production vertex
  if (_debug)
    printInfo << "reading production kinematics particles from "
	      << "'" << prodVertKey->getPath() << "'" << endl;
  return mapProductionVertexType(vertType, *isPartKeys, X);
}


bool
keyFileParser::mapProductionVertexType(const string&      vertType,
				       const Setting&     isPartKeys,
				       const particlePtr& X)
{
  if (vertType == "diffractiveDissociation") {
    if (isPartKeys.getLength() > 1) {
      printWarn << "list of production kinematics particles in '" << isPartKeys.getPath() << "' "
		<< "has " << isPartKeys.getLength()  << " entries. using only first one." << endl;
    }
    particlePtr beam;
    if (not constructParticle(isPartKeys[0], beam))
      return false;
    _prodVert = createDiffractiveDissVertex(beam, X);

  } else {
    printWarn << "unknown production vertex type '" << vertType << "'" << endl;
    return false;
  }
  if (_debug)
    printInfo << "constructed " << *_prodVert << endl;
  return true;
}
