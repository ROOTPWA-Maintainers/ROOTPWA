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
bool                         keyFileParser::_debug = false;


isobarDecayTopologyPtr
keyFileParser::parse(const string& keyFileName)
{
  printInfo << "parsing key file '" << keyFileName << "'" << endl;

  _prodVert = interactionVertexPtr();
  _decayVertices.clear();
  _fsParticles.clear();
  isobarDecayTopologyPtr topoNullPtr;

  // parse key file
  Config key;
  try {
    key.readFile(keyFileName.c_str());
  } catch(const FileIOException& ioEx) {
    printWarn << "I/O error while reading key file '" << keyFileName << "'. "
	      << "cannot construct decay topology." << endl;
    return topoNullPtr;
  } catch(const ParseException&  parseEx) {
    printWarn << "parse error in '" << parseEx.getFile() << "' line " << parseEx.getLine()
	      << ": " << parseEx.getError() << ". cannot construct decay topology." << endl;
    return topoNullPtr;
  }

  const Setting& rootKey = key.getRoot();

  // find 'wave' group
  const Setting* waveKey = findGroup(rootKey, "wave");
  if (!waveKey) {
    printWarn << "cannot find 'wave' group. cannot construct decay topology." << endl;
    return topoNullPtr;
  }

  // find 'options' group
  const Setting* waveOptionsKey = findGroup(*waveKey, "options");
  // default option values
  bool boseSymmetrize       = true;
  bool useReflectivityBasis = false;
  if (waveOptionsKey) {
    if (waveOptionsKey->lookupValue("boseSymmetrize",       boseSymmetrize) && _debug)
      printInfo << "setting wave option 'boseSymmetrize' to "
		<< ((boseSymmetrize) ? "true" : "false") << endl;
    if (waveOptionsKey->lookupValue("useReflectivityBasis", useReflectivityBasis) && _debug)
      printInfo << "setting wave option 'useReflectivityBasis' to "
		<< ((useReflectivityBasis) ? "true" : "false") << endl;
  }

  // find 'XDecay' group
  const Setting* XDecayKey = findGroup(*waveKey, "XDecay");
  if (!XDecayKey) {
    printWarn << "cannot find 'XDecay' group. cannot construct decay topology." << endl;
    return topoNullPtr;
  }

  // find 'XQuantumNumbers' group
  const Setting* XQnKey = findGroup(*waveKey, "XQuantumNumbers");
  if (!XQnKey) {
    printWarn << "cannot find 'XQuantumNumbers' group. cannot construct decay topology." << endl;
    return topoNullPtr;
  }
  // get X quantum numbers
  map<string, int> XQn;
  XQn["isospin"];
  XQn["G"];
  XQn["J"];
  XQn["P"];
  XQn["C"];
  XQn["M"];
  bool XQnOkay = true;
  for (map<string, int>::iterator i = XQn.begin(); i != XQn.end(); ++i)
    if (!XQnKey->lookupValue(i->first, i->second)) {
      printWarn << "cannot find integer field '" << i->first << "in "
		<< "'" << XQnKey->getPath() << "'." << endl;
      XQnOkay = false;
    }
  if (!XQnOkay) {
    printWarn << "cannot find X quantum numbers. cannot construct decay topology." << endl;
    return topoNullPtr;
  }
  // create X particle
  particlePtr X = createParticle("X", XQn["isospin"], XQn["G"],
				 XQn["J"], XQn["P"], XQn["C"], XQn["M"]);
  if (_debug)
    printInfo << "constructed X particle: " << X->qnSummary() << endl;

  // create production vertex
  if (!constructProductionVertex(rootKey, X)) {
    printWarn << "problems constructing production vertex. cannot construct decay topology." << endl;
    return topoNullPtr;
  }

  // traverse decay chain and create final state particles and isobar decay vertices
  if (!constructDecayVertex(*XDecayKey, X)) {
    printWarn << "problems constructing decay chain. cannot construct decay topology." << endl;
    return topoNullPtr;
  }

  // create decay isobar decay topology
  isobarDecayTopologyPtr topo(new isobarDecayTopology(_prodVert, _decayVertices, _fsParticles));
  topo->calcIsobarCharges();
  
  printInfo << "successfully constructed decay topology from key file "
	    << "'" << keyFileName << "'." << endl;
  return topo;
}


const Setting*
keyFileParser::findGroup(const Setting&     parent,
			 const std::string& groupName,
			 const bool         mustExist)
{
  // find field
  if (!parent.exists(groupName)) {
    if (mustExist)
      printWarn << "cannot find '" << groupName << "' field in '" << parent.getPath() << "'"
		<< "of key file '" << parent.getSourceFile() << "'." << endl;
    return 0;
  }
  const Setting* groupKey = &parent[groupName];
  if (!groupKey->isGroup()) {
    printWarn << "'" << groupName << "' field in '" << parent.getPath() << "'"
	      << "of key file '" << parent.getSourceFile() << "' is not a group." << endl;
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
  if (!parent.exists(listName)) {
    if (mustExist)
      printWarn << "cannot find '" << listName << "' field in '" << parent.getPath() << "'"
		<< "of key file '" << parent.getSourceFile() << "'." << endl;
    return 0;
  }
  const Setting* listKey = &parent[listName];
  if (!listKey->isList()) {
    printWarn << "'" << listName << "' field in '" << parent.getPath() << "'"
	      << "of key file '" << parent.getSourceFile() << "' is not a list. "
	      << "check that braces are correct." << endl;
    return 0;
  }
  if (listKey->getLength() < 1) {
    printWarn << "list '" << listName << "' in '" << parent.getPath() << "'"
	      << "of key file '" << parent.getSourceFile() << "' is empty." << endl;
    return 0;
  }
  return listKey;
}


bool
keyFileParser::constructParticle(const Setting& particleKey,
				 particlePtr&   particle)
{
  string name;
  if (particleKey.lookupValue("name", name)) {
    particle = createParticle(name);
    if (_debug)
      printInfo << "created particle " << particle->qnSummary() << endl;
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
  
  if (!success)
    return false;

  // get isobar vertex parameters
  int L = 0, S = 0;
  if (   !parentKey.lookupValue("L", L)
      || !parentKey.lookupValue("S", S))
    printWarn << "Either L or S are not specified in '" << parentKey.getPath() << "'. "
	      << "using zero." << endl;
  string massDepString = "";
  parentKey.lookupValue("massDep", massDepString);
  massDependencePtr massDep;
  if (parentParticle->bareName() != "X") {
    if (   (massDepString == "BreitWigner")
	|| (massDepString == ""))
      massDep = createRelativisticBreitWigner();
    else {
      printWarn << "unknown mass dependence '" << massDepString << "'. using Breit-Wigner." << endl;
      massDep = createRelativisticBreitWigner();
    }
  }
  
  // construct isobar decay vertex
  _decayVertices.push_back(createIsobarDecayVertex(parentParticle, daughters[0],
						   daughters[1], L, S, massDep));
  if (_debug)
    printInfo << "constructed " << *_decayVertices.back() << endl;
  return true;
}


bool
keyFileParser::constructProductionVertex(const Setting&     rootKey,
					 const particlePtr& X)
{
  // find 'productionVertex' group
  const Setting* prodVertKey = findGroup(rootKey, "productionVertex");
  if (!prodVertKey) {
    printWarn << "cannot find 'productionVertex' group" << endl;
    return false;
  }

  bool success = true;

  // get vertex type
  string vertType;
  if (!prodVertKey->lookupValue("type", vertType)) {
    printWarn << "cannot find 'type' entry in '" << prodVertKey->getPath() << "'" << endl;
    success = false;
  }

  // get beam particles
  if (!prodVertKey->exists("beamParticles")) {
    printWarn << "cannot find 'beamParticles' entry in '" << prodVertKey->getPath() << "'" << endl;
    success = false;
  }
  if (_debug)
    printInfo << "reading beam particles from '" << prodVertKey->getPath() << "':" << endl;
  const Setting* beamPartKeys = findList(*prodVertKey, "beamParticles");
  if (!beamPartKeys)
    success = false;
  
  if (!success)
    return false;

  // create production vertex
  if (vertType == "diffractiveDissociation") {
    if (beamPartKeys->getLength() > 1) {
      printWarn << "list of beam particles in '" << beamPartKeys->getPath() << "' "
		<< "has " << beamPartKeys->getLength()  << " entries. using only first one." << endl;
    }
    particlePtr beam;
    if (!constructParticle((*beamPartKeys)[0], beam))
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
