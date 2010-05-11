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
vector<isobarDecayVertexPtr> keyFileParser::_decayVertices;
vector<particlePtr>          keyFileParser::_fsParticles;
bool                         keyFileParser::_debug = false;


isobarDecayTopologyPtr
keyFileParser::parse(const string& keyFileName)
{
  printInfo << "parsing key file '" << keyFileName << "'" << endl;

  isobarDecayTopologyPtr topoNullPtr;
  _decayVertices.clear();
  _fsParticles.clear();

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
    printWarn << "problems finding 'wave' group. cannot construct decay topology." << endl;
    return topoNullPtr;
  }

  // find 'options' group
  const Setting* optionsKey = findGroup(*waveKey, "options");
  if (optionsKey) {
  }

  // find 'XDecay' group
  const Setting* XDecay = findGroup(*waveKey, "XDecay");
  if (!XDecay) {
    printWarn << "problems finding 'XDecay' group. cannot construct decay topology." << endl;
    return topoNullPtr;
  }

  // find 'XQuantumNumbers' group
  const Setting* XQnKey = findGroup(*waveKey, "XQuantumNumbers");
  if (!XQnKey) {
    printWarn << "problems finding 'XQuantumNumbers' group. cannot construct decay topology." << endl;
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
    printWarn << "problems finding X quantum numbers. cannot construct decay topology." << endl;
    return topoNullPtr;
  }
  // create X particle
  particlePtr X = createParticle("X", XQn["isospin"], XQn["G"],
				 XQn["J"], XQn["P"], XQn["C"], XQn["M"]);
  if (_debug)
    printInfo << "constructed X particle: " << X->qnSummary() << endl;

  // traverse decay chain and create final state particles and isobar decay vertices
  constructVertex(*XDecay, X);

//!!! still needs to be implemented
  // create production vertex
  particlePtr              beam     = createParticle("pi-");
  diffractiveDissVertexPtr prodVert = createDiffractiveDissVertex(beam, X);
  

  // create decay isobar decay topology
  isobarDecayTopologyPtr topo(new isobarDecayTopology(prodVert, _decayVertices, _fsParticles));
  topo->calcIsobarCharges();
  
  printInfo << "successfully constructed decay topology from key file "
	    << "'" << keyFileName << "'." << endl;
  return topo;
}


const Setting*
keyFileParser::findGroup(const Setting&     parent,
			 const std::string& groupName)
{
  // find field
  if (!parent.exists(groupName)) {
    printWarn << "cannot find '" << groupName << "' field in '" << parent.getPath() << "'"
	      << "of key file '" << parent.getSourceFile() << "'." << endl;
    return 0;
  }
  const Setting* group = &parent[groupName];
  if (!group->isGroup()) {
    printWarn << "'" << groupName << "' field in '" << parent.getPath() << "'"
	      << "of key file '" << parent.getSourceFile() << "' is not a group." << endl;
    return 0;
  }
  return group;
}


bool
keyFileParser::constructVertex(const Setting&     parentKey,
			       const particlePtr& parentParticle)
{
  bool                success = true;
  vector<particlePtr> daughters;
  if (parentKey.exists("fsParticles")) {
    if (_debug)
      printInfo << "reading final state particles from '" << parentKey.getPath() << "':" << endl;
    const Setting& fsPartKeys = parentKey["fsParticles"];
    for (int i = 0; i < fsPartKeys.getLength(); ++i) {
      const Setting& fsPartKey = fsPartKeys[i];
      // mandatory fields
      string name;
      if (fsPartKey.lookupValue("name", name)) {
	particlePtr p = createParticle(name);
	daughters.push_back   (p);
	_fsParticles.push_back(p);
      }	else {
	printWarn << "final state particle group '" << fsPartKey.getPath() << "' "
		  << "has no name field." << endl;
	success = false;
      }
      if (_debug)
	cout << "        created final state particle[" << i << "]: "
	     << daughters.back()->qnSummary() << endl;
    }
  }
  if (parentKey.exists("isobars")) {
    if (_debug)
      printInfo << "reading isobars from '" << parentKey.getPath() << "':" << endl;
    const Setting& isobarKeys = parentKey["isobars"];
    for (int i = 0; i < isobarKeys.getLength(); ++i) {
      const Setting& isobarKey = isobarKeys[i];
      // mandatory fields
      string name;
      if (isobarKey.lookupValue("name", name))
	daughters.push_back(createParticle(name));
      else {
	printWarn << "isobar group '" << isobarKey.getPath() << "' "
		  << "has no name field." << endl;
	success = false;
      }
      if (_debug)
	cout << "        created isobar[" << i << "]: " << daughters.back()->qnSummary() << endl;
      success &= constructVertex(isobarKeys[i], daughters.back());
    }
  }
  if (!success) {
    printWarn << "problems constructing decay chain from key "
	      << "file '" << parentKey.getSourceFile() << "" << endl;
    return false;
  }
  if (daughters.size() != 2) {
    printWarn << "number of daughters specified in '" << parentKey.getPath() << "'"
	      << "is " << daughters.size() << " != 2. " << flush;
    if (daughters.size() > 2)
      cout << "constructing isobar vertex using " << daughters[0]->qnSummary() << " "
	   << "and " << daughters[2]->qnSummary() << endl;
    else if (daughters.size() < 2) {
      cout << "cannot construct isobar vertex." << endl;
      return false;
    }
  }
  
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
  _decayVertices.push_back(createIsobarDecayVertex(parentParticle, daughters[0], daughters[1],
						   L, S, massDep));
  if (_debug)
    printInfo << "constructed " << _decayVertices.back();

  return true;
}
