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
//      container class that holds all external information for
//      amplitude calculation
//      internally the decay process is represented as a graph
//      the graph is constraint to contain exactly one production
//      vertex and at least one interaction vertex; in addtion for
//      each final-state particle a corresponding final-state vertex
//      is created for internal use
//
//      "final-state" particles are the measured decay daughters;
//      additional final-state particles that belong to the production
//      process are handled by the production vertex
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <vector>
#include <map>
#include <algorithm>

#include "TClonesArray.h"
#include "TClass.h"
#include "TObjString.h"
#include "TVector3.h"

#include "reportingUtilsRoot.hpp"
#include "conversionUtils.hpp"
#include "decayTopology.h"


using namespace std;
using namespace boost;
using namespace rpwa;


bool decayTopology::_debug = false;


decayTopology::decayTopology()
	: decayTopologyGraphType(),
	  _prodVertex           (),
	  _decayVertices        (),
	  _fsParticles          (),
	  _fsDataPartIndexMap   (),
	  _fsDataPartMomCache   ()
{
}


decayTopology::decayTopology(const productionVertexPtr&          productionVertex,
                             const vector<interactionVertexPtr>& decayVertices,
                             const vector<particlePtr>&          fsParticles,
                             const bool                          performTopologyCheck)
{
	constructDecay(productionVertex, decayVertices, fsParticles, performTopologyCheck);
}


decayTopology::decayTopology(const decayTopology& topo)
	: decayTopologyGraphType()
{
	*this = topo;
}


decayTopology::decayTopology(const decayTopologyGraphType& graph)
{
	*this = graph;
}


decayTopology::~decayTopology()
{ }


decayTopology&
decayTopology::operator =(const decayTopology& topo)
{
	if (this != &topo) {
		decayTopologyGraphType::operator =(topo);
		_prodVertex         = topo._prodVertex;
		_decayVertices      = topo._decayVertices;
		_fsParticles        = topo._fsParticles;
		_fsDataPartIndexMap = topo._fsDataPartIndexMap;
		_fsDataPartMomCache = topo._fsDataPartMomCache;
	}
	return *this;
}


decayTopology&
decayTopology::operator =(const decayTopologyGraphType& graph)
{
	if (this != &graph) {
		decayTopologyGraphType::operator =(graph);
		buildInternalData();
	}
	return *this;
}


decayTopology*
decayTopology::doClone(const bool cloneFsParticles,
                       const bool cloneProdKinematics) const
{
	if (_debug)
		printDebug << "cloning decay topology '" << name() << "' "
		           << ((cloneFsParticles   ) ? "in" : "ex") << "cluding final-state particles, "
		           << ((cloneProdKinematics) ? "in" : "ex") << "cluding production kinematics particles"
		           << endl;
	// copy graph data structure
	decayTopology* topoClone = new decayTopology(*this);
	// clone nodes
	nodeIterator iNd, iNdEnd;
	vector<interactionVertexPtr> newVerts;
	for (tie(iNd, iNdEnd) = topoClone->nodes(); iNd != iNdEnd; ++iNd) {
		const interactionVertexPtr& vert = topoClone->vertex(*iNd);
		// clone vertex
		if (   isDecayVertex(vert)
		    or isProductionVertex(vert)
		    or (cloneFsParticles and isFsVertex(vert)))
			topoClone->isDecayVertex(vert);
			newVerts.push_back(topoClone->cloneNode(*iNd));
	}
	// looping over incoming particles of respective vertices and clone them
	// this also clones dangling particles that have no associated edge
	//!!! what about outgoing particles in production vertex other than X?
	for (unsigned int i = 0; i < newVerts.size(); ++i) {
		if (not topoClone->isProductionVertex(newVerts[i]) or cloneProdKinematics)
			for (unsigned int j = 0; j < newVerts[i]->nmbInParticles(); ++j) {
				const particlePtr part = newVerts[i]->inParticles()[j];
				if (not topoClone->isFsParticle(part) or cloneFsParticles) {
					if (topoClone->isEdge(part))
						topoClone->cloneEdge(topoClone->edge(part));
					else {
						// dangling particle
						if (_debug)
							printDebug << "cloning dangling " << *part << endl;
						// clone particle and replace it in array of incoming particles in vertex
						particlePtr newPart(part->clone());
						newVerts[i]->inParticles()[j] = newPart;
					}
				}
			}
	}
	topoClone->name() = "topoClone";
	return topoClone;
}


void
decayTopology::clear()
{
	decayTopologyGraphType::clear();
	_prodVertex.reset();
	_decayVertices.clear();
	_fsParticles.clear();
	_fsDataPartIndexMap.clear();
	_fsDataPartMomCache.clear();
}


map<string, unsigned int>
decayTopology::nmbIndistFsParticles() const
{
	map<string, unsigned int> partMult;
	for (unsigned int i = 0; i < nmbFsParticles(); ++i) {
		const string partName = fsParticles()[i]->name();
		map<string, unsigned int>::iterator entry = partMult.find(partName);
		if (entry != partMult.end())
			++(entry->second);
		else
			partMult[partName] = 1;
	}
	return partMult;
}


int
decayTopology::fsParticlesIntrinsicParity() const
{
	// calculate product of final-state particles
	int intrinsicParity = 1;
	for (unsigned int i = 0; i < nmbFsParticles(); ++i)
		intrinsicParity *= fsParticles()[i]->P();
	return intrinsicParity;
}


int
decayTopology::spaceInvEigenValue() const
{
	// parity of X is P = P_spatial * P_intrinsic; we want to know P_spatial
	return XParticle()->P() / fsParticlesIntrinsicParity();
}


int
decayTopology::reflectionEigenValue() const
{
	// eigenvalue of reflection through production plane is r_spatial / P_intrinsic
	// where r = -1 / refl is the reflectivity eiganvalue
	return -1 / (XParticle()->reflectivity() * fsParticlesIntrinsicParity());
}


void
decayTopology::transformFsParticles(const TLorentzRotation& L)
{
	for (unsigned int i = 0; i < nmbFsParticles(); ++i)
		fsParticles()[i]->transform(L);
}


bool
decayTopology::isVertex(const interactionVertexPtr& vert) const
{
	if (not vert)
		return false;
	nodeIterator iNd, iNdEnd;
	for (tie(iNd, iNdEnd) = nodes(); iNd != iNdEnd; ++iNd)
		if (vert == vertex(*iNd))
			return true;
	return false;
}


bool
decayTopology::isParticle(const particlePtr& part) const
{
	if (not part)
		return false;
	edgeIterator iEd, iEdEnd;
	for (tie(iEd, iEdEnd) = edges(); iEd != iEdEnd; ++iEd)
		if (part == particle(*iEd))
			return true;
	return false;
}


bool
decayTopology::isProductionVertex(const interactionVertexPtr& vert) const
{
	if (not vert or not isVertex(vert))
		return false;
	return (vert == _prodVertex);
}


bool
decayTopology::isDecayVertex(const interactionVertexPtr& vert) const
{
	if (not vert)
		return false;
	if (not isVertex(vert))
		return false;
	if (isProductionVertex(vert) or isFsVertex(vert))
		return false;
	return true;
}


int
decayTopology::decayVertexIndex(const interactionVertexPtr& vert) const
{
	if (not vert or not isDecayVertex(vert))
		return -1;
	for (unsigned int i = 0; i < nmbDecayVertices(); ++i)
		if (vert == decayVertices()[i])
			return i;
	return -1;
}


bool
decayTopology::isFsVertex(const interactionVertexPtr& vert) const
{
	if (not vert)
		return false;
	if (dynamic_pointer_cast<fsVertex>(vert))
		return isVertex(vert);
	return false;
}


bool
decayTopology::isFsParticle(const particlePtr& part) const
{
	if (not part or not isParticle(part))
		return false;
	for (unsigned int i = 0; i < nmbFsParticles(); ++i)
		if (part == fsParticles()[i])
			return true;
	return false;
}


int
decayTopology::fsParticlesIndex(const particlePtr& part) const
{
	if (not part or not isParticle(part))
		return -1;
	for (unsigned int i = 0; i < nmbFsParticles(); ++i)
		if (part == fsParticles()[i])
			return i;
	return -1;
}


bool
decayTopology::checkTopology() const
{
	bool topologyIsOkay = true;
	// make sure there are no dangling particles
	nodeIterator iNd, iNdEnd;
	for (tie(iNd, iNdEnd) = nodes(); iNd != iNdEnd; ++iNd) {
		const interactionVertexPtr& vert = vertex(*iNd);
		if (not vert) {
			printWarn << "node[" << *iNd << "] has vertex null pointer" << endl;
			topologyIsOkay = false;
			continue;
		}
		if (isProductionVertex(vert)) {
			// for production vertex only X particle has associated edge
			if (not isEdge(static_pointer_cast<rpwa::productionVertex>(vert)->XParticle())) {
				printWarn << "X particle of " << *vert << " has no associated edge" << endl;
				topologyIsOkay = false;
			} else if (_debug)
				printDebug << "success: X particle of " << *vert << " has associated edge" << endl;
		} else {
			for (unsigned int i = 0; i < vert->nmbInParticles(); ++i)
				if (not isEdge(vert->inParticles()[i])) {
					printWarn << "incoming particle[" << i << "] of " << *vert
					          << " has no associated edge" << endl;
					topologyIsOkay = false;
				} else if (_debug)
					printDebug << "success: incoming particle[" << i << "] of " << *vert
					           << " has associated edge" << endl;
			for (unsigned int i = 0; i < vert->nmbOutParticles(); ++i)
				if (not isEdge(vert->outParticles()[i])) {
					printWarn << "outgoing particle[" << i << "] of " << *vert
					          << " has no associated edge" << endl;
					topologyIsOkay = false;
				} else if (_debug)
					printDebug << "success: outgoing particle[" << i << "] of " << *vert
					           << " has associated edge" << endl;
		}
	}
	// check production vertex
	if (not productionVertex()) {
		printWarn << "production vertex is null pointer" << endl;
		topologyIsOkay = false;
	} else {
		if (not isNode(_prodVertex) or (vertex(node(_prodVertex)) != _prodVertex)) {
			printWarn << *_prodVertex << " has no associated node" << endl;
			topologyIsOkay = false;
		} else if (_debug)
			printDebug << "success: " << *_prodVertex << " has associated node" << endl;
		if (nmbOutEdges(_prodVertex) != 1) {
			printWarn << *_prodVertex << " has " << nmbOutEdges(_prodVertex) << " != 1 "
			          << "outgoing edges" << endl;
			topologyIsOkay = false;
		} else if (nmbInEdges(_prodVertex) != 0) {
			printWarn << *_prodVertex << " has " << nmbInEdges(_prodVertex) << " != 0 "
			          << "incoming edges" << endl;
			topologyIsOkay = false;
		} else if (_debug)
			printDebug << "success: " << *_prodVertex
			           << " has exactly 1 outgoing and no incoming edges" << endl;
	}
	// make sure that all nodes are reachable from top node
	vector<nodeDesc> sortedNds;
	unsigned int     nmbVertices = nmbDecayVertices() + nmbFsParticles();
	if (productionVertex()) {
		sortedNds = sortNodesDfs(node(productionVertex()));
		++nmbVertices;
	} else {
		nodeDesc topNd = topNode();
		if (_debug)
			printDebug << "using node[" << topNd << "] as top node" << endl;
		sortedNds = sortNodesDfs(topNd);
	}
	if (sortedNds.size() != nmbVertices) {
		printWarn << "number of nodes reachable from top node is "
		          << sortedNds.size() << " != " << nmbVertices << endl;
		topologyIsOkay = false;
	} else if (_debug)
		printDebug << "success: all nodes are reachable from top node" << endl;
	// make sure all interaction vertices have corresponding nodes
	for (unsigned int i = 0; i < nmbDecayVertices(); ++i) {
		const interactionVertexPtr& vert = decayVertices()[i];
		if (not vert) {
			printWarn << "interaction vertex[" << i << "] is null pointer" << endl;
			topologyIsOkay = false;
			continue;
		}
		if (not isNode(vert) or (vertex(node(vert)) != vert)) {
			printWarn << *vert << " at " << vert << " has no associated node" << endl;
			topologyIsOkay = false;
		} else if (_debug)
			printDebug << "success: " << *vert << " has associated node" << endl;
		if (nmbOutEdges(vert) < 1) {
			printWarn << *vert << " has " << nmbOutEdges(vert) << " < 1 "
			          << "outgoing edges" << endl;
			topologyIsOkay = false;
		} else if (nmbInEdges(vert) < 1) {
			printWarn << *vert << " has " << nmbInEdges(vert) << " < 1 "
			          << "incoming edges" << endl;
			topologyIsOkay = false;
		} else if (_debug)
			printDebug << "success: " << *vert
			           << " has at least 1 outgoing and 1 incoming edge" << endl;
	}
	// make sure all final-state particles have corresponding edges
	for (unsigned int i = 0; i < nmbFsParticles(); ++i) {
		const particlePtr& part = fsParticles()[i];
		if (not part) {
			printWarn << "final-state particle[" << i << "] is null pointer" << endl;
			topologyIsOkay = false;
			continue;
		}
		if (not isEdge(part) or (particle(this->edge(part)) != part)) {
			printWarn << "final state " << *part << " has no associated edge" << endl;
			topologyIsOkay = false;
		} else if (_debug)
			printDebug << "success: final state " << *part << " has associated edge" << endl;
		if (isEdge(part) and not isFsVertex(toVertex(part))) {
			printWarn << "vertex associated to final-state particle "
			          << "'" << part->name() << "' is not a final-state vertex" << endl;
			topologyIsOkay = false;
		} else if (_debug)
			printDebug << "success: vertex associated to final-state particle "
			           << "'" << part->name() << "' is a final-state vertex" << endl;
	}
	// make sure edges connect the right vertices
	edgeIterator iEd, iEdEnd;
	for (tie(iEd, iEdEnd) = edges(); iEd != iEdEnd; ++iEd) {
		const particlePtr& part = particle(*iEd);
		// check outgoing particles of the vertex the edge is coming from
		const interactionVertexPtr& fromVert     = fromVertex(*iEd);
		bool                        isInFromVert = false;
		for (unsigned int i = 0; i < fromVert->nmbOutParticles(); ++i)
			if (part == fromVert->outParticles()[i]) {
				isInFromVert = true;
				break;
			}
		if (not isInFromVert) {
			printWarn << "particle '" << part->name() << "' = " << part << " associated to edge " << *iEd
			          << " is not in list of outgoing particles of " << *fromVert << " at " << fromVert
			          << ": ";
			for (unsigned int i = 0; i < fromVert->nmbOutParticles(); ++i)
				cout << "[" << i << "] = " << fromVert->outParticles()[i] << "  ";
			cout << endl;
			topologyIsOkay = false;
		} else if (_debug)
			printDebug << "success: found particle '" << part->name() << "' "
			           << "in list of outgoing particles of vertex " << *fromVert << endl;
		// check outgoing particles of the vertex the edge is coming from
		const interactionVertexPtr& toVert     = toVertex  (*iEd);
		bool                        isInToVert = false;
		for (unsigned int i = 0; i < toVert->nmbInParticles(); ++i)
			if (part == toVert->inParticles()[i]) {
				isInToVert = true;
				break;
			}
		if (not isInToVert) {
			printWarn << "particle '" << part->name() << "' = " << part << " associated to edge " << *iEd
			          << " is not in list of incoming particles of " << *toVert << " at " << toVert
			          << ": ";
			for (unsigned int i = 0; i < toVert->nmbOutParticles(); ++i)
				cout << "[" << i << "] = " << toVert->outParticles()[i] << "  ";
			cout << endl;
			topologyIsOkay = false;
		} else if (_debug)
			printDebug << "success: found particle '" << part->name() << "' "
			           << "in list of incoming particles of vertex " << *toVert << endl;
	}
	// make sure final-state indices make sense
	// two cases are allowed:
	// i)  all final-state particles have indices at default value -1
	// ii) all final-state particles have an index != -1; the indices
	// range from 0 to (number of final-state particles - 1)
	// sort final-state particles w.r.t. name
	map<string, vector<particlePtr> > fsPartSorted;
	bool                              allIndicesAtDefault = true;
	for (unsigned int i = 0; i < nmbFsParticles(); ++i) {
		const particlePtr& part = fsParticles()[i];
		fsPartSorted[part->name()].push_back(part);
		if (part->index() >= 0)
			allIndicesAtDefault = false;
	}
	if (allIndicesAtDefault) {
		if (_debug)
			printDebug << "success: all final-state particles have default indices" << endl;
	} else
		for (map<string, vector<particlePtr> >::iterator i = fsPartSorted.begin();
		     i != fsPartSorted.end(); ++i) {
			// sort particles with same name by their index
			vector<particlePtr>& parts = i->second;
			sort(parts.begin(), parts.end(), rpwa::compareIndicesAsc);
			// check that indices are consecutive and include all unmbers
			// from 0 to (number of final-state particles - 1)
			for (unsigned int j = 0; j < parts.size(); ++j)
				if (parts[j]->index() != (int)j) {
					printWarn << "indices for final-state particles '" << i->first << "' are not consecutive. "
					          << "expected index [" << j << "], found [" << parts[j]->index() << "]." << endl;
					topologyIsOkay = false;
				} else if (_debug)
					printDebug << "success: final-state particle '" << i->first << "' "
					           << "has expected index [" << j << "]" << endl;
		}
	if (_debug)
		printDebug << "decay topology " << ((topologyIsOkay) ? "passed" : "did not pass")
		           << " all tests" << endl;
	return topologyIsOkay;
}


decayTopology
decayTopology::subDecay(const nodeDesc& startNd,
                        const bool      linkToMotherTopo)
{
	decayTopology subTopo(dfsSubGraph(startNd, linkToMotherTopo));
	subTopo.name() = "subdecay";
	return subTopo;
}


void
decayTopology::addDecay(const decayTopology& topo)
{
	decayTopologyGraphType::addGraph(topo);
	buildInternalData();
}


void decayTopology::setProductionVertex(const productionVertexPtr& productionVertex)
{
	if (not productionVertex) {
		printErr << "null pointer for production vertex. Aborting..." << endl;
		throw;
	}
	if (not productionVertex->XParticle()) {
		printErr << "null pointer for particle[0] coming out of production vertex. Aborting..." << endl;
		throw;
	}
	name() = "\"" + productionVertex->XParticle()->qnSummary() + "\"";
	if (not _prodVertex) {
		// topology does not have production vertex -> create graph node
		addVertex(productionVertex);
	} else {
		// delete existing production vertex -> update graph node
		vertex(node(_prodVertex)) = productionVertex;
		_prodVertex.reset();
	}
	_prodVertex = productionVertex;
}


bool
decayTopology::initKinematicsData(const vector<string>& prodKinPartNamesVec,
                                  const vector<string>& decayKinPartNamesVec)
{
	TClonesArray prodKinPartNames("TObjString", prodKinPartNamesVec.size());
	TClonesArray decayKinPartNames("TObjString", decayKinPartNamesVec.size());
	for(unsigned int i = 0; i < prodKinPartNamesVec.size(); ++i) {
		new (prodKinPartNames[i]) TObjString(prodKinPartNamesVec[i].c_str());
	}
	for(unsigned int i = 0; i < decayKinPartNamesVec.size(); ++i) {
		new (decayKinPartNames[i]) TObjString(decayKinPartNamesVec[i].c_str());
	}
	return initKinematicsData(prodKinPartNames, decayKinPartNames);
}


bool
decayTopology::initKinematicsData(const TClonesArray& prodKinPartNames,
                                  const TClonesArray& decayKinPartNames)
{
	// two modes are supported:
	// i) all final-state particles have indices at default value -1: in
	//    case there are more than one final-state particles of the same
	//    type, the momenta are assigned in the same order the particles
	//    are stored in the final-state particle array
	// ii) all final-state particles have an index != -1: the data are
	//     assigned according to the indices of the particles
	// in case only part of the final-state particles have non-default
	// indices, the result is undefined

	// init production kinematics
	bool success = productionVertex()->initKinematicsData(prodKinPartNames);

	// check decay kinematics data
	const string partClassName = decayKinPartNames.GetClass()->GetName();
	if (partClassName != "TObjString") {
		printWarn << "decay kinematics particle names are of type '" << partClassName
		          << "' and not TObjString. cannot read decay kinematics." << endl;
		success = false;
	}
	const int nmbFsPart = decayKinPartNames.GetEntriesFast();
	if ((nmbFsPart < 0) or ((unsigned int)nmbFsPart != nmbFsParticles())) {
		printWarn << "array of decay kinematics particle names has wrong size: "
		          << nmbFsPart << " (expected " << nmbFsParticles() << ")" << endl;
		success = false;
	}
	if (not success)
		return false;
	// sort decay kinematics particles w.r.t. particle name
	map<string, vector<unsigned int> > fsDataPartIndices;
	map<string, unsigned int>          fsDataPartCount;
	for (unsigned int i = 0; i < nmbFsParticles(); ++i) {
		const string name = ((TObjString*)decayKinPartNames[i])->GetString().Data();
		fsDataPartIndices[name].push_back(i);
		fsDataPartCount[name] = 0;
	}

	// create index map: index of final-state particle in decay graph -> index in data array
	_fsDataPartIndexMap.clear();
	for (unsigned int i = 0; i < nmbFsParticles(); ++i) {
		const particlePtr& part     = fsParticles()[i];
		const string       partName = part->name();
		map<string, vector<unsigned int> >::const_iterator entry = fsDataPartIndices.find(partName);
		if (entry != fsDataPartIndices.end()) {
			int partIndex = part->index();
			if (partIndex < 0) {
				partIndex = fsDataPartCount[partName];
				++fsDataPartCount[partName];
			}
			if ((unsigned int)partIndex < entry->second.size()) {
				if (_debug)
					printDebug << "assigning decay kinematics data index [" << entry->second[partIndex] << "] "
					           << "to final-state particle '" << partName << "'[" << partIndex << "] "
					           << "at index [" << i << "]" << endl;
				_fsDataPartIndexMap[i] = entry->second[partIndex];
			} else {
				printWarn << "index [" << partIndex << "] for final-state particle "
				          << "'" << partName << "' out of range. input data contain only "
				          << entry->second.size() << " entries for this particle type." << endl;
				success = false;
			}
		} else {
			printWarn << "cannot find entry for final-state particle '" << part->name() << "' "
			          << "in input data." << endl;
			success = false;
		}
	}
	if (_fsDataPartIndexMap.size() != nmbFsParticles()) {
		printWarn << "could not find all final-state particles in input data." << endl;
		success = false;
	}

	return success;
}


bool
decayTopology::readKinematicsData(const TClonesArray& prodKinMomenta,
                                  const TClonesArray& decayKinMomenta)
{
	// set production kinematics
	bool success = productionVertex()->readKinematicsData(prodKinMomenta);

	// check momentum array
	const int nmbFsPart = decayKinMomenta.GetEntriesFast();
	if ((nmbFsPart < 0) or ((unsigned int)nmbFsPart != nmbFsParticles())) {
		printWarn << "array of decay kinematics particle momenta has wrong size: "
		          << nmbFsPart << " (expected " << nmbFsParticles() << "). "
		          << "cannot read decay kinematics." << endl;
		success = false;
	}
	if (not success)
		return false;

	// set decay kinematics
	for (unsigned int i = 0; i < nmbFsParticles(); ++i) {
		const particlePtr& part      = fsParticles()[i];
		const unsigned int partIndex = _fsDataPartIndexMap[i];
		const TVector3*    mom       = dynamic_cast<TVector3*>(decayKinMomenta[partIndex]);
		if (not mom) {
			printWarn << "decay kinematics data entry [" << partIndex << "] is not of type TVector3. "
			          << "cannot read decay kinematics momentum for particle '" << part->name() << "'. "
			          << "skipping." << endl;
			success = false;
			continue;
		}
		if (_debug)
			printDebug << "setting momentum of final-state particle '" << part->name() << "' "
			           << "at index [" << i << "] to " << *mom << " GeV "
			           << "at input data index [" << partIndex << "]" << endl;
		part->setMomentum(*mom);
	}
	fillKinematicsDataCache();

	return success;
}


void
decayTopology::fillKinematicsDataCache()
{
	// adjust cache size
	_fsDataPartMomCache.resize(nmbFsParticles(), TVector3());
	// set decay kinematics
	for (unsigned int i = 0; i < nmbFsParticles(); ++i) {
		const particlePtr& part = fsParticles()[i];
		_fsDataPartMomCache[i] = part->momentum();
	}
}


bool
decayTopology::revertMomenta()
{
	// revert production kinematics
	bool success = productionVertex()->revertMomenta();
	// revert decay kinematics
	if (_fsDataPartMomCache.size() != nmbFsParticles()) {
		if (_debug)
			printWarn << "cache for final-state particle momenta has size "
			          << _fsDataPartMomCache.size() << " != # of final-state particles = "
			          << nmbFsParticles() << endl;
		return false;
	}
	for (unsigned int i = 0; i < nmbFsParticles(); ++i) {
		const particlePtr& part = fsParticles()[i];
		part->setMomentum(_fsDataPartMomCache[i]);
		if (_debug)
			printDebug << "resetting momentum of final-state particle '" << part->name() << "'"
			           << "[" << i << "] to " << _fsDataPartMomCache[i] << " GeV" << endl;
	}
	return success;
}


bool
decayTopology::revertMomenta(const vector<unsigned int>& fsPartPermMap)  // final-state permutation map
{
	// revert production kinematics
	bool success = productionVertex()->revertMomenta();
	// revert decay kinematics
	if (   (_fsDataPartMomCache.size() != nmbFsParticles())
	    or (fsPartPermMap.size      () != nmbFsParticles())) {
		if (_debug)
			printWarn << "cache size for final-state particle momenta (= "
			          << _fsDataPartMomCache.size() << ") or size of permutation map (= "
			          << fsPartPermMap.size() << ") do not match # of final-state particles (= "
			          << nmbFsParticles() << ")" << endl;
		return false;
	}
	for (unsigned int i = 0; i < nmbFsParticles(); ++i) {
		const unsigned int newIndex = fsPartPermMap[i];
		const particlePtr& part     = fsParticles()[i];
		part->setMomentum(_fsDataPartMomCache[newIndex]);
		if (_debug)
			printDebug << "(re)setting momentum of final-state particle "
			           << "'" << part->name() << "'[" << i << "] "
			           << "to that of '" << fsParticles()[newIndex]->name()
			           << "'[" << newIndex << "] = " << _fsDataPartMomCache[newIndex] << " GeV" << endl;
	}
	return success;
}


void
decayTopology::saveDecayToVertices(const decayTopologyPtr decay)
{
	productionVertex()->setDecay(decay);
	const vector<interactionVertexPtr>& decayVertices = this->decayVertices();
	for(unsigned int i = 0; i < decayVertices.size(); ++i) {
		decayVertices[i]->setDecay(decay);
	}
}


ostream&
decayTopology::print(ostream& out) const
{
	// print nodes
	out << "decay topology '" << name() << "' has " << nmbNodes() << " node(s):" << endl;
	if (productionVertex())
		out << "    production  node[" << node(productionVertex()) << "] = "
		    << productionVertex() << ": " << *productionVertex() << endl;
	else
		out << "    topology has no production node." << endl;
	for (unsigned int i = 0; i < nmbDecayVertices(); ++i)
		out << "    interaction node[" << node(decayVertices()[i]) << "] = "
		    << decayVertices()[i] << ": " << *decayVertices()[i] << endl;
	nodeIterator iNd, iNdEnd;
	for (tie(iNd, iNdEnd) = nodes(); iNd != iNdEnd; ++iNd)
		if (isFsVertex(vertex(*iNd)))
			out << "    final-state node[" << *iNd << "] = " << vertex(*iNd) << ": "
			    << *vertex(*iNd) << endl;
	// print edges
	out << "decay topology '" << name() << "' has " << nmbEdges() << " edge(s):" << endl;
	edgeIterator iEd, iEdEnd;
	for (tie(iEd, iEdEnd) = edges(); iEd != iEdEnd; ++iEd) {
		const particlePtr& part = particle(*iEd);
		out << "    edge[" << *iEd << "] = [" << fromNode(*iEd) << ", "
		    << toNode(*iEd) << "] = " << part << ": '" << part->name() << "'";
		if (isFsParticle(part))
			out << " (final state)";
		out << endl;
	}
	return out;
}


ostream&
decayTopology::printProdKinParticles(ostream& out) const
{
	out << "production kinematics particles:" << endl;
	for (unsigned int i = 0; i < productionVertex()->nmbInParticles(); ++i) {
		const particlePtr& part = productionVertex()->inParticles()[i];
		out << "    incoming particle[" << i << "]: " << part->qnSummary() << ", "
		    << "index = " << part->index() << ", p = " << part->lzVec() << " GeV" << endl;
	}
	for (unsigned int i = 0; i < productionVertex()->nmbOutParticles(); ++i) {
		const particlePtr& part = productionVertex()->outParticles()[i];
		out << "    outgoing particle[" << i << "]: " << part->qnSummary() << ", "
		    << "index = " << part->index() << ", p = " << part->lzVec() << " GeV" << endl;
	}
	return out;
}


ostream&
decayTopology::printDecayKinParticles(ostream& out) const
{
	out << "decay kinematics particles:" << endl;
	for (unsigned int i = 0; i < nmbFsParticles(); ++i) {
		const particlePtr& part = fsParticles()[i];
		out << "    particle[" << i << "]: " << part->qnSummary() << ", "
		    << "index = " << part->index() << ", p = " << part->lzVec() << " GeV" << endl;
	}
	return out;
}


decayTopology&
decayTopology::constructDecay(const productionVertexPtr&          productionVertex,
                              const vector<interactionVertexPtr>& decayVertices,
                              const vector<particlePtr>&          fsParticles,
                              const bool                          performTopologyCheck)
{
	clear();
	const unsigned int nmbDecayVert = decayVertices.size();
	const unsigned int nmbFsPart    = fsParticles.size();
	if (_debug)
		printDebug << "constructing decay topology with " << nmbDecayVert   << " decay vertices and "
		           << nmbFsPart << " final-state particles" << endl;
	if (nmbFsPart < 1) {
		printErr << "cannot construct decay topology without final-state particles. Aborting..." << endl;
		throw;
	}
	if (nmbDecayVert < 1) {
		printWarn << "need at least production and X-decay vertex "
		          << "to construct decay topology. Aborting..." << endl;
		throw;
	}

	// create graph node for production vertex and store pointer
	if (not productionVertex) {
		printErr << "null pointer for production vertex. Aborting..." << endl;
		throw;
	}
	if (not productionVertex->XParticle()) {
		printErr << "null pointer for X particle coming out of production vertex. Aborting..." << endl;
		throw;
	}
	name() = "\"" + productionVertex->XParticle()->qnSummary() + "\"";
	_prodVertex = productionVertex;
	addVertex(productionVertex);
	// create graph nodes for interaction vertices and store pointers
	for (unsigned int i = 0; i < nmbDecayVert; ++i) {
		if (not decayVertices[i]) {
			printErr << "null pointer for decay vertex[" << i << "]. Aborting..." << endl;
			throw;
		}
		addVertex(decayVertices[i]);
	}
	// create final-state nodes and vertices
	for (unsigned int i = 0; i < nmbFsPart; ++i) {
		if (not fsParticles[i]) {
			printErr << "null pointer for final-state particle[" << i << "]. Aborting..." << endl;
			throw;
		}
		const fsVertexPtr& fsVert = createFsVertex(fsParticles[i]);
		addVertex(fsVert);
	}
	// memorize depth-first sorted interaction vertices and final-state particles
	vector<nodeDesc> sortedNds = sortNodesDfs(node(_prodVertex));
	for (unsigned int i = 0; i < sortedNds.size(); ++i) {
		const interactionVertexPtr& vert   = vertex(sortedNds[i]);
		const fsVertexPtr&          fsVert = dynamic_pointer_cast<fsVertex>(vert);
		if (fsVert) {
			const particlePtr& part = fsVert->inParticles()[0];
			for (unsigned int i = 0; i < nmbFsPart; ++i)
				if (part == fsParticles[i]) {
					_fsParticles.push_back(part);
					break;
				}
		} else if (vert != _prodVertex) {
			for (unsigned int i = 0; i < nmbDecayVert; ++i)
				if (vert == decayVertices[i]) {
					_decayVertices.push_back(vert);
					break;
				}
		}
	}
	// check that topology makes sense
	if (performTopologyCheck and not checkTopology()) {
		printErr << "topology has problems that need to be fixed. Aborting..." << endl;
		throw;
	}
	return *this;
}


void
decayTopology::buildInternalData()
{
	bool success = true;
	// find production vertex
	// assumes that production vertex is the only vertex in graph that
	// has no incoming edges and exactly one outgoing edge
	_prodVertex.reset();
	unsigned int nmbProdVertCandidates = 0;
	nodeIterator iNd, iNdEnd;
	for (tie(iNd, iNdEnd) = nodes(); iNd != iNdEnd; ++iNd)
		if ((nmbInEdges(*iNd) == 0) and (nmbOutEdges(*iNd) == 1)) {
			_prodVertex = dynamic_pointer_cast<rpwa::productionVertex>(vertex(*iNd));
			if (_prodVertex)
				++nmbProdVertCandidates;
		}
	// set final-state particles
	_fsParticles.clear();
	vector<nodeDesc> sortedNds;
	if (_prodVertex)
		sortedNds = sortNodesDfs(node(_prodVertex));
	else {
		// take non-FS vertex with no incoming particle as top node
		nodeDesc     topNode;
		unsigned int nmbTopVertCandidates = 0;
		for (tie(iNd, iNdEnd) = nodes(); iNd != iNdEnd; ++iNd)
			if ((nmbInEdges(*iNd) == 0) and not isFsVertex(vertex(*iNd))) {
				topNode = *iNd;
				++nmbTopVertCandidates;
			}
		if (nmbTopVertCandidates == 1)
			sortedNds = sortNodesDfs(topNode);
		else {
			if (_debug)
				printWarn << "ill-formed graph. vertex order will be undefined." << endl;
			for (tie(iNd, iNdEnd) = nodes(); iNd != iNdEnd; ++iNd)
				sortedNds.push_back(*iNd);
		}
	}
	for (unsigned int i = 0; i < sortedNds.size(); ++i) {
		const interactionVertexPtr& vert = vertex(sortedNds[i]);
		if (isFsVertex(vert))
			_fsParticles.push_back(vert->inParticles()[0]);
	}
	if (nmbProdVertCandidates > 1) {
		printWarn << "found " << nmbProdVertCandidates << " instead of 0 or 1 candidate "
		          << "for production vertex in graph '" << name() << "'" << endl;
		success = false;
	}
	if (_fsParticles.size() < 1) {
		printWarn << "cannot find final-state particles in graph '" << name() << "'" << endl;
		success = false;
	}
	if (not success) {
		printErr << "cannot construct decay topology from graph '" << name() << "'. Aborting..." << endl;
		throw;
	}
	// find interaction vertices
	_decayVertices.clear();
	for (unsigned int i = 0; i < sortedNds.size(); ++i) {
		const interactionVertexPtr& vert = vertex(sortedNds[i]);
		if (isDecayVertex(vert))
			_decayVertices.push_back(vert);
	}
}


interactionVertexPtr
decayTopology::cloneNode(const nodeDesc& nd,
                         const bool      cloneInParticles,
                         const bool      cloneOutParticles)
{
	const interactionVertexPtr  vert        = vertex(nd);  // this must not be a reference
	const bool                  isProdVert  = isProductionVertex(vert);
	const bool                  isDecayVert = isDecayVertex(vert);
	const interactionVertexPtr& newVert = decayTopologyGraphType::cloneNode(nd, cloneInParticles, cloneOutParticles);  // created new vertex for this node
	// update member variables
	if (isProdVert)
		_prodVertex = static_pointer_cast<rpwa::productionVertex>(newVert);
	else if (isDecayVert)
		for (unsigned int i = 0; i < nmbDecayVertices(); ++i)
			if (vert == _decayVertices[i])
				_decayVertices[i] = newVert;
	return newVert;
}


particlePtr
decayTopology::cloneEdge(const edgeDesc& ed)
{
	const particlePtr  part    = particle(ed);  // this must not be a reference
	const particlePtr& newPart = decayTopologyGraphType::cloneEdge(ed);
	// update member variable
	if (isFsParticle(part))
		for (unsigned int i = 0; i < nmbFsParticles(); ++i)
			if (part == fsParticles()[i])
				_fsParticles[i] = newPart;
	return newPart;
}
