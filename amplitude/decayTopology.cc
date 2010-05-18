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
//      container class that holds all external information for
//      amplitude calculation
//      internally the decay process is represented as a graph
//      the graph is constraint to contain exactly one production
//      vertex and at least one interaction vertex; in addtion for
//      each final state particle a corresponding final state vertex
//      is created for internal use
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

#include "utilities.h"
#include "decayTopology.h"

	
using namespace std;
using namespace boost;
using namespace rpwa;


bool decayTopology::_debug = false;


decayTopology::decayTopology()
  : decayTopologyGraphType(),
    _prodVertex           (),
    _intVertices          (),
    _fsParticles          ()
{
}


decayTopology::decayTopology(const interactionVertexPtr&              productionVertex,
			     const std::vector<interactionVertexPtr>& interactionVertices,
			     const std::vector<particlePtr>&          fsParticles)
{
  constructDecay(productionVertex, interactionVertices, fsParticles);
}


decayTopology::decayTopology(const decayTopology& topo)
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
    _prodVertex  = topo._prodVertex;
    _intVertices = topo._intVertices;
    _fsParticles = topo._fsParticles;
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
decayTopology::clone(const bool cloneFsParticles,
		      const bool cloneProductionVertex) const
{
  if (_debug)
    printInfo << "cloning decay topology '" << name() << "'; "
	      << "cloneFsParticles = "      << cloneFsParticles << ", "
	      << "cloneProductionVertex = " << cloneProductionVertex << endl;
  // copy graph data structure
  decayTopology* topoClone = new decayTopology(*this);
  nodeIterator iNd, iNdEnd;
  for (tie(iNd, iNdEnd) = topoClone->nodes(); iNd != iNdEnd; ++iNd) {
    const interactionVertexPtr v = topoClone->vertex(*iNd);
    // clone vertex
    interactionVertexPtr newV;
    if (   isInteractionVertex(v)
	|| (cloneFsParticles      && isFsVertex        (v))
        || (cloneProductionVertex && isProductionVertex(v)))
      newV = topoClone->cloneNode(*iNd);
    // clone particles by looping over the incoming particles of the respective vertices
    // this also copies dangling particles that have no associated edge
    //!!! what about beam particle in production vertex?
    if (newV && !isProductionVertex(v))
      for (unsigned int i = 0; i < newV->nmbInParticles(); ++i) {
	const particlePtr p = newV->inParticles()[i];
	if (!isFsParticle(p) || cloneFsParticles) {
	  if (topoClone->isEdge(p))
	    particlePtr newP = topoClone->cloneEdge(topoClone->edge(p));
	  else {
	    // dangling particle
	    if (_debug)
	      printInfo << "cloning dangling " << *p << endl;
	    // clone particle and store it in array of incoming particles in vertex
	    particlePtr newP(p->clone());
	    newV->inParticles()[i] = newP;
	  }
	}
      }
  }
  topoClone->buildReverseMaps();
  return topoClone;
}


void
decayTopology::clear()
{
  decayTopologyGraphType::clear();
  _prodVertex.reset();
  _intVertices.clear();
  _fsParticles.clear();
}


bool
decayTopology::isInteractionVertex(const interactionVertexPtr& vert) const
{
  if (isProductionVertex(vert) || isFsVertex(vert))
    return false;
  return true;
}


bool
decayTopology::isFsVertex(const interactionVertexPtr& vert) const
{
  if (dynamic_pointer_cast<fsVertex>(vert))
    return true;
  return false;
}


bool
decayTopology::isFsParticle(const particlePtr& part) const
{
  for (unsigned int i = 0; i < _fsParticles.size(); ++i)
    if (part == _fsParticles[i])
      return true;
  return false;
}


bool
decayTopology::checkTopology() const
{
  bool topologyIsOkay = true;
  // make sure there are no dangling particles
  nodeIterator iNd, iNdEnd;
  for (tie(iNd, iNdEnd) = nodes(); iNd != iNdEnd; ++iNd) {
    const interactionVertexPtr& vert = vertex(*iNd);
    if (!vert) {
      printWarn << "node[" << *iNd << "] has vertex null pointer" << endl;
      topologyIsOkay = false;
      continue;
    }
    if (vert != _prodVertex)  // incoming particles of production vertex have no edges
      for (unsigned int i = 0; i < vert->nmbInParticles(); ++i)
	if (!isEdge(vert->inParticles()[i])) {
	  printWarn << "incoming particle[" << i << "] of " << *vert
		    << " has no associated edge" << endl;
	  topologyIsOkay = false;
	} else if (_debug)
	  printInfo << "success: incoming particle[" << i << "] of " << *vert
		    << " has associated edge" << endl;
    for (unsigned int i = 0; i < vert->nmbOutParticles(); ++i)
      if (!isEdge(vert->outParticles()[i])) {
	printWarn << "outgoing particle[" << i << "] of " << *vert
		  << " has no associated edge" << endl;
	topologyIsOkay = false;
      } else if (_debug)
	printInfo << "success: outgoing particle[" << i << "] of " << *vert
		  << " has associated edge" << endl;
  }
  // check production vertex
  if (!_prodVertex) {
    printWarn << "production vertex is null pointer" << endl;
    topologyIsOkay = false;
  } else {
    if (!isNode(_prodVertex) || (vertex(node(_prodVertex)) != _prodVertex)) {
      printWarn << *_prodVertex << " has no associated node" << endl;
      topologyIsOkay = false;
    } else if (_debug)
      printInfo << "success: " << *_prodVertex << " has associated node" << endl;
    if (nmbOutEdges(_prodVertex) != 1) {
      printWarn << *_prodVertex << " has " << nmbOutEdges(_prodVertex) << " != 1 "
		<< "outgoing edges" << endl;
      topologyIsOkay = false;
    } else if (nmbInEdges(_prodVertex) != 0) {
      printWarn << *_prodVertex << " has " << nmbInEdges(_prodVertex) << " != 0 "
		<< "incoming edges" << endl;
      topologyIsOkay = false;
    } else if (_debug)
      printInfo << "success: " << *_prodVertex
		<< " has exactly 1 outgoing and no incoming edges" << endl;
  }
  // make sure that all nodes are reachable from top node
  vector<nodeDesc> sortedNds;
  unsigned int     nmbVertices = _intVertices.size() + _fsParticles.size();
  if (_prodVertex) {
    sortedNds = sortNodesDfs(node(_prodVertex));
    ++nmbVertices;
  } else {
    nodeDesc topNd = topNode();
    if (_debug)
      printInfo << "using node[" << topNd << "] as top node" << endl;
    sortedNds = sortNodesDfs(topNd);
  }
  if (sortedNds.size() != nmbVertices) {
    printWarn << "number of nodes reachable from top node is "
	      << sortedNds.size() << " != " << nmbVertices << endl;
    topologyIsOkay = false;
  } else if (_debug)
    printInfo << "success: all nodes are reachable from top node" << endl;
  // make sure all interaction vertices have corresponding nodes
  for (unsigned int i = 0; i < _intVertices.size(); ++i) {
    const interactionVertexPtr& vert = _intVertices[i];
    if (!vert) {
      printWarn << "interaction vertex[" << i << "] is null pointer" << endl;
      topologyIsOkay = false;
      continue;
    }
    if (!isNode(vert) || (vertex(node(vert)) != vert)) {
      printWarn << *vert << " has no associated node" << endl;
      topologyIsOkay = false;
    } else if (_debug)
      printInfo << "success: " << *vert << " has associated node" << endl;
    if (nmbOutEdges(vert) < 1) {
      printWarn << *vert << " has " << nmbOutEdges(vert) << " < 1 "
		<< "outgoing edges" << endl;
      topologyIsOkay = false;
    } else if (nmbInEdges(vert) < 1) {
      printWarn << *vert << " has " << nmbInEdges(vert) << " < 1 "
		<< "incoming edges" << endl;
      topologyIsOkay = false;
    } else if (_debug)
      printInfo << "success: " << *vert
		<< " has at least 1 outgoing and 1 incoming edge" << endl;
  }
  // make sure all final state particles have corresponding edges
  for (unsigned int i = 0; i < _fsParticles.size(); ++i) {
    const particlePtr& part = _fsParticles[i];
    if (!part) {
      printWarn << "final state particle[" << i << "] is null pointer" << endl;
      topologyIsOkay = false;
      continue;
    }
    if (!isEdge(part) || (particle(this->edge(part)) != part)) {
      printWarn << "final state " << *part << " has no associated edge" << endl;
      topologyIsOkay = false;
    } else if (_debug)
      printInfo << "success: final state " << *part << " has associated edge" << endl;
    if (isEdge(part) && !isFsVertex(toVertex(part))) {
      printWarn << "vertex associated to final state particle "
    		<< "'" << part->name() << "' is not a final state Vertex" << endl;
      topologyIsOkay = false;
    } else if (_debug)
      printInfo << "success: vertex associated to final state particle "
    		<< "'" << part->name() << "' is a final state vertex" << endl;
  }
  // make sure final state indices make sense
  // sort final state particles w.r.t. name
  map<string, vector<particlePtr> > fsPartSorted;
  for (unsigned int i = 0; i < _fsParticles.size(); ++i) {
    const particlePtr& part = _fsParticles[i];
    fsPartSorted[part->name()].push_back(part);
  }
  for (map<string, vector<particlePtr> >::iterator i = fsPartSorted.begin();
       i != fsPartSorted.end(); ++i) {
    // sort particles of the same name by their index
    vector<particlePtr>& parts = i->second;
    sort(parts.begin(), parts.end(), rpwa::compareIndicesAsc);
    // check indices
    for (unsigned int j = 0; j < parts.size(); ++j)
      if (parts[j]->index() != (int)j) {
  	printWarn << "indices for final state particles '" << i->first << "' are not consecutive. "
  		  << "expected index [" << j << "], found [" << parts[j]->index() << "]." << endl;
  	topologyIsOkay = false;
      } else if (_debug)
	printInfo << "success: final state particle '" << i->first << "' "
		  << "has expected index [" << j << "]" << endl;
  }
  if (_debug)
    printInfo << "decay topology " << ((topologyIsOkay) ? "passed" : "did not pass")
	      << " all tests" << endl;
  return topologyIsOkay;
}


decayTopology
decayTopology::subDecay(const nodeDesc& startNd)
{
  decayTopology subTopo(dfsSubGraph(startNd));
  subTopo.name() = "subdecay";
  return subTopo;
}


void
decayTopology::addDecay(const decayTopology& topo)
{
  decayTopologyGraphType::addGraph(topo);
  buildInternalData();
}


void decayTopology::setProductionVertex(const interactionVertexPtr& productionVertex)
{
  if (!productionVertex) {
    printErr << "null pointer for production vertex. aborting." << endl;
    throw;
  }
  if (!productionVertex->outParticles()[0]) {
    printErr << "null pointer for particle[0] coming out of production vertex. aborting." << endl;
    throw;
  }
  name() = productionVertex->outParticles()[0]->qnSummary();
  if (!_prodVertex) {
    // topology does not have production vertex -> create graph node
    addVertex(productionVertex);
  } else {
    // delete existsing production vertex -> update graph node
    vertex(node(_prodVertex)) = productionVertex;
    _prodVertex.reset();
  }
  _prodVertex = productionVertex;
}


bool
decayTopology::readData(const TClonesArray& initialStateNames,
			const TClonesArray& initialStateMomenta,
			const TClonesArray& finalStateNames,
			const TClonesArray& finalStateMomenta)
{
  bool success = true;
  // set initial state
  // production vertex class has to implement readData()
  if(!productionVertex()->readData(initialStateNames, initialStateMomenta))
    success = false;
  // check final state data
  const string nameClassName = finalStateNames.GetClass()->GetName();
  if (nameClassName != "TObjString") {
    printWarn << "final state names are not of type TObjString." << endl;
    success = false;
  }
  const string momClassName = finalStateMomenta.GetClass()->GetName();
  if (momClassName != "TVector3") {
    printWarn << "final state momenta are not of type TVector3." << endl;
    success = false;
  }
  if (finalStateNames.GetEntriesFast() != finalStateMomenta.GetEntriesFast()) {
    printWarn << "arrays for final state names and momenta have different sizes: "
	      << finalStateNames.GetEntriesFast() << " vs. "
	      << finalStateMomenta.GetEntriesFast() << endl;
    success = false;
  }
  if (!success)
    return false;
  // sort final state momenta w.r.t. particle name
  const unsigned int nmbFsPart = finalStateNames.GetEntriesFast();
  map<string, vector<const TVector3*> > fsMomenta;
  for (unsigned int i = 0; i < nmbFsPart; ++i) {
    const string    name = ((TObjString*)finalStateNames[i])->GetString().Data();
    const TVector3* mom  = (TVector3*)finalStateMomenta[i];
    fsMomenta[name].push_back(mom);
  }
  // set final state
  for (unsigned int i = 0; i < nmbFsParticles(); ++i) {
    const particlePtr& part     = fsParticles()[i];
    const string       partName = part->name();
    map<string, vector<const TVector3*> >::const_iterator entry = fsMomenta.find(partName);
    if (entry != fsMomenta.end()) {
      unsigned int partIndex = part->index();
      if (partIndex < 0)
      	partIndex = 0;
      if (partIndex < entry->second.size()) {
      	if (_debug)
      	  printInfo << "setting momentum of final state particle " << partName << "["
		    << partIndex << "] to " << *(entry->second[partIndex]) << " GeV" << endl;
      	part->setMomentum(*(entry->second[partIndex]));
      } else {
      	printWarn << "index [" << partIndex << "] for final state particle "
      		  << "'" << part->name() << "' out of range. data contain only "
      		  << entry->second.size() << " entries for this particle type." << endl;
      	success = false;
      }
    } else {
      printWarn << "cannot find entry for final state particle '" << part->name() << "' "
		<< "in data." << endl;
      success = false;
    }
  }
  return success;
}


ostream&
decayTopology::print(ostream& out) const
{
  // print nodes
  out << "decay topology '" << name() << "' has " << nmbNodes() << " node(s):" << endl;
  if (_prodVertex)
    out << "    production  node[" << node(_prodVertex) << "] = " << *_prodVertex << endl;
  else
    out << "    topology has no production node." << endl;
  for (unsigned int i = 0; i < nmbInteractionVertices(); ++i)
    out << "    interaction node[" << node(_intVertices[i]) << "] = " << *_intVertices[i] << endl;
  nodeIterator iNd, iNdEnd;
  for (tie(iNd, iNdEnd) = nodes(); iNd != iNdEnd; ++iNd)
    if (isFsVertex(vertex(*iNd)))
      out << "    final state node[" << *iNd << "] = " << *vertex(*iNd) << endl;
  // print edges
  out << "decay topology '" << name() << "' has " << nmbEdges() << " edge(s):" << endl;
  edgeIterator iEd, iEdEnd;
  for (tie(iEd, iEdEnd) = edges(); iEd != iEdEnd; ++iEd) {
    const particlePtr part = particle(*iEd);
    out << "    edge[" << *iEd << "] = [" << fromNode(*iEd) << ", "
	<< toNode(*iEd) << "] = '" << part->name() << "'";
    if (isFsParticle(part))
      out << " (final state)";
    out << endl;
  }
  return out;
}


ostream&
decayTopology::printIsParticles(ostream& out) const
{
  out << "initial state particles:" << endl;
  for (unsigned int i = 0; i < productionVertex()->nmbInParticles(); ++i) {
    const particlePtr& part = productionVertex()->inParticles()[i];
    out << "    particle[" << i << "]: " << part->qnSummary() << ", "
	<< "index = " << part->index() << ", p = " << part->lzVec() << " GeV" << endl;
  }
  return out;
}


ostream&
decayTopology::printFsParticles(ostream& out) const
{
  out << "final state particles:" << endl;
  for (unsigned int i = 0; i < nmbFsParticles(); ++i) {
    const particlePtr& part = fsParticles()[i];
    out << "    particle[" << i << "]: " << part->qnSummary() << ", "
	<< "index = " << part->index() << ", p = " << part->lzVec() << " GeV" << endl;
  }
  return out;
}


decayTopology&
decayTopology::constructDecay(const interactionVertexPtr&              productionVertex,
			      const std::vector<interactionVertexPtr>& interactionVertices,
			      const std::vector<particlePtr>&          fsParticles)
{
  clear();
  const unsigned int nmbVert   = interactionVertices.size();
  const unsigned int nmbFsPart = fsParticles.size();
  if (_debug)
    printInfo << "constructing decay topology with "
	      << nmbVert   << " interaction vertices and "
	      << nmbFsPart << " final state particles" << endl;
  if (nmbFsPart < 1) {
    printErr << "cannot construct decay topology without final state particles. aborting." << endl;
    throw;
  }
  if (nmbVert < 1) {
    printWarn << "need at least production and X-decay vertex "
	      << "to construct decay topology. aborting." << endl;
    throw;
  }

  // create graph node for production vertex and store pointer
  if (!productionVertex) {
    printErr << "null pointer for production vertex. aborting." << endl;
    throw;
  }
  if (!productionVertex->outParticles()[0]) {
    printErr << "null pointer for particle[0] coming out of production vertex. aborting." << endl;
    throw;
  }
  name() = productionVertex->outParticles()[0]->qnSummary();
  _prodVertex = productionVertex;
  addVertex(productionVertex);
  // create graph nodes for interaction vertices and store pointers
  for (unsigned int i = 0; i < nmbVert; ++i) {
    if (!interactionVertices[i]) {
      printErr << "null pointer for decay vertex[" << i << "]. aborting." << endl;
      throw;
    }
    addVertex(interactionVertices[i]);
  }
  // create final state nodes and vertices
  for (unsigned int i = 0; i < nmbFsPart; ++i) {
    if (!fsParticles[i]) {
      printErr << "null pointer for final state particle[" << i << "]. aborting." << endl;
      throw;
    }
    fsVertexPtr fsVert = createFsVertex(fsParticles[i]);
    addVertex(fsVert);
  }
  // memorize depth-first sorted interaction vertices and final state particles
  vector<nodeDesc> sortedNds = sortNodesDfs(node(_prodVertex));
  for (unsigned int i = 0; i < sortedNds.size(); ++i) {
    const interactionVertexPtr& vert   = vertex(sortedNds[i]);
    const fsVertexPtr           fsVert = dynamic_pointer_cast<fsVertex>(vert);
    if (fsVert) {
      const particlePtr& part = fsVert->inParticles()[0];
      for (unsigned int i = 0; i < nmbFsPart; ++i)
	if (part == fsParticles[i]) {
	  _fsParticles.push_back(part);
	  break;
	}
    } else if (vert != _prodVertex) {
      for (unsigned int i = 0; i < nmbVert; ++i)
	if (vert == interactionVertices[i]) {
	  _intVertices.push_back(vert);
	  break;
	}
    }
  }
  // check that topology makes sense
  if (!checkTopology()) {
    printErr << "topology has problems that need to be fixed. aborting." << endl;
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
    if ((nmbInEdges(*iNd) == 0) && (nmbOutEdges(*iNd) == 1)) {
      _prodVertex = vertex(*iNd);
      ++nmbProdVertCandidates;
    }
  // set final state particles
  _fsParticles.clear();
  vector<nodeDesc> sortedNds;
  if (_prodVertex)
    sortedNds = sortNodesDfs(node(_prodVertex));
  else {
    // take non-FS vertex with no incoming particle as top node
    nodeDesc     topNode;
    unsigned int nmbTopVertCandidates = 0;
    for (tie(iNd, iNdEnd) = nodes(); iNd != iNdEnd; ++iNd)
      if ((nmbInEdges(*iNd) == 0) && !isFsVertex(vertex(*iNd))) {
	topNode = *iNd;
	++nmbTopVertCandidates;
      }
    if (nmbTopVertCandidates == 1)
      sortedNds = sortNodesDfs(topNode);
    else {
      if (_debug)
	printWarn << "ill-formed graph. vertex order will be un-defined." << endl;
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
    printWarn << "cannot find final state particles in graph '" << name() << "'" << endl;
    success = false;
  }
  if (!success) {
    printErr << "cannot construct decay topology from graph '" << name() << "'. aborting." << endl;
    throw;
  }
  // find interaction vertices
  _intVertices.clear();
  for (unsigned int i = 0; i < sortedNds.size(); ++i) {
    const interactionVertexPtr& vert = vertex(sortedNds[i]);
    if (isInteractionVertex(vert))
      _intVertices.push_back(vert);
  }
}


interactionVertexPtr
decayTopology::cloneNode(const nodeDesc& nd)
{
  const interactionVertexPtr v    = vertex(nd);
  interactionVertexPtr       newV = decayTopologyGraphType::cloneNode(nd);
  // update member variables
  if (isProductionVertex(v))
    _prodVertex = newV;
  else if (isInteractionVertex(v))
    for (unsigned int i = 0; i < _intVertices.size(); ++i)
      if (v == _intVertices[i])
	_intVertices[i] = newV;
  return newV;
}


particlePtr
decayTopology::cloneEdge(const edgeDesc& ed)
{
  const particlePtr p    = particle(ed);
  particlePtr       newP = decayTopologyGraphType::cloneEdge(ed);
  // update member variable
  if (isFsParticle(p))
    for (unsigned int i = 0; i < _fsParticles.size(); ++i)
      if (p == _fsParticles[i])
  	_fsParticles[i] = newP;
  edgeIterator iEd, iEdEnd;
  return newP;
}
