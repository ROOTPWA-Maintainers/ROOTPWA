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


// #include <list>

#include "utilities.h"
#include "decayTopology2.h"

	
using namespace std;
using namespace boost;
using namespace rpwa;


bool decayTopology2::_debug = false;


decayTopology2::decayTopology2()
  : decayGraphType(),
    _prodVertex   (),
    _intVertices  (),
    _fsParticles  ()
{
  _prodVertex.reset();
}


decayTopology2::decayTopology2(const decayTopology2& topo)
{
  *this = topo;
}


decayTopology2::decayTopology2(const interactionVertexPtr&              productionVertex,
			       const std::vector<interactionVertexPtr>& interactionVertices,
			       const std::vector<particlePtr>&          fsParticles)
{
  constructDecay(productionVertex, interactionVertices, fsParticles);
}


decayTopology2::decayTopology2(const decayGraphType& graph)
{
  *this = graph;
}


decayTopology2::~decayTopology2()
{ }


decayTopology2&
decayTopology2::operator =(const decayTopology2& topo)
{
  if (this != &topo) {
    decayGraphType::operator =(topo);
    _prodVertex  = topo._prodVertex;
    _intVertices = topo._intVertices;
    _fsParticles = topo._fsParticles;
  }
  return *this;
}


decayTopology2&
decayTopology2::operator =(const decayGraphType& graph)
{
  if (this != &graph) {
    decayGraphType::operator =(graph);
    bool success = true;
    // find production vertex and final state particles
    vector<nodeDesc> sortedNds = sortNodesDfs(node(_prodVertex));
    for (unsigned int i = 0; i < sortedNds.size(); ++i) {
      if (nmbInEdges(sortedNds[i]) == 0)
	_prodVertex = vertex(sortedNds[i]);
      if (nmbOutEdges(sortedNds[i]) == 0) {
	const interactionVertexPtr& vert = vertex(sortedNds[i]);
	if (dynamic_pointer_cast<fsVertex>(vert))
	  _fsParticles.push_back(vert->inParticles()[0]);
      }
    }
    if (!_prodVertex) {
      printWarn << "cannot find production vertex in graph '" << graph.name() << "'" << endl;
      success = false;
    }
    if (_fsParticles.size() < 1) {
      printWarn << "cannot find final state particles in graph '" << graph.name() << "'" << endl;
      success = false;
    }
    if (!success) {
      printErr << "cannot construct decay topology from graph '" << graph.name() << "'. aborting." << endl;
      throw;
    }
    // find interaction vertices
    for (unsigned int i = 0; i < sortedNds.size(); ++i) {
      const interactionVertexPtr& vert = vertex(sortedNds[i]);
      if (isInteractionVertex(vert))
	_intVertices.push_back(vert);
    }
  }
  return *this;
}


void
decayTopology2::clear()
{
  decayGraphType::clear();
  _prodVertex.reset();
  _intVertices.clear();
  _fsParticles.clear();
}


// decayTopology2::decayGraph
// decayTopology2::deepCopyGraph(const decayGraph& srcGraph,
// 			     const bool        copyFsParticles)
// {
//   // copy graph data structure
//   decayGraph newGraph = srcGraph;
//   // copy vertices
//   nodeVertexType nodeVertexMap = get(vertex_vertexPointer, newGraph);
//   for (nodeIterator iNode = vertices(newGraph).first; iNode != vertices(newGraph).second; ++iNode) {
//     interactionVertex* srcVertex = nodeVertexMap[*iNode];
//     if (srcVertex) {
//       if (_debug)
//   	printInfo << "copying vertex " << *iNode << endl;
//       nodeVertexMap[*iNode] = &srcVertex->clone();
//     }
//   }  
//   // copy particles
//   //!!! what about beam particle?
//   edgeParticleType edgeParticleMap = get(edge_particlePointer, newGraph);
//   for (edgeIterator iEdge = edges(newGraph).first; iEdge != edges(newGraph).second; ++iEdge) {
//     particle* srcParticle = edgeParticleMap[*iEdge];
//     particle* newParticle = srcParticle;
//     if (!srcParticle) {
//       printErr << "null pointer to particle for edge " << *iEdge << ". aborting." << endl;
//       throw;
//     }
//     if (_debug) {
//       if (nodeVertexMap[target(*iEdge, newGraph)])
//   	printInfo << "copying particle ";
//       else if (copyFsParticles)
//   	printInfo << "copying final state particle ";
//     }
//     if (nodeVertexMap[target(*iEdge, newGraph)] || copyFsParticles) {
//       if (_debug)
//   	cout << *iEdge << " '" << srcParticle->name() << "'" << endl;
//       newParticle             = &srcParticle->clone();
//       edgeParticleMap[*iEdge] = newParticle;
//     }
//     // modify the vertices that the particle connects
//     nodeDesc           srcNode   = source(*iEdge, newGraph);
//     interactionVertex* srcVertex = nodeVertexMap[srcNode];
//     if (!srcVertex) {
//       printErr << "null pointer to source vertex. aborting." << endl;
//       throw;
//     }
//     for (unsigned int i = 0; i < srcVertex->nmbOutParticles(); ++i)
//       if (srcVertex->outParticles()[i] == srcParticle) {
// 	srcVertex->outParticles()[i] = newParticle;
// 	break;
//       }
//     nodeDesc           targetNode   = target(*iEdge, newGraph);
//     interactionVertex* targetVertex = nodeVertexMap[targetNode];
//     if (targetVertex) {
//       for (unsigned int i = 0; i < targetVertex->nmbInParticles(); ++i)
// 	if (targetVertex->inParticles()[i] == srcParticle) {
// 	  targetVertex->inParticles()[i] = newParticle;
// 	  break;
// 	}
//     }
//   }
//   return newGraph;
// }


bool
decayTopology2::isInteractionVertex(const interactionVertexPtr& vert) const
{
  if (isProductionVertex(vert) || isFsVertex(vert))
    return false;
  return true;
}


bool
decayTopology2::isFsVertex(const interactionVertexPtr& vert) const
{
  for (unsigned int i = 0; i < _fsParticles.size(); ++i)
    if (_fsParticles[i] == vert->inParticles()[0])
      return true;
  return false;
}


bool
decayTopology2::isFsParticle(const particlePtr& part) const
{
  for (unsigned int i = 0; i < _fsParticles.size(); ++i)
    if (part == _fsParticles[i])
      return true;
  return false;
}


ostream&
decayTopology2::print(ostream& out) const
{
  out << "decay topology nodes:" << endl;
  out << "    production node[" << node(_prodVertex) << "] = " << *_prodVertex << endl;
  for (unsigned int i = 0; i < nmbInteractionVertices(); ++i)
    out << "    interaction node[" << node(_intVertices[i]) << "] = " << *_intVertices[i] << endl;
  out << "decay topology particles:" << endl;
  edgeIterator iEdge, iEdgeEnd;
  for (tie(iEdge, iEdgeEnd) = edges(); iEdge != iEdgeEnd; ++iEdge) {
    const particlePtr part = particle(*iEdge);
    out << "    edge[" << *iEdge << "] = [" << fromNode(*iEdge) << ", " << toNode(*iEdge) << "] = '" << part->name() << "'";
    if (isFsParticle(part))
      out << " (final state)";
    out << endl;
  }
  return out;
}


bool
decayTopology2::checkTopology() const
{
  bool topologyIsOkay = true;
  // make sure there are no dangling particles
  nodeIterator iNode, iNodeEnd;
  for (tie(iNode, iNodeEnd) = nodes(); iNode != iNodeEnd; ++iNode) {
    const interactionVertexPtr& vert = vertex(*iNode);
    if (!vert) {
      printWarn << "node[" << *iNode << "] has vertex null pointer" << endl;
      topologyIsOkay = false;
      continue;
    }
    if (vert != _prodVertex)  // incoming particles of production vertex have no edges
      for (unsigned int i = 0; i < vert->nmbInParticles(); ++i)
	if (!isEdge(vert->inParticles()[i])) {
	  printWarn << "incoming particle[" << i << "] of vertex " << *vert << " has no associated graph edge" << endl;
	  topologyIsOkay = false;
	} else if (_debug)
	  printInfo << "success: incoming particle[" << i << "] of vertex " << *vert << " has associated graph edge" << endl;
    for (unsigned int i = 0; i < vert->nmbOutParticles(); ++i)
      if (!isEdge(vert->outParticles()[i])) {
	printWarn << "outgoing particle[" << i << "] of vertex " << *vert << " has no associated graph edge" << endl;
	topologyIsOkay = false;
      } else if (_debug)
	printInfo << "success: outgoing particle[" << i << "] of vertex " << *vert << " has associated graph edge" << endl;
  }
  // check production vertex
  if (!_prodVertex) {
    printWarn << "production vertex is null pointer" << endl;
    topologyIsOkay = false;
  }
  if (!isNode(_prodVertex) || (vertex(node(_prodVertex)) != _prodVertex)) {
    printWarn << *_prodVertex << " has no associated graph node" << endl;
    topologyIsOkay = false;
  } else if (_debug)
    printInfo << "success: " << *_prodVertex << " has node in graph" << endl;
  if (nmbOutEdges(_prodVertex) != 1) {
    printWarn << *_prodVertex << " has " << nmbOutEdges(_prodVertex) << " != 1 "
	      << "outgoing edges" << endl;
    topologyIsOkay = false;
  } else if (nmbInEdges(_prodVertex) != 0) {
    printWarn << *_prodVertex << " has " << nmbInEdges(_prodVertex) << " != 0 "
	      << "incoming edges" << endl;
    topologyIsOkay = false;
  } else if (_debug)
    printInfo << "success: " << *_prodVertex << " has exactly 1 outgoing and no incoming edges" << endl;
  // make sure that all nodes are reachable from production node
  //!!! this might not work; check!
  const vector<nodeDesc> sortedNds   = sortNodesDfs(node(_prodVertex));
  const unsigned int     nmbVertices = 1 + _intVertices.size() + _fsParticles.size();
  if (sortedNds.size() != nmbVertices) {
    printWarn << "number of nodes reachable from production node is " << sortedNds.size() << " != "
	      << nmbVertices << endl;
    topologyIsOkay = false;
  } else if (_debug)
    printInfo << "success: all nodes are reachable from poduction node" << endl;
  // make sure all interaction vertices have corresponding nodes
  for (unsigned int i = 0; i < _intVertices.size(); ++i) {
    const interactionVertexPtr& vert = _intVertices[i];
    if (!vert) {
      printWarn << "interaction vertex[" << i << "] is null pointer" << endl;
      topologyIsOkay = false;
      continue;
    }
    if (!isNode(vert) || (vertex(node(vert)) != vert)) {
      printWarn << *vert << " has no associated graph node" << endl;
      topologyIsOkay = false;
    } else if (_debug)
      printInfo << "success: " << *vert << " has node in graph" << endl;
    if (nmbOutEdges(vert) < 1) {
      printWarn << *vert << " has " << nmbOutEdges(vert) << " < 1 "
		<< "outgoing edges" << endl;
      topologyIsOkay = false;
    } else if (nmbInEdges(vert) < 1) {
      printWarn << *vert << " has " << nmbInEdges(vert) << " < 1 "
		<< "incoming edges" << endl;
      topologyIsOkay = false;
    } else if (_debug)
      printInfo << "success: " << *vert << " has at least 1 outgoing and 1 incoming edge" << endl;
  }
  // make sure all final state particles have corresponding edges
  for (unsigned int i = 0; i < _fsParticles.size(); ++i) {
    const particlePtr& part = _fsParticles[i];
    if (!part) {
      printWarn << "final state particle[" << i << "] is null pointer" << endl;
      topologyIsOkay = false;
      continue;
    }
    if (!isEdge(part) || (particle(edge(part)) != part)) {
      printWarn << *part << " has no associated graph edge" << endl;
      topologyIsOkay = false;
    } else if (_debug)
      printInfo << "success: " << *part << " has edge in graph" << endl;
    const fsVertexPtr fsVert = dynamic_pointer_cast<fsVertex>(toVertex(part));
    if (!fsVert) {
      printWarn << "vertex associated to " << *part << " is not a fsVertex" << endl;
      topologyIsOkay = false;
    } else if (_debug)
      printInfo << "success: vertex associated to " << *part << " is an fsVertex" << endl;
  }
  return topologyIsOkay;
}


// decayTopology2::decayGraph
// decayTopology2::subGraph(interactionVertex& startVertex)
// {
//   // find all nodes below start node
//   nodeDesc         startNode = _vertexNodeMap[&startVertex];
//   vector<nodeDesc> subGraphNodes;
//   {
//     vector<default_color_type> colors(num_vertices(_graph));
//     depth_first_visit(_graph, startNode,
// 		      makeDfsRecorder(back_inserter(subGraphNodes)),
// 		      make_iterator_property_map(colors.begin(),
// 						 get(vertex_index, _graph), colors[0]));
//   }
//   if (_debug)
//     printInfo << "creating subgraph starting at node " << startNode
// 	      << " using nodes: " << flush;
//   decayGraph subGraph = _graph.create_subgraph();
//   for (unsigned int j = 0; j < subGraphNodes.size(); ++j) {
//     add_vertex(subGraphNodes[j], subGraph);
//     if (_debug)
//       cout << subGraphNodes[j] << "    ";
//   }
//   if (_debug)
//     cout << endl;
//   return subGraph;
// }


decayTopology2&
decayTopology2::constructDecay(const interactionVertexPtr&              productionVertex,
			       const std::vector<interactionVertexPtr>& interactionVertices,
			       const std::vector<particlePtr>& fsParticles)
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
  _prodVertex     = productionVertex;
  nodeDesc prodNd = addVertex(productionVertex);
  if (_debug)
    printInfo << "added " << *vertex(prodNd) << " as tree node [" << prodNd << "]" << endl;
  // create graph nodes for interaction vertices and store pointers
  for (unsigned int i = 0; i < nmbVert; ++i) {
    if (!interactionVertices[i]) {
      printErr << "null pointer for decay vertex[" << i << "]. aborting." << endl;
      throw;
    }
    nodeDesc intNd = addVertex(interactionVertices[i]);
    if (_debug)
      printInfo << "added " << *vertex(intNd) << " as tree node [" << intNd << "]" << endl;
  }
  // create final state nodes and vertices
  for (unsigned int i = 0; i < nmbFsPart; ++i) {
    if (!fsParticles[i]) {
      printErr << "null pointer for final state particle[" << i << "]. aborting." << endl;
      throw;
    }
    fsVertexPtr fsVert = createFsVertex(fsParticles[i]);
    nodeDesc    fsNd   = addVertex(fsVert);
    if (_debug)
      printInfo << "added " << *vertex(fsNd) << " as tree node [" << fsNd << "]" << endl;
  }
  // check that topolgy makes sense
  if (!checkTopology()) {
    printErr << "topology has problems that need to be fixed. aborting." << endl;
    throw;
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
  // set graph name to particle coming out of production vertex
  name() = _prodVertex->outParticles()[0]->name();
  return *this;
}
