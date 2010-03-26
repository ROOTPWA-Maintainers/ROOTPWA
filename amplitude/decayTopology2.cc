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
//      vertex and at least one interaction veretx; in addtion for
//      each final state particle a corresponding final state vertex
//      will be created for internal use
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
    _prodVertex   ()
{
  _prodVertex.reset();
}


// decayTopology2::decayTopology2(const decayTopology2& topo)
// {
//   *this = topo;
// }


decayTopology2::decayTopology2(const VPtr&                         productionVertex,
			       const std::vector<VPtr>&            interactionVertices,
			       const std::vector<rpwa::particle*>& fsParticles)
{
  constructDecay(productionVertex, interactionVertices, fsParticles);
}


decayTopology2::~decayTopology2()
{ }


// decayTopology2&
// decayTopology2::operator = (const decayTopology2& topo)
// {
//   if (this != &topo) {
//     _graph            = topo._graph;
//     _nodeDataMap      = topo._nodeDataMap;
//     _nodeVertexMap    = topo._nodeVertexMap;
//     _nodeIndexMap     = topo._nodeIndexMap;
//     _edgeDataMap      = topo._edgeDataMap;
//     _edgeParticleMap  = topo._edgeParticleMap;
//     _edgeIndexMap     = topo._edgeIndexMap;
//     _vertexNodeMap    = topo._vertexNodeMap;
//     _particleEdgeMap  = topo._particleEdgeMap;
//     _prodVertex = topo._prodVertex;
//     _vertices         = topo._vertices;
//     _fsParticles      = topo._fsParticles;
//   }
//   return *this;
// }


// ostream&
// decayTopology2::print(ostream& out) const
// {
//   out << "decay topology nodes:" << endl;
//   cout << "    production node[" << _vertexNodeMap.find(_prodVertex)->second << "] = ";
//   out << _prodVertex->inParticles() [0]->summary() << "  --->  "
//       << _prodVertex->outParticles()[0]->summary() << endl;
//   for (unsigned int i = 0; i < nmbVertices(); ++i) {
//     cout << "    decay node[" << _vertexNodeMap.find(_vertices[i])->second << "] = ";
//     for (unsigned int j = 0; j < _vertices[i]->nmbInParticles(); ++j) {
//       out << _vertices[i]->inParticles()[j]->summary();
//       if (j < _vertices[i]->nmbInParticles() - 1)
//   	out << "  +  ";
//     }
//     out << "  --->  ";
//     for (unsigned int j = 0; j < _vertices[i]->nmbOutParticles(); ++j) {
//       out << _vertices[i]->outParticles()[j]->summary();
//       if (j < _vertices[i]->nmbOutParticles() - 1)
//   	out << "  +  ";
//     }
//     out << endl;
//   }
//   out << "decay topology particles:" << endl;
//   //cout << "### node[7] " << ((_nodeVertexMap[7]) ? "OOOPS!" : "okay") << endl;
//   edgeIterator iEdge, iEdgeEnd;
//   for (tie(iEdge, iEdgeEnd) = edges(_graph); iEdge != iEdgeEnd; ++iEdge) {
//     const particle* p = _edgeParticleMap[*iEdge];
//     out << "    particle" << *iEdge << " ";
//     if (p) {
//       out << p->name() << " " << target(*iEdge, _graph) << " ";
//       if (!_nodeVertexMap[target(*iEdge, _graph)])
// 	out << "(final state)";
//       out << endl;
//     } else
//       out << "error! null pointer to particle for edge" << *iEdge << endl;
//   }
//   return out;
// }


// ostream&
// decayTopology2::writeGraphViz(ostream& out) const
// {
//   write_graphviz(out, _graph);
//   return out;
// }


void
decayTopology2::clear()
{
  decayGraphType::clear();  //!!! does not exist for subgraph
  _prodVertex.reset();
  _intVertices.clear();
  _fsParticleVertexMap.clear();
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


// interactionVertex*
// decayTopology2::vertex(particle& part)
// {
//   edgeDesc edge = _particleEdgeMap[&part];
//   nodeDesc node = target(edge, _graph);
//   return _nodeVertexMap[node];
// }


decayTopology2&
decayTopology2::constructDecay(const VPtr&                         productionVertex,
			       const std::vector<VPtr>&            interactionVertices,
			       const std::vector<rpwa::particle*>& fsParticles)
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
    fsVertexPtr fsVert = createFsVertex(*fsParticles[i]);
    nodeDesc    fsNd   = addVertex(fsVert);
    if (_debug)
      printInfo << "added " << *vertex(fsNd) << " as tree node [" << fsNd << "]" << endl;
  }
  // memorize depth-first sorted nodes
  vector<nodeDesc> sortedNds = sortNodesDfs(node(_prodVertex));
  for (unsigned int i = 0; i < sortedNds.size(); ++i) {
    const VPtr& v = vertex(sortedNds[i]);
    for (unsigned int i = 0; i < nmbVert; ++i)
      if (v == interactionVertices[i]) {
	_intVertices.push_back(v);
	break;
      }
    const fsVertexPtr fsV = dynamic_pointer_cast<fsVertex>(v);
    if (fsV)
      for (unsigned int i = 0; i < nmbFsPart; ++i) {
	rpwa::particle* p = fsV->inParticles()[0];
	if (p == fsParticles[i]) {
	  _fsParticleVertexMap[p] = fsV;
	  break;
	}
      }			   
  }
  // set graph name to particle coming out of production vertex
  if (nmbOutParticles(_prodVertex) != 1) {
    printErr << "production vertex has " << nmbOutParticles(_prodVertex) << " != 1 "
	     << "outgoing particles. aborting." << endl;
    throw;
  }
  name() = _prodVertex->outParticles()[0]->name();
  return *this;
}
