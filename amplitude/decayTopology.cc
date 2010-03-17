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
//      internally the decay process is represented as a graph using
//      the Boost Graph Library
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <list>

#include <boost/graph/depth_first_search.hpp>
#if BOOST_GRAPH_LIB
#include <boost/graph/graphviz.hpp>
#endif

#include "utilities.h"
#include "decayTopology.h"

	
using namespace std;
using namespace boost;
using namespace rpwa;


bool decayTopology::_debug = false;


decayTopology::decayTopology()
  : _productionVertex(0)
{ }


decayTopology::decayTopology(const decayTopology& topo)
{
  *this = topo;
}


decayTopology::decayTopology(const vector<particle*>&          fsParticles,
			     interactionVertex&                productionVertex,
			     const vector<interactionVertex*>& vertices)
{
  constructDecay(fsParticles, productionVertex, vertices);
}


decayTopology::~decayTopology()
{ }


decayTopology&
decayTopology::operator = (const decayTopology& topo)
{
  if (this != &topo) {
    _graph            = topo._graph;
    _nodeDataMap      = topo._nodeDataMap;
    _nodeVertexMap    = topo._nodeVertexMap;
    _nodeIndexMap     = topo._nodeIndexMap;
    _edgeDataMap      = topo._edgeDataMap;
    _edgeParticleMap  = topo._edgeParticleMap;
    _edgeIndexMap     = topo._edgeIndexMap;
    _vertexNodeMap    = topo._vertexNodeMap;
    _particleEdgeMap  = topo._particleEdgeMap;
    _productionVertex = topo._productionVertex;
    _vertices         = topo._vertices;
    _fsParticles      = topo._fsParticles;
  }
  return *this;
}


decayTopology&
decayTopology::constructDecay(const vector<particle*>&          fsParticles,
			      interactionVertex&                productionVertex,
			      const vector<interactionVertex*>& interactionVertices)
{
  const unsigned int nmbVert   = interactionVertices.size();
  const unsigned int nmbFsPart = fsParticles.size();
  if (_debug)
    printInfo << "constructing decay topology with "
	      << nmbFsPart << " final state particles and "
	      << nmbVert   << " interaction vertices" << endl;
  clear();
  _productionVertex = &productionVertex;
  if (nmbFsPart < 1) {
    printWarn << "cannot construct decay topology without final state particles" << endl;
    clear();
    return *this;
  }
  if (nmbVert < 1) {
    printWarn << "need at least production and X-decay vertex to construct decay topology" << endl;
    clear();
    return *this;
  }
  get_property(_graph, graph_name) = "root graph";
  // get node and edge maps
  _nodeDataMap     = get(vertex_bundle,        _graph);
  _nodeVertexMap   = get(vertex_vertexPointer, _graph);
  _nodeIndexMap    = get(vertex_index,         _graph);
  _edgeDataMap     = get(edge_bundle,          _graph);
  _edgeParticleMap = get(edge_particlePointer, _graph);
  _edgeIndexMap    = get(edge_index,           _graph);
  // create graph node for production vertex and store pointer
  {
    nodeDesc node = add_vertex(_graph);
    _nodeVertexMap[node]              = &productionVertex;
    _vertexNodeMap[&productionVertex] = node;
    if (_debug)
      printInfo << "added " << productionVertex << " as tree node [" << node << "]" << endl;
  }
  // create graph nodes for interaction vertices and store pointers
  for (unsigned int i = 0; i < nmbVert; ++i) {
    nodeDesc node = add_vertex(_graph);
    if (!interactionVertices[i]) {
      printErr << "null pointer for decay vertex[" << i << "]. aborting." << endl;
      throw;
    }
    _nodeVertexMap[node]                   = interactionVertices[i];
    _vertexNodeMap[interactionVertices[i]] = node;
    if (_debug)
      printInfo << "added " << *interactionVertices[i] << " as tree node [" << node << "]" << endl;
  }
  // create final state nodes
  for (unsigned int i = 0; i < nmbFsPart; ++i) {
    nodeDesc node        = add_vertex(_graph);
    _nodeVertexMap[node] = 0;
    if (_debug)
      printInfo << "created tree node [" << node << "] for final state "
		<< *fsParticles[i] << endl;
  }
  // create edges that connect the intercation nodes and store pointers to particles
  for (vertexNodeMapIt iFromVert = _vertexNodeMap.begin();
       iFromVert != _vertexNodeMap.end(); ++iFromVert)
    for (vertexNodeMapIt iToVert = _vertexNodeMap.begin();
	 iToVert != _vertexNodeMap.end(); ++iToVert) {
      interactionVertex* fromVertex = iFromVert->first;
      for (unsigned int iOutPart = 0; iOutPart < fromVertex->nmbOutParticles(); ++iOutPart) {
	interactionVertex* toVertex = iToVert->first;
	for (unsigned int iInPart = 0; iInPart < toVertex->nmbInParticles(); ++iInPart) {
	  particle* p = fromVertex->outParticles()[iOutPart];
	  if (p == toVertex->inParticles()[iInPart]) {
	    bool     inserted;
	    edgeDesc edge;
	    tie(edge, inserted) = add_edge(iFromVert->second, iToVert->second, _graph);
	    if (inserted) {
	      if (!p) {
		printErr << "null pointer to isobar particle. aborting." << endl;
		throw;
	      }
	      _edgeParticleMap[edge] = p;
	      _particleEdgeMap[p]    = edge;
	      if (_debug)
		printInfo << "adding edge from node [" << iFromVert->second << "] "
			  << "to node [" << iToVert->second << "]" << endl;
	    }
	  }
	}
      }
    }
  // create edges for final state particles and store pointers
  for (vertexNodeMapIt iFromVert = _vertexNodeMap.begin();
       iFromVert != _vertexNodeMap.end(); ++iFromVert) {
    interactionVertex* fromVertex = iFromVert->first;
    for (unsigned int iOutPart = 0; iOutPart < fromVertex->nmbOutParticles(); ++iOutPart)
      for (unsigned int iFsPart = 0; iFsPart < nmbFsPart; ++iFsPart) {
	particle* p = fromVertex->outParticles()[iOutPart];
	if (p == fsParticles[iFsPart]) {
	  bool     inserted;
	  edgeDesc edge;
	  tie(edge, inserted) = add_edge(iFromVert->second, nmbVert + 1 + iFsPart, _graph);
	  if (inserted) {
	    if (!p) {
	      printErr << "null pointer for final state particle[" << iFsPart << "]. "
		       << "aborting." << endl;
	      throw;
	    }
	    _edgeParticleMap[edge] = p;
	    _particleEdgeMap[p]    = edge;
	    if (_debug)
	      printInfo << "adding edge from node [" << iFromVert->second << "] "
			<< "to final state node [" << nmbVert + 1 + iFsPart << "]" << endl;
	  }
	}
      }
  }
  // sorting nodes depth-first
  vector<nodeDesc> nodes;
  vector<default_color_type> colors(num_vertices(_graph));
  depth_first_visit(_graph, _vertexNodeMap[&productionVertex],
		    makeDfsRecorder(back_inserter(nodes)),
		    make_iterator_property_map(colors.begin(),
						 get(vertex_index, _graph), colors[0]));
  if (_debug) {
    printInfo << "depth-first node order: ";
    for (unsigned int i = 0; i < nodes.size(); ++i)
      cout << nodes[i] << ((i < nodes.size() - 1) ? ", " : "");
    cout << endl;
  }
  for (unsigned int i = 0; i < nodes.size(); ++i) {
    interactionVertex* vertex = _nodeVertexMap[nodes[i]];
    if (vertex && vertex != _productionVertex)  // skip production vertex
      _vertices.push_back(vertex);
  }
  // copy final state particles
  _fsParticles = fsParticles;
  return *this;
}


ostream&
decayTopology::print(ostream& out) const
{
  out << "decay topology nodes:" << endl;
  cout << "    production node[" << _vertexNodeMap.find(_productionVertex)->second << "] = ";
  out << _productionVertex->inParticles() [0]->name() << "  --->  "
      << _productionVertex->outParticles()[0]->name() << endl;
  for (unsigned int i = 0; i < nmbVertices(); ++i) {
    cout << "    decay node[" << _vertexNodeMap.find(_vertices[i])->second << "] = ";
    for (unsigned int j = 0; j < _vertices[i]->nmbInParticles(); ++j) {
      out << _vertices[i]->inParticles()[j]->name();
      if (j < _vertices[i]->nmbInParticles() - 1)
  	out << "  +  ";
    }
    out << "  --->  ";
    for (unsigned int j = 0; j < _vertices[i]->nmbOutParticles(); ++j) {
      out << _vertices[i]->outParticles()[j]->name();
      if (j < _vertices[i]->nmbOutParticles() - 1)
  	out << "  +  ";
    }
    out << endl;
  }
  out << "decay topology particles:" << endl;
  //cout << "### node[7] " << ((_nodeVertexMap[7]) ? "OOOPS!" : "okay") << endl;
  edgeIterator iEdge, iEdgeEnd;
  for (tie(iEdge, iEdgeEnd) = edges(_graph); iEdge != iEdgeEnd; ++iEdge) {
    const particle* p = _edgeParticleMap[*iEdge];
    out << "    particle" << *iEdge << " ";
    if (p) {
      out << p->name() << " " << target(*iEdge, _graph) << " ";
      if (!_nodeVertexMap[target(*iEdge, _graph)])
	out << "(final state)";
      out << endl;
    } else
      out << "error! null pointer to particle for edge" << *iEdge << endl;
  }
  return out;
}


ostream&
decayTopology::writeGraphViz(ostream& out) const
{
#if BOOST_GRAPH_LIB
  write_graphviz(out, _graph);
#else
  printWarn << "cannot write GraphViz file. libboost_graph.so needed." << endl;
#endif
  return out;
}


void
decayTopology::clear()
{
  //_graph.clear();  //!!! does not exist for subgraph
  decayGraph empty;
  _graph           = empty;  // workaround
  _nodeDataMap     = get(vertex_bundle,        _graph);
  _nodeVertexMap   = get(vertex_vertexPointer, _graph);
  _nodeIndexMap    = get(vertex_index,         _graph);
  _edgeDataMap     = get(edge_bundle,          _graph);
  _edgeParticleMap = get(edge_particlePointer, _graph);
  _edgeIndexMap    = get(edge_index,           _graph);
  _vertexNodeMap.clear();
  _particleEdgeMap.clear();
  _productionVertex = 0;
  _vertices.clear();
  _fsParticles.clear();
}


decayTopology::decayGraph
decayTopology::deepCopyGraph(const decayGraph& srcGraph,
			     const bool        copyFsParticles)
{
  // copy graph data structure
  decayGraph newGraph = srcGraph;
  // copy vertices
  nodeVertexType nodeVertexMap = get(vertex_vertexPointer, newGraph);
  for (nodeIterator iNode = vertices(newGraph).first; iNode != vertices(newGraph).second; ++iNode) {
    interactionVertex* srcVertex = nodeVertexMap[*iNode];
    if (srcVertex) {
      if (_debug)
  	printInfo << "copying vertex " << *iNode << endl;
      nodeVertexMap[*iNode] = &srcVertex->clone();
    }
  }  
  // copy particles
  //!!! what about beam particle?
  edgeParticleType edgeParticleMap = get(edge_particlePointer, newGraph);
  for (edgeIterator iEdge = edges(newGraph).first; iEdge != edges(newGraph).second; ++iEdge) {
    particle* srcParticle = edgeParticleMap[*iEdge];
    particle* newParticle = srcParticle;
    if (!srcParticle) {
      printErr << "null pointer to particle for edge " << *iEdge << ". aborting." << endl;
      throw;
    }
    if (_debug) {
      if (nodeVertexMap[target(*iEdge, newGraph)])
  	printInfo << "copying particle ";
      else if (copyFsParticles)
  	printInfo << "copying final state particle ";
    }
    if (nodeVertexMap[target(*iEdge, newGraph)] || copyFsParticles) {
      if (_debug)
  	cout << *iEdge << " '" << srcParticle->name() << "'" << endl;
      newParticle             = &srcParticle->clone();
      edgeParticleMap[*iEdge] = newParticle;
    }
    // modify the vertices that the particle connects
    nodeDesc           srcNode   = source(*iEdge, newGraph);
    interactionVertex* srcVertex = nodeVertexMap[srcNode];
    if (!srcVertex) {
      printErr << "null pointer to source vertex. aborting." << endl;
      throw;
    }
    for (unsigned int i = 0; i < srcVertex->nmbOutParticles(); ++i)
      if (srcVertex->outParticles()[i] == srcParticle) {
	srcVertex->outParticles()[i] = newParticle;
	break;
      }
    nodeDesc           targetNode   = target(*iEdge, newGraph);
    interactionVertex* targetVertex = nodeVertexMap[targetNode];
    if (targetVertex) {
      for (unsigned int i = 0; i < targetVertex->nmbInParticles(); ++i)
	if (targetVertex->inParticles()[i] == srcParticle) {
	  targetVertex->inParticles()[i] = newParticle;
	  break;
	}
    }
  }
  return newGraph;
}
