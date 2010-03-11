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

#include <boost/graph/topological_sort.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/breadth_first_search.hpp>

#include "utilities.h"
#include "isobarDecayVertex.h"
#include "decayTopology.h"

	
using namespace std;
using namespace boost;
using namespace rpwa;


bool decayTopology::_debug = false;


decayTopology::decayTopology()
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
    _graph           = topo._graph;
    _nodeProp        = topo._nodeProp;
    _edgeProp        = topo._edgeProp;
    _vertices        = topo._vertices;
    _fsParticles     = topo._fsParticles;
    _vertexNodeMap   = topo._vertexNodeMap;
    _particleEdgeMap = topo._particleEdgeMap;
  }
  return *this;
}


decayTopology&
decayTopology::constructDecay(const vector<particle*>&          fsParticles,
			      interactionVertex&                productionVertex,
			      const vector<interactionVertex*>& decayVertices)
{
  clear();
  const unsigned int nmbVert   = decayVertices.size();
  const unsigned int nmbFsPart = fsParticles.size();
  if (nmbFsPart == 0) {
    printWarn << "cannot construct decay topology without final state particles" << endl;
    clear();
    return *this;
  }
  if (nmbVert == 0) {
    printWarn << "cannot construct decay topology without vertices" << endl;
    clear();
    return *this;
  }
  // create graph node for production vertex and store pointer
  _nodeProp = get(vertex_name, _graph);
  {
    nodeDesc node = add_vertex(_graph);
    _nodeProp[node]                   = &productionVertex;
    _vertexNodeMap[&productionVertex] = node;
    if (_debug)
      printInfo << "added " << productionVertex << " as tree node [" << node << "]" << endl;
  }
  // create graph nodes for interaction vertices and store pointers
  for (unsigned int i = 0; i < nmbVert; ++i) {
    nodeDesc node = add_vertex(_graph);
    if (!decayVertices[i])
      printWarn << "zero pointer for decay vertex[" << i << "]" << endl;
    _nodeProp[node]                  = decayVertices[i];
    _vertexNodeMap[decayVertices[i]] = node;
    if (_debug)
      printInfo << "added " << *decayVertices[i] << " as tree node [" << node << "]" << endl;
  }
  // create final state nodes
  for (unsigned int i = 0; i < nmbFsPart; ++i) {
    nodeDesc node   = add_vertex(_graph);
    _nodeProp[node] = 0;
    if (_debug)
      printInfo << "created tree node [" << node << "] for final state "
		<< *fsParticles[i] << endl;
  }
  // create edges that connect the intercation nodes and store pointers to particles
  _edgeProp = get(edge_name, _graph);
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
	      if (!p)
		printWarn << "zero pointer to isobar particle" << endl;
	      _edgeProp[edge]     = p;
	      _particleEdgeMap[p] = edge;
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
	    if (!p)
	      printWarn << "zero pointer for final state particle[" << iFsPart << "]" << endl;
	    _edgeProp[edge]     = p;
	    _particleEdgeMap[p] = edge;
	    if (_debug)
	      printInfo << "adding edge from node [" << iFromVert->second << "] "
			<< "to final state node [" << nmbVert + 1 + iFsPart << "]" << endl;
	  }
	}
      }
  }
  // sorting nodes depth first
  vector<nodeDesc> nodes;
  depth_first_visit(_graph, _vertexNodeMap[&productionVertex],
		    makeDfsRecorder(back_inserter(nodes)),
		    get(vertex_color, _graph));
  if (_debug)
    printInfo << "depth first node order: ";
  for (unsigned int i = 0; i < nodes.size(); ++i) {
    cout << nodes[i] << ((i < nodes.size() - 1) ? ", " : "");
    interactionVertex* vertex = _nodeProp[nodes[i]];
    if (vertex)
      _vertices.push_back(vertex);
  }
  cout << endl;
  // copy final state particles
  _fsParticles = fsParticles;
  return *this;
}


bool
decayTopology::verifyTopology() const
{
  // check that decay topology is a tree of isobar decays
  bool          topologyOkay = true;
  nodeIndexType nodeIndex    = get(vertex_index, _graph);
  nodeIterator  iNode, iNodeEnd;
  unsigned int  countNodesWithNoInputEdge = 0;
  for (tie(iNode, iNodeEnd) = boost::vertices(_graph); iNode != iNodeEnd; ++iNode) {
    const unsigned int i = get(nodeIndex, *iNode);
    // check that each node has exactly 1 incoming edge (isobar)
    unsigned int nmbInEdges = in_degree(*iNode, _graph);
    if (nmbInEdges == 0) {
      ++countNodesWithNoInputEdge;
      if (_debug)
	printInfo << "assuming node[" << i << "] is production node" << endl;
    } else if (nmbInEdges != 1) {
      printWarn << "number of input edges of node[" << i << "] is "
    		<< nmbInEdges << " != 1" << endl;
      topologyOkay = false;
    } else if (_debug)
      printInfo << "number of input edges of node[" << i << "] is correct" << endl;
    if (countNodesWithNoInputEdge > 1) {
      printWarn << "there are " << countNodesWithNoInputEdge
		<< " nodes with no no input edges." << endl;
      topologyOkay = false;
    }
    // check that for each node the number of outgoing edges is either 2 (decay node) or 0 (final state node)
    unsigned int nmbOutEdges = out_degree(*iNode, _graph);
    if (nmbOutEdges == 0) {
      if (_nodeProp[*iNode]) {
	printWarn << "node[" << i << "] has no outgoing edges, "
		  << "but has a interactionVertex pointer assigned." << endl;
	topologyOkay = false;
      } else if (_debug)
	printInfo << "final state node[" << i << "] is correct" << endl;
    } else if (nmbOutEdges == 2) {
      if (!_nodeProp[*iNode]) {
	printWarn << "node[" << i << "] has 2 outgoing edges, "
		  << "but no interactionVertex pointer assigned." << endl;
	topologyOkay = false;
      } else if (_debug)
	printInfo << "interaction node[" << i << "] is correct" << endl;
    } else {
      printWarn << "number of output edges of node[" << i << "] is "
		<< nmbOutEdges << " != 0 or 2" << endl;
      topologyOkay = false;
    }
  }
  return topologyOkay;
}


const TLorentzVector&
decayTopology::updateIsobarLzVec()
{
  // loop over decay vertices (leaving out production vertex)
  // propagate changes from final state particles up to X-system
  for (unsigned int i = nmbVertices() - 1; i > 0; --i) {
    if (_debug)
      printInfo << "updating Lorentz-vector of mother isobar of "
		<< "node[" << _vertexNodeMap.find(_vertices[i])->second << "]" << endl;
    static_cast<isobarDecayVertex*>(_vertices[i])->updateMotherLzVec();
  }
  return static_cast<isobarDecayVertex*>(&xDecayVertex())->mother().lzVec();
}


ostream&
decayTopology::print(ostream& out) const
{
  out << "decay topology nodes:" << endl;
  for (unsigned int i = 0; i < nmbVertices(); ++i) {
    cout << "    " << ((i == 0) ? "production" : "decay")
	 << " node[" << _vertexNodeMap.find(_vertices[i])->second << "] = ";
    for (unsigned int j = 0; j < _vertices[i]->nmbInParticles(); ++j) {
      const particle* p = _vertices[i]->inParticles()[j];
      out << p->name() << sign(p->charge());
      if (j < _vertices[i]->nmbInParticles() - 1)
	out << " + ";
    }
    out << "  --->  ";
    for (unsigned int j = 0; j < _vertices[i]->nmbOutParticles(); ++j) {
      const particle* p = _vertices[i]->outParticles()[j];
      out << p->name() << sign(p->charge());
      if (j < _vertices[i]->nmbOutParticles() - 1)
	out << " + ";
    }
    out << endl;
  }
  out << "decay topology particles:" << endl;
  edgeIterator iEdge, iEdgeEnd;
  for (tie(iEdge, iEdgeEnd) = edges(_graph); iEdge != iEdgeEnd; ++iEdge) {
    const particle* p = _edgeProp[*iEdge];
    out << "    particle" << *iEdge << " ";
    if (p) {
      out << p->name() << sign(p->charge()) << " ";
      if (!_nodeProp[target(*iEdge, _graph)])
	out << "(final state)";
      out << endl;
    } else
      out << "error! zero pointer to particle" << endl;
  }
  return out;
}


ostream&
decayTopology::writeGraphViz(ostream& out) const
{
  write_graphviz(out, _graph);
  return out;
}


void
decayTopology::clear()
{
  _graph.clear();
  _nodeProp = get(vertex_name, _graph);
  _edgeProp = get(edge_name, _graph);
  _vertices.clear();
  _fsParticles.clear();
  _vertexNodeMap.clear();
  _particleEdgeMap.clear();
}
