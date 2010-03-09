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

#include "utilities.h"
#include "decayTopology.h"

	
using namespace std;
using namespace boost;
using namespace rpwa;


bool decayTopology::_debug = false;


decayTopology::decayTopology()
  : _nmbVertices   (0),
    _nmbFsParticles(0)
{ }


decayTopology::decayTopology(const decayTopology& topo)
{
  *this = topo;
}


decayTopology::~decayTopology()
{ }


decayTopology&
decayTopology::operator = (const decayTopology& topo)
{
  if (this != &topo) {
    _graph          = topo._graph;
    _vertexProp     = topo._vertexProp;
    _edgeProp       = topo._edgeProp;
    _nmbVertices    = topo._nmbVertices;
    _nmbFsParticles = topo._nmbFsParticles;
  }
  return *this;
}


decayTopology&
decayTopology::constructDecay(const vector<particle*>&          fsParticles,
			      const vector<interactionVertex*>& vertices)
{
  _nmbVertices    = vertices.size();
  _nmbFsParticles = fsParticles.size();
  if (_nmbFsParticles == 0) {
    printWarn << "cannot construct decay topology without final state particles" << endl;
    clear();
    return *this;
  }
  if (_nmbVertices == 0) {
    printWarn << "cannot construct decay topology without vertices" << endl;
    clear();
    return *this;
  }
  // create graph vertices for interaction vertices and store pointers
  _vertexProp = get(vertex_name, _graph);
  for (unsigned int i = 0; i < _nmbVertices; ++i) {
    vertexDesc v  = add_vertex(_graph);
    _vertexProp[v] = vertices[i];
  }
  // create final state vertices
  for (unsigned int i = 0; i < _nmbFsParticles; ++i) {
    vertexDesc v  = add_vertex(_graph);
    _vertexProp[v] = 0;
  }
  // create edges that connect the intercation vertices and store pointers to particles
  _edgeProp = get(edge_name, _graph);
  for (unsigned int iFromVert = 0; iFromVert < _nmbVertices; ++iFromVert)
    for (unsigned int iToVert = 0; iToVert < _nmbVertices; ++iToVert) {
      interactionVertex* fromVertex = vertices[iFromVert];
      for (unsigned int iOutPart = 0; iOutPart < fromVertex->nmbOutParticles(); ++iOutPart) {
	interactionVertex* toVertex = vertices[iToVert];
	for (unsigned int iInPart = 0; iInPart < toVertex->nmbInParticles(); ++iInPart)
	  if (fromVertex->outParticles()[iOutPart] == toVertex->inParticles()[iInPart]) {
	    bool     inserted;
	    edgeDesc edge;
	    tie(edge, inserted) = add_edge(iFromVert, iToVert, _graph);
	    if (inserted)
	      _edgeProp[edge] = fromVertex->outParticles()[iOutPart];
	  }
      }
    }
  // create edges for final state particles and store pointers
  for (unsigned int iFromVert = 0; iFromVert < _nmbVertices; ++iFromVert) {
    interactionVertex* fromVertex = vertices[iFromVert];
    for (unsigned int iOutPart = 0; iOutPart < fromVertex->nmbOutParticles(); ++iOutPart)
      for (unsigned int iFsPart = 0; iFsPart < _nmbFsParticles; ++iFsPart)
	if (fromVertex->outParticles()[iOutPart] == fsParticles[iFsPart]) {
	  bool     inserted;
	  edgeDesc edge;
	  tie(edge, inserted) = add_edge(iFromVert, _nmbVertices + iFsPart, _graph);
	  if (inserted)
	    _edgeProp[edge] = fromVertex->outParticles()[iOutPart];
	}
  }
  return *this;
}


ostream&
decayTopology::print(ostream& out) const
{
  out << "decay topology vertices:" << endl;
  list<vertexDesc> vertices;
  topological_sort(_graph, front_inserter(vertices));
  for (list<vertexDesc>::const_iterator iVert = vertices.begin(); iVert != vertices.end(); ++iVert) {
    const interactionVertex* v = _vertexProp[*iVert];
    if (v)
      out << "    decay ";
    else
      out << "    final state ";
    out << "vertex[" << *iVert << "] ";
    if (v) {
      for (unsigned int i = 0; i < v->nmbInParticles(); ++i) {
	const particle* p = v->inParticles()[i];
	out << p->name() << sign(p->charge());
	if (i < v->nmbInParticles() - 1)
	  out << " + ";
      }
      out << " -> ";
      for (unsigned int i = 0; i < v->nmbOutParticles(); ++i) {
	const particle* p = v->outParticles()[i];
	out << p->name() << sign(p->charge());
	if (i < v->nmbOutParticles() - 1)
	  out << " + ";
      }
      out << endl;
    }
  }
  out << "decay topology particles:" << endl;
  edgeIterator iEdge, iEdgeEnd;
  for (tie(iEdge, iEdgeEnd) = edges(_graph); iEdge != iEdgeEnd; ++iEdge) {
    const particle* p = _edgeProp[*iEdge];
    out << "    particle" << *iEdge << " ";
    if (p) {
      out << p->name() << sign(p->charge()) << " ";
      if (!_vertexProp[target(*iEdge, _graph)])
	out << "final state particle";
      out << endl;
    } else
      out << "error! zero pointer to particle" << endl;
  }
  return out;
}


void
decayTopology::clear()
{
  _graph.clear();
  _vertexProp     = get(vertex_name, _graph);
  _edgeProp       = get(edge_name, _graph);
  _nmbVertices    = 0;
  _nmbFsParticles = 0;
}
