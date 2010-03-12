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
//      amplitude calculation of isobar decay
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "utilities.h"
#include "isobarDecayTopology.h"

	
using namespace std;
using namespace boost;
using namespace rpwa;


bool isobarDecayTopology::_debug = false;


isobarDecayTopology::isobarDecayTopology()
  : decayTopology()
{ }


isobarDecayTopology::isobarDecayTopology(const isobarDecayTopology& topo)
  : decayTopology()
{
  *this = topo;
}


isobarDecayTopology::isobarDecayTopology(const vector<particle*>&          fsParticles,
					 interactionVertex&                productionVertex,
					 const vector<isobarDecayVertex*>& decayVertices)
  : decayTopology()
{
  constructDecay(fsParticles, productionVertex, decayVertices);
}


isobarDecayTopology::~isobarDecayTopology()
{ }


isobarDecayTopology&
isobarDecayTopology::operator = (const isobarDecayTopology& topo)
{
  if (this != &topo) {
    decayTopology::operator = (topo);
    _vertices = topo._vertices;
  }
  return *this;
}


isobarDecayTopology&
isobarDecayTopology::constructDecay(const vector<particle*>&          fsParticles,
				    interactionVertex&                productionVertex,
				    const vector<isobarDecayVertex*>& isobarDecayVertices)
{
  const unsigned int nmbVert = isobarDecayVertices.size();
  if (_debug)
    printInfo << "constructing isobar decay topology with "
	      << fsParticles.size() << " final state particles and "
	      << nmbVert            << " isobar decay vertices" << endl;
  clear();
  vector<interactionVertex*> vertices(nmbVert, 0);
  for (unsigned int i = 0; i < nmbVert; ++i)
    vertices[i] = isobarDecayVertices[i];
  decayTopology::constructDecay(fsParticles, productionVertex, vertices);
  _vertices.resize(nmbVertices(), 0);
  for (unsigned int i = 0; i < nmbVertices(); ++i)
    _vertices[i] = static_cast<isobarDecayVertex*>(interactionVertices()[i]);
  return *this;
}


bool
isobarDecayTopology::verifyTopology() const
{
  // check that decay topology is a tree of isobar decays
  bool          topologyOkay = true;
  nodeIndexType nodeIndex    = get(vertex_index, _graph);
  nodeIterator  iNode, iNodeEnd;
  unsigned int  countNodesWithNoInputEdge   = 0;
  unsigned int  countNodesWithOneOutputEdge = 0;
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
		<< " nodes with no input edges." << endl;
      topologyOkay = false;
    }
    // check that for each node the number of outgoing edges is either 2 (decay node) or 0 (final state node)
    unsigned int nmbOutEdges = out_degree(*iNode, _graph);
    if (nmbOutEdges == 0) {
      if (_nodeProp[*iNode]) {
	printWarn << "node[" << i << "] has no outgoing edges, "
		  << "but has a isobar decay vertex pointer assigned." << endl;
	topologyOkay = false;
      } else if (_debug)
	printInfo << "final state node[" << i << "] is correct" << endl;
    } else if (nmbOutEdges == 2) {
      if (!_nodeProp[*iNode]) {
	printWarn << "node[" << i << "] has 2 outgoing edges, "
		  << "but no isobar decay vertex pointer assigned." << endl;
	topologyOkay = false;
      } else if (_debug)
	printInfo << "interaction node[" << i << "] is correct" << endl;
    } else if (nmbOutEdges == 1) {
      ++countNodesWithOneOutputEdge;
      if (_debug)
	printInfo << "assuming node[" << i << "] is production node" << endl;
    } else {
      printWarn << "number of output edges of node[" << i << "] is "
		<< nmbOutEdges << " != 0 or 2" << endl;
      topologyOkay = false;
    }
    if (countNodesWithOneOutputEdge > 1) {
      printWarn << "there are " << countNodesWithOneOutputEdge
		<< " nodes with one output edge." << endl;
      topologyOkay = false;
    }
  }
  return topologyOkay;
}


const TLorentzVector&
isobarDecayTopology::updateIsobarLzVec()
{
  // loop over decay vertices (leaving out production vertex)
  // propagate changes from final state particles up to X-system
  for (int i = nmbVertices() - 1; i >= 0; --i) {
    if (_debug)
      printInfo << "updating Lorentz-vector of mother isobar of "
		<< "node[" << _vertexNodeMap.find(_vertices[i])->second << "]" << endl;
   _vertices[i]->updateMotherLzVec();
  }
  return xIsobarDecayVertex().mother().lzVec();
}


void
isobarDecayTopology::clear()
{
  decayTopology::clear();
}
