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


#include <algorithm>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/copy.hpp>

#include "utilities.h"
#include "diffractiveDissVertex.h"
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
					 const vector<interactionVertex*>& interactionVertices)
  : decayTopology()
{
  constructDecay(fsParticles, productionVertex, interactionVertices);
}


isobarDecayTopology::isobarDecayTopology(const vector<particle*>&          fsParticles,
					 interactionVertex&                productionVertex,
					 const vector<isobarDecayVertex*>& isobarDecayVertices)
  : decayTopology()
{
  constructDecay(fsParticles, productionVertex, isobarDecayVertices);
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
				    const vector<interactionVertex*>& interactionVertices)
{
  const unsigned int nmbVert = interactionVertices.size();
  if (_debug)
    printInfo << "constructing isobar decay topology with "
	      << fsParticles.size() << " final state particles and "
	      << nmbVert            << " isobar decay vertices" << endl;
  clear();
  for (unsigned int i = 0; i < nmbVert; ++i)
    if (!dynamic_cast<isobarDecayVertex*>(interactionVertices[i])) {
      printErr << "interaction vertex[" << i << "] is not an isobarDecayVertex. aborting." << endl;
      throw;
    }
  decayTopology::constructDecay(fsParticles, productionVertex, interactionVertices);
  _vertices.resize(nmbVertices(), 0);
  for (unsigned int i = 0; i < nmbVertices(); ++i)
    _vertices[i] = static_cast<isobarDecayVertex*>(this->interactionVertices()[i]);
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



isobarDecayTopology& 
isobarDecayTopology::addSubDecay(const isobarDecayTopology& subDecay) {

  cerr << "Before copying vertices" << endl;
  
  std::vector<interactionVertex*> vertices=this->interactionVertices();
  std::vector<particle*>          fsPart=this->fsParticles(); 

  cerr << "Before adding vertices" << endl;
  cerr << vertices.size() << "   " << this->interactionVertices().size()<<endl;
  cerr << fsPart.size() << endl;
  vertices.insert(vertices.end(),subDecay._vertices.begin(),subDecay._vertices.end());
 
  cerr << "Before adding final state particles" << endl;
 
  fsPart.insert(fsPart.end(),subDecay._fsParticles.begin(),subDecay._fsParticles.end());

  interactionVertex* prod=_productionVertex;

 cerr << vertices.size() << endl;
  cerr << fsPart.size() << endl;
  cerr << "Before ConstructDecay" << endl;

  return constructDecay(fsPart,*prod,vertices);
 }



bool
isobarDecayTopology::verifyTopology() const
{
  // check that decay topology is a tree of isobar decays
  bool          topologyOkay = true;
  nodeIterator  iNode, iNodeEnd;
  unsigned int  countNodesWithNoInputEdge   = 0;
  unsigned int  countNodesWithOneOutputEdge = 0;
  for (tie(iNode, iNodeEnd) = boost::vertices(_graph); iNode != iNodeEnd; ++iNode) {
    const unsigned int i = get(_nodeIndexMap, *iNode);
    // check that each node has exactly 1 incoming edge (isobar)
    const unsigned int nmbInEdges = in_degree(*iNode, _graph);
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
    const unsigned int nmbOutEdges = out_degree(*iNode, _graph);
    if (nmbOutEdges == 0) {
      if (_nodeVertexMap[*iNode]) {
	printWarn << "node[" << i << "] has no outgoing edges, "
		  << "but has an isobar decay vertex pointer assigned." << endl;
	topologyOkay = false;
      } else if (_debug)
	printInfo << "final state node[" << i << "] is correct" << endl;
    } else if (nmbOutEdges == 2) {
      if (!_nodeVertexMap[*iNode]) {
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


bool 
isobarDecayTopology::checkConsistency() const {
  bool allVertConsistent = true;
  for (unsigned int i = 0; i < _vertices.size(); ++i) {
    if (_debug)
      printInfo << "checking consistency of isobar decay vertex "
		<< _vertices[i]->mother().name() << "  --->  "
		<< _vertices[i]->daughter1().name() << "  +  "
		<< _vertices[i]->daughter2().name() << endl;
    if (!_vertices[i]->checkConsistency())
      allVertConsistent = false;
  }
  return allVertConsistent;
}


// !NOTE! the following is all very ugly
// in order to make such kind of operations easier and more
// transparent we need a better class structure for the topologies
// and the graph
vector<isobarDecayTopology*>
isobarDecayTopology::possibleDecays()
{
  // order nodes pedth-first
  const unsigned int nmbVert = num_vertices(_graph) - 1;
  vector<nodeDesc>   startNodes;
  {
    vector<default_color_type> colors(num_vertices(_graph));
    depth_first_visit(_graph, _vertexNodeMap[&xInteractionVertex()],
		      makeDfsRecorder(back_inserter(startNodes)),
		      make_iterator_property_map(colors.begin(),
						 get(vertex_index, _graph), colors[0]));
  }
  // create array of possible subraphs
  vector<decayGraph> subGraphs(nmbVert);
  for (unsigned int i = 0; i < nmbVert; ++i) {
    // find all nodes below current node
    vector<nodeDesc> subGraphNodes;
    {
      vector<default_color_type> colors(num_vertices(_graph));
      depth_first_visit(_graph, startNodes[i],
			makeDfsRecorder(back_inserter(subGraphNodes)),
			make_iterator_property_map(colors.begin(),
						   get(vertex_index, _graph), colors[0]));
    }
    if (_debug)
      printInfo << "creating subgraph[" << i << "] starting at node " << startNodes[i]
		<< " using nodes: " << flush;
    subGraphs[i] = _graph.create_subgraph();
    for (unsigned int j = 0; j < subGraphNodes.size(); ++j) {
      add_vertex(subGraphNodes[j], subGraphs[i]);
      if (_debug)
	cout << subGraphNodes[j] << "    ";
    }
    if (_debug)
      cout << endl;
  }
  reverse(subGraphs.begin(),  subGraphs.end());
  reverse(startNodes.begin(), startNodes.end());
  if (_debug)
    for (unsigned int i = 0; i < subGraphs.size(); ++i) {
      // printInfo << "created subgraph[" << i << "]:" << endl;
      // print_graph (subGraphs[i], get(vertex_index, subGraphs[i]));
      // print_edges2(subGraphs[i], get(vertex_index, subGraphs[i]), get(edge_index, subGraphs[i]));
    }
  // create all possible decays
  map<nodeDesc, vector<decayGraph> > decayPossibilities;  // memorizes all possible decay graphs starting at the respective node
  for (unsigned int i = 0; i < nmbVert; ++i) {
    if (!_nodeVertexMap[startNodes[i]]) {
      // final state node
      decayPossibilities[startNodes[i]].push_back(subGraphs[i]);
      continue;
    }
    // get daughter nodes
    vector<nodeDesc> daughterNodes;
    for (adjIterator iNode = adjacent_vertices(startNodes[i], _graph).first;
	 iNode != adjacent_vertices(startNodes[i], _graph).second; ++iNode)
      daughterNodes.push_back(*iNode);
    if (daughterNodes.size() != 2) {
      printErr << "number of daughter vertices of subgraph top node is "
	       << daughterNodes.size() << " != 2. aborting." << endl;
      throw;
    }
    // sort daughter nodes
    {
      isobarDecayVertex* dv[2] = {
	static_cast<isobarDecayVertex*>(_nodeVertexMap[daughterNodes[0]]),
	static_cast<isobarDecayVertex*>(_nodeVertexMap[daughterNodes[1]])
      };
      isobarDecayVertex* mv = static_cast<isobarDecayVertex*>(_nodeVertexMap[startNodes[i]]);
      unsigned int index = 0;
      if (!dv[0] && dv[1])
	index = 1;
      if (dv[0] || dv[1])
	if (&dv[index]->mother() != mv->outParticles()[index])
	  swap(daughterNodes[0], daughterNodes[1]);
    }
    // get daughter nodes
    vector<decayGraph>& decays          = decayPossibilities[startNodes[i]];
    vector<decayGraph>& daughter1Decays = decayPossibilities[daughterNodes[0]];
    vector<decayGraph>& daughter2Decays = decayPossibilities[daughterNodes[1]];
    for (unsigned int i1 = 0; i1 < daughter1Decays.size(); ++i1)
      for(unsigned int i2 = 0; i2 < daughter2Decays.size(); ++i2) {
	//cout << "!!! " << i1 << ", " << i2 << endl;
	// join subgraphs
	//!!! if not taken care of by the calling code this is a potential memory leak
	nodeDesc   topNode;
	decayGraph joinedGraph = joinDaughterGraphs(*static_cast<isobarDecayVertex*>(_nodeVertexMap[startNodes[i]]),
						    daughter1Decays[i1],
						    daughter2Decays[i2],
						    topNode);
	isobarDecayVertex* vertex = static_cast<isobarDecayVertex*>(get(vertex_vertexPointer, joinedGraph, topNode));
	//cout << "!!! " << topNode << ", " << *vertex << endl;
	particle&          d1     = vertex->daughter1();
	particle&          d2     = vertex->daughter2();
	//cout << "!!! " << d2 << endl;
	// create all possible quantum number combinations for this node
	// quantum numbers fixed by daughter quantum numbers
	const int baryonNmb   = d1.baryonNmb()   + d2.baryonNmb();
	const int charge      = d1.charge()      + d2.charge();
	const int strangeness = d1.strangeness() + d2.strangeness();
	const int charm       = d1.charm()       + d2.charm();
	const int beauty      = d1.beauty()      + d2.beauty();
	const int G           = d1.G()           * d2.G();
	const int C           = d1.C()           * d2.C();
	// daughter quantum numbers that define ranges of possible mother quantum numbers
	const int s1 = d1.J();
	const int s2 = d2.J();
	const int I1 = d1.isospin();
	const int I2 = d2.isospin();

	const int  minS = 0;
	const int  maxS = 2;
	const int  minL = 0;
	const int  maxL = 2;
	const int  minJ = 0;
	const int  maxJ = 2;
	const int  minI = 2;
	const int  maxI = 2;
	const bool allowJpcExotic = true;
  
	for (int S = max(abs(s1 - s2), minS); S <= min(s1 + s2, maxS); S += 2) {        // loop over all allowed total spins
	  for (int L = max(0, minL); L <= maxL; L += 2) {                               // loop over all allowed relative orbital angular momenta
	    const int P = d1.P() * d2.P() * (L % 4 == 0 ? 1 : -1);  // parity
	    for (int J = max(abs(L - S), minJ); J <= min(L + S, maxJ); J += 2) {        // L-S coupling loop
	      for (int I = max(abs(I1 - I2), minI); I <= min(I1 + I2, maxI); I += 2) {  // isospin loop
		//int C = G * (I % 4 == 0 ? 1 : -1);  // C-Parity???
		if (_debug)
		  printInfo << vertex->mother().name() << "  --->  "
			    << d1.name() << "  +   " << d2.name() << "    "
			    << "2IG = "  << I << sign(G) << ", "
			    << "2JPC = " << J << sign(P) << sign(C) << ", "
			    << "2S = "   << S <<", 2L = " << L
			    << endl;
		if (!allowJpcExotic) {
		  // quark model boundary conditions:
		  // for P == C everything ok; check P = (-1)^(J + 1)
		  if ((P != C) && (   (C != (J % 4 == 0     ? 1 : -1))
				   || (P != (J + 2 % 4 == 0 ? 1 : -1)))) {
		    if (_debug)
		      printInfo << "disregarding spin-exotic isobar above" << endl;
		    continue;
		  }
		}
		//!!! if not taken care of by the calling code this is a potential memory leak
		decayGraph graphCopy = deepCopyGraph(joinedGraph, false);
		// set mother vertex quantum numbers
		isobarDecayVertex* v = static_cast<isobarDecayVertex*>(get(vertex_vertexPointer, graphCopy, topNode));
		v->setL(L);
		v->setS(S);
		// copy mother particle
		particle* m = &v->mother().clone();
		v->inParticles()[0] = m;
		// m->setName("foo");
		m->setCharge     (charge);
		m->setBaryonNmb  (baryonNmb);
		m->setIsospin    (I);
		m->setStrangeness(strangeness);
		m->setCharm      (charm);
		m->setBeauty     (beauty);
		m->setG          (G);
		m->setJ          (J);
		m->setP          (P);
		m->setC          (C);
		//cout << "!!! " << topNode << ", " << *v << endl;
		decays.push_back(graphCopy);
	      }  // isospin loop
	    }  // L-S coupling loop
	  }  // L loop
	}  // S loop
      }  // loop over daughter decays
  }  // loop over all start nodes
  // construct isobar decay topologies for all decays
  vector<decayGraph>&          decays = decayPossibilities[startNodes[nmbVert - 1]];
  vector<isobarDecayTopology*> decayTopos;
  for (unsigned int i = 0; i < decays.size(); ++i) {
    // get decay vertices
    vector<isobarDecayVertex*> decayVertices;
    for (nodeIterator iNode = vertices(decays[i]).first;
	 iNode != vertices(decays[i]).second; ++iNode) {
      isobarDecayVertex* v;
      v = static_cast<isobarDecayVertex*>(get(vertex_vertexPointer, decays[i], *iNode));
      if (v)
	decayVertices.push_back(v);
    }
    // construct production vertex
    particle& beam = static_cast<diffractiveDissVertex*>(&productionVertex())->beam();
    particle& X    = static_cast<isobarDecayVertex*>(get(vertex_vertexPointer, decays[i], 0))->mother();
    diffractiveDissVertex* prodVert = new diffractiveDissVertex(beam, X);
    // construct isobar decay topology
    decayTopos.push_back(new isobarDecayTopology(fsParticles(), *prodVert, decayVertices));
  }
  return decayTopos;
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


decayTopology::decayGraph
isobarDecayTopology::joinDaughterGraphs(const isobarDecayVertex& motherVertex,
					const decayGraph&        daughterGraph1,
					const decayGraph&        daughterGraph2,
					nodeDesc&                newMotherNode)
{
  //cout << "!!! old mother vertex: " << motherVertex << endl;
  decayGraph newGraph;
  newMotherNode = add_vertex(newGraph);
  copy_graph(daughterGraph1, newGraph);
  copy_graph(daughterGraph2, newGraph);
  nodeDesc daughterNodes[2] = {1, 1 + num_vertices(daughterGraph1)};
  // copy isobar decay vertex at topNode
  isobarDecayVertex* newMotherVertex = &motherVertex.clone();
  // connect subgraphs to new top node
  for (unsigned int i = 0; i < 2; ++i) {
    bool     inserted;
    edgeDesc edge;
    tie(edge, inserted) = add_edge(newMotherNode, daughterNodes[i], newGraph);
    isobarDecayVertex* daughterVertex;
    daughterVertex = static_cast<isobarDecayVertex*>(get(vertex_vertexPointer, newGraph, daughterNodes[i]));
    if (!daughterVertex) {
      // final state vertex
      put(edge_particlePointer, newGraph, edge, motherVertex.outParticles()[i]);
    } else {
      // isobar decay vertex
      put(edge_particlePointer, newGraph, edge, &daughterVertex->mother());
      newMotherVertex->outParticles()[i] = &daughterVertex->mother();
      //cout << "!!! daughter vertex " << i << ": " << *daughterVertex << endl;
    }
  }
  put(vertex_vertexPointer, newGraph, newMotherNode, newMotherVertex);
  //cout << "!!! new mother vertex: " << *newMotherVertex << endl;
  return newGraph;
}
