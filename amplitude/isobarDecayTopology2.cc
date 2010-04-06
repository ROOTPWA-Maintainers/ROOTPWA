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
#include "diffractiveDissVertex2.h"
#include "isobarDecayTopology2.h"
#include "particleDataTable.h"

	
using namespace std;
using namespace boost;
using namespace rpwa;


bool isobarDecayTopology2::_debug = false;


isobarDecayTopology2::isobarDecayTopology2()
  : decayTopology2()
{ }


isobarDecayTopology2::isobarDecayTopology2(const isobarDecayTopology2& topo)
  : decayTopology2()
{
  *this = topo;
}


isobarDecayTopology2::isobarDecayTopology2(const interactionVertexPtr&         productionVertex,
					   const vector<isobarDecayVertexPtr>& isobarDecayVertices,
					   const vector<particlePtr>&          fsParticles)
  : decayTopology2()
{
  constructDecay(productionVertex, isobarDecayVertices, fsParticles);
}


isobarDecayTopology2::isobarDecayTopology2(const interactionVertexPtr&         productionVertex,
					   const vector<interactionVertexPtr>& isobarDecayVertices,
					   const vector<particlePtr>&          fsParticles)
  : decayTopology2()
{
  constructDecay(productionVertex, isobarDecayVertices, fsParticles);
}


isobarDecayTopology2::~isobarDecayTopology2()
{ }


isobarDecayTopology2&
isobarDecayTopology2::operator =(const isobarDecayTopology2& topo)
{
  if (this != &topo) {
    decayTopology2::operator =(topo);
    _isobarVertices = topo._isobarVertices;
  }
  return *this;
}


void
isobarDecayTopology2::clear()
{
  decayTopology2::clear();
  _isobarVertices.clear();
}


isobarDecayTopology2&
isobarDecayTopology2::constructDecay(const interactionVertexPtr&         productionVertex,
				     const vector<isobarDecayVertexPtr>& isobarDecayVertices,
				     const vector<particlePtr>&          fsParticles)
{
  const unsigned int nmbVert = isobarDecayVertices.size();
  if (_debug)
    printInfo << "constructing isobar decay topology with "
	      << fsParticles.size() << " final state particles and "
	      << nmbVert            << " isobar decay vertices" << endl;
  vector<interactionVertexPtr> intVertices(nmbVert);
  for (unsigned int i = 0; i < nmbVert; ++i)
    intVertices[i] = isobarDecayVertices[i];
  decayTopology2::constructDecay(productionVertex, intVertices, fsParticles);
  // copy sorted vertices
  _isobarVertices.resize(nmbInteractionVertices());
  for (unsigned int i = 0; i < nmbInteractionVertices(); ++i)
    _isobarVertices[i] = static_pointer_cast<isobarDecayVertex2>(interactionVertices()[i]);
  if (nmbInteractionVertices() != nmbVert) {
    printErr << "number of interaction vertices  = " << nmbInteractionVertices()
	     << " does not match number of vertices given in parameter array = " << nmbVert
	     << ". aborting." << endl;
    throw;
  }
  return *this;
}


isobarDecayTopology2&
isobarDecayTopology2::constructDecay(const interactionVertexPtr&         productionVertex,
				     const vector<interactionVertexPtr>& isobarDecayVertices,
				     const vector<particlePtr>&          fsParticles)
{
  const unsigned int nmbVert = isobarDecayVertices.size();
  vector<isobarDecayVertexPtr> decayVertices(nmbVert);
  for (unsigned int i = 0; i < nmbVert; ++i) {
    decayVertices[i] = dynamic_pointer_cast<isobarDecayVertex2>(isobarDecayVertices[i]);
    if (!decayVertices[i]) {
      printErr << "interaction vertex[" << i << "] is not an isobarDecayVertex2. aborting." << endl;
      throw;
    }
  }
  return constructDecay(productionVertex, decayVertices, fsParticles);
}


bool
isobarDecayTopology2::checkTopology() const
{
  // perform basic checks of topology
  bool topologyOkay = decayTopology2::checkTopology();
  // check that decay topology is a tree of isobar decays
  for (unsigned int i = 0; i < nmbInteractionVertices(); ++i) {
    const int index = node(_isobarVertices[i]);
    // check that each isobar decay node has exactly 1 incoming edge (isobar)
    const unsigned int nmbIn = nmbInEdges(_isobarVertices[i]);
    if (nmbIn != 1) {
      printWarn << "number of incoming edges of node[" << index << "] is "
		<< nmbIn << " != 1" << endl;
      topologyOkay = false;
    } else if (_debug)
      printInfo << "number of incoming edges of node[" << index << "] = " << nmbIn
		<< " is correct" << endl;
    // check that for each isobar decay node the number of outgoing edges is 2
    const unsigned int nmbOut = nmbOutEdges(_isobarVertices[i]);
    if (nmbOut != 2) {
      printWarn << "number of outgoing edges of node[" << index << "] is "
    		<< nmbOut << " != 2" << endl;
      topologyOkay = false;
    } else if (_debug)
      printInfo << "number of outgoing edges of node[" << index << "] = " << nmbOut
		<< " is correct" << endl;
  }
  return topologyOkay;
}


bool 
isobarDecayTopology2::checkConsistency() const
{
  bool allVertConsistent = true;
  for (unsigned int i = 0; i < _isobarVertices.size(); ++i) {
    const bool vertexConsistent = _isobarVertices[i]->checkConsistency();
    if (!vertexConsistent) {
      allVertConsistent = false;
      if (_debug)
	printInfo << "isobar decay vertex " << *_isobarVertices[i] << " is "
		  << ((vertexConsistent) ? "" : "NOT ") << "consistent" << endl;
    }
  }
  return allVertConsistent;
}


// isobarDecayTopology2& 
// isobarDecayTopology2::addSubDecay(const isobarDecayTopology2& subDecay) {

//   cerr << "Before copying vertices" << endl;
  
//   std::vector<interactionVertex*> vertices=this->interactionVertices();
//   std::vector<particle*>          fsPart=this->fsParticles(); 

//   cerr << "Before adding vertices" << endl;
//   cerr << vertices.size() << "   " << this->interactionVertices().size()<<endl;
//   cerr << fsPart.size() << endl;
//   vertices.insert(vertices.end(),subDecay._isobarVertices.begin(),subDecay._isobarVertices.end());
 
//   cerr << "Before adding final state particles" << endl;
 
//   fsPart.insert(fsPart.end(),subDecay._fsParticles.begin(),subDecay._fsParticles.end());

//   interactionVertex* prod=_productionVertex;

//   cerr << vertices.size() << endl;
//   cerr << fsPart.size() << endl;
//   cerr << "Before ConstructDecay" << endl;

//   return constructDecay(fsPart,*prod,vertices);
// }



// // !NOTE! the following is all very ugly
// // in order to make such kind of operations easier and more
// // transparent we need a better class structure for the topologies
// // and the graph
// vector<isobarDecayTopology2*>
// isobarDecayTopology2::possibleDecays()
// {
//   // order nodes depth-first
//   const unsigned int nmbVert = num_vertices(_graph) - 1;
//   vector<nodeDesc>   startNodes;
//   {
//     vector<default_color_type> colors(num_vertices(_graph));
//     depth_first_visit(_graph, _vertexNodeMap[&xInteractionVertex()],
// 		      makeDfsRecorder(back_inserter(startNodes)),
// 		      make_iterator_property_map(colors.begin(),
// 						 get(vertex_index, _graph), colors[0]));
//   }
//   // create array of possible subraphs
//   vector<decayGraph> subGraphs(nmbVert);
//   for (unsigned int i = 0; i < nmbVert; ++i) {
//     // find all nodes below current node
//     vector<nodeDesc> subGraphNodes;
//     {
//       vector<default_color_type> colors(num_vertices(_graph));
//       depth_first_visit(_graph, startNodes[i],
// 			makeDfsRecorder(back_inserter(subGraphNodes)),
// 			make_iterator_property_map(colors.begin(),
// 						   get(vertex_index, _graph), colors[0]));
//     }
//     if (_debug)
//       printInfo << "creating subgraph[" << i << "] starting at node " << startNodes[i]
// 		<< " using nodes: " << flush;
//     subGraphs[i] = _graph.create_subgraph();
//     for (unsigned int j = 0; j < subGraphNodes.size(); ++j) {
//       add_vertex(subGraphNodes[j], subGraphs[i]);
//       if (_debug)
// 	cout << subGraphNodes[j] << "    ";
//     }
//     if (_debug)
//       cout << endl;
//   }
//   reverse(subGraphs.begin(),  subGraphs.end());
//   reverse(startNodes.begin(), startNodes.end());
//   if (_debug)
//     for (unsigned int i = 0; i < subGraphs.size(); ++i) {
//       // printInfo << "created subgraph[" << i << "]:" << endl;
//       // print_graph (subGraphs[i], get(vertex_index, subGraphs[i]));
//       // print_edges2(subGraphs[i], get(vertex_index, subGraphs[i]), get(edge_index, subGraphs[i]));
//     }
//   // create all possible decays
//   map<nodeDesc, vector<decayGraph> > decayPossibilities;  // memorizes all possible decay graphs starting at the respective node
//   for (unsigned int i = 0; i < nmbVert; ++i) {
//     if (!_nodeVertexMap[startNodes[i]]) {
//       // final state node
//       decayPossibilities[startNodes[i]].push_back(subGraphs[i]);
//       continue;
//     }
//     // get daughter nodes
//     vector<nodeDesc> daughterNodes;
//     for (adjIterator iNode = adjacent_isobarVertices(startNodes[i], _graph).first;
// 	 iNode != adjacent_isobarVertices(startNodes[i], _graph).second; ++iNode)
//       daughterNodes.push_back(*iNode);
//     if (daughterNodes.size() != 2) {
//       printErr << "number of daughter vertices of subgraph top node is "
// 	       << daughterNodes.size() << " != 2. aborting." << endl;
//       throw;
//     }
//     // sort daughter nodes
//     {
//       isobarDecayVertex2* dv[2] = {
// 	static_cast<isobarDecayVertex2*>(_nodeVertexMap[daughterNodes[0]]),
// 	static_cast<isobarDecayVertex2*>(_nodeVertexMap[daughterNodes[1]])
//       };
//       isobarDecayVertex2* mv = static_cast<isobarDecayVertex2*>(_nodeVertexMap[startNodes[i]]);
//       unsigned int index = 0;
//       if (!dv[0] && dv[1])
// 	index = 1;
//       if (dv[0] || dv[1])
// 	if (&dv[index]->mother() != mv->outParticles()[index])
// 	  swap(daughterNodes[0], daughterNodes[1]);
//     }
//     // get daughter nodes
//     vector<decayGraph>& decays          = decayPossibilities[startNodes[i]];
//     vector<decayGraph>& daughter1Decays = decayPossibilities[daughterNodes[0]];
//     vector<decayGraph>& daughter2Decays = decayPossibilities[daughterNodes[1]];
//     for (unsigned int i1 = 0; i1 < daughter1Decays.size(); ++i1)
//       for(unsigned int i2 = 0; i2 < daughter2Decays.size(); ++i2) {
// 	//cout << "!!! " << i1 << ", " << i2 << endl;
// 	// join subgraphs
// 	//!!! if not taken care of by the calling code this is a potential memory leak
// 	nodeDesc   topNode;
// 	decayGraph joinedGraph = joinDaughterGraphs(*static_cast<isobarDecayVertex2*>(_nodeVertexMap[startNodes[i]]),
// 						    daughter1Decays[i1],
// 						    daughter2Decays[i2],
// 						    topNode);
// 	isobarDecayVertex2* vertex = static_cast<isobarDecayVertex2*>(get(vertex_vertexPointer, joinedGraph, topNode));
// 	//cout << "!!! " << topNode << ", " << *vertex << endl;
// 	particle&          d1     = vertex->daughter1();
// 	particle&          d2     = vertex->daughter2();
// 	//cout << "!!! " << d2 << endl;
// 	// create all possible quantum number combinations for this node
// 	// quantum numbers fixed by daughter quantum numbers
// 	const int baryonNmb   = d1.baryonNmb()   + d2.baryonNmb();
// 	const int charge      = d1.charge()      + d2.charge();
// 	const int strangeness = d1.strangeness() + d2.strangeness();
// 	const int charm       = d1.charm()       + d2.charm();
// 	const int beauty      = d1.beauty()      + d2.beauty();
// 	const int G           = d1.G()           * d2.G();
// 	//const int C           = d1.C()           * d2.C();
// 	// daughter quantum numbers that define ranges of possible mother quantum numbers
// 	const int s1 = d1.J();
// 	const int s2 = d2.J();
// 	const int I1 = d1.isospin();
// 	const int I2 = d2.isospin();

// 	const int  minS = 0;
// 	const int  maxS = 6;
// 	const int  minL = 0;
// 	const int  maxL = 6;
// 	const int  minJ = 0;

// 	const int  maxJ = 8;

// 	const int  minI = 0;
// 	const int  maxI = 2;
// 	const bool allowJpcExotic = false;
  
// 	for (int S = max(abs(s1 - s2), minS); S <= min(s1 + s2, maxS); S += 2) {        // loop over all allowed total spins
// 	  for (int L = max(0, minL); L <= maxL; L += 2) {                               // loop over all allowed relative orbital angular momenta
// 	    const int P = d1.P() * d2.P() * (L % 4 == 0 ? 1 : -1);  // parity
// 	    for (int J = max(abs(L - S), minJ); J <= min(L + S, maxJ); J += 2) {        // L-S coupling loop
// 	      for (int I = max(abs(I1 - I2), minI); I <= min(I1 + I2, maxI); I += 2) {  // isospin loop
// 		// check if charged state is allowed... 
// 		if(abs(charge-0.5*(baryonNmb+strangeness+charm+beauty))!=I*0.5)continue;
// 		int C = G * (I % 4 == 0 ? 1 : -1);  // C-Parity???
// 		if (_debug)
// 		  printInfo << vertex->mother().name() << "  --->  "
// 			    << d1.name() << "  +   " << d2.name() << "    "
// 			    << "2IG = "  << I << sign(G) << ", "
// 			    << "2JPC = " << J << sign(P) << sign(C) << ", "
// 			    << "2S = "   << S <<", 2L = " << L
// 			    << endl;
// 		if (!allowJpcExotic) {
// 		  // quark model boundary conditions:
// 		  // for P == C everything ok; check P = (-1)^(J + 1)
// 		  if ((P != C) && (   (C != (J % 4 == 0     ? 1 : -1))
// 				   || (P != (J + 2 % 4 == 0 ? 1 : -1)))) {
// 		    if (_debug)
// 		      printInfo << "disregarding spin-exotic isobar above" << endl;
// 		    continue;
// 		  }
// 		}
// 		//
// 		//Get fitting candidates from PDTable
		
// 		stringstream isobarname("isobar");
// 		isobarname << (abs(charge)>0 ? (string)sign(charge) : "");
// 		if(abs(charge)>1)isobarname << abs(charge);
// 		particleProperties prop(isobarname.str(),I,G,J,P,C);
// 		particleDataTable& pdt = particleDataTable::instance();
// 		string opt="IGJPC";
// 		vector<const particleProperties*> selection;
// 		if(vertex->mother().name()!="X-")selection=pdt.entrylist(prop,opt);
// 		else selection.push_back(&vertex->mother());
// 		unsigned int niso=selection.size();
// 		cout << "Found "<<niso<<" isobar candidates"<<endl;
// 		for(unsigned int iiso=0;iiso<niso;++iiso){
// 		  // check proposed mass for this isobar
// 		  // mother mass has to be larger than both daughters
// 		  if(selection[iiso]->mass()+selection[iiso]->width()<d1.mass()+d2.mass()){
// 		    cout << "Out of mass range m=" << selection[iiso]->mass() 
// 			 << "  md1="<<d1.mass()
// 			 << "  md2="<<d2.mass()<< endl << "Skipping candidate!" << endl;
// 		    continue;
// 		  }


// 		  //!!! if not taken care of by the calling code this is a potential memory leak
// 		  decayGraph graphCopy = deepCopyGraph(joinedGraph, false);
// 		  // set mother vertex quantum numbers
// 		  isobarDecayVertex2* v = static_cast<isobarDecayVertex2*>(get(vertex_vertexPointer, graphCopy, topNode));
// 		  v->setL(L);
// 		  v->setS(S);
// 		  // copy mother particle
// 		  particle* m = &v->mother().clone();
// 		  v->inParticles()[0] = m;
// 		  m->setName(selection[iiso]->name());
// 		  m->setMass(selection[iiso]->mass());
// 		  m->setWidth(selection[iiso]->width());
// 		  m->setCharge     (charge);
// 		  m->setBaryonNmb  (baryonNmb);
// 		  m->setIsospin    (I);
// 		  m->setStrangeness(strangeness);
// 		  m->setCharm      (charm);
// 		  m->setBeauty     (beauty);
// 		  m->setG          (G);
// 		  m->setJ          (J);
// 		  m->setP          (P);
// 		  m->setC          (C);
// 		  //cout << "!!! " << topNode << ", " << *v << endl;
// 		  decays.push_back(graphCopy);
// 		} // isobar candidate
// 	      }  // isospin loop
// 	    }  // L-S coupling loop
// 	  }  // L loop
// 	}  // S loop
//       }  // loop over daughter decays
//   }  // loop over all start nodes
//   // construct isobar decay topologies for all decays
//   vector<decayGraph>&          decays = decayPossibilities[startNodes[nmbVert - 1]];
//   vector<isobarDecayTopology2*> decayTopos;
//   for (unsigned int i = 0; i < decays.size(); ++i) {
//     // get decay vertices
//     vector<isobarDecayVertex2*> decayVertices;
//     for (nodeIterator iNode = vertices(decays[i]).first;
// 	 iNode != vertices(decays[i]).second; ++iNode) {
//       isobarDecayVertex2* v;
//       v = static_cast<isobarDecayVertex2*>(get(vertex_vertexPointer, decays[i], *iNode));
//       if (v)
// 	decayVertices.push_back(v);
//     }
//     // construct production vertex
//     particle& beam = static_cast<diffractiveDissVertex*>(&productionVertex())->beam();
//     particle& X    = static_cast<isobarDecayVertex2*>(get(vertex_vertexPointer, decays[i], 0))->mother();
//     diffractiveDissVertex* prodVert = new diffractiveDissVertex(beam, X);
//     // construct isobar decay topology
//     decayTopos.push_back(new isobarDecayTopology2(fsParticles(), *prodVert, decayVertices));
//   }
//   return decayTopos;
// }


// decayTopology::decayGraph
// isobarDecayTopology2::joinDaughterGraphs(const isobarDecayVertex2& motherVertex,
// 					const decayGraph&        daughterGraph1,
// 					const decayGraph&        daughterGraph2,
// 					nodeDesc&                newMotherNode)
// {
//   //cout << "!!! old mother vertex: " << motherVertex << endl;
//   decayGraph newGraph;
//   newMotherNode = add_vertex(newGraph);
//   copy_graph(daughterGraph1, newGraph);
//   copy_graph(daughterGraph2, newGraph);
//   nodeDesc daughterNodes[2] = {1, 1 + num_vertices(daughterGraph1)};
//   // copy isobar decay vertex at topNode
//   isobarDecayVertex2* newMotherVertex = &motherVertex.clone();
//   // connect subgraphs to new top node
//   for (unsigned int i = 0; i < 2; ++i) {
//     bool     inserted;
//     edgeDesc edge;
//     tie(edge, inserted) = add_edge(newMotherNode, daughterNodes[i], newGraph);
//     isobarDecayVertex2* daughterVertex;
//     daughterVertex = static_cast<isobarDecayVertex2*>(get(vertex_vertexPointer, newGraph, daughterNodes[i]));
//     if (!daughterVertex) {
//       // final state vertex
//       put(edge_particlePointer, newGraph, edge, motherVertex.outParticles()[i]);
//     } else {
//       // isobar decay vertex
//       put(edge_particlePointer, newGraph, edge, &daughterVertex->mother());
//       newMotherVertex->outParticles()[i] = &daughterVertex->mother();
//       //cout << "!!! daughter vertex " << i << ": " << *daughterVertex << endl;
//     }
//   }
//   put(vertex_vertexPointer, newGraph, newMotherNode, newMotherVertex);
//   //cout << "!!! new mother vertex: " << *newMotherVertex << endl;
//   return newGraph;
// }


const TLorentzVector&
isobarDecayTopology2::calcIsobarLzVec()
{
  // loop over isobar decay vertices and propagate changes from final state particles up to X-system
  for (int i = nmbInteractionVertices() - 1; i >= 0; --i) {
    if (_debug)
      printInfo << "calculating Lorentz-vector of mother isobar '"
		<< _isobarVertices[i]->mother()->name() << "' "
		<< "of node[" << node(_isobarVertices[i]) << "]" << endl;
   _isobarVertices[i]->calcMotherLzVec();
  }
  return xIsobarDecayVertex()->mother()->lzVec();
}


ostream&
isobarDecayTopology2::print(ostream& out) const
{
  out << "isobar decay topology '" << name() << "' has " << nmbNodes() << " nodes:" << endl;
  out << "    production   node[" << node(productionVertex()) << "] = "
      << *productionVertex() << endl;
  for (unsigned int i = 0; i < nmbInteractionVertices(); ++i)
    out << "    isobar decay node[" << node(interactionVertices()[i]) << "] = "
	<< *interactionVertices()[i] << endl;
  out << "isobar decay topology '" << name() << "' has " << nmbEdges() << " edges:" << endl;
  edgeIterator iEdge, iEdgeEnd;
  for (tie(iEdge, iEdgeEnd) = edges(); iEdge != iEdgeEnd; ++iEdge) {
    const particlePtr part = particle(*iEdge);
    out << "    edge[" << *iEdge << "] = [" << fromNode(*iEdge) << ", "
	<< toNode(*iEdge) << "] = '" << part->name() << "'";
    if (isFsParticle(part))
      out << " (final state)";
    out << endl;
  }
  return out;
}
