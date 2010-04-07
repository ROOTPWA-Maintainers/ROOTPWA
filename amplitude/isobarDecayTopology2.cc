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


isobarDecayTopology2::isobarDecayTopology2(const isobarDecayTopology2& topo)
  : decayTopology2()
{
  *this = topo;
}


isobarDecayTopology2::isobarDecayTopology2(const decayTopology2& topo)
  : decayTopology2()
{
  *this = topo;
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


isobarDecayTopology2&
isobarDecayTopology2::operator =(const decayTopology2& topo)
{
  if (this != &topo) {
    decayTopology2::operator =(topo);
    buildIsobarVertexArray();
  }
  return *this;
}


isobarDecayTopology2*
isobarDecayTopology2::clone(const bool cloneFsParticles,
			    const bool cloneProductionVertex) const
{
  decayTopology2        topoClone       = *decayTopology2::clone(cloneFsParticles, cloneProductionVertex);
  isobarDecayTopology2* isobarTopoClone = new isobarDecayTopology2(topoClone);
  isobarTopoClone->buildIsobarVertexArray();
  return isobarTopoClone;
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
  buildIsobarVertexArray();
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
  bool topologyIsOkay = decayTopology2::checkTopology();
  // check that decay topology is a tree of isobar decays
  for (unsigned int i = 0; i < nmbInteractionVertices(); ++i) {
    const int index = node(_isobarVertices[i]);
    // check that each isobar decay node has exactly 1 incoming edge (isobar)
    const unsigned int nmbIn = nmbInEdges(_isobarVertices[i]);
    if (nmbIn != 1) {
      printWarn << "number of incoming edges of node[" << index << "] is "
		<< nmbIn << " != 1" << endl;
      topologyIsOkay = false;
    } else if (_debug)
      printInfo << "number of incoming edges of node[" << index << "] = " << nmbIn
		<< " is correct" << endl;
    // check that for each isobar decay node the number of outgoing edges is 2
    const unsigned int nmbOut = nmbOutEdges(_isobarVertices[i]);
    if (nmbOut != 2) {
      printWarn << "number of outgoing edges of node[" << index << "] is "
    		<< nmbOut << " != 2" << endl;
      topologyIsOkay = false;
    } else if (_debug)
      printInfo << "number of outgoing edges of node[" << index << "] = " << nmbOut
		<< " is correct" << endl;
  }
  if (_debug)
    printInfo << "isobar decay topology " << ((topologyIsOkay) ? "passed" : "did not pass")
	      << " all tests" << endl;
  return topologyIsOkay;
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
  if (_debug)
    printInfo << "information in isobar decay vertices is " << ((allVertConsistent) ? "" : "NOT ")
	      << "consistent" << endl;
  return allVertConsistent;
}


vector<isobarDecayTopology2>
isobarDecayTopology2::possibleDecays()
{
  if (!checkTopology()) {
    printErr << "wrong graph topology. aborting." << endl;
    throw;
  }
  if (!checkConsistency()) {
    printErr << "isobar decay vertex data are not consistent. aborting." << endl;
    throw;
  }
  // order nodes depth-first
  vector<nodeDesc> startNodes = sortNodesDfs(xIsobarDecayVertex());
  // create subgraphs for all subdecays
  vector<isobarDecayTopology2> subDecays(startNodes.size());
  for (unsigned int i = 0; i < startNodes.size(); ++i)
    subDecays[i] = subDecay(startNodes[i]);
  reverse(startNodes.begin(), startNodes.end());
  reverse(subDecays.begin(),  subDecays.end());
  if (_debug)
    for (unsigned int i = 0; i < subDecays.size(); ++i)
      printInfo << "created subdecay[" << i << "]: " << subDecays[i] << endl;
  // // create all possible subdecays
  // map<nodeDesc, vector<decayGraph> > decayPossibilities;  // possible decay graphs starting at node
  // for (unsigned int i = 0; i < startNodes.size(); ++i) {

  //   if (!_nodeVertexMap[startNodes[i]]) {
  //     // final state node
  //     decayPossibilities[startNodes[i]].push_back(subDecays[i]);
  //     continue;
  //   }

  //   // get daughter nodes
  //   vector<nodeDesc> daughterNodes;
  //   for (adjIterator iNode = adjacent_isobarVertices(startNodes[i], _graph).first;
  // 	 iNode != adjacent_isobarVertices(startNodes[i], _graph).second; ++iNode)
  //     daughterNodes.push_back(*iNode);

  //   // sort daughter nodes
  //   {
  //     isobarDecayVertex2* dv[2] = {
  // 	static_cast<isobarDecayVertex2*>(_nodeVertexMap[daughterNodes[0]]),
  // 	static_cast<isobarDecayVertex2*>(_nodeVertexMap[daughterNodes[1]])
  //     };
  //     isobarDecayVertex2* mv = static_cast<isobarDecayVertex2*>(_nodeVertexMap[startNodes[i]]);
  //     unsigned int index = 0;
  //     if (!dv[0] && dv[1])
  // 	index = 1;
  //     if (dv[0] || dv[1])
  // 	if (&dv[index]->mother() != mv->outParticles()[index])
  // 	  swap(daughterNodes[0], daughterNodes[1]);
  //   }

  //   // get daughter nodes
  //   vector<decayGraph>& decays          = decayPossibilities[startNodes[i]];
  //   vector<decayGraph>& daughter1Decays = decayPossibilities[daughterNodes[0]];
  //   vector<decayGraph>& daughter2Decays = decayPossibilities[daughterNodes[1]];
  //   for (unsigned int i1 = 0; i1 < daughter1Decays.size(); ++i1)
  //     for(unsigned int i2 = 0; i2 < daughter2Decays.size(); ++i2) {
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
  //     }  // loop over daughter decays
  // }  // loop over all start nodes

  // // construct isobar decay topologies for all decays
  // vector<decayGraph>&          decays = decayPossibilities[startNodes[nmbVert - 1]];
  vector<isobarDecayTopology2> decayTopos;
  // for (unsigned int i = 0; i < decays.size(); ++i) {
  //   // get decay vertices
  //   vector<isobarDecayVertex2*> decayVertices;
  //   for (nodeIterator iNode = vertices(decays[i]).first;
  // 	 iNode != vertices(decays[i]).second; ++iNode) {
  //     isobarDecayVertex2* v;
  //     v = static_cast<isobarDecayVertex2*>(get(vertex_vertexPointer, decays[i], *iNode));
  //     if (v)
  // 	decayVertices.push_back(v);
  //   }
  //   // construct production vertex
  //   particle& beam = static_cast<diffractiveDissVertex*>(&productionVertex())->beam();
  //   particle& X    = static_cast<isobarDecayVertex2*>(get(vertex_vertexPointer, decays[i], 0))->mother();
  //   diffractiveDissVertex* prodVert = new diffractiveDissVertex(beam, X);
  //   // construct isobar decay topology
  //   decayTopos.push_back(new isobarDecayTopology2(fsParticles(), *prodVert, decayVertices));
  // }
  return decayTopos;
}


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
  // print nodes
  out << "isobar decay topology '" << name() << "' has " << nmbNodes() << " node(s):" << endl;
  if (productionVertex())
    out << "    production   node[" << node(productionVertex()) << "] = "
	<< *productionVertex() << endl;
  else
    out << "    topology has no production node." << endl;
  for (unsigned int i = 0; i < nmbInteractionVertices(); ++i)
    out << "    isobar decay node[" << node(interactionVertices()[i]) << "] = "
	<< *interactionVertices()[i] << endl;
  nodeIterator iNode, iNodeEnd;
  for (tie(iNode, iNodeEnd) = nodes(); iNode != iNodeEnd; ++iNode)
    if (isFsVertex(vertex(*iNode)))
      out << "    final state node[" << *iNode << "] = " << *vertex(*iNode) << endl;
  // print edges
  out << "isobar decay topology '" << name() << "' has " << nmbEdges() << " edge(s):" << endl;
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


void
isobarDecayTopology2::buildIsobarVertexArray()
{
  bool success = true;
  _isobarVertices.resize(nmbInteractionVertices());
  for (unsigned int i = 0; i < nmbInteractionVertices(); ++i) {
    const interactionVertexPtr& v = interactionVertices()[i];
    _isobarVertices[i] = dynamic_pointer_cast<isobarDecayVertex2>(v);
    if (!_isobarVertices[i]) {
      printWarn << *v << " is not of type isobarDecayVertex." << endl;
      success = false;
    }
  }
  if (!success) {
    printErr << "incompatible topology to copy from. some interaction vertices are "
	     << "not of type isobarDecayVertex. aborting." << endl;
    throw;
  }
}


isobarDecayTopology2
isobarDecayTopology2::subDecay(const nodeDesc& startNd)
{
  isobarDecayTopology2 subTopo(decayTopology2::subDecay(startNd));
  return subTopo;
}


void
isobarDecayTopology2::addDecay(const isobarDecayTopology2& topo)
{
  decayTopology2::addDecay(topo);
  buildIsobarVertexArray();
}


isobarDecayTopology2
isobarDecayTopology2::joinDaughterGraphs(const isobarDecayVertexPtr&         motherVertex,
					 const vector<isobarDecayTopology2>& daughterDecays)  ///< joins daughter decay graphs and connects them to a common mother vertex
{
  if (_debug) {
    printInfo << "joining " << daughterDecays.size() << " daughter graphs with mother vertex "
	 << *motherVertex << endl;
  }
  isobarDecayTopology2 newTopo;
  newTopo.addVertex(motherVertex);
  for (unsigned int i = 0; i < daughterDecays.size(); ++i)
    newTopo.addDecay(daughterDecays[i]);
  return newTopo;
}


isobarDecayTopology2
isobarDecayTopology2::joinDaughterGraphs(const isobarDecayVertexPtr& motherVertex,
					 const isobarDecayTopology2& daughter1Decay,
					 const isobarDecayTopology2& daughter2Decay)  ///< joins daughter decay graphs and connects them to a common mother vertex
{
  std::vector<isobarDecayTopology2> daughterDecays(2);
  daughterDecays[0] = daughter1Decay;
  daughterDecays[1] = daughter2Decay;
  return joinDaughterGraphs(motherVertex, daughterDecays);
}
