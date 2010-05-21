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
#include "diffractiveDissVertex.h"
#include "isobarDecayTopology.h"
#include "particleDataTable.h"

	
using namespace std;
using namespace boost;
using namespace rpwa;


bool isobarDecayTopology::_debug = false;


isobarDecayTopology::isobarDecayTopology()
  : decayTopology()
{ }


isobarDecayTopology::isobarDecayTopology(const interactionVertexPtr&         productionVertex,
					 const vector<isobarDecayVertexPtr>& isobarDecayVertices,
					 const vector<particlePtr>&          fsParticles)
  : decayTopology()
{
  constructDecay(productionVertex, isobarDecayVertices, fsParticles);
}


isobarDecayTopology::isobarDecayTopology(const interactionVertexPtr&         productionVertex,
					 const vector<interactionVertexPtr>& isobarDecayVertices,
					 const vector<particlePtr>&          fsParticles)
  : decayTopology()
{
  constructDecay(productionVertex, isobarDecayVertices, fsParticles);
}


isobarDecayTopology::isobarDecayTopology(const isobarDecayTopology& topo)
  : decayTopology()
{
  *this = topo;
}


isobarDecayTopology::isobarDecayTopology(const decayTopology& topo)
  : decayTopology()
{
  *this = topo;
}


isobarDecayTopology::~isobarDecayTopology()
{ }


isobarDecayTopology&
isobarDecayTopology::operator =(const isobarDecayTopology& topo)
{
  if (this != &topo) {
    decayTopology::operator =(topo);
    _isobarVertices = topo._isobarVertices;
  }
  return *this;
}


isobarDecayTopology&
isobarDecayTopology::operator =(const decayTopology& topo)
{
  if (this != &topo) {
    decayTopology::operator =(topo);
    buildIsobarVertexArray();
  }
  return *this;
}


isobarDecayTopology*
isobarDecayTopology::clone(const bool cloneFsParticles,
			   const bool cloneProductionVertex) const
{
  if (_debug)
    printInfo << "cloning isobar decay topology '" << name() << "'; "
	      << "cloneFsParticles = "      << cloneFsParticles << ", "
	      << "cloneProductionVertex = " << cloneProductionVertex << endl;
  decayTopology         topoClone       = *decayTopology::clone(cloneFsParticles, cloneProductionVertex);
  isobarDecayTopology* isobarTopoClone = new isobarDecayTopology(topoClone);
  isobarTopoClone->buildIsobarVertexArray();
  return isobarTopoClone;
}


void
isobarDecayTopology::clear()
{
  decayTopology::clear();
  _isobarVertices.clear();
}


isobarDecayTopology&
isobarDecayTopology::constructDecay(const interactionVertexPtr&         productionVertex,
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
  decayTopology::constructDecay(productionVertex, intVertices, fsParticles);
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


isobarDecayTopology&
isobarDecayTopology::constructDecay(const interactionVertexPtr&         productionVertex,
				    const vector<interactionVertexPtr>& isobarDecayVertices,
				    const vector<particlePtr>&          fsParticles)
{
  const unsigned int nmbVert = isobarDecayVertices.size();
  vector<isobarDecayVertexPtr> decayVertices(nmbVert);
  for (unsigned int i = 0; i < nmbVert; ++i) {
    decayVertices[i] = dynamic_pointer_cast<isobarDecayVertex>(isobarDecayVertices[i]);
    if (!decayVertices[i]) {
      printErr << "interaction vertex[" << i << "] is not an isobarDecayVertex. aborting." << endl;
      throw;
    }
  }
  return constructDecay(productionVertex, decayVertices, fsParticles);
}


bool
isobarDecayTopology::checkTopology() const
{
  // perform basic checks of topology
  bool topologyIsOkay = decayTopology::checkTopology();
  // check that decay topology is a tree of isobar decays
  for (unsigned int i = 0; i < nmbInteractionVertices(); ++i) {
    const nodeDesc& nd = node(_isobarVertices[i]);
    // check that each isobar decay node has exactly 1 incoming edge (isobar)
    const unsigned int nmbIn = nmbInEdges(_isobarVertices[i]);
    if (nmbIn != 1) {
      printWarn << "number of incoming edges of node[" << nd << "] is "
		<< nmbIn << " != 1" << endl;
      topologyIsOkay = false;
    } else if (_debug)
      printInfo << "number of incoming edges of node[" << nd << "] = " << nmbIn
		<< " is correct" << endl;
    // check that for each isobar decay node the number of outgoing edges is 2
    const unsigned int nmbOut = nmbOutEdges(_isobarVertices[i]);
    if (nmbOut != 2) {
      printWarn << "number of outgoing edges of node[" << nd << "] is "
    		<< nmbOut << " != 2" << endl;
      topologyIsOkay = false;
    } else if (_debug)
      printInfo << "number of outgoing edges of node[" << nd << "] = " << nmbOut
		<< " is correct" << endl;
  }
  if (_debug)
    printInfo << "isobar decay topology " << ((topologyIsOkay) ? "passed" : "did not pass")
	      << " all tests" << endl;
  return topologyIsOkay;
}


bool 
isobarDecayTopology::checkConsistency() const
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


vector<isobarDecayTopology>
isobarDecayTopology::possibleDecays(const int  minI,
				    const int  maxI,
				    const int  minJ,
				    const int  maxJ,
				    const int  minL,
				    const int  maxL,
				    const int  minS,
				    const int  maxS,
				    const bool allowJpcExotic)
{
  if (!checkTopology()) {
    printErr << "ill-formed graph topology. aborting." << endl;
    throw;
  }
  // order nodes depth-first
  vector<nodeDesc> startNds = sortNodesDfs(xIsobarDecayVertex());
  // create decay topologies of all subdecays
  vector<isobarDecayTopology> subDecays(startNds.size());
  for (unsigned int i = 0; i < startNds.size(); ++i)
    subDecays[i] = subDecay(startNds[i]);
  // reverse order, because decay trees are build starting from the final state nodes
  reverse(startNds.begin(), startNds.end());
  reverse(subDecays.begin (), subDecays.end ());
  if (_debug)
    for (unsigned int i = 0; i < subDecays.size(); ++i)
      printInfo << "created subdecay[" << i << "]: " << subDecays[i];

  // create all possible subdecays
  map<nodeDesc, vector<isobarDecayTopology> > decayPossibilities;  // all possible decay graphs starting at certain node
  for (unsigned int iStart = 0; iStart < startNds.size(); ++iStart) {
    
    // special case for final state nodes
    if (isFsVertex(vertex(startNds[iStart]))) {
      // final state vertices have no daughter tracks; subdecay
      // topology consist of just the final state vertex
      decayPossibilities[startNds[iStart]].push_back(subDecays[iStart]);
      continue;
    }

    // get daughter decay topologies
    vector<vector<isobarDecayTopology>* > daughterDecays;
    adjIterator iNd, iNdEnd;
    for (tie(iNd, iNdEnd) = adjacentVertices(startNds[iStart]); iNd != iNdEnd; ++iNd)
      daughterDecays.push_back(&decayPossibilities[*iNd]);
    if (daughterDecays.size() != 2) {
      printErr << "node[" << startNds[iStart] << "]: " << *vertex(startNds[iStart]) << " has "
	       << daughterDecays.size() << " daughters. Exactly two are required."
	       << "aborting." << endl;
      throw;
    }

    // loop over all combinations of daughter decays
    unsigned int iDaughter[2];
    for (iDaughter[0] = 0; iDaughter[0] < daughterDecays[0]->size(); ++iDaughter[0])
      for(iDaughter[1] = 0; iDaughter[1] < daughterDecays[1]->size(); ++iDaughter[1]) {

	// copy mother vertex
	isobarDecayVertexPtr motherVertex(new isobarDecayVertex(*static_pointer_cast<isobarDecayVertex>(vertex(startNds[iStart]))));
	particlePtr          mother = motherVertex->mother();
	// get daughter particles from the respective decay topologies
	const nodeDesc topNodes[2] = 
	  {(*daughterDecays[0])[iDaughter[0]].topNode(),
	   (*daughterDecays[1])[iDaughter[1]].topNode()};
	const particlePtr daughters[2] = 
	  {(*daughterDecays[0])[iDaughter[0]].vertex(topNodes[0])->inParticles()[0],
	   (*daughterDecays[1])[iDaughter[1]].vertex(topNodes[1])->inParticles()[0]};
	// set daughters of mother vertex
	motherVertex->daughter1() = daughters[0];
	motherVertex->daughter2() = daughters[1];
  	// join daughter subdecays and mother vertex
  	isobarDecayTopology motherDecay = joinDaughterDecays(motherVertex,
							      (*daughterDecays[0])[iDaughter[0]],
							      (*daughterDecays[1])[iDaughter[1]]);
	
	// calculate mother quantum numbers fixed by daughter quantum numbers
	const int baryonNmb   = daughters[0]->baryonNmb()   + daughters[1]->baryonNmb();
	const int charge      = daughters[0]->charge()      + daughters[1]->charge();
	const int strangeness = daughters[0]->strangeness() + daughters[1]->strangeness();
	const int charm       = daughters[0]->charm()       + daughters[1]->charm();
	const int beauty      = daughters[0]->beauty()      + daughters[1]->beauty();
	const int G           = daughters[0]->G()           * daughters[1]->G();
	// daughter quantum numbers that define ranges of possible mother quantum numbers
	const int spins   [2] = {daughters[0]->J(),       daughters[1]->J()      };
	const int isospins[2] = {daughters[0]->isospin(), daughters[1]->isospin()};

	// loop over all allowed combinations of daughter quantum numbers
	// loop over allowed total spins
	for (int S = max(abs(spins[0] - spins[1]), minS);
	     S <= min(spins[0] + spins[1], maxS); S += 2) {
	  // loop over allowed relative orbital angular momenta
  	  for (int L = max(0, minL); L <= maxL; L += 2) {
  	    const int P = daughters[0]->P() * daughters[1]->P() * (L % 4 == 0 ? 1 : -1);  // parity
	    // loop over allowed total angular momenta
  	    for (int J = max(abs(L - S), minJ); J <= min(L + S, maxJ); J += 2) {
	      // loop over allowed isospins
  	      for (int I = max(abs(isospins[0] - isospins[1]), minI);
		   I <= min(isospins[0] + isospins[1], maxI); I += 2) {
  		// check whether charge state is allowed
		// !!! the Gell-Mann-Nishijima formula cannot be used
		//     here, because it employs the the z-component of
		//     the isospin, I_z (not I); see PDG 2008 eq. 14.1
  		// if (abs(charge - 0.5 * (baryonNmb + strangeness + charm + beauty)) != 0.5 * I)
		//   continue;
  		const int C = G * (I % 4 == 0 ? 1 : -1);  // C-parity
  		if (!allowJpcExotic)
  		  // quark model restrictions: P == C is always allowed
		  // check that P = (-1)^(J + 1)
  		  if ((P != C) && (   (C != (J % 4     == 0 ? 1 : -1))
				   || (P != (J + 2 % 4 == 0 ? 1 : -1)))) {
  		    if (_debug)
  		      printInfo << "disregarding spin-exotic isobar with IG(JPC) = "
				<< 0.5 * I << sign(G)
				<< "(" << 0.5 * J << sign(P) << sign(C) << ")" << endl;
  		    continue;
  		  }
		
  		// find candidates for mother isobar in particle data table
  		particleProperties isobarProp(*mother);
		isobarProp.setBaryonNmb(baryonNmb);
		isobarProp.setSCB      (strangeness, charm, beauty);
		isobarProp.setIGJPC    (I, G, J, P, C);
  		particleDataTable& pdt = particleDataTable::instance();
  		vector<const particleProperties*> isobarCandidates;
		const string motherName    = mother->name();
		const double minIsobarMass = daughters[0]->mass() + daughters[1]->mass();
  		if (   (motherName == "X-") || (motherName == "X0")
		    || (motherName == "X+") || (motherName == "X")) {
		  if (mother->mass() + mother->width() < minIsobarMass) {
		    cout << isobarProp.qnSummary() << " mass too low. skipping candidate." << endl;
		    continue;
		  }
		  isobarCandidates.push_back(&isobarProp);
  		} else {
		  isobarCandidates = pdt.entriesMatching(isobarProp, "allQn", minIsobarMass);
		  if (_debug)
		    printInfo << "found " << isobarCandidates.size() << " isobar candidates for "
			      << isobarProp.qnSummary() << " in particle data table" << endl;
		}

		// loop over isobar candidates
  		for (unsigned int iIsobar = 0; iIsobar < isobarCandidates.size(); ++iIsobar) {
		  // clone topology
  		  isobarDecayTopology* decayCopy = motherDecay.clone(false, false);
  		  // set mother vertex quantum numbers
  		  const isobarDecayVertexPtr& motherVertexCopy =
		    static_pointer_cast<isobarDecayVertex>(decayCopy->vertex(motherDecay.node(motherVertex)));
  		  motherVertexCopy->setL(L);
  		  motherVertexCopy->setS(S);
		  // set mother particle
  		  const particlePtr& motherCopy = motherVertexCopy->mother();
		  motherCopy->setProperties(*isobarCandidates[iIsobar]);
  		  motherCopy->setCharge(charge);
		  if (_debug)
		    printInfo << "created decay topology for " << *motherVertexCopy << endl;
  		  decayPossibilities[startNds[iStart]].push_back(*decayCopy);
   		} // isobar candidate loop
  	      }  // isospin loop
  	    }  // L-S coupling loop
  	  }  // L loop
  	}  // S loop
      }  // loop over daughter decays
  }  // loop over all start nodes

  // extract decays for X-decay vertex and add production vertex
  vector<isobarDecayTopology> decays = decayPossibilities[startNds.back()];
  for (unsigned int i = 0; i < decays.size(); ++i) {
    // clone production vertex and set X-particle
    const interactionVertexPtr newProdVert(productionVertex()->clone(false, false));
    const particlePtr&         newX = decays[i].xIsobarDecayVertex()->mother();
    newProdVert->outParticles()[0] = newX;
    // add production vertex
    decays[i].setProductionVertex(newProdVert);
  }
  return decays;
}


const TLorentzVector&
isobarDecayTopology::calcIsobarLzVec()
{
  // loop over isobar decay vertices and propagate Lorentz-vectors from final state particles up to X-system
  for (int i = nmbInteractionVertices() - 1; i >= 0; --i) {
    if (_debug)
      printInfo << "calculating Lorentz-vector of mother isobar '"
		<< _isobarVertices[i]->mother()->name() << "' "
		<< "of node[" << node(_isobarVertices[i]) << "]" << endl;
    _isobarVertices[i]->calcMotherLzVec();
  }
  return xIsobarDecayVertex()->mother()->lzVec();
}


void
isobarDecayTopology::calcIsobarCharges()
{
  // loop over isobar decay vertices and propagate charges from final state particles up to X-system
  for (int i = nmbInteractionVertices() - 1; i >= 0; --i) {
    if (_debug)
      printInfo << "calculating charge of mother isobar '"
		<< _isobarVertices[i]->mother()->name() << "' "
		<< "of node[" << node(_isobarVertices[i]) << "]" << endl;
    _isobarVertices[i]->calcMotherCharge();
  }
  // update graph name
  name() = productionVertex()->outParticles()[0]->qnSummary();
}


ostream&
isobarDecayTopology::print(ostream& out) const
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
  nodeIterator iNd, iNdEnd;
  for (tie(iNd, iNdEnd) = nodes(); iNd != iNdEnd; ++iNd)
    if (isFsVertex(vertex(*iNd)))
      out << "    final state node[" << *iNd << "] = " << *vertex(*iNd) << endl;
  // print edges
  out << "isobar decay topology '" << name() << "' has " << nmbEdges() << " edge(s):" << endl;
  edgeIterator iEd, iEdEnd;
  for (tie(iEd, iEdEnd) = edges(); iEd != iEdEnd; ++iEd) {
    const particlePtr part = particle(*iEd);
    out << "    edge[" << *iEd << "] = " << part->name();
    const int fsIndex = fsParticleIndex(part);
    if (fsIndex >= 0)
      out << "[" << fsIndex << "] (final state)";
    out << endl;
  }
  return out;
}


void
isobarDecayTopology::buildIsobarVertexArray()
{
  bool success = true;
  _isobarVertices.resize(nmbInteractionVertices());
  for (unsigned int i = 0; i < nmbInteractionVertices(); ++i) {
    const interactionVertexPtr& v = interactionVertices()[i];
    _isobarVertices[i] = dynamic_pointer_cast<isobarDecayVertex>(v);
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


isobarDecayTopology
isobarDecayTopology::subDecay(const nodeDesc& startNd)
{
  isobarDecayTopology subTopo(decayTopology::subDecay(startNd));
  subTopo.name() = vertex(startNd)->inParticles()[0]->qnSummary();
  return subTopo;
}


void
isobarDecayTopology::addDecay(const isobarDecayTopology& topo)
{
  decayTopology::addDecay(topo);
  buildIsobarVertexArray();
}


isobarDecayTopology
isobarDecayTopology::joinDaughterDecays(const isobarDecayVertexPtr&        motherVertex,
					const vector<isobarDecayTopology>& daughterDecays)  ///< joins daughter decay graphs and connects them to a common mother vertex
{
  if (_debug) {
    printInfo << "joining " << daughterDecays.size() << " daughter graphs with mother "
	 << *motherVertex << endl;
  }
  isobarDecayTopology newTopo;
  newTopo.addVertex(motherVertex);
  for (unsigned int i = 0; i < daughterDecays.size(); ++i)
    newTopo.addDecay(daughterDecays[i]);
  return newTopo;
}


isobarDecayTopology
isobarDecayTopology::joinDaughterDecays(const isobarDecayVertexPtr& motherVertex,
					const isobarDecayTopology&  daughter1Decay,
					const isobarDecayTopology&  daughter2Decay)  ///< joins daughter decay graphs and connects them to a common mother vertex
{
  std::vector<isobarDecayTopology> daughterDecays(2);
  daughterDecays[0] = daughter1Decay;
  daughterDecays[1] = daughter2Decay;
  return joinDaughterDecays(motherVertex, daughterDecays);
}
