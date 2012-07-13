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


isobarDecayTopology::isobarDecayTopology(const productionVertexPtr&          productionVertex,
                                         const vector<isobarDecayVertexPtr>& isobarDecayVertices,
                                         const vector<particlePtr>&          fsParticles)
	: decayTopology()
{
	constructDecay(productionVertex, isobarDecayVertices, fsParticles);
}


isobarDecayTopology::isobarDecayTopology(const productionVertexPtr&          productionVertex,
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
isobarDecayTopology::doClone(const bool cloneFsParticles,
                             const bool cloneProdKinematics) const
{
	if (_debug)
		printDebug << "cloning isobar decay topology '" << name() << "'; "
		           << ((cloneFsParticles   ) ? "in" : "ex") << "cluding final state particles, "
		           << ((cloneProdKinematics) ? "in" : "ex") << "cluding production kinematics particles"
		           << endl;
	decayTopology*       topoClone       = decayTopology::doClone(cloneFsParticles, cloneProdKinematics);
	isobarDecayTopology* isobarTopoClone = new isobarDecayTopology(*topoClone);
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
isobarDecayTopology::constructDecay(const productionVertexPtr&          productionVertex,
                                    const vector<isobarDecayVertexPtr>& isobarDecayVertices,
                                    const vector<particlePtr>&          fsParticles)
{
	const unsigned int nmbVert = isobarDecayVertices.size();
	if (_debug)
		printDebug << "constructing isobar decay topology with "
		           << fsParticles.size() << " final state particles and "
		           << nmbVert            << " isobar decay vertices" << endl;
	vector<interactionVertexPtr> intVertices(nmbVert);
	for (unsigned int i = 0; i < nmbVert; ++i)
		intVertices[i] = isobarDecayVertices[i];
	decayTopology::constructDecay(productionVertex, intVertices, fsParticles);
	// copy sorted vertices
	buildIsobarVertexArray();
	if (nmbDecayVertices() != nmbVert) {
		printErr << "number of interaction vertices  = " << nmbDecayVertices()
		         << " does not match number of vertices given in parameter array = " << nmbVert
		         << ". aborting." << endl;
		throw;
	}
	return *this;
}


isobarDecayTopology&
isobarDecayTopology::constructDecay(const productionVertexPtr&          productionVertex,
                                    const vector<interactionVertexPtr>& isobarDecayVertices,
                                    const vector<particlePtr>&          fsParticles)
{
	const unsigned int nmbVert = isobarDecayVertices.size();
	vector<isobarDecayVertexPtr> decayVertices(nmbVert);
	for (unsigned int i = 0; i < nmbVert; ++i) {
		decayVertices[i] = dynamic_pointer_cast<isobarDecayVertex>(isobarDecayVertices[i]);
		if (not decayVertices[i]) {
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
	for (unsigned int i = 0; i < nmbDecayVertices(); ++i) {
		const nodeDesc& nd = node(_isobarVertices[i]);
		// check that each isobar decay node has exactly 1 incoming edge (isobar)
		const unsigned int nmbIn = nmbInEdges(_isobarVertices[i]);
		if (nmbIn != 1) {
			printWarn << "number of incoming edges of node[" << nd << "] is "
			          << nmbIn << " != 1" << endl;
			topologyIsOkay = false;
		} else if (_debug)
			printDebug << "number of incoming edges of node[" << nd << "] = " << nmbIn
			           << " is correct" << endl;
		// check that for each isobar decay node the number of outgoing edges is 2
		const unsigned int nmbOut = nmbOutEdges(_isobarVertices[i]);
		if (nmbOut != 2) {
			printWarn << "number of outgoing edges of node[" << nd << "] is "
			          << nmbOut << " != 2" << endl;
			topologyIsOkay = false;
		} else if (_debug)
			printDebug << "number of outgoing edges of node[" << nd << "] = " << nmbOut
			           << " is correct" << endl;
	}
	if (_debug) {
		if (topologyIsOkay)
			printDebug << "isobar decay topology passed all tests" << endl;
		else
			printWarn << "isobar decay topology did not pass all tests" << endl;
	}
	return topologyIsOkay;
}


bool 
isobarDecayTopology::checkConsistency() const
{
	bool allVertConsistent = true;
	for (unsigned int i = 0; i < _isobarVertices.size(); ++i) {
		const bool vertexConsistent = _isobarVertices[i]->checkConsistency();
		if (not vertexConsistent) {
			allVertConsistent = false;
			if (_debug) {
				if (vertexConsistent)
					printDebug << "success: isobar decay vertex " << *_isobarVertices[i]
					           << " is consistent" << endl;
				else
					printWarn << "isobar decay vertex " << *_isobarVertices[i] << " is not consistent" << endl;
			}
		}
	}
	if (_debug) {
		if (allVertConsistent)
			printDebug << "success: information in isobar decay vertices is consistent" << endl;
		else
			printWarn << "information in isobar decay vertices is not consistent" << endl;
	}
	return allVertConsistent;
}


const TLorentzVector&
isobarDecayTopology::calcIsobarLzVec()
{
	// loop over isobar decay vertices and propagate Lorentz-vectors from final state particles up to X-system
	for (int i = nmbDecayVertices() - 1; i >= 0; --i) {
		if (_debug)
			printDebug << "calculating Lorentz-vector of parent isobar '"
			           << _isobarVertices[i]->parent()->name() << "' "
			           << "of node[" << node(_isobarVertices[i]) << "]" << endl;
		_isobarVertices[i]->calcParentLzVec();
	}
	return XIsobarDecayVertex()->parent()->lzVec();
}


void
isobarDecayTopology::calcIsobarCharges()
{
	// loop over isobar decay vertices and propagate charges from final state particles up to X-system
	for (int i = nmbDecayVertices() - 1; i >= 0; --i) {
		const particlePtr& isobar = _isobarVertices[i]->parent();
		if (_debug)
			printDebug << "calculating charge of parent isobar '"
			           << isobar->name() << "' "
			           << "of node[" << node(_isobarVertices[i]) << "]" << endl;
		const int isobarCharge = isobar->charge();
		if (isobarCharge != _isobarVertices[i]->calcParentCharge())
			if (isobar != XParticle()) {
				printWarn << "fixed charge of isobar '" << isobar->name() << "' "
				          << "from " << isobarCharge << " to " << isobar->charge() << ". "
				          << "please fix wave definition." << endl;
				isobar->fillFromDataTable(isobar->name());
			}
	}
	// correct C-parity of X, if necessary
	if ((abs(XParticle()->charge()) > 0) and (XParticle()->C() != 0)) {
		printWarn << "X is charged, but has C-parity = " << XParticle()->C()
		          << ". setting C-parity to zero. please fix wave definition." << endl;
		XParticle()->setC(0);
	}
	// update graph name
	name() = "\"" + XParticle()->qnSummary() + "\"";
}


void
isobarDecayTopology::calcIsobarBaryonNmbs()
{
	// loop over isobar decay vertices and propagate baryon numbers from final state particles up to X-system
	for (int i = nmbDecayVertices() - 1; i >= 0; --i) {
		if (_debug)
			printDebug << "calculating baryon number of parent isobar '"
			           << _isobarVertices[i]->parent()->name() << "' "
			           << "of node[" << node(_isobarVertices[i]) << "]" << endl;
		_isobarVertices[i]->calcParentBaryonNmb();
	}
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
	for (unsigned int i = 0; i < nmbDecayVertices(); ++i)
		out << "    isobar decay node[" << node(decayVertices()[i]) << "] = "
		    << *decayVertices()[i] << endl;
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
		const int fsIndex = fsParticlesIndex(part);
		if (fsIndex >= 0)
			out << "[" << fsIndex << "] (final state)";
		out << endl;
	}
	return out;
}


ostream&
isobarDecayTopology::writeGraphViz(ostream& out)
{
	if (_debug)
		printDebug << "generating graphViz file for graph '" << name() << "'" << endl;
	// set global properties
	graphNodeAttribute()["style"    ] = "filled";
	graphNodeAttribute()["fillcolor"] = "white";
	// set node names
	nodeIterator iNd, iNdEnd;
	for (tie(iNd, iNdEnd) = nodes(); iNd != iNdEnd; ++iNd) {
		stringstream label;
		label << vertex(*iNd)->name();
		const isobarDecayVertexPtr& vert = dynamic_pointer_cast<isobarDecayVertex>(vertex(*iNd));
		if (vert)
			label << ": L = " << spinQn(vert->L()) << ", S = " << spinQn(vert->S());
		nodeAttribute(*iNd)["label"] = label.str();
	}
	// set node shapes
	nodeAttribute(productionVertex())["shape"] = "box";
	for (unsigned int i = 0; i < nmbFsParticles(); ++i) {
		nodeAttribute(toNode(fsParticles()[i]))["shape"] = "diamond";
		nodeAttribute(toNode(fsParticles()[i]))["style"] = "dashed,filled";
	}
	// set edge names
	edgeIterator iEd, iEdEnd;
	for (tie(iEd, iEdEnd) = edges(); iEd != iEdEnd; ++iEd)
		edgeAttribute(*iEd)["label"] = particle(*iEd)->label();
	// set X edge name
	edgeAttribute(edge(XParticle()))["label"] = XParticle()->qnSummary();
	decayTopologyGraphType::writeGraphViz(out);
	return out;
}


bool
isobarDecayTopology::writeGraphViz(const string& outFileName)
{
	ofstream graphVizFile(outFileName.c_str());
	if (not graphVizFile) {
		printWarn << "cannot create file '" << outFileName << "'. graph is not written." << endl;
		return false;
	}
	writeGraphViz(graphVizFile);
	return true;
}


void
isobarDecayTopology::buildIsobarVertexArray()
{
	bool success = true;
	_isobarVertices.resize(nmbDecayVertices());
	for (unsigned int i = 0; i < nmbDecayVertices(); ++i) {
		const interactionVertexPtr& v = decayVertices()[i];
		_isobarVertices[i] = dynamic_pointer_cast<isobarDecayVertex>(v);
		if (not _isobarVertices[i]) {
			printWarn << *v << " is not of type isobarDecayVertex." << endl;
			success = false;
		}
	}
	if (not success) {
		printErr << "incompatible topology to copy from. some interaction vertices are "
		         << "not of type isobarDecayVertex. aborting." << endl;
		throw;
	}
}


isobarDecayTopology
isobarDecayTopology::subDecay(const nodeDesc& startNd,
                              const bool      linkToParentTopo)
{
	isobarDecayTopology subTopo(decayTopology::subDecay(startNd, linkToParentTopo));
	subTopo.name() = "\"" + vertex(startNd)->inParticles()[0]->qnSummary() + "\"";
	return subTopo;
}


void
isobarDecayTopology::addDecay(const isobarDecayTopology& topo)
{
	decayTopology::addDecay(topo);
	buildIsobarVertexArray();
}


isobarDecayTopology
isobarDecayTopology::joinDaughterDecays(const isobarDecayVertexPtr&        parentVertex,
                                        const vector<isobarDecayTopology>& daughterDecays)  ///< joins daughter decay graphs and connects them to a common parent vertex
{
	if (_debug) {
		printDebug << "joining " << daughterDecays.size() << " daughter graphs with parent "
		           << *parentVertex << endl;
	}
	isobarDecayTopology newTopo;
	newTopo.addVertex(parentVertex);
	for (unsigned int i = 0; i < daughterDecays.size(); ++i)
		newTopo.addDecay(daughterDecays[i]);
	return newTopo;
}


isobarDecayTopology
isobarDecayTopology::joinDaughterDecays(const isobarDecayVertexPtr& parentVertex,
                                        const isobarDecayTopology&  daughter1Decay,
                                        const isobarDecayTopology&  daughter2Decay)  ///< joins daughter decay graphs and connects them to a common parent vertex
{
	std::vector<isobarDecayTopology> daughterDecays(2);
	daughterDecays[0] = daughter1Decay;
	daughterDecays[1] = daughter2Decay;
	return joinDaughterDecays(parentVertex, daughterDecays);
}

void isobarDecayTopology::doStuff() {

	const std::vector<particlePtr> fsParts = fsParticles();

	//Get a vector of groups of particles which have to be permutated.
	std::vector< std::vector<unsigned int> > groups;
	for(unsigned int i = 0; i < fsParts.size(); ++i) {
		const particlePtr& particle = fsParts.at(i);
		bool inserted = false;
		for(unsigned int j = 0; j < groups.size(); ++j) {
			const particlePtr& compPart = fsParts.at(groups.at(j).at(0));
			if(particle->isospin() == compPart->isospin() &&
			   particle->J() == compPart->J() &&
			   particle->P() == compPart->P() &&
//			   particle->isospinProj() == compPart->isospinProj() &&
			   particle->C() == compPart->C())
			{
				groups.at(j).push_back(i);
				inserted = true;
			}
		}
		if(!inserted) {
			groups.push_back(std::vector<unsigned int>(1,i));
		}
	}

	printDebug<<"Doing stuff"<<std::endl;
	for(unsigned int i = 0; i < groups.size(); ++i) {
		printDebug<<"Group "<<i<<std::endl;
		for(unsigned int j = 0; j < groups.at(i).size(); ++j) {
			printDebug<<j<<": "<<groups.at(i).at(j)<<" ("<<fsParts.at(groups.at(i).at(j))->name()<<")"<<std::endl;
		}
	}

	// Saving all the isospins
	std::vector<int> isospinProjs;
	for(unsigned int j = 0; j < fsParts.size(); ++j) {
		isospinProjs.push_back(fsParts.at(j)->isospinProj());
	}

	// Permutating all the particles and checking if the permutation makes sense	
	for(unsigned int i = 0; i < groups.size(); ++i) {
		std::vector<unsigned int> group = groups.at(i);

		printDebug<<"Group permutations "<<i<<std::endl;

		// The unity permutation of the right lenght
		std::vector<unsigned int> permutations;
		for(unsigned int j = 0; j < group.size(); ++j) {
			permutations.push_back(j);
		}
		// Permutate the permutations
		do {
//			printDebug<<"New Permutation"<<std::endl;
			// Create an virgin map
			std::vector<unsigned int> map;
			for(unsigned int j = 0; j < fsParts.size(); ++j) {
				map.push_back(j);
			}

			// Swap stuff in the map
			for(unsigned int j = 0; j < group.size(); ++j) {
				map.at(group.at(j)) = group.at(permutations.at(j));
			}

			// Make sure we are not swapping two indistinguishable particles.
			bool breaking = false;
			for(unsigned int j = 0; j < map.size(); ++j) {
				if(j == map.at(j)) {
					continue;
				}
				if(fsParts.at(j)->name() == fsParts.at(map.at(j))->name()) {
					breaking = true;
//					printDebug<<"Swapping same particles ("<<j<<"<->"<<map.at(j)<<" | "<<fsParts.at(j)->name()<<"<->"<<fsParts.at(map.at(j))->name()<<")"<<std::endl;
					break;
				}
			}
			if(breaking) {
				continue;
			}
			printDebug<<"Found valid permutation: ";
			for(unsigned int j = 0; j < map.size(); ++j) {
				std::cout<<map.at(j);
			}
			std::cout<<std::endl;

			// Set the isospin in the final state particles
			for(unsigned int j = 0; j < map.size(); ++j) {
				fsParts.at(j)->setIsospinProj(isospinProjs.at(j));
			}



		} while(next_permutation(permutations.begin(), permutations.end()));



	}


























}
