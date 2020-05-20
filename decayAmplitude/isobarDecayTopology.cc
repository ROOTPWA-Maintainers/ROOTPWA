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
#include "nonInteractionVertex.h"
#include "particleDataTable.h"
#include "spinUtils.hpp"


using namespace std;
using namespace boost;
using namespace rpwa;


namespace {

	vector< vector<unsigned int> >	getFsParticleGroupsByJPI(const vector<particlePtr>& fsParts);

	vector< vector< vector<unsigned int> > > getAllGroupPermutations(const vector< vector<unsigned int> >& groupsExt);

	std::string	permMapString(const vector<unsigned int>& permMap);

}


bool isobarDecayTopology::_debug = false;


isobarDecayTopology::isobarDecayTopology()
	: decayTopology()
{ }


isobarDecayTopology::isobarDecayTopology(const productionVertexPtr&          productionVertex,
                                         const vector<isobarDecayVertexPtr>& isobarDecayVertices,
                                         const vector<particlePtr>&          fsParticles,
                                         const bool                          performTopologyCheck)
	: decayTopology()
{
	constructDecay(productionVertex, isobarDecayVertices, fsParticles, performTopologyCheck);
}


isobarDecayTopology::isobarDecayTopology(const productionVertexPtr&          productionVertex,
                                         const vector<interactionVertexPtr>& isobarDecayVertices,
                                         const vector<particlePtr>&          fsParticles,
                                         const bool                          performTopologyCheck)
	: decayTopology()
{
	constructDecay(productionVertex, isobarDecayVertices, fsParticles, performTopologyCheck);
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
		           << ((cloneFsParticles   ) ? "in" : "ex") << "cluding final-state particles, "
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
                                    const vector<particlePtr>&          fsParticles,
                                    const bool                          performTopologyCheck)
{
	const unsigned int nmbVert = isobarDecayVertices.size();
	if (_debug)
		printDebug << "constructing isobar decay topology with "
		           << fsParticles.size() << " final-state particles and "
		           << nmbVert            << " isobar decay vertices" << endl;
	vector<interactionVertexPtr> intVertices(nmbVert);
	for (unsigned int i = 0; i < nmbVert; ++i)
		intVertices[i] = isobarDecayVertices[i];
	decayTopology::constructDecay(productionVertex, intVertices, fsParticles, performTopologyCheck);
	// copy sorted vertices
	buildIsobarVertexArray();
	if (performTopologyCheck and (nmbDecayVertices() != nmbVert)) {
		printErr << "number of interaction vertices  = " << nmbDecayVertices()
		         << " does not match number of vertices given in parameter array = " << nmbVert
		         << ". Aborting..." << endl;
		throw;
	}
	return *this;
}


isobarDecayTopology&
isobarDecayTopology::constructDecay(const productionVertexPtr&          productionVertex,
                                    const vector<interactionVertexPtr>& isobarDecayVertices,
                                    const vector<particlePtr>&          fsParticles,
                                    const bool                          performTopologyCheck)
{
	const unsigned int nmbVert = isobarDecayVertices.size();
	vector<isobarDecayVertexPtr> decayVertices(nmbVert);
	for (unsigned int i = 0; i < nmbVert; ++i) {
		decayVertices[i] = dynamic_pointer_cast<isobarDecayVertex>(isobarDecayVertices[i]);
		if (not decayVertices[i]) {
			printErr << "interaction vertex[" << i << "] is not an isobarDecayVertex. Aborting..." << endl;
			throw;
		}
	}
	return constructDecay(productionVertex, decayVertices, fsParticles, performTopologyCheck);
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
	// loop over isobar decay vertices and propagate Lorentz-vectors from final-state particles up to X-system
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
isobarDecayTopology::calcIsobarCharges(const bool warnIfNotExistent)
{
	// loop over isobar decay vertices and propagate charges from final-state particles up to X-system
	for (int i = nmbDecayVertices() - 1; i >= 0; --i) {
		const particlePtr& isobar = _isobarVertices[i]->parent();
		if (_debug and warnIfNotExistent)
			printDebug << "calculating charge of parent isobar '"
			           << isobar->name() << "' "
			           << "of node[" << node(_isobarVertices[i]) << "]" << endl;
		const int isobarCharge = isobar->charge();
		if (isobarCharge != _isobarVertices[i]->calcParentCharge())
			if (isobar != XParticle()) {
				if(warnIfNotExistent) {
					printWarn << "fixed charge of isobar '" << isobar->name() << "' "
					          << "from " << isobarCharge << " to " << isobar->charge() << ". "
					          << "please fix wave definition." << endl;
				}
				isobar->fillFromDataTable(isobar->name(), warnIfNotExistent);
			}
	}
	// correct C-parity of X, if necessary
	if ((abs(XParticle()->charge()) > 0) and (XParticle()->C() != 0)) {
		if(warnIfNotExistent) {
			printWarn << "X is charged, but has C-parity = " << XParticle()->C()
			          << ". setting C-parity to zero. please fix wave definition." << endl;
			XParticle()->setC(0);
		}
	}
	// update graph name
	name() = "\"" + XParticle()->qnSummary() + "\"";
}


void
isobarDecayTopology::calcIsobarBaryonNmbs()
{
	// loop over isobar decay vertices and propagate baryon numbers from final-state particles up to X-system
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
			out << "    final-state node[" << *iNd << "] = " << *vertex(*iNd) << endl;
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
		printErr << "incompatible topology. some interaction vertices are "
		         << "not of type isobarDecayVertex. Aborting..." << endl;
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


isobarDecayTopology
isobarDecayTopology::subDecayConsistent(const isobarDecayVertexPtr& startVertex)
{
	isobarDecayTopology subTopo = subDecay(startVertex);
	nonInteractionVertexPtr prodVertex = createNonInteractionVertex(startVertex->parent());
	subTopo.setProductionVertex(prodVertex);
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
	vector<isobarDecayTopology> daughterDecays(2);
	daughterDecays[0] = daughter1Decay;
	daughterDecays[1] = daughter2Decay;
	return joinDaughterDecays(parentVertex, daughterDecays);
}


double
isobarDecayTopology::getIsospinClebschGordanProduct(isobarDecayVertexPtr vertex) const
{
	if (not vertex)
		vertex = static_pointer_cast<isobarDecayVertex>(XDecayVertex());

	const particlePtr daughter1 = vertex->daughter1();
	const particlePtr daughter2 = vertex->daughter2();
	const particlePtr parent    = vertex->parent   ();

	const double clebsch = clebschGordanCoeff<double>(daughter1->isospin(), daughter1->isospinProj(),
	                                                  daughter2->isospin(), daughter2->isospinProj(),
	                                                  parent->isospin   (), parent->isospinProj   ());
	double clebschDaughter1 = 1.;
	double clebschDaughter2 = 1.;
	if (not isFsParticle(daughter1))
		clebschDaughter1
			= getIsospinClebschGordanProduct(static_pointer_cast<isobarDecayVertex>(toVertex(daughter1)));
	if (not isFsParticle(daughter2))
		clebschDaughter2
			= getIsospinClebschGordanProduct(static_pointer_cast<isobarDecayVertex>(toVertex(daughter2)));
	return clebsch * clebschDaughter1 * clebschDaughter2;
}


vector<symTermMap>
isobarDecayTopology::getIsospinSymmetrization()
{
	const vector<particlePtr> fsParts = fsParticles();

	// Construct a vector of groups of indices of final-state particles that have the same J, P, and I.
	const vector< vector<unsigned int> > groups = getFsParticleGroupsByJPI(fsParts);

	// Saving the z-projections of all isospins .
	// Needed to reset everything as it was at the end.
	vector<int> isospinProjs;
	for(unsigned int i = 0; i < fsParts.size(); ++i) {
		isospinProjs.push_back(fsParts[i]->isospinProj());
	}

	// A vector to save the found permutations and their Clebsch-Gordans
	vector<symTermMap> symAmplitudes;

	// Permuting all the particles and checking if the permutation makes sense
	const vector< vector< vector<unsigned int> > > permutations = getAllGroupPermutations(groups);
	for(unsigned int i_permutation = 0; i_permutation < permutations.size(); ++i_permutation) {
		const vector< vector<unsigned int> >& permutation = permutations[i_permutation];
		vector<unsigned int> permMap(fsParts.size());
		for(unsigned int i_group = 0; i_group < groups.size(); ++i_group) {
			for(unsigned int i_groupEntry = 0; i_groupEntry < groups[i_group].size(); ++i_groupEntry) {
				permMap[ groups[i_group][i_groupEntry] ] = permutation[i_group][i_groupEntry];
			}
		}
		// check that the to-be-added permutation permMap is not equivalent to identity
		// and not identical to any of the already added permutations
		// assumes that first permMap is the identity
		if(_debug) {
			printDebug << "processing permutation " << permMapString(permMap) << endl;
		}
		bool discardPermutation = false;
		for(vector<symTermMap>::const_iterator symAmp = symAmplitudes.begin(); symAmp != symAmplitudes.end(); ++symAmp) {

			if(_debug) {
				cout << "    ... checking against permutation " << permMapString(symAmp->fsPartPermMap) << endl;
			}

			for(unsigned int i = 0; i < permMap.size(); ++i) {

				if (symAmp->fsPartPermMap[i] == permMap[i]) {
					continue; // particle i was not swapped -> go on with check of next particle
				}

				// Make sure we are not swapping two indistinguishable particles.
				if(fsParts.at(symAmp->fsPartPermMap[i])->name() == fsParts.at(permMap[i])->name()) {
					discardPermutation = true;
					break;
				}

				// Check if two particles from the same isobar are being swapped.
				for(unsigned int j = i + 1; j < permMap.size(); ++j) {
					if (symAmp->fsPartPermMap[j] == permMap[j]) {
						continue;
					}
					if(    fromVertex(fsParts.at(symAmp->fsPartPermMap[i])) == fromVertex(fsParts.at(symAmp->fsPartPermMap[j]))
					   and fromVertex(fsParts.at(symAmp->fsPartPermMap[i])) == fromVertex(fsParts.at(permMap[i]))
					   and fromVertex(fsParts.at(permMap[i]))               == fromVertex(fsParts.at(permMap[j]))) {
						discardPermutation = true;
						break;
					}
				}
				if(discardPermutation) {
					break;
				}
			}
			if(discardPermutation) {
				break;
			}
		}

		if(discardPermutation) { // check next permutation
			continue;
		}

		// Set the isospin of the final-state particles as given by the permutation.
		for(unsigned int i = 0; i < permMap.size(); ++i) {
			fsParts.at(i)->setIsospinProj(isospinProjs.at(permMap[i]));
		}
		calcIsobarCharges(false);

		// Check for isospin consistency in all vertices.
		for(unsigned int i = 0; i < _isobarVertices.size(); ++i) {
			const particlePtr d1 = _isobarVertices[i]->daughter1();
			const particlePtr d2 = _isobarVertices[i]->daughter2();
			const particlePtr p  = _isobarVertices[i]->parent();
			if(not spinStatesCanCouple(d1->isospin(), d1->isospinProj(),
			                           d2->isospin(), d2->isospinProj(),
			                            p->isospin(),  p->isospinProj())) {
				discardPermutation = true;
			}
		}

		// If something is amiss, reset isospin projections and proceed to the next permutation.
		if(discardPermutation) {
			for(unsigned int i = 0; i < fsParts.size(); ++i) {
				fsParts[i]->setIsospinProj(isospinProjs.at(i));
			}
			calcIsobarCharges();
			continue;
		}

		// Save permutation that survived all criteria
		symTermMap symTerm(getIsospinClebschGordanProduct(), permMap);
		symAmplitudes.push_back(symTerm);

		if(_debug) {
			cout << "    ... found valid permutation: " << permMapString(permMap) << " = (";
			for(unsigned int i = 0; i < permMap.size(); ++i) {
				cout << fsParts.at(i)->name();
			}
			cout << "), Clebsch-Gordan = " << getIsospinClebschGordanProduct() << endl;
		}

		// Resetisospin projections for the next permutation
		for(unsigned int i = 0; i < fsParts.size(); ++i) {
			fsParts[i]->setIsospinProj(isospinProjs.at(i));
		}
		calcIsobarCharges();

	}  // End of the loop over the permutations

	for(unsigned int i = 0; i < symAmplitudes.size(); ++i) {
		symAmplitudes[i].factor /= sqrt(symAmplitudes.size());
	}

	return symAmplitudes;
}


vector<symTermMap>
isobarDecayTopology::getBoseSymmetrization() const
{
	// get final state indistinguishable particles
	typedef map<string, unsigned int>::const_iterator indistFsPartIt;
	const map<string, unsigned int> indistFsPart = nmbIndistFsParticles();
	printInfo << "indistinguishable final-state multiplicities "
	          << "(marked FS particles will be Bose symmetrized): ";
	for (indistFsPartIt i = indistFsPart.begin(); i != indistFsPart.end(); ++i)
		cout << i->first << " = " << i->second << ((i->second) >= 2 ? " <<<  " : "  ");
	cout << endl;

	// calculate normalization factor
	double nmbCombinations = 1;
	for (indistFsPartIt i = indistFsPart.begin(); i != indistFsPart.end(); ++i)
		nmbCombinations *= factorial<double>(i->second);
	const double normFactor = 1 / sqrt(nmbCombinations);

	// initialize indices used to generate final-state permutation maps
	// in order to get all permutations with std::next_permutation indices have to be sorted ascending
	map<string, vector<unsigned int> > origFsPartIndices;
	for (unsigned int i = 0; i < nmbFsParticles(); ++i) {
		const string partName = fsParticles()[i]->name();
		origFsPartIndices[partName].push_back(i);
	}
	map<string, vector<unsigned int> >           newFsPartIndices = origFsPartIndices;
	map<string, vector<unsigned int> >::iterator firstEntry       = newFsPartIndices.begin();
	vector<symTermMap> symTermMaps;
	genBoseSymTermMaps(origFsPartIndices, newFsPartIndices, firstEntry, symTermMaps);
	for(unsigned int i = 0; i < symTermMaps.size(); ++i) {
		symTermMaps[i].factor = normFactor;
	}
	return symTermMaps;
}


bool
isobarDecayTopology::isobarIsAffectedByPermutation(const isobarDecayVertexPtr& vertex,
                                                   const vector<unsigned int>& permutation) const
{
	vector<unsigned int> fsPartIndices = getFsPartIndicesConnectedToVertex(vertex);
	for(unsigned int i = 0; i < permutation.size(); ++i) {
		if(permutation[i] == i) {
			continue;
		}
		if((find(fsPartIndices.begin(), fsPartIndices.end(), permutation[i]) == fsPartIndices.end() and
			find(fsPartIndices.begin(), fsPartIndices.end(), i) != fsPartIndices.end()) or
		   (find(fsPartIndices.begin(), fsPartIndices.end(), i) == fsPartIndices.end() and
			find(fsPartIndices.begin(), fsPartIndices.end(), permutation[i]) != fsPartIndices.end()))
		{
			return true;
		}
	}
	return false;
}


bool
isobarDecayTopology::daughtersAreAffectedByPermutation(const isobarDecayVertexPtr& vertex,
                                                       const vector<unsigned int>& permutation) const
{
	vector<isobarDecayVertexPtr> verticesToTest(1, vertex);
	isobarDecayVertexPtr daughter;
	if(not isFsParticle(vertex->daughter1())) {
		daughter = dynamic_pointer_cast<isobarDecayVertex>(toVertex(vertex->daughter1()));
		if(not daughter) {
			printErr << "Got NULL pointer while getting toVertex. Aborting..." << endl;
			throw;
		}
		verticesToTest.push_back(daughter);
	}
	if(not isFsParticle(vertex->daughter2())) {
		daughter = dynamic_pointer_cast<isobarDecayVertex>(toVertex(vertex->daughter2()));
		if(not daughter) {
			printErr << "Got NULL pointer while getting toVertex. Aborting..." << endl;
			throw;
		}
		verticesToTest.push_back(daughter);
	}
	for(unsigned int i = 0; i < verticesToTest.size(); ++i) {
		if(isobarIsAffectedByPermutation(verticesToTest[i], permutation)) {
			return true;
		}
	}
	return false;
}


vector<unsigned int>
isobarDecayTopology::getFsPartIndicesConnectedToVertex(const isobarDecayVertexPtr& vertex) const
{
	vector<unsigned int> indices;
	int index1 = fsParticlesIndex(vertex->daughter1());
	int index2 = fsParticlesIndex(vertex->daughter2());
	if(index1 >= 0) {
		indices.push_back(static_cast<unsigned int>(index1));
	} else {
		const particlePtr& daughter1 = vertex->daughter1();
		if(not isFsParticle(daughter1)) {
			isobarDecayVertexPtr daughterVertex1 = dynamic_pointer_cast<isobarDecayVertex>(toVertex(daughter1));
			if(not daughterVertex1) {
				printErr << "Got NULL pointer while getting toVertex. Aborting..." << endl;
				throw;
			}
			vector<unsigned int> indices1 = getFsPartIndicesConnectedToVertex(daughterVertex1);
			indices.insert(indices.end(), indices1.begin(), indices1.end());
		}
	}
	if(index2 >= 0) {
		indices.push_back(static_cast<unsigned int>(index2));
	} else {
		const particlePtr& daughter2 = vertex->daughter2();
		if(not isFsParticle(daughter2)) {
			isobarDecayVertexPtr daughterVertex2 = dynamic_pointer_cast<isobarDecayVertex>(toVertex(daughter2));
			if(not daughterVertex2) {
				printErr << "Got NULL pointer while getting toVertex. Aborting..." << endl;
				throw;
			}
			vector<unsigned int> indices2 = getFsPartIndicesConnectedToVertex(daughterVertex2);
			indices.insert(indices.end(), indices2.begin(), indices2.end());
		}
	}
	return indices;
}


void
isobarDecayTopology::genBoseSymTermMaps
(const map<string, vector<unsigned int> >&     origFsPartIndices,      // original final-state particle ordering sorted by species
 const map<string, vector<unsigned int> >&     newFsPartIndices,       // new final-state particle ordering sorted by species
 map<string, vector<unsigned int> >::iterator& newFsPartIndicesEntry,  // entry of particle species that is currently being symmetrized
 vector<symTermMap>&                           newSymTermMaps) const   // generated permutation maps
{
	// loop over all permutations for current final-state particle species
	do {
		map<string, vector<unsigned int> >::iterator nextFsPartIndicesEntry = newFsPartIndicesEntry;
		if (++nextFsPartIndicesEntry != newFsPartIndices.end())
			// recurse to permutations of other final-state particle species
			genBoseSymTermMaps(origFsPartIndices, newFsPartIndices, nextFsPartIndicesEntry,
			                   newSymTermMaps);
		else {
			// build map for current permutation
			vector<unsigned int> bosePermMap(nmbFsParticles(), 0);
			if (_debug)
				printDebug << "Bose-symmetrization final-state permutation: ";
			for (map<string, vector<unsigned int> >::const_iterator i = origFsPartIndices.begin();
			     i != origFsPartIndices.end(); ++i) {
				const string partName = i->first;
				map<string, vector<unsigned int> >::const_iterator entry = newFsPartIndices.find(partName);
				assert(entry != newFsPartIndices.end());
				for (unsigned int j = 0; j < i->second.size(); ++j) {
					assert(entry->second.size() == i->second.size());
					const unsigned int origFsPartIndex = i->second[j];
					const unsigned int newFsPartIndex  = entry->second[j];
					if (_debug)
						cout << partName << "[" << origFsPartIndex << " -> " << newFsPartIndex << "]  ";
					bosePermMap[origFsPartIndex] = newFsPartIndex;
				}
			}
			if (_debug)
				cout << endl;
			// compute effective permutation map by reshuffling Bose
			// symmetrization permuation map according to base permutation map
			vector<unsigned int> fsPartPermMap(bosePermMap.size(), 0);
			for (unsigned int i = 0; i < fsPartPermMap.size(); ++i) {
				fsPartPermMap[i] = bosePermMap[i];
			}
			if (_debug)
				printDebug << "effective Bose-symmetrization is: final-state permutation map = "
				           << fsPartPermMap << endl;
			symTermMap symTerm(1, fsPartPermMap);
			newSymTermMaps.push_back(symTerm);
		}
	} while (next_permutation(newFsPartIndicesEntry->second.begin(),
	                          newFsPartIndicesEntry->second.end()));
}


vector<unsigned int>
isobarDecayTopology::findIsobarBoseSymVertices() const
{
	vector<unsigned int>   boseSymVertices;
	isobarDecayTopologyPtr topo = this->clone();  // needed, because subDecay below is not const
	for (unsigned int i = 0; i < topo->nmbDecayVertices(); ++i) {
		// find decay vertex with two isobars as daughters
		isobarDecayVertexPtr vert         = topo->isobarDecayVertices()[i];
		const particlePtr    daughters[2] = {vert->daughter1(), vert->daughter2()};
		if (topo->isFsParticle(daughters[0]) or topo->isFsParticle(daughters[1]))
			continue;
		// make sure the two daughters have the same final state
		const decayTopology isobarTopos[2] = {topo->subDecay(topo->toNode(daughters[0])),
		                                      topo->subDecay(topo->toNode(daughters[1]))};
		if (isobarTopos[0].nmbIndistFsParticles() == isobarTopos[1].nmbIndistFsParticles()) {
			if (_debug)
				printDebug << "found isobar decay vertex " << *vert << " "
				           << "that has two isobar daughters which decay into the same final state: " << endl
				           << "daughter 1 " << isobarTopos[0]
				           << "daughter 2 " << isobarTopos[1];
			boseSymVertices.push_back(i);
		}
	}
	return boseSymVertices;
}


namespace {

	// Constructs a vector of groups of indices of particles which have to be permutated.
	// A group is a vector of indices of those final-state particles that have the same J, P, and I.
	vector< vector<unsigned int> >
	getFsParticleGroupsByJPI(const vector<particlePtr>& fsParts)
	{
		vector< vector<unsigned int> > groups;
		for(unsigned int i = 0; i < fsParts.size(); ++i) {
			const particlePtr& particle = fsParts[i];
			bool addedToGroup = false;
			for(unsigned int j = 0; j < groups.size(); ++j) {
				const particlePtr& compPart = fsParts.at(groups[j][0]);
				if(    particle->isospin() == compPart->isospin()
				   and particle->J()       == compPart->J()
				   and particle->P()       == compPart->P())
				{
					groups[j].push_back(i);
					addedToGroup = true;
				}
			}
			if(not addedToGroup) {
				// create a new group
				groups.push_back(vector<unsigned int>(1, i));
			}
		}

		if(isobarDecayTopology::debug()) {
			printDebug << "found isospin symmetrization groups:" << endl;
			for(unsigned int i = 0; i < groups.size(); ++i) {
				cout << "    ... group " << i << ":" << endl;
				for(unsigned int j = 0; j < groups[i].size(); ++j) {
					cout << "        ... " << j << ": " << groups[i][j]
						<< " (" << fsParts.at(groups[i][j])->name() << ")" << endl;
				}
			}
		}

		return groups;
	}


	// Returns a list with all possible permutations of all groups of particles
	// groupsExt is a vector of groups where each group is a vector of particle ids that belong to that group ( [<group>][<particleInGroup>] )
	// The return value is a list of all permutations where each permutation is a list of groups
	// where each group is a vector of particle ids that belong to that group ( [<permutation>][<group>][<particleInGroup>] )
	vector< vector< vector<unsigned int> > >
	getAllGroupPermutations(const vector< vector<unsigned int> >& groupsExt)
	{
		vector< vector< vector<unsigned int> > > permutations;

		vector< vector<unsigned int> > subGroups(groupsExt.size() - 1);
		for (size_t i = 0; i < groupsExt.size() - 1; ++i) {
			subGroups[i] = groupsExt[i];
		}

		vector<unsigned int> curGroup = groupsExt[groupsExt.size()-1];
		// make sure that the cur_group is sorted such that next_permutation works
		std::sort(curGroup.begin(), curGroup.end());

		if(subGroups.empty()){ // this is the lowest group
			do{
				vector< vector<unsigned int> > lastElement(1, curGroup);
				permutations.push_back(lastElement);
			} while(next_permutation(curGroup.begin(), curGroup.end()));

		} else { // build all permutation, take into account all permutations of all subgroups

			// get all permutations of all subgroups
			const vector< vector< vector<unsigned int> > > subPermutations = getAllGroupPermutations(subGroups);

			do{
				for(unsigned int i_subPermutation = 0; i_subPermutation < subPermutations.size(); ++i_subPermutation){
					vector< vector< unsigned int> > permutation(subPermutations[i_subPermutation]);
					permutation.push_back(curGroup);
					permutations.push_back(permutation);
				}
			} while(next_permutation(curGroup.begin(), curGroup.end()));
		}

		return permutations;
	}


	// converts permutation index map into string
	std::string
	permMapString(const vector<unsigned int>& permMap)
	{
		stringstream s;
		for(unsigned int i = 0; i < permMap.size(); ++i) {
			s << permMap[i];
		}
		return s.str();
	}

}
