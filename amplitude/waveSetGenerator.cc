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
//      class that generates all possible waves given the final state
//      particles, a list of isobars and constraints on I, J, L, and S
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "particleDataTable.h"
#include "waveSetGenerator.h"

	
using namespace std;
using namespace boost;
using namespace rpwa;


bool waveSetGenerator::_debug = false;


waveSetGenerator::waveSetGenerator()
{
	reset();
}


waveSetGenerator::~waveSetGenerator()
{ }


void
waveSetGenerator::reset()
{
	// _isospinRange    = make_pair(0, 0);
	// _JRange          = make_pair(0, 0);
	// _LRange          = make_pair(0, 0);
	// _SRange          = make_pair(0, 0);
	_isospinRange    = make_pair(0, 2);
	_JRange          = make_pair(0, 8);
	_LRange          = make_pair(0, 6);
	_SRange          = make_pair(0, 6);
	_allowJpcExotics = false;
	_isobarBlackList.clear();
	_isobarWhiteList.clear();
	_waveSet.clear();
}


size_t
waveSetGenerator::generateWaveSet(const isobarDecayTopologyPtr& templateTopo)
{
	if (not templateTopo) {
		printErr << "null pointer to template graph topology. aborting." << endl;
		throw;
	}
	if (not templateTopo->checkTopology()) {
		printErr << "ill-formed template graph topology. aborting." << endl;
		throw;
	}
	// order nodes depth-first
	vector<nodeDesc> startNds
		= templateTopo->sortNodesDfs(templateTopo->XIsobarDecayVertex());
	// create decay topologies of all subdecays
	vector<isobarDecayTopology> subDecays(startNds.size());
	for (size_t i = 0; i < startNds.size(); ++i)
		subDecays[i] = templateTopo->subDecay(startNds[i]);
	// reverse order, because decay trees are build starting from the final state nodes
	reverse(startNds.begin (), startNds.end ());
	reverse(subDecays.begin(), subDecays.end());
	if (_debug)
		for (size_t i = 0; i < subDecays.size(); ++i)
			printInfo << "created subdecay[" << i << "]: " << subDecays[i];

	// create all possible subdecays
	map<nodeDesc, vector<isobarDecayTopology> > decayPossibilities;  // all possible decay graphs starting at certain node
	for (size_t iStart = 0; iStart < startNds.size(); ++iStart) {
    
		// special case for final state nodes
		if (templateTopo->isFsVertex(templateTopo->vertex(startNds[iStart]))) {
			// final state vertices have no daughter tracks; subdecay
			// topology consist of just the final state vertex
			decayPossibilities[startNds[iStart]].push_back(subDecays[iStart]);
			continue;
		}

		// get daughter decay topologies
		vector<vector<isobarDecayTopology>* > daughterDecays;
		adjIterator iNd, iNdEnd;
		for (tie(iNd, iNdEnd) = templateTopo->adjacentVertices(startNds[iStart]); iNd != iNdEnd; ++iNd)
			daughterDecays.push_back(&decayPossibilities[*iNd]);

		// loop over all combinations of daughter decays
		size_t iDaughter[2];
		for (iDaughter[0] = 0; iDaughter[0] < daughterDecays[0]->size(); ++iDaughter[0])
			for (iDaughter[1] = 0; iDaughter[1] < daughterDecays[1]->size(); ++iDaughter[1]) {

				// copy parent vertex
				isobarDecayVertexPtr parentVertex
					(new isobarDecayVertex(*static_pointer_cast<isobarDecayVertex>
					                       (templateTopo->vertex(startNds[iStart]))));
				particlePtr parent = parentVertex->parent();
				// get daughter particles from the respective decay topologies
				const nodeDesc topNodes[2] = 
					{(*daughterDecays[0])[iDaughter[0]].topNode(),
					 (*daughterDecays[1])[iDaughter[1]].topNode()};
				const particlePtr daughters[2] = 
					{(*daughterDecays[0])[iDaughter[0]].vertex(topNodes[0])->inParticles()[0],
					 (*daughterDecays[1])[iDaughter[1]].vertex(topNodes[1])->inParticles()[0]};
				// set daughters of parent vertex
				parentVertex->daughter1() = daughters[0];
				parentVertex->daughter2() = daughters[1];
				// join daughter subdecays and parent vertex
				isobarDecayTopology parentDecay
					= templateTopo->joinDaughterDecays(parentVertex,
					                                   (*daughterDecays[0])[iDaughter[0]],
					                                   (*daughterDecays[1])[iDaughter[1]]);
	
				// calculate parent quantum numbers fixed by daughter quantum numbers
				const int baryonNmb   = daughters[0]->baryonNmb()   + daughters[1]->baryonNmb();
				const int charge      = daughters[0]->charge()      + daughters[1]->charge();
				const int strangeness = daughters[0]->strangeness() + daughters[1]->strangeness();
				const int charm       = daughters[0]->charm()       + daughters[1]->charm();
				const int beauty      = daughters[0]->beauty()      + daughters[1]->beauty();
				const int G           = daughters[0]->G()           * daughters[1]->G();
				// daughter quantum numbers that define ranges of possible parent quantum numbers
				const int spins   [2] = {daughters[0]->J(),       daughters[1]->J()      };
				const int isospins[2] = {daughters[0]->isospin(), daughters[1]->isospin()};

				// loop over all allowed combinations of daughter quantum numbers
				// loop over allowed total spins
				for (int S = max(abs(spins[0] - spins[1]), _SRange.first);
				     S <= min(spins[0] + spins[1], _SRange.second); S += 2) {
					// loop over allowed relative orbital angular momenta
					for (int L = max(0, _LRange.first); L <= _LRange.second; L += 2) {
						const int P = daughters[0]->P() * daughters[1]->P() * (L % 4 == 0 ? 1 : -1);  // parity
						// loop over allowed total angular momenta
						for (int J = max(abs(L - S), _JRange.first); J <= min(L + S, _JRange.second); J += 2) {
							// loop over allowed isospins
							for (int I = max(abs(isospins[0] - isospins[1]), _isospinRange.first);
							     I <= min(isospins[0] + isospins[1], _isospinRange.second); I += 2) {
								// check whether charge state is allowed
								// !!! the Gell-Mann-Nishijima formula cannot be used
								//     here, because it employs the the z-component of
								//     the isospin, I_z (not I); see PDG 2008 eq. 14.1
								// if (abs(charge - 0.5 * (baryonNmb + strangeness + charm + beauty)) != 0.5 * I)
								//   continue;
								const int C = G * (I % 4 == 0 ? 1 : -1);  // C-parity
								if (not _allowJpcExotics)
									// quark model restrictions: P == C is always allowed
									// check that P = (-1)^(J + 1)
									if ((P != C) and (   (C != (J % 4     == 0 ? 1 : -1))
									                  or (P != (J + 2 % 4 == 0 ? 1 : -1)))) {
										if (_debug)
											printInfo << "disregarding spin-exotic isobar with IG(JPC) = "
											          << 0.5 * I << sign(G)
											          << "(" << 0.5 * J << sign(P) << sign(C) << ")" << endl;
										continue;
									}
		
								// find candidates for parent isobar in particle data table
								particleProperties isobarProp(*parent);
								isobarProp.setBaryonNmb(baryonNmb);
								isobarProp.setSCB      (strangeness, charm, beauty);
								isobarProp.setIGJPC    (I, G, J, P, C);
								particleDataTable& pdt = particleDataTable::instance();
								vector<const particleProperties*> isobarCandidates;
								const string parentName    = parent->name();
								const double minIsobarMass = daughters[0]->mass() + daughters[1]->mass();
								if (   (parentName == "X-") or (parentName == "X0")
								    or (parentName == "X+") or (parentName == "X")) {
									if (parent->mass() + parent->width() < minIsobarMass) {
										cout << isobarProp.qnSummary() << " mass too low. skipping candidate." << endl;
										continue;
									}
									isobarCandidates.push_back(&isobarProp);
								} else {
									isobarCandidates = pdt.entriesMatching(isobarProp, "allQn", minIsobarMass,
									                                       _isobarWhiteList, _isobarBlackList);
									if (_debug)
										printInfo << "found " << isobarCandidates.size() << " isobar candidates for "
										          << isobarProp.qnSummary() << " in particle data table" << endl;
								}

								// loop over isobar candidates
								for (size_t iIsobar = 0; iIsobar < isobarCandidates.size(); ++iIsobar) {
									// clone topology
									const isobarDecayTopologyPtr& decayCopy = parentDecay.clone(false, false);
									// set parent vertex quantum numbers
									const isobarDecayVertexPtr& parentVertexCopy =
										static_pointer_cast<isobarDecayVertex>
										(decayCopy->vertex(parentDecay.node(parentVertex)));
									parentVertexCopy->setL(L);
									parentVertexCopy->setS(S);
									// set parent particle
									const particlePtr& parentCopy = parentVertexCopy->parent();
									parentCopy->setProperties(*isobarCandidates[iIsobar]);
									parentCopy->setCharge(charge);
									if (_debug)
										printInfo << "created decay topology for " << *parentVertexCopy << endl;
									decayPossibilities[startNds[iStart]].push_back(*decayCopy);
								} // isobar candidate loop
							}  // isospin loop
						}  // L-S coupling loop
					}  // L loop
				}  // S loop
			}  // loop over daughter decays
	}  // loop over all start nodes
	
	// extract decays for X-decay vertex and add production vertex
	_waveSet = decayPossibilities[startNds.back()];
	for (size_t i = 0; i < _waveSet.size(); ++i) {
		// clone production vertex and set X-particle
		const productionVertexPtr newProdVert(static_pointer_cast<rpwa::productionVertex>
		                                      (templateTopo->productionVertex()->clone(false, false)));
		const particlePtr&        newX = _waveSet[i].XIsobarDecayVertex()->parent();
		newProdVert->outParticles()[0] = newX;
		// add production vertex
		_waveSet[i].setProductionVertex(newProdVert);
	}

	return _waveSet.size();
}


ostream&
waveSetGenerator::print(ostream& out) const
{
	out << "wave set generator:" << endl
	    << "    isospin range ....... [" << 0.5 * _isospinRange.first << ", "
	    << 0.5 * _isospinRange.second << "]" << endl
	    << "    J range ............. [" << 0.5 * _JRange.first       << ", "
	    << 0.5 * _JRange.second       << "]" << endl
	    << "    L range ............. [" << 0.5 * _LRange.first       << ", "
	    << 0.5 * _LRange.second       << "]" << endl
	    << "    S range ............. [" << 0.5 * _SRange.first       << ", "
	    << 0.5 * _SRange.second       << "]" << endl
	    << "    allow JPC exotics ... " << ((_allowJpcExotics) ? "true" : "false")  << endl;
	out << "    isobar black list:";
	if (_isobarBlackList.size() == 0)
		out << " empty" << endl;
	else
		for (size_t i = 0; i < _isobarBlackList.size(); ++i)
			out << "        " << _isobarBlackList[i] << endl;
	out << "    isobar white list:";
	if (_isobarWhiteList.size() == 0)
		out << " all isobars" << endl;
	else
		for (size_t i = 0; i < _isobarWhiteList.size(); ++i)
			out << "        " << _isobarWhiteList[i] << endl;
	return out;
}
