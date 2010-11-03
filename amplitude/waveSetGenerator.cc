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


#include "mathUtils.hpp"
#include "clebschGordanCoeff.hpp"
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
	_isospinRange          = make_pair(0, 2);
	_JRange                = make_pair(0, 8);
	_LRange                = make_pair(0, 6);
	_SRange                = make_pair(0, 6);
	_allowJpcExotics       = false;
	_requireMinIsobarMass  = true;
	_isobarMassWindowSigma = 1;
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
				const int parentBaryonNmb   = daughters[0]->baryonNmb()   + daughters[1]->baryonNmb();
				const int parentCharge      = daughters[0]->charge()      + daughters[1]->charge();
				const int parentStrangeness = daughters[0]->strangeness() + daughters[1]->strangeness();
				const int parentCharm       = daughters[0]->charm()       + daughters[1]->charm();
				const int parentBeauty      = daughters[0]->beauty()      + daughters[1]->beauty();
				const int parentG           = daughters[0]->G()           * daughters[1]->G();
				// daughter quantum numbers that define ranges of possible parent quantum numbers
				const int daughterI  [2] = {daughters[0]->isospin(),     daughters[1]->isospin()    };
				const int daughterI_z[2] = {daughters[0]->isospinProj(), daughters[1]->isospinProj()};
				const int daughterJ  [2] = {daughters[0]->J(),           daughters[1]->J()          };

				// loop over all allowed combinations of daughter quantum numbers
				// loop over allowed total spins
				for (int S = max(abs(daughterJ[0] - daughterJ[1]), _SRange.first);
				     S <= min(daughterJ[0] + daughterJ[1], _SRange.second); S += 2) {
					// loop over allowed relative orbital angular momenta
					for (int L = max(0, _LRange.first); L <= _LRange.second; L += 2) {
						const int parentP = daughters[0]->P() * daughters[1]->P() * (L % 4 == 0 ? 1 : -1);
						// loop over allowed total angular momenta
						for (int parentJ = max(abs(L - S), _JRange.first);
						     parentJ <= min(L + S, _JRange.second); parentJ += 2) {
							// loop over allowed isospins
							for (int parentI = max(abs(daughterI[0] - daughterI[1]), _isospinRange.first);
							     parentI <= min(daughterI[0] + daughterI[1], _isospinRange.second); parentI += 2) {
								// check whether charge state is allowed using
								// Gell-Mann-Nishijima formula (see PDG 2008 eq. 14.1)
								const int parentI_z = 2 * parentCharge - (parentBaryonNmb + parentStrangeness
								                                          + parentCharm + parentBeauty);
								if ((abs(parentI_z) > parentI) or isOdd(parentI - parentI_z))
									continue;
								// make sure that isospin Clebsch-Gordan is not zero
								if (clebschGordanCoeff<double>(daughterI[0], daughterI_z[0],
								                               daughterI[1], daughterI_z[1],
								                               parentI, parentI_z) == 0)
									continue;
								const int parentC = parentG * (parentI % 4 == 0 ? 1 : -1);  // C-parity
								if (not _allowJpcExotics)
									// quark model restrictions: P == C is always allowed
									// check that P = (-1)^(J + 1)
									if ((parentP != parentC)
									    and (   (parentC != (parentJ % 4     == 0 ? 1 : -1))
									         or (parentP != (parentJ + 2 % 4 == 0 ? 1 : -1)))) {
										if (_debug)
											printInfo << "disregarding spin-exotic isobar with IG(JPC) = "
											          << 0.5 * parentI << sign(parentG)
											          << "(" << 0.5 * parentJ << sign(parentP) << sign(parentC) << ")"
											          << endl;
										continue;
									}
		
								// find candidates for parent isobar in particle data table
								particleProperties isobarProp(*parent);
								isobarProp.setBaryonNmb(parentBaryonNmb);
								isobarProp.setSCB      (parentStrangeness, parentCharm, parentBeauty);
								isobarProp.setIGJPC    (parentI, parentG, parentJ, parentP, parentC);
								particleDataTable& pdt = particleDataTable::instance();
								vector<const particleProperties*> isobarCandidates;
								const string parentName = parent->name();
								if (   (parentName == "X-") or (parentName == "X0")
								    or (parentName == "X+") or (parentName == "X"))
									isobarCandidates.push_back(&isobarProp);
								else {
									const double minIsobarMass = (_requireMinIsobarMass) ?
										  (daughters[0]->mass() - _isobarMassWindowSigma * daughters[0]->width())
										+ (daughters[1]->mass() - _isobarMassWindowSigma * daughters[1]->width())
										: 0;
									isobarCandidates
										= pdt.entriesMatching(isobarProp, "allQn", minIsobarMass, _isobarMassWindowSigma,
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
									parentCopy->setCharge(parentCharge);
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
	    << "    isospin range .............. [" << 0.5 * _isospinRange.first << ", "
	    << 0.5 * _isospinRange.second << "]" << endl
	    << "    J range .................... [" << 0.5 * _JRange.first       << ", "
	    << 0.5 * _JRange.second       << "]" << endl
	    << "    L range .................... [" << 0.5 * _LRange.first       << ", "
	    << 0.5 * _LRange.second       << "]" << endl
	    << "    S range .................... [" << 0.5 * _SRange.first       << ", "
	    << 0.5 * _SRange.second       << "]" << endl
	    << "    allow JPC exotics .......... " << ((_allowJpcExotics)      ? "true" : "false") << endl
	    << "    require min. isobar mass ... " << ((_requireMinIsobarMass) ? "true" : "false") << endl
	    << "    isobar mass window par. .... " << _isobarMassWindowSigma << " [Gamma]" << endl;
	out << "    isobar black list:";
	if (_isobarBlackList.size() == 0)
		out << " empty" << endl;
	else {
		out << endl;
		for (size_t i = 0; i < _isobarBlackList.size(); ++i)
			out << "        " << _isobarBlackList[i] << endl;
	}
	out << "    isobar white list:";
	if (_isobarWhiteList.size() == 0)
		out << " all isobars" << endl;
	else {
		out << endl;
		for (size_t i = 0; i < _isobarWhiteList.size(); ++i)
			out << "        " << _isobarWhiteList[i] << endl;
	}
	return out;
}
