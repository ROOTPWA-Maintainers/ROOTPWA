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
//      class that generates all possible waves given the final state
//      particles, a list of isobars and constraints on I, J, L, and S
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <boost/assign.hpp>

#include "libconfig.h++"

#include "physUtils.hpp"
#include "spinUtils.hpp"
#include "conversionUtils.hpp"
#include "libConfigUtils.hpp"
#include "arrayUtils.hpp"
#include "particleDataTable.h"
#include "waveDescription.h"
#include "waveSetGenerator.h"


using namespace std;
using namespace boost;
namespace bt = boost::tuples;
using namespace libconfig;
using namespace rpwa;


bool waveSetGenerator::_debug = false;


waveSetGenerator::waveSetGenerator()
{
	reset();
}


waveSetGenerator::~waveSetGenerator()
{ }


bool
waveSetGenerator::setWaveSetParameters(const string& templateKeyFileName)
{
	// construct template decay topology
	_templateTopo.reset();
	waveDescription waveDesc;
	if (   not waveDesc.parseKeyFile(templateKeyFileName)
	    or not waveDesc.constructDecayTopology(_templateTopo, true)) {
		printWarn << "problems constructing template decay topology from key file "
		          << "'" << templateKeyFileName << "'. cannot generate wave set." << endl;
		return false;
	}
	// read wave set parameters from template key file
	printInfo << "reading wave set parameters from key file '" << templateKeyFileName << "'" << endl;
	Config key;
	if (not parseLibConfigFile(templateKeyFileName, key, _debug)) {
		printWarn << "problems reading wave set parameters from '" << templateKeyFileName << "'. "
		          << "cannot generate wave set." << endl;
		return false;
	}
	// find and parse group with wave set parameters
	const Setting* waveSetParKey = findLibConfigGroup(key.getRoot(), "waveSetParameters", false);
	if (waveSetParKey) {
		const Setting* isoSpinRangeKey = findLibConfigArray(*waveSetParKey, "isospinRange", false);
		if (isoSpinRangeKey)
			_isospinRange = bt::make_tuple<int, int>((*isoSpinRangeKey)[0], (*isoSpinRangeKey)[1]);
		const Setting* JRangeKey = findLibConfigArray(*waveSetParKey, "JRange", false);
		if (JRangeKey)
			_JRange = bt::make_tuple<int, int>((*JRangeKey)[0], (*JRangeKey)[1]);
		const Setting* MRangeKey = findLibConfigArray(*waveSetParKey, "MRange", false);
		waveSetParKey->lookupValue("reflectivity",     _reflectivity    );
		waveSetParKey->lookupValue("useReflectivity",  _useReflectivity );
		waveSetParKey->lookupValue("allowSpinExotics", _allowSpinExotics);
		if (MRangeKey)
			_spinProjRange = bt::make_tuple<int, int>((*MRangeKey)[0], (*MRangeKey)[1]);
		const Setting* LRangeKey = findLibConfigArray(*waveSetParKey, "LRange", false);
		if (LRangeKey)
			_LRange = bt::make_tuple<int, int>((*LRangeKey)[0], (*LRangeKey)[1]);
		const Setting* SRangeKey = findLibConfigArray(*waveSetParKey, "SRange", false);
		if (SRangeKey)
			_SRange = bt::make_tuple<int, int>((*SRangeKey)[0], (*SRangeKey)[1]);
		const Setting* isobarBlackListKey = findLibConfigArray(*waveSetParKey,
		                                                       "isobarBlackList", false);
		if (isobarBlackListKey) {
			_isobarBlackList.clear();
			for (int i = 0; i < isobarBlackListKey->getLength(); ++i)
				_isobarBlackList.push_back((*isobarBlackListKey)[i]);
		}
		const Setting* isobarWhiteListKey = findLibConfigArray(*waveSetParKey,
		                                                       "isobarWhiteList", false);
		if (isobarWhiteListKey) {
			_isobarWhiteList.clear();
			for (int i = 0; i < isobarWhiteListKey->getLength(); ++i)
				_isobarWhiteList.push_back((*isobarWhiteListKey)[i]);
		}
		waveSetParKey->lookupValue("requireMinIsobarMass",  _requireMinIsobarMass );
		waveSetParKey->lookupValue("isobarMassWindowSigma", _isobarMassWindowSigma);
	}
	if (_debug)
		printDebug << "parameters of " << *this;
	return true;
}


size_t
waveSetGenerator::generateWaveSet()
{
	if (not _templateTopo) {
		printWarn << "template decay topology was not yet constructed. "
		          << "did you call waveSetGenerator::setWaveSetParameters()?" << endl;
		return 0;
	}
	if (not _templateTopo->checkTopology()) {
		printWarn << "ill-formed template decay topology. cannot generate wave set." << endl;
		return 0;
	}
	// order nodes depth-first
	vector<nodeDesc> startNds
		= _templateTopo->sortNodesDfs(_templateTopo->XIsobarDecayVertex());
	// create decay topologies of all subdecays
	vector<isobarDecayTopology> subDecays(startNds.size());
	for (size_t i = 0; i < startNds.size(); ++i)
		subDecays[i] = _templateTopo->subDecay(startNds[i]);
	// reverse order, because decay trees are build starting from the final state nodes
	reverse(startNds.begin (), startNds.end ());
	reverse(subDecays.begin(), subDecays.end());

	// create all possible subdecays
	map<nodeDesc, vector<isobarDecayTopology> > decayPossibilities;  // all possible decay graphs starting at certain node
	for (size_t iStart = 0; iStart < startNds.size(); ++iStart) {

		// special case for final state nodes
		if (_templateTopo->isFsVertex(_templateTopo->vertex(startNds[iStart]))) {
			// final state vertices have no daughter tracks; subdecay
			// topology consist of just the final state vertex
			decayPossibilities[startNds[iStart]].push_back(subDecays[iStart]);
			continue;
		}

		if (_debug)
			printDebug << "generating decay topologies for subdecay[" << iStart << "]: "
			           << subDecays[iStart];

		// get daughter decay topologies
		vector<vector<isobarDecayTopology>* > daughterDecays;
		adjIterator iNd, iNdEnd;
		for (bt::tie(iNd, iNdEnd) = _templateTopo->adjacentVertices(startNds[iStart]); iNd != iNdEnd; ++iNd)
			daughterDecays.push_back(&decayPossibilities[*iNd]);

		// loop over all combinations of daughter decays
		size_t iDaughter[2];
		for (iDaughter[0] = 0; iDaughter[0] < daughterDecays[0]->size(); ++iDaughter[0])
			for (iDaughter[1] = 0; iDaughter[1] < daughterDecays[1]->size(); ++iDaughter[1]) {

				// get daughter particles from the respective decay topologies
				const nodeDesc topNodes[2] =
					{(*daughterDecays[0])[iDaughter[0]].topNode(),
					 (*daughterDecays[1])[iDaughter[1]].topNode()};
				const particlePtr daughters[2] =
					{(*daughterDecays[0])[iDaughter[0]].vertex(topNodes[0])->inParticles()[0],
					 (*daughterDecays[1])[iDaughter[1]].vertex(topNodes[1])->inParticles()[0]};
				if (_debug)
					printDebug << "combining decay[" << iDaughter[0] << "] of first daughter "
					           << "'" << daughters[0]->name() << "' "
					           << "with decay[" << iDaughter[1] << "] of second daughter "
					           << "'" << daughters[1]->name() << "'" << endl
					           << (*daughterDecays[0])[iDaughter[0]]
					           << (*daughterDecays[1])[iDaughter[1]];

				// copy parent vertex
				isobarDecayVertexPtr parentVertex
					(new isobarDecayVertex(*static_pointer_cast<isobarDecayVertex>
					                       (_templateTopo->vertex(startNds[iStart]))));
				particlePtr parent = parentVertex->parent();
				// set daughters of parent vertex
				parentVertex->daughter1() = daughters[0];
				parentVertex->daughter2() = daughters[1];
				// join daughter subdecays and parent vertex
				isobarDecayTopology parentDecay
					= _templateTopo->joinDaughterDecays(parentVertex,
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
				int S, SMax;
				for (bt::tie(S, SMax) = getSpinRange(daughterJ[0], daughterJ[1], _SRange); S <= SMax; S += 2) {
					if (_debug)
						printDebug << "setting total spin = " << spinQn(S) << " "
						           << "(max spin = " << spinQn(SMax) << ")" << endl;

					// loop over allowed relative orbital angular momenta
					for (int L = max(0, bt::get<0>(_LRange)); L <= bt::get<1>(_LRange); L += 2) {
						if (_debug)
							printDebug << "setting relative orbital angular momentum = " << spinQn(L) << " "
							           << "(max L = " << spinQn(bt::get<1>(_LRange)) << ")" << endl;
						const int parentP = daughters[0]->P() * daughters[1]->P() * powMinusOne(L / 2);

						// loop over allowed total angular momenta
						int parentJ, parentJMax;
						for (bt::tie(parentJ, parentJMax) = getSpinRange(L, S, _JRange);
						     parentJ <= parentJMax; parentJ += 2) {
							if (_debug)
								printDebug << "setting parent spin = " << spinQn(parentJ) << " "
								           << "(max spin = " << spinQn(parentJMax) << ")" << endl;

							// loop over allowed isospins
							int parentI, parentIMax;
							for (bt::tie(parentI, parentIMax) = getSpinRange(daughterI[0], daughterI[1], _isospinRange);
							     parentI <= parentIMax; parentI += 2) {
								if (_debug)
									printDebug << "setting parent isospin = " << spinQn(parentI) << " "
									           << "(max isospin = " << spinQn(parentIMax) << ")" << endl;

								// get I_z from Gell-Mann-Nishijima formula (see PDG 2008 eq. 14.1)
								const int parentI_z = 2 * parentCharge - (parentBaryonNmb + parentStrangeness
								                                          + parentCharm + parentBeauty);

								// C-parity
								// defined only for mesons with zero flavor quantum numbers
								int parentC = 0;
								if ((parentBaryonNmb == 0) and (parentI_z == 0) and (parentStrangeness == 0)
								    and (parentCharm == 0) and (parentBeauty == 0))
									parentC = parentG * powMinusOne(parentI / 2);
								if (_debug)
									printDebug << "trying isobar quantum numbers IG(JPC) = "
									           << spinQn(parentI) << sign(parentG)
									           << "(" << spinQn(parentJ) << sign(parentP) << sign(parentC) << ")"
									           << ", with L = " << spinQn(L) << ", S = " << spinQn(S) << endl;

								// check whether charge state is allowed using
								if (not spinStatesCanCouple(daughterI, daughterI_z, parentI, parentI_z)) {
									if (_debug)
										cout << "        rejected, because (I, I_z)[0] = (" << spinQn(daughterI[0])
										     << ", " << spinQn(daughterI_z[0]) << ") and (I, I_z)[1] = ("
										     << spinQn(daughterI[1]) << ", " << spinQn(daughterI_z[1])
										     << ") cannot couple to parent (I, I_z) = ("
										     << spinQn(parentI) << ", " << spinQn(parentI_z) << ")" << endl;
									continue;
								}
								if (not _allowSpinExotics)
									// check whether quantum number combination is spin-exotic
									if (igjpIsExotic(parentI, parentG, parentJ, parentP)) {
										if (_debug)
											cout << "        rejected, because quantum number combination IG(JPC) = "
											     << spinQn(parentI) << sign(parentG) << "(" << spinQn(parentJ)
											     << sign(parentP) << sign(parentC) << ") is spin-exotic" << endl;
										continue;
									}

								// enforce rules for particle-antiparticle systems
								if (   daughters[1]->antiPartProperties()
								    == *(static_pointer_cast<particleProperties>(daughters[0]))) {
									// !NOTE! since equality test function is virtual,
									// order in the above comparison is essential
									// (see S.U. Chung: "C- and G-Parity: a New Definition and Applications, Version V"
									//  eq. 32)
									if (parentG != powMinusOne((parentI + L + S) / 2)) {
										if (_debug)
											cout << "        rejected, because I[= " << spinQn(parentI) << "], "
											     << "L[= " << spinQn(L) << "], and S[= " << spinQn(S) << "] "
											     << "are not consistent with G-parity[= " << sign(parentG) << "]. "
											     << "for particle-antiparticle system ["
											     << daughters[0]->name() << ", " << daughters[1]->name() << "] "
											     << "(-)^(I + L + S) = " << sign(powMinusOne((parentI + L + S) / 2))
											     << " must be equal to G-parity." << endl;
										continue;
									}
									// I_z of the parent and C = (-)^(L + S) should be fulfilled automatically
									// (see S.U. Chung: "C- and G-Parity: a New Definition and Applications, Version V"
									//  eq. 31)
									assert(parentI_z == 0);
									assert(parentC == powMinusOne((L + S) / 2));
								}

								// enforce rules for particle-particle systems
								// (see S.U. Chung: "C- and G-Parity: a New Definition and Applications, Version V"
								//  eq. 36)
								if (*(daughters[0]) == *(daughters[1]))
									if (powMinusOne((parentI + L + S) / 2) != powMinusOne(daughterI[0])) {
										if (_debug)
											cout << "        rejected, because system with two identical particles ["
											     << daughters[0]->name() << ", " << daughters[1]->name() << "] must have "
											     << ((powMinusOne(daughterI[0]) == +1) ? "even" : "odd")
											     << " (I[= " << spinQn(parentI) << "] + L[= " << spinQn(L) << "] "
											     << "+ S[= " << spinQn(S) << "]) = "
											     << spinQn(parentI) + spinQn(L) + spinQn(S) << endl;
										continue;
									}

								// set isobar properties
								particleProperties isobarProp(*parent);
								isobarProp.setBaryonNmb(parentBaryonNmb);
								isobarProp.setSCB      (parentStrangeness, parentCharm, parentBeauty);
								isobarProp.setIGJPC    (parentI, parentG, parentJ, parentP, parentC);

								// create all allowed decay topologies
								if (parent->isXParticle()) {
									// for X loop over allowed M and reflectivity states
									bt::tuple<int, int> parentMRange =
										bt::make_tuple(spinQnLargerEqual((_useReflectivity) ?
										                                 spinQnSmallerEqual(parentJ, 0) : -parentJ,
										                                 bt::get<0>(_spinProjRange)),
										               spinQnSmallerEqual(parentJ, bt::get<1>(_spinProjRange)));
									bt::tuple<int, int> parentReflRange(0, 0);
									if (_useReflectivity) {
										if (_reflectivity == 0)
											parentReflRange = bt::make_tuple(-1, +1);
										else
											parentReflRange = bt::make_tuple(signum(_reflectivity), signum(_reflectivity));
									}
									for (int parentM = bt::get<0>(parentMRange); parentM <= bt::get<1>(parentMRange);
									     parentM += 2) {
										if (_debug)
											printDebug << "setting X spin projection to " << spinQn(parentM) << " "
											           << "(max spin projection = " << spinQn(bt::get<1>(parentMRange))
											           << ")" << endl;
										for (int parentRefl = bt::get<0>(parentReflRange);
										     parentRefl <= bt::get<1>(parentReflRange); parentRefl += 2) {
											if (_debug and (parentRefl != 0))
												printDebug << "setting X reflectivity to "
												           << parityQn(parentRefl) << endl;

											// check whether amplitude would be zero for given reflectivity
											if (    _useReflectivity and (parentM == 0)
											    and (reflectivityFactor(parentJ, parentP, parentM, parentRefl) == +1)) {
												if (_debug)
													printDebug << "rejected, because for quantum number combination "
													           << "JP(M refl) = " << spinQn(parentJ) << sign(parentP) << "("
													           << spinQn(parentM) << sign(parentRefl) << ") "
													           << "amplitude is always zero"<< endl;
												continue;
											}

											// clone topology
											const isobarDecayTopologyPtr& decayCopy
												= createNewDecayTopology(parentDecay, parentVertex, L, S,
												                         isobarProp, parentCharge);
											decayCopy->XDecayVertex()->inParticles()[0]->setSpinProj(parentM);
											decayCopy->XDecayVertex()->inParticles()[0]->setReflectivity(parentRefl);
											decayPossibilities[startNds[iStart]].push_back(*decayCopy);
										}
									}
								} else {
									// for intermediate decay vertices find isobar candidates in particle data table
									const double minIsobarMass = (_requireMinIsobarMass) ?
										  (daughters[0]->mass() - _isobarMassWindowSigma * daughters[0]->width())
										+ (daughters[1]->mass() - _isobarMassWindowSigma * daughters[1]->width())
										: 0;

									const particleProperties::decayMode
										decay(assign::list_of(daughters[0]->name())(daughters[1]->name()), L, S);
									vector<const particleProperties*> possibleIsobars
									    = particleDataTable::entriesMatching(isobarProp, "allQn",
									                                         minIsobarMass, _isobarMassWindowSigma,
									                                         _isobarWhiteList, _isobarBlackList,
									                                         decay, _forceDecayCheck);
									if (_debug) {
										if (possibleIsobars.empty())
											printDebug << "rejected, because no matching isobar could be found for "
											           << isobarProp.qnSummary() << " in particle data table" << endl;
										else
											printDebug << "found " << possibleIsobars.size() << " isobar candidate(s) for "
											           << isobarProp.qnSummary() << " in particle data table" << endl;
									}

									// loop over isobar candidates
									for (size_t iIsobar = 0; iIsobar < possibleIsobars.size(); ++iIsobar) {
										if (_debug)
											printDebug << "setting isobar = '"
											           << possibleIsobars[iIsobar]->name() << "'" << endl;
										// clone topology
										const isobarDecayTopologyPtr& decayCopy
											= createNewDecayTopology(parentDecay, parentVertex, L, S,
											                         *possibleIsobars[iIsobar], parentCharge);
										decayPossibilities[startNds[iStart]].push_back(*decayCopy);
									}
								}
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
		                                      (_templateTopo->productionVertex()->clone(false, false)));
		const particlePtr& newX = _waveSet[i].XIsobarDecayVertex()->parent();
		newProdVert->outParticles()[0] = newX;
		// add production vertex
		_waveSet[i].setProductionVertex(newProdVert);
		// correct quantum numbers
		_waveSet[i].calcIsobarCharges();
		_waveSet[i].productionVertex()->setXFlavorQN();  // sets baryon nmb, S, C, and B of X
	}

	// post process wave set
	set<size_t> boseSymWaveIndices = findBoseSymDecays();
	if (boseSymWaveIndices.size() > 0) {
		_waveSet.erase(remove_if(_waveSet.begin(), _waveSet.end(),
		                         rpwa::isInListOfIndices<isobarDecayTopology>(_waveSet.begin(),
			boseSymWaveIndices)),
		               _waveSet.end());
		printSucc << "removed " << boseSymWaveIndices.size() << " Bose duplicates "
		          << "from wave set" << endl;
	}
	return _waveSet.size();
}


bool
waveSetGenerator::writeKeyFiles(const string& dirName,
                                const bool    newKeyFileNameConvention)
{
	size_t countSuccess = 0;
	for (size_t i = 0; i < _waveSet.size(); ++i) {
		const string keyFileName = dirName + "/"
			+ waveDescription::waveNameFromTopology(_waveSet[i], newKeyFileNameConvention) + ".key";
		if (waveDescription::writeKeyFile(keyFileName, _waveSet[i]))
			++countSuccess;
	}
	printInfo << "wrote " << countSuccess << " out of " << _waveSet.size() << " key files" << endl;
	if (countSuccess != _waveSet.size()) {
	  printWarn << "writing of " << _waveSet.size() - countSuccess << " key files failed" << endl;
	  return false;
	}
	return true;
}


void
waveSetGenerator::reset()
{
	_isospinRange          = bt::make_tuple(0, 2);
	_JRange                = bt::make_tuple(0, 0);
	_spinProjRange         = bt::make_tuple(0, 0);
	_reflectivity          = 0;
	_useReflectivity       = false;
	_allowSpinExotics      = false;
	_LRange                = bt::make_tuple(0, 0);
	_SRange                = bt::make_tuple(0, 0);
	_requireMinIsobarMass  = false;
	_forceDecayCheck       = false;
	_isobarMassWindowSigma = 0;
	_isobarBlackList.clear();
	_isobarWhiteList.clear();
	_templateTopo.reset();
	_waveSet.clear();
}


ostream&
waveSetGenerator::print(ostream& out) const
{
	out << "wave set generator:" << endl
	    << "    isospin range .............. [" << spinQn(bt::get<0>(_isospinRange))  << ", "
	    << spinQn(bt::get<1>(_isospinRange))  << "]" << endl
	    << "    J range .................... [" << spinQn(bt::get<0>(_JRange))        << ", "
	    << spinQn(bt::get<1>(_JRange))        << "]" << endl
	    << "    M range .................... [" << spinQn(bt::get<0>(_spinProjRange)) << ", "
	    << spinQn(bt::get<1>(_spinProjRange)) << "]" << endl
	    << "    allow spin-exotics ......... " << yesNo(_allowSpinExotics) << endl;
	if (_useReflectivity)
		out << "    generate reflectivity ...... " << ((_reflectivity == 0) ? "both"
		                                               : sign(_reflectivity)) << endl;
	out << "    L range .................... [" << spinQn(bt::get<0>(_LRange))        << ", "
	    << spinQn(bt::get<1>(_LRange))        << "]" << endl
	    << "    S range .................... [" << spinQn(bt::get<0>(_SRange))        << ", "
	    << spinQn(bt::get<1>(_SRange))        << "]" << endl
	    << "    require min. isobar mass ... " << yesNo(_requireMinIsobarMass) << endl
	    << "    isobar mass window par. .... " << _isobarMassWindowSigma << " [Gamma]" << endl
	    << "    force decay checks ......... " << yesNo(_forceDecayCheck) << endl;
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
	if (_templateTopo)
		out << "template " << *_templateTopo;
	return out;
}


const isobarDecayTopologyPtr
waveSetGenerator::createNewDecayTopology(const isobarDecayTopology&  parentDecay,
                                         const isobarDecayVertexPtr& parentVertex,
                                         const int                   L,
                                         const int                   S,
                                         const particleProperties&   isobar,
                                         const int                   parentCharge)
{
	// clone topology
	const isobarDecayTopologyPtr& decayCopy = parentDecay.clone(false, false);
	// set parent vertex quantum numbers
	const isobarDecayVertexPtr& parentVertexCopy
		= static_pointer_cast<isobarDecayVertex>(decayCopy->vertex(parentDecay.node(parentVertex)));
	parentVertexCopy->setL(L);
	parentVertexCopy->setS(S);
	// set parent particle
	const particlePtr& parentCopy = parentVertexCopy->parent();
	parentCopy->setProperties(isobar);
	parentCopy->setCharge(parentCharge);
	if (_debug)
		printDebug << "created decay topology for " << *parentVertexCopy << endl;
	return decayCopy;
}


set<size_t>
waveSetGenerator::findBoseSymDecays() const
{
	printInfo << "looking for Bose duplicates in wave set..." << endl;
	set<size_t> waveIndicesToRemove;
	set<size_t> waveIndicesToKeep;
	for (size_t waveIndex = 0; waveIndex < _waveSet.size(); ++waveIndex) {

		// don't consider entries that are already tagged
		set<size_t>::const_iterator indexEntry = waveIndicesToRemove.find(waveIndex);
		if (indexEntry != waveIndicesToRemove.end())
			continue;
		indexEntry = waveIndicesToKeep.find(waveIndex);
		if (indexEntry != waveIndicesToKeep.end())
			continue;

		// get Bose-symmetric vertices
		isobarDecayTopologyPtr     topo              = _waveSet[waveIndex].clone();
		const vector<unsigned int> boseSymVertIds    = topo->findIsobarBoseSymVertices();
		const size_t               nmbBoseSymVertIds = boseSymVertIds.size();
		if (nmbBoseSymVertIds > 1) {
			printErr << "the code has not yet been tested for decay topologies with "
			         << nmbBoseSymVertIds << " Bose-symmetric decay vertices. Aborting..." << endl;
			throw;
		}
		if (nmbBoseSymVertIds != 1)
			continue;

		// make sure the two daughters are not the same particle
		isobarDecayVertexPtr vert         = topo->isobarDecayVertices()[boseSymVertIds[0]];
		const particlePtr    daughters[2] = {vert->daughter1(), vert->daughter2()};
		if (*(daughters[0]) == *(daughters[1]))
			continue;
		const decayTopology daughterTopos[2] = {topo->subDecay(topo->toNode(daughters[0])),
		                                        topo->subDecay(topo->toNode(daughters[1]))};
		if (_debug)
			printDebug << "wave[" << waveIndex << "] "
			           << waveDescription::waveNameFromTopology(*topo, true)
			           << " has Bose-symmetric vertex " << *vert << ". "
			           << "looking for partner waves in wave set." << endl;

		// find indices of all decay vertices that are not below the Bose-symmetric vertex
		vector<unsigned int> nonSymVertIds;
		{
			// get graph nodes below Bose-symmetric vertex
			const vector<nodeDesc> nodes = topo->sortNodesDfs(vert);
			for (unsigned int i = 0; i < topo->nmbDecayVertices(); ++i) {
				const nodeDesc nd        = topo->node(topo->isobarDecayVertices()[i]);
				bool           isSymNode = false;
				for (unsigned int j = 0; j < nodes.size(); ++j)
					if (nd == nodes[j]) {
						isSymNode = true;
						break;
					}
				if (not isSymNode)
					nonSymVertIds.push_back(i);
			}
		}

		// look for Bose-symmetric topology
		// this uses the fact that all waves have the same decay toplogy,
		// because they were generated from the same template
		for (size_t i = 0; i < _waveSet.size(); ++i) {
			if (i == waveIndex)
				continue;
			isobarDecayTopologyPtr symTopo = _waveSet[i].clone();

			// compare vertices that are not below the Bose-symmetric vertex
			bool foundSymTopo = true;
			for (size_t j = 0; j < nonSymVertIds.size(); ++j)
				if (   *(topo->isobarDecayVertices   ()[nonSymVertIds[j]])
				    != *(symTopo->isobarDecayVertices()[nonSymVertIds[j]])) {
					foundSymTopo = false;
					break;
				}
			if (not foundSymTopo)
				continue;
			// compare Bose-symmetric vertex
			isobarDecayVertexPtr symVert = symTopo->isobarDecayVertices()[boseSymVertIds[0]];
			// swap daughters 1 and 2
			swap(symVert->outParticles()[0], symVert->outParticles()[1]);
			if (*symVert != *vert)
				continue;
			swap(symVert->outParticles()[1], symVert->outParticles()[0]);  // swap back, because otherwise printout is confusing
			// compare daughter topologies
			const particlePtr   symDaughters    [2] = {symVert->daughter2(), symVert->daughter1()};  // !NOTE! reversed order
			const decayTopology symDaughterTopos[2]
				= {symTopo->subDecay(symTopo->toNode(symDaughters[0])),
				   symTopo->subDecay(symTopo->toNode(symDaughters[1]))};
			for (unsigned int k = 0; k < 2; ++k)
				for (unsigned int l = 0; l < daughterTopos[k].nmbDecayVertices(); ++l)
					if (   *(daughterTopos   [k].decayVertices()[l])
					    != *(symDaughterTopos[k].decayVertices()[l])) {
						foundSymTopo = false;
						break;
					}
			if (not foundSymTopo)
				continue;
			if (_debug)
				printDebug << "found Bose-partner wave[" << i << "] "
				           << waveDescription::waveNameFromTopology(*symTopo, true) << endl;
			// keep topology with the heavier first particle
			if (daughters[0]->mass() > daughters[1]->mass()) {
				if (   not waveIndicesToRemove.insert(i).second
				    or not waveIndicesToKeep.insert(waveIndex).second) {
					printErr << "wave[" << i << "]([" << waveIndex << "]) was tagged for "
					         << "removal(for keeping) twice. this should never happen. Aborting..." << endl;
					throw;
				}
				printInfo << "tagged wave " << waveDescription::waveNameFromTopology(*symTopo, true)
				          << " as superfluous" << endl;
			} else {
				if (   not waveIndicesToRemove.insert(waveIndex).second
				    or not waveIndicesToKeep.insert(i).second) {
					printErr << "wave[" << waveIndex << "]([" << i << "]) was tagged for "
					         << "removal(for keeping) twice. this should never happen. Aborting..." << endl;
					throw;
				}
				printInfo << "tagged wave " << waveDescription::waveNameFromTopology(*topo, true)
				          << " as superfluous" << endl;
			}
		}  // inner loop over wave set

	}  // outer loop over wave set

	printSucc << "found " << waveIndicesToRemove.size() << " Bose duplicates "
	          << "in wave set" << endl;
	return waveIndicesToRemove;
}
