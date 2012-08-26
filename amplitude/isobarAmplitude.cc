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
//      virtual base class for isobar decay amplitude independent of formalism
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <algorithm>
#include <cassert>

#include "TLorentzRotation.h"
#include "TMath.h"

#include "conversionUtils.hpp"
#include "factorial.hpp"
#include "isobarAmplitude.h"

  
using namespace std;
using namespace boost;
using namespace rpwa;


bool isobarAmplitude::_debug = false;


isobarAmplitude::isobarAmplitude()
	: _decay               (),
	  _useReflectivityBasis(true),
	  _boseSymmetrize      (true),
	  _isospinSymmetrize   (true),
	  _doSpaceInversion    (false),
	  _doReflection        (false)
{ }


isobarAmplitude::isobarAmplitude(const isobarDecayTopologyPtr& decay)
	: _useReflectivityBasis(true),
	  _boseSymmetrize      (true),
	  _isospinSymmetrize   (true),
	  _doSpaceInversion    (false),
	  _doReflection        (false)
{
	setDecayTopology(decay);
}


isobarAmplitude::~isobarAmplitude()
{ }


void
isobarAmplitude::setDecayTopology(const isobarDecayTopologyPtr& decay)
{
	if (not decay) {
		printErr << "null pointer to decay topology. aborting." << endl;
		throw;
	}
	if (not decay->checkTopology()) {
		printErr << "decay does not have the correct topology. aborting." << endl;
		throw;
	}
	if (not decay->checkConsistency()) {
		printErr << "decay is not consistent. aborting." << endl;
		throw;
	}
	_decay = decay;
}


void
isobarAmplitude::init()
{
	_symTermMaps.clear();
	// create first symmetrization entry with identity permutation map
	vector<unsigned int> identityPermMap;
	for (unsigned int i = 0; i < _decay->nmbFsParticles() ; ++i)
		identityPermMap.push_back(i);
	_symTermMaps.push_back(symTermMap(1, identityPermMap));
	// create final-state symmetriztion terms
	if (_boseSymmetrize)
		initBoseSymTermMaps();
	if (_isospinSymmetrize)
		initIsospinSymTermMaps();
}


// at the moment it is not possible to fix the helicity state of
// final state particles for cases where the amplitudes for the
// different helicity states have to be added incoherently
complex<double>
isobarAmplitude::amplitude() const
{
	const unsigned int nmbSymTerms = _symTermMaps.size();
	if (nmbSymTerms < 1) {
		printErr << "array of symmetrization terms is empty. make sure isobarAmplitude::init() "
		         << "was called. cannot calculate amplitude. returning 0." << endl;
		return 0;
	}
	// loop over all symmetrization terms; assumes that init() was called before
	complex<double> amp = 0;
	for (unsigned int i = 0; i < nmbSymTerms; ++i)
		amp += _symTermMaps[i].factor * symTermAmp(_symTermMaps[i].fsPartPermMap);
	return amp;
}


TLorentzRotation
isobarAmplitude::gjTransform(const TLorentzVector& beamLv,  // beam Lorentz-vector
                             const TLorentzVector& XLv)     // X  Lorentz-vector
{
	TLorentzVector beam    = beamLv;
	TLorentzVector X       = XLv;
	const TVector3 yGjAxis = beam.Vect().Cross(X.Vect());  // y-axis of Gottfried-Jackson frame
	// rotate so that yGjAxis becomes parallel to y-axis and beam momentum ends up in (x, z)-plane
	TRotation rot1;
	rot1.RotateZ(-yGjAxis.Phi());
	rot1.RotateY(piHalf - yGjAxis.Theta());
	rot1.RotateZ(piHalf);
	beam *= rot1;
	X    *= rot1;
	// boost to X RF
	TLorentzRotation boost;
	boost.Boost(-X.BoostVector());
	beam *= boost;
	// rotate about yGjAxis so that beam momentum is along z-axis
	TRotation rot2;
	rot2.RotateY(-signum(beam.X()) * beam.Theta());
	// construct total transformation
	TLorentzRotation gjTransform(rot1);
	gjTransform.Transform(boost);
	gjTransform.Transform(rot2);
	return gjTransform;
}


void
isobarAmplitude::spaceInvertDecay() const
{
	if (_debug)
		printDebug << "space inverting final state momenta." << endl;
	// transform final state particles into X rest frame
	const TLorentzVector&  beamLv  = _decay->productionVertex()->referenceLzVec();
	const TLorentzVector&  XLv     = _decay->XParticle()->lzVec();
	const TLorentzRotation gjTrans = gjTransform(beamLv, XLv);
	_decay->transformFsParticles(gjTrans);
	// perform parity transformation on final state particles in X rest frame
	for (unsigned int i = 0; i < _decay->nmbFsParticles(); ++i) {
		const particlePtr&    part     = _decay->fsParticles()[i];
		const TLorentzVector& fsPartLv = part->lzVec();
		part->setLzVec(TLorentzVector(-fsPartLv.Vect(), fsPartLv.E()));
	}
	// transform final state particles back to lab frame
	_decay->transformFsParticles(gjTrans.Inverse());
}


void
isobarAmplitude::reflectDecay() const
{
	if (_debug)
		printDebug << "reflecting final state momenta through production plane." << endl;
	// transform final state particles into X rest frame
	const TLorentzVector&  beamLv  = _decay->productionVertex()->referenceLzVec();
	const TLorentzVector&  XLv     = _decay->XParticle()->lzVec();
	const TLorentzRotation gjTrans = gjTransform(beamLv, XLv);
	_decay->transformFsParticles(gjTrans);
	// reflect final state particles through production plane
	for (unsigned int i = 0; i < _decay->nmbFsParticles(); ++i) {
		const particlePtr&    part     = _decay->fsParticles()[i];
		const TLorentzVector& fsPartLv = part->lzVec();
		part->setLzVec(TLorentzVector(fsPartLv.X(), -fsPartLv.Y(), fsPartLv.Z(), fsPartLv.E()));
	}
	// transform final state particles back to lab frame
	_decay->transformFsParticles(gjTrans.Inverse());
}


// assumes that all particles in the decay are mesons
// !!! in this primitive recursion scheme daughter amplitudes might be
//     called multiple times; for amplitudes with high-spin isobars a
//     large speed gain might be achievable by precalculating the
//     amplitudes at each vertex for all possible helicity values;
//     additional speed can be gained by checking that daughter
//     amplitudes are not zero before calculating the remaining parts
//     of the amplitude
complex<double>
isobarAmplitude::twoBodyDecayAmplitudeSum(const isobarDecayVertexPtr& vertex,           // current vertex
                                          const bool                  topVertex) const  // switches special treatment of X decay vertex; needed for reflectivity basis
{
	const particlePtr& parent    = vertex->parent();
	const particlePtr& daughter1 = vertex->daughter1();
	const particlePtr& daughter2 = vertex->daughter2();
	complex<double>    ampSum    = 0;
	for (int lambda1 = -daughter1->J(); lambda1 <= +daughter1->J(); lambda1 += 2) {
		// calculate decay amplitude for daughter 1
		daughter1->setSpinProj(lambda1);
		const isobarDecayVertexPtr& daughter1Vertex =
			dynamic_pointer_cast<isobarDecayVertex>(_decay->toVertex(daughter1));
		complex<double> daughter1Amp = 0;
		if (daughter1Vertex) {
			daughter1Amp = twoBodyDecayAmplitudeSum(daughter1Vertex, false);
		} else
			daughter1Amp = 1;
		for (int lambda2 = -daughter2->J(); lambda2 <= +daughter2->J(); lambda2 += 2) {
			// calculate decay amplitude for daughter 2
			daughter2->setSpinProj(lambda2);
			const isobarDecayVertexPtr& daughter2Vertex =
				dynamic_pointer_cast<isobarDecayVertex>(_decay->toVertex(daughter2));
			complex<double> daughter2Amp = 0;
			if (daughter2Vertex) {
				daughter2Amp = twoBodyDecayAmplitudeSum(daughter2Vertex, false);
			} else
				daughter2Amp = 1;
			complex<double> parentAmp = twoBodyDecayAmplitude(vertex, topVertex);
			complex<double> amp       = parentAmp * daughter1Amp * daughter2Amp;
			if (_debug)
				printDebug << "amplitude term for : "
				           << parent->name()    << " [lambda = " << spinQn(parent->spinProj()) << "] -> "
				           << daughter1->name() << " [lambda = " << spinQn(lambda1           ) << "] + "
				           << daughter2->name() << " [lambda = " << spinQn(lambda2           ) << "] = "
				           << "parent amp. = " << maxPrecisionDouble(parentAmp)
				           << " * daughter_1 amp = " << maxPrecisionDouble(daughter1Amp)
				           << " * daughter_2 amp = " << maxPrecisionDouble(daughter2Amp) << " = "
				           << maxPrecisionDouble(amp) << endl;
			ampSum += amp;
		}
	}
	if (_debug)
		printDebug << "decay amplitude for " << *vertex << " = " << maxPrecisionDouble(ampSum) << endl;
	return ampSum;
}


complex<double>
isobarAmplitude::symTermAmp(const std::vector<unsigned int>& fsPartPermMap) const
{
	// (re)set final state momenta
	if (not _decay->revertMomenta(fsPartPermMap)) {
		printErr << "problems reverting momenta in decay topology. cannot calculate amplitude. "
		         << "returning 0." << endl;
		return 0;
	}
	// transform daughters into their respective RFs
	transformDaughters();
	// calculate amplitude
	return twoBodyDecayAmplitudeSum(_decay->XIsobarDecayVertex(), true);
}


void
isobarAmplitude::genBoseSymTermMaps
(const map<string, vector<unsigned int> >&     origFsPartIndices,      // original final-state particle ordering sorted by species
 const map<string, vector<unsigned int> >&     newFsPartIndices,       // new final-state particle ordering sorted by species
 map<string, vector<unsigned int> >::iterator& newFsPartIndicesEntry,  // entry of particle species that is currently being symmetrized
 const vector<unsigned int>&                   baseFsPartPermMap,      // final-state permutation map that is permutated
 vector<symTermMap>&                           symTermMaps) const      // generated permutation maps
{
	// loop over all permutations for current final-state particle species
	do {
		map<string, vector<unsigned int> >::iterator nextFsPartIndicesEntry = newFsPartIndicesEntry;
		if (++nextFsPartIndicesEntry != newFsPartIndices.end())
			// recurse to permutations of other final-state particle species
			genBoseSymTermMaps(origFsPartIndices, newFsPartIndices, nextFsPartIndicesEntry,
			                   baseFsPartPermMap, symTermMaps);
		else {
			// build map for current permutation
			vector<unsigned int> permMap(_decay->nmbFsParticles(), 0);
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
					permMap[origFsPartIndex] = newFsPartIndex;
				}
			}
			if (_debug)
				cout << endl;
			// compute effective permutation map taking into account base permutation map
			assert(baseFsPartPermMap.size() > 0);
			if (_debug) {
				printDebug << "effective Bose-symmetrization final-state permutation map for "
				           << "base permutation map " << baseFsPartPermMap << " is ";
			}
			vector<unsigned int> fsPartPermMap(permMap.size(), 0);
			assert(fsPartPermMap.size() == baseFsPartPermMap.size());
			for (unsigned int i = 0; i < fsPartPermMap.size(); ++i) {
				const unsigned int newIndex = baseFsPartPermMap[i];
				fsPartPermMap[i] = permMap[newIndex];
			}
			if (_debug)
				cout << fsPartPermMap << endl;
			symTermMap symTerm(1, fsPartPermMap);
			symTermMaps.push_back(symTerm);
		}
	} while (next_permutation(newFsPartIndicesEntry->second.begin(),
	                          newFsPartIndicesEntry->second.end()));
}


void
isobarAmplitude::initBoseSymTermMaps()
{
	// get final state indistinguishable particles
	typedef map<string, unsigned int>::const_iterator indistFsPartIt;
	const map<string, unsigned int> indistFsPart = _decay->nmbIndistFsParticles();
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
	for (unsigned int i = 0; i < _decay->nmbFsParticles(); ++i) {
		const string partName = _decay->fsParticles()[i]->name();
		origFsPartIndices[partName].push_back(i);
	}
	map<string, vector<unsigned int> > newFsPartIndices = origFsPartIndices;

	// generate new Bose-permutation maps for all existing final-state maps
	const unsigned int nmbSymTerms = _symTermMaps.size();
	vector<symTermMap> newSymTermMaps[nmbSymTerms];
	for (unsigned int i = 0; i < nmbSymTerms; ++i) {
		// recursively create permutation maps for all Bose-symmetrization terms
		printInfo << "generating permutation maps for " << nmbCombinations << " "
		          << "Bose-symmetrization terms with base symmetrization term [" << i << " / "
		          << nmbSymTerms << "]: factor = " << maxPrecisionDouble(_symTermMaps[i].factor) << "; "
		          << "final-state permutation map = " << _symTermMaps[i].fsPartPermMap << endl;
		
		map<string, vector<unsigned int> >::iterator firstEntry        = newFsPartIndices.begin();
		const vector<unsigned int>&                  baseFsPartPermMap = _symTermMaps[i].fsPartPermMap;
		genBoseSymTermMaps(origFsPartIndices, newFsPartIndices, firstEntry,
		                   baseFsPartPermMap, newSymTermMaps[i]);
		// set factor for all terms
		for (unsigned int j = 0; j < newSymTermMaps[i].size(); ++j)
			newSymTermMaps[i][j].factor = normFactor * _symTermMaps[i].factor;
	}

	// rebuild array of permutation maps for final-state symmetrization with new terms
	_symTermMaps.clear();
	for (unsigned int i = 0; i < nmbSymTerms; ++i)
		for (unsigned int j = 0; j < newSymTermMaps[i].size(); ++j)
			_symTermMaps.push_back(newSymTermMaps[i][j]);
}


void
isobarAmplitude::initIsospinSymTermMaps()
{
	// vector<symTermMap> symMaps = _decay->getIsospinSymmetrization();
	// const unsigned int nmbSymTerms = symMaps.size();
	// if (nmbSymTerms < 1) {
	// 	printErr << "array of isospin symmetrization terms is empty. cannot calculate amplitude. "
	// 	         << "returning 0." << endl;
	// 	return 0;
	// }
	// if (nmbSymTerms == 1) {
	// 	printDebug << "no isospin symmetrization needed for this amplitude." << endl;
	// 	return unSymmetrizedAmp();
	// }
	// // if (_debug)
	// 	printDebug << "symmetrizing amplitude w.r.t. isospin using " << nmbSymTerms << " terms." << endl;
	// complex<double> amp = 0;
	// for(unsigned int i = 0; i < nmbSymTerms; ++i) {
	// 	const double                clebsch       = symMaps[i].factor.real();
	// 	const vector<unsigned int>& fsPartPermMap = symMaps[i].fsPartPermMap;
	// 	// if (_debug) {
	// 		printDebug << "calculating amplitude for isospin-symmetrization term "
	// 		           << "with Clebsch-Gordan coefficient = " << clebsch << " "
	// 		           << "and final state permutation ";
	// 		for (unsigned int j = 0; j < fsPartPermMap.size(); ++j)
	// 			cout << _decay->fsParticles()[j]->name() << "[" << j << " -> " << fsPartPermMap[j] << "]  ";
	// 		cout << endl;
	// 	// }
	// 	complex<double> foo = unSymmetrizedAmp(fsPartPermMap);
	// 	cout << "    term[" << i << "] = " << maxPrecisionDouble(foo) << ", clebsch = " << clebsch
	// 	     << ", " << signum(clebsch) << endl;
	// 	amp += signum(clebsch) * foo;
	// 	// amp += signum(clebsch) * unSymmetrizedAmp(fsPartPermMap);
	// }
	// return amp / sqrt(nmbSymTerms);
}


ostream&
isobarAmplitude::printParameters(ostream& out) const
{
	out << name() << ": " << endl
	    << "    reflectivity basis ............... " << enDisabled(_useReflectivityBasis) << endl
	    << "    Bose symmetrization .............. " << enDisabled(_boseSymmetrize      ) << endl
	    << "    isospin symmetrization ........... " << enDisabled(_isospinSymmetrize   ) << endl
	    << "    space inversion of FS momenta .... " << enDisabled(_doSpaceInversion    ) << endl
	    << "    reflection through prod. plane ... " << enDisabled(_doReflection        ) << endl;
	return out;
}


ostream&
isobarAmplitude::print(ostream& out) const
{
	printParameters(out);
	out << *_decay;
	return out;
}
