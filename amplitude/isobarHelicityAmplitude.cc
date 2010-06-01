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
//      prototype implementation for general isobar decay amplitude in
//      the helicity framework
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

#include "utilities.h"
#include "isobarHelicityAmplitude.h"

	
using namespace std;
using namespace boost;
using namespace rpwa;


bool isobarHelicityAmplitude::_debug = false;


isobarHelicityAmplitude::isobarHelicityAmplitude()
  : _decay               (),
    _useReflectivityBasis(true),
    _boseSymmetrize      (true),
    _doSpaceInversion    (false),
    _doReflection        (false)
{ }


isobarHelicityAmplitude::isobarHelicityAmplitude(const isobarDecayTopologyPtr& decay)
  : _useReflectivityBasis(true),
    _boseSymmetrize      (true),
    _doSpaceInversion    (false),
    _doReflection        (false)
{
  setDecayTopology(decay);
}


isobarHelicityAmplitude::~isobarHelicityAmplitude()
{ }


void
isobarHelicityAmplitude::setDecayTopology(const isobarDecayTopologyPtr& decay)
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


// at the moment it is not possible to fix the helicity state of
// final state particles for cases where the amplitudes for the
// different helicity states have to be added incoherently
complex<double>
isobarHelicityAmplitude::amplitude() const
{
  // recursively sum over all possible helicities of the decay particles
  complex<double> amp = 0;
  if (_boseSymmetrize)
    amp = boseSymmetrizedAmp();
  else {
    // transform daughters into their respective RFs
    transformDaughters();
    // calculate amplitude
    amp = twoBodyDecayAmplitudeSum(_decay->XIsobarDecayVertex(), true);
  }
  return amp;
}


TLorentzRotation
isobarHelicityAmplitude::hfTransform(const TLorentzVector& daughterLv)
{
  TLorentzVector daughter = daughterLv;
  const TVector3 zAxisMother(0, 0, 1);  // take z-axis as defined in mother frame
  const TVector3 yHfAxis = zAxisMother.Cross(daughter.Vect());  // y-axis of helicity frame
  // rotate so that yHfAxis becomes parallel to y-axis and zHfAxis ends up in (x, z)-plane
  TRotation rot1;
  rot1.RotateZ(-yHfAxis.Phi());
  rot1.RotateY(piHalf - yHfAxis.Theta());
  rot1.RotateZ(piHalf);
  daughter *= rot1;
  // rotate about yHfAxis so that daughter momentum is along z-axis
  TRotation rot2;
  rot2.RotateY(-signum(daughter.X()) * daughter.Theta());
  daughter *= rot2;
  // boost to daughter RF
  rot1.Transform(rot2);
  TLorentzRotation hfTransform(rot1);
  hfTransform.Boost(-daughter.BoostVector());
  return hfTransform;
}


TLorentzRotation
isobarHelicityAmplitude::gjTransform(const TLorentzVector& beamLv,
				     const TLorentzVector& XLv)
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
isobarHelicityAmplitude::spaceInvertDecay() const
{
  if (_debug)
    printInfo << "space inverting final state momenta." << endl;
  // transform final state particles into X rest frame
  const TLorentzVector&  beamLv  = _decay->productionVertex  ()->inParticles()[0]->lzVec();
  const TLorentzVector&  XLv     = _decay->XIsobarDecayVertex()->mother()->lzVec();
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
isobarHelicityAmplitude::reflectDecay() const
{
  if (_debug)
    printInfo << "reflecting final state momenta through production plane." << endl;
  // transform final state particles into X rest frame
  const TLorentzVector&  beamLv  = _decay->productionVertex  ()->inParticles()[0]->lzVec();
  const TLorentzVector&  XLv     = _decay->XIsobarDecayVertex()->mother()->lzVec();
  const TLorentzRotation gjTrans = gjTransform(beamLv, XLv);
  _decay->transformFsParticles(gjTrans);
  // reflect final state particles through production plane
  for (unsigned int i = 0; i < _decay->nmbFsParticles(); ++i) {
    const particlePtr&    part     = _decay->fsParticles()[i];
    const TLorentzVector& fsPartLv = part->lzVec();
    part->setLzVec(TLorentzVector(fsPartLv.X(), -fsPartLv.Y(), -fsPartLv.Z(), fsPartLv.E()));
  }
  // transform final state particles back to lab frame
  _decay->transformFsParticles(gjTrans.Inverse());
}


void
isobarHelicityAmplitude::transformDaughters() const
{
  // calculate Lorentz-vectors of all isobars
  _decay->calcIsobarLzVec();
  // modify event for testing purposes
  if (_doSpaceInversion) {
    spaceInvertDecay();
    // recalculate Lorentz-vectors of all isobars
    _decay->calcIsobarLzVec();
  }
  if (_doReflection) {
    reflectDecay();
    // recalculate Lorentz-vectors of all isobars
    _decay->calcIsobarLzVec();
  }
  // calculate Lorentz-transformations into the correct frames for the
  // daughters in the decay vertices
  // 1) transform daughters of all decay vertices into Gottfried-Jackson frame
  //!!! this assumes that beam particle is first incoming particle in production vertex
  //    this should be solved in a more general way so that the production vertex is asked
  //    for the beam
  const TLorentzVector&  beamLv  = _decay->productionVertex  ()->inParticles()[0]->lzVec();
  const TLorentzVector&  XLv     = _decay->XIsobarDecayVertex()->mother()->lzVec();
  const TLorentzRotation gjTrans = gjTransform(beamLv, XLv);
  for (unsigned int i = 0; i < _decay->nmbInteractionVertices(); ++i) {
    const isobarDecayVertexPtr& vertex = _decay->isobarDecayVertices()[i];
    if (_debug)
      printInfo << "transforming outgoing particles of vertex " << *vertex
		<< " into " << vertex->mother()->name() << " Gottfried-Jackson RF" << endl;
    vertex->transformOutParticles(gjTrans);
  }
  // 2) transform daughters of isobar decay vertices to the respective helicity frames
  for (unsigned int i = 1; i < _decay->nmbInteractionVertices(); ++i) {  // exclude X-decay vertex
    const isobarDecayVertexPtr& vertex  = _decay->isobarDecayVertices()[i];
    if (_debug)
      printInfo << "transforming all child particles of vertex " << *vertex
  		<< " into " << vertex->mother()->name() << " helicity RF" << endl;
    const TLorentzRotation hfTrans = hfTransform(vertex->mother()->lzVec());
    // get all particles downstream of this vertex
    decayTopologyGraphType subGraph = _decay->dfsSubGraph(vertex);
    decayTopologyGraphType::edgeIterator iEd, iEdEnd;
    for (tie(iEd, iEdEnd) = subGraph.edges(); iEd != iEdEnd; ++iEd) {
      const particlePtr& part = subGraph.particle(*iEd);
      if (_debug)
    	cout << "    transforming " << part->name() << " into "
    	     << vertex->mother()->name() << " helicity RF" << endl;
      part->transform(hfTrans);
    }
  }
}


// assumes that daughters were transformed into mother RF
complex<double>
isobarHelicityAmplitude::twoBodyDecayAmplitude(const isobarDecayVertexPtr& vertex,
					       const bool                  topVertex) const
{
  const particlePtr& mother    = vertex->mother();
  const particlePtr& daughter1 = vertex->daughter1();
  const particlePtr& daughter2 = vertex->daughter2();

  // calculate Clebsch-Gordan coefficient for L-S coupling
  const int    L         = vertex->L();
  const int    S         = vertex->S();
  const int    J         = mother->J();
  const int    lambda1   = daughter1->spinProj();
  const int    lambda2   = daughter2->spinProj();
  const int    lambda    = lambda1 - lambda2;
  const double lsClebsch = cgCoeff(L, 0, S, lambda, J, lambda, _debug);
  if (lsClebsch == 0)
    return 0;

  // calculate Clebsch-Gordan coefficient for S-S coupling
  const int    s1        = daughter1->J();
  const int    s2        = daughter2->J();
  const double ssClebsch = cgCoeff(s1, lambda1, s2, -lambda2, S, lambda, _debug);
  if (ssClebsch == 0)
    return 0;

  // calculate D-function
  const int       Lambda = mother->spinProj();
  const int       P      = mother->P();
  const int       refl   = mother->reflectivity();
  const double    phi    = daughter1->lzVec().Phi();  // use daughter1 as analyzer
  const double    theta  = daughter1->lzVec().Theta();
  complex<double> DFunc;
  if (topVertex and _useReflectivityBasis)
    DFunc = DFuncConjRefl(J, Lambda, lambda, P, refl, phi, theta, _debug);
  else
    DFunc = DFuncConj(J, Lambda, lambda, phi, theta, _debug);

  // calulate barrier factor
  const double q  = daughter1->lzVec().Vect().Mag();
  const double bf = barrierFactor(L, q, _debug);

  // calculate Breit-Wigner
  const complex<double> bw = vertex->massDepAmplitude();

  // calculate normalization factor
  const double norm = normFactor(L, _debug);

  // calculate decay amplitude
  complex<double> amp = norm * DFunc * lsClebsch * ssClebsch * bf * bw;
  
  if (_debug)
    printInfo << "two-body decay amplitude = " << amp << endl;
  return amp;
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
isobarHelicityAmplitude::twoBodyDecayAmplitudeSum(const isobarDecayVertexPtr& vertex,
						  const bool                  topVertex) const
{
  if (_debug)
    printInfo << "calculating decay amplitude for " << *vertex << endl;
  const particlePtr& mother    = vertex->mother();
  const particlePtr& daughter1 = vertex->daughter1();
  const particlePtr& daughter2 = vertex->daughter2();
  complex<double>    ampSum    = 0;
  for (int lambda1 = -daughter1->J(); lambda1 <= +daughter1->J(); lambda1 += 2) {
    // calculate decay amplitude for daughter 1
    const isobarDecayVertexPtr& daughter1Vertex =
      dynamic_pointer_cast<isobarDecayVertex>(_decay->toVertex(daughter1));
    complex<double> daughter1Amp = 0;
    if (daughter1Vertex) {
      daughter1->setSpinProj(lambda1);
      daughter1Amp = twoBodyDecayAmplitudeSum(daughter1Vertex, false);
    } else
      daughter1Amp = 1;
    for (int lambda2 = -daughter2->J(); lambda2 <= +daughter2->J(); lambda2 += 2) {
      // calculate decay amplitude for daughter 2
      const isobarDecayVertexPtr& daughter2Vertex =
	dynamic_pointer_cast<isobarDecayVertex>(_decay->toVertex(daughter2));
      complex<double> daughter2Amp = 0;
      if (daughter2Vertex) {
	daughter2->setSpinProj(lambda2);
	daughter2Amp = twoBodyDecayAmplitudeSum(daughter2Vertex, false);
      } else
	daughter2Amp = 1;
      complex<double> motherAmp = twoBodyDecayAmplitude(vertex, topVertex);
      complex<double> amp       = motherAmp * daughter1Amp * daughter2Amp;
      if (_debug)
	printInfo << "amplitude term for : "
		  << mother->name()    << " [lambda = " << 0.5 * mother->spinProj() << "] -> "
		  << daughter1->name() << " [lambda = " << 0.5 * lambda1 << "] + "
		  << daughter2->name() << " [lambda = " << 0.5 * lambda2 << "] = "
		  << "mother amp. = " << motherAmp << " * daughter_1 amp = "
		  << daughter1Amp << " * daughter_2 amp = "<< daughter2Amp << " = " << amp << endl;
      ampSum += amp;
    }
  }
  if (_debug)
    printInfo << "decay amplitude for " << *vertex << " = " << ampSum << endl;
  return ampSum;
}


complex<double>
isobarHelicityAmplitude::sumBoseSymTerms(const map<string, vector<unsigned int> >&     origFsPartIndices,
					 const map<string, vector<unsigned int> >&     newFsPartIndices,
					 map<string, vector<unsigned int> >::iterator& newFsPartIndicesEntry) const
{
  
  complex<double> amp = 0;
  do {
    map<string, vector<unsigned int> >::iterator nextFsPartIndicesEntry = newFsPartIndicesEntry;
    if (++nextFsPartIndicesEntry != newFsPartIndices.end())
      // recurse to other permutations
      amp += sumBoseSymTerms(origFsPartIndices, newFsPartIndices, nextFsPartIndicesEntry);
    else {
      // build final state index map for this permutation
      vector<unsigned int> fsPartIndexMap(_decay->nmbFsParticles(), 0);
      if (_debug)
      	printInfo << "calculating amplitude for Bose term final state permutation ";
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
	  fsPartIndexMap[origFsPartIndex] = newFsPartIndex;
	}
      }
      if (_debug)
       	cout << endl;
      // (re)set final state momenta
      if (not _decay->revertMomenta(fsPartIndexMap)) {
      	printWarn << "problems reverting momenta in decay topology. returning 0." << endl;
      	return 0;
      }
      // transform daughters into their respective RFs
      transformDaughters();
      // calculate amplitude
      amp += twoBodyDecayAmplitudeSum(_decay->XIsobarDecayVertex(), true);
    }
  } while (next_permutation(newFsPartIndicesEntry->second.begin(),
			    newFsPartIndicesEntry->second.end()));
  return amp;
}


complex<double>
isobarHelicityAmplitude::boseSymmetrizedAmp() const
{
  // get final state indistinguishable particles
  typedef map<string, unsigned int>::const_iterator indistFsPartIt;
  const map<string, unsigned int> indistFsPart = _decay->nmbIndistFsParticles();
  if (_debug) {
    printInfo << "Bose-symmetrizing amplitude: indistinguishable final state particle "
	      << "multiplicities (marked FS particles will be Bose-symmetrized): ";
    for (indistFsPartIt i = indistFsPart.begin(); i != indistFsPart.end(); ++i)
      cout << i->first << " = " << i->second << ((i->second) >= 2 ? " <<<  " : "  ");
    cout << endl;
  }

  // calculate normalization factor
  double nmbCombinations = 1;
  for (indistFsPartIt i = indistFsPart.begin(); i != indistFsPart.end(); ++i)
    nmbCombinations *= TMath::Factorial(i->second);
  const double normFactor = 1 / sqrt(nmbCombinations);

  // initialize indices used to generate final state permutations
  // in order to get all permutations with std::next_permutation
  // indices have to be sorted ascending
  map<string, vector<unsigned int> > origFsPartIndices, newFsPartIndices;
  for (unsigned int i = 0; i < _decay->nmbFsParticles(); ++i) {
    const string partName = _decay->fsParticles()[i]->name();
    origFsPartIndices[partName].push_back(i);
    newFsPartIndices [partName].push_back(i);
  }

  // Bose-symmetrize amplitudes
  if (_debug)
    printInfo << "Bose-symmetrizing amplitude using " << nmbCombinations << " terms" << endl;
  map<string, vector<unsigned int> >::iterator firstEntry = newFsPartIndices.begin();
  return normFactor * sumBoseSymTerms(origFsPartIndices, newFsPartIndices, firstEntry);
}


ostream&
isobarHelicityAmplitude::print(ostream& out) const
{
  out << "isobar helicity amplitude: "
      << *_decay
      << "reflectivity basis : "            << ((_useReflectivityBasis) ? "en" : "dis") << "abled" << endl
      << "Bose-symmetrization: "            << ((_boseSymmetrize)       ? "en" : "dis") << "abled" << endl
      << "space inversion of FS momenta: "  << ((_doSpaceInversion)     ? "en" : "dis") << "abled" << endl
      << "reflection through prod. plane: " << ((_doReflection)         ? "en" : "dis") << "abled" << endl;
  return out;
}


