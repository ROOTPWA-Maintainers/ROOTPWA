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


#include "TLorentzRotation.h"

#include "utilities.h"
#include "isobarHelicityAmplitude.h"
#include "isobarHelicityAmplitude2.h"

	
using namespace std;
using namespace boost;
using namespace rpwa;


bool isobarHelicityAmplitude2::_debug = false;


isobarHelicityAmplitude2::isobarHelicityAmplitude2()
  : _decay(0)
{ }


isobarHelicityAmplitude2::isobarHelicityAmplitude2(isobarDecayTopology2& decay)
{
  setDecayTopology(decay);
}


isobarHelicityAmplitude2::~isobarHelicityAmplitude2()
{ }


void
isobarHelicityAmplitude2::setDecayTopology(isobarDecayTopology2& decay)
{
  if (!decay.checkTopology()) {
    printErr << "decay does not have the correct topology. aborting." << endl;
    throw;
  }
  if (!decay.checkConsistency()) {
    printErr << "decay is not consistent. aborting." << endl;
    throw;
  }
  _decay = &decay;
}


// at the moment it is not possible to fix the helicity state of
// final state particles for cases where the amplitudes for the
// different helicity states have to be added incoherently
complex<double>
isobarHelicityAmplitude2::amplitude()
{
  // transform daughters into their respective RFs
  transformDaughters();
  // recursively sum over all possible helicities of the decay particles
  return twoBodyDecayAmplitudeSum(_decay->xIsobarDecayVertex(), true);
}


TLorentzRotation
isobarHelicityAmplitude2::hfTransform(const TLorentzVector& daughterLv)
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
isobarHelicityAmplitude2::gjTransform(const TLorentzVector& beamLv,
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
isobarHelicityAmplitude2::transformDaughters()
{
  // calculate Lorentz-vectors of all isobars
  _decay->calcIsobarLzVec();
  // calculate Lorentz-transformations into the correct frames for the
  // daughters in the decay vertices
  // 1) transform daughters of all decay vertices into Gottfried-Jackson frame
  //!!! this assumes that beam particle is first incoming particle in production vertex
  //    this should be solved in a more general way so that the production vertex is asked
  //    for the beam
  const TLorentzVector&  beamLv  = _decay->productionVertex  ()->inParticles()[0]->lzVec();
  const TLorentzVector&  XLv     = _decay->xIsobarDecayVertex()->mother()->lzVec();
  const TLorentzRotation gjTrans = gjTransform(beamLv, XLv);
  for (unsigned int i = 0; i < _decay->nmbInteractionVertices(); ++i) {
    const isobarDecayVertexPtr& vertex = _decay->isobarDecayVertices()[i];
    if (_debug)
      printInfo << "transforming outgoing particles of vertex " << *vertex
		<< "into " << vertex->mother()->name() << " Gottfried-Jackson RF" << endl;
    vertex->transformOutParticles(gjTrans);
  }
  // 2) transform daughters of isobar decay vertices to the respective helicity frames
  for (unsigned int i = 1; i < _decay->nmbInteractionVertices(); ++i) {  // exclude X-decay vertex
    const isobarDecayVertexPtr& vertex  = _decay->isobarDecayVertices()[i];
    if (_debug)
      printInfo << "transforming all child particles of vertex " << *vertex
		<< "into " << vertex->mother()->name() << " helicity RF" << endl;
    const TLorentzRotation hfTrans = hfTransform(vertex->mother()->lzVec());
    // get all particles downstream of this vertex
    decayTopologyGraphType subGraph = _decay->dfsSubGraph(vertex);
    decayTopologyGraphType::edgeIterator iEd, iEdEnd;
    for (tie(iEd, iEdEnd) = subGraph.edges(); iEd != iEdEnd; ++iEd) {
      const particlePtr& p = subGraph.particle(*iEd);
      if (_debug)
	cout << "    transforming " << p->name() << " into "
	     << vertex->mother()->name() << " helicity RF" << endl;
      p->transform(hfTrans);
    }
  }
}


// assumes that daughters were transformed into mother RF
complex<double>
isobarHelicityAmplitude2::twoBodyDecayAmplitude(const isobarDecayVertexPtr& vertex,
						const bool                  topVertex) const
{
  const particlePtr& mother    = vertex->mother();
  const particlePtr& daughter1 = vertex->daughter1();
  const particlePtr& daughter2 = vertex->daughter2();
  ostringstream      debug;

  // calculate normalization factor
  const int    L    = vertex->L();
  const double norm = normFactor(L);
  if (_debug)
    debug << "< sqrt(" << L << " + 1) = " << norm << " >  ";

  // calculate D-function
  const int             J       = mother->J();
  const int             Lambda  = mother->spinProj();
  const int             lambda1 = daughter1->spinProj();
  const int             lambda2 = daughter2->spinProj();
  const int             lambda  = lambda1 - lambda2;
  const double          phi     = daughter1->lzVec().Phi();  // use daughter1 as analyzer
  const double          theta   = daughter1->lzVec().Theta();
  const complex<double> DFunc   = DFuncConj(J, Lambda, lambda, phi, theta);
  if (_debug)
    debug << "< D^{" << J << " *}" << "_{" << Lambda << ", " << lambda << "}"
	  << "(" << phi << ", " << theta << ", 0) = " << DFunc << " >  ";

  // calculate Clebsch-Gordan coefficient for L-S coupling
  const int    S         = vertex->S();
  const double lsClebsch = cgCoeff(L, 0, S, lambda, J, lambda);
  if (_debug)
    debug << "< (" << L << ", 0; " << S << ", " << lambda
	  << " | " << J << ", " << lambda << ") = " << lsClebsch << " >  ";
  // if (lsClebsch == 0)
  //   return 0;

  // calculate Clebsch-Gordan coefficient for S-S coupling
  const int    s1        = daughter1->J();
  const int    s2        = daughter2->J();
  const double ssClebsch = cgCoeff(s1, lambda1, s2, -lambda2, S, lambda);
  if (_debug)
    debug << "< (" << s1 << ", " << lambda1 << "; " << s2 << ", " << -lambda2
	  << " | " << S << ", " << lambda << ") = " << ssClebsch << " >  ";
  // if (ssClebsch == 0)
  //   return 0;

  // calulate barrier factor
  const double q  = daughter1->lzVec().Vect().Mag();
  const double bf = barrierFactor(L, q);
  if (_debug)
    debug << "< Blatt-Weisskopf(" << L << ", " << q << ") = " << bf << " >  ";

  // calculate Breit-Wigner
  const double    M      = mother->lzVec().M();
  const double    M0     = mother->mass();
  const double    Gamma0 = mother->width();
  const double    m1     = daughter1->lzVec().M();
  const double    m2     = daughter2->lzVec().M();
  const double    q0     = breakupMomentum(M0, m1, m2);
  complex<double> bw     = 1;
  if (!topVertex)
    bw = breitWigner(M, M0, Gamma0, L, q, q0);
  if (_debug) {
    if (topVertex)
      debug << "< no mass dep. = " << bw << " >";
    else
      debug << "< Breit-Wigner(" << M << ", " << M0 << ", " << Gamma0
	    << ", " << L << ", " << q << ", " << q0 << ") = " << bw << " >";
  }

  // calculate decay amplitude
  complex<double> amp = norm * DFunc * lsClebsch * ssClebsch * bf * bw;
  
  if (_debug)
    printInfo << "calculating two-body decay amplitude for " << *vertex << ": "
	      << debug.str() << " = " << amp << endl;
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
isobarHelicityAmplitude2::twoBodyDecayAmplitudeSum(const isobarDecayVertexPtr& vertex,
						   const bool                  topVertex)
{
  const particlePtr& daughter1 = vertex->daughter1();
  const particlePtr& daughter2 = vertex->daughter2();
  complex<double>    ampSum    = 0;
  for (int lambda1 = -daughter1->J(); lambda1 <= +daughter1->J(); lambda1 += 2) {
    // calculate decay amplitude for daughter 1
    const isobarDecayVertexPtr& daughter1Vertex =
      dynamic_pointer_cast<isobarDecayVertex2>(_decay->toVertex(daughter1));
    complex<double> amp1 = 1;
    if (daughter1Vertex) {
      daughter1->setSpinProj(lambda1);
      amp1 *= twoBodyDecayAmplitudeSum(daughter1Vertex, false);
    }
    for (int lambda2 = -daughter2->J(); lambda2 <= +daughter2->J(); lambda2 += 2) {
      // calculate decay amplitude for daughter 2
      const isobarDecayVertexPtr& daughter2Vertex =
	dynamic_pointer_cast<isobarDecayVertex2>(_decay->toVertex(daughter2));
      complex<double> amp2 = 1;
      if (daughter2Vertex) {
	daughter2->setSpinProj(lambda2);
	amp2 *= twoBodyDecayAmplitudeSum(daughter2Vertex, false);
      }
      complex<double> amp = twoBodyDecayAmplitude(vertex, topVertex) * amp1 * amp2;
      if (_debug)
	printInfo << "amplitude for " << *vertex << ": "
		  << "[lambda = " << vertex->mother()->spinProj() << "]  --->  "
		  << daughter1->name() << " [lambda = " << lambda1 << "]  +  "
		  << daughter2->name() << " [lambda = " << lambda2 << "]  =  " << amp << endl;
      ampSum += amp;
    }
  }
  return ampSum;
}
