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
//      container class for all external particle related information
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "reportingUtils.hpp"
#include "reportingUtilsRoot.hpp"
#include "conversionUtils.hpp"
#include "particle.h"


using namespace std;
using namespace rpwa;


bool particle::_debug = false;


particle::particle()
	: particleProperties(),
	  _spinProj         (0),
	  _lzVecs           (),
	  _index            (-1),
	  _refl             (0)
{ }


particle::particle(const particle& part)
{
	*this = part;
}


particle::particle(const particleProperties&    partProp,
                   const int                    index,
                   const int                    spinProj,
                   const int                    refl,
                   const std::vector<TVector3>& momenta)
	: particleProperties(partProp),
	  _spinProj         (spinProj),
	  _lzVecs           (),
	  _index            (index),
	  _refl             (refl)
{
	setMomenta(momenta);
}


particle::particle(const string&                partName,
                   const bool                   requirePartInTable,
                   const int                    index,
                   const int                    spinProj,
                   const int                    refl,
                   const std::vector<TVector3>& momenta)
	: particleProperties(),
	  _spinProj(spinProj),
	  _lzVecs  (),
	  _index   (index),
	  _refl    (refl)
{
	if (not fillFromDataTable(partName, requirePartInTable))
		// set at least name
		setName(partName);
	setMomenta(momenta);
}


particle::particle(const string& partName,
                   const int     isospin,
                   const int     G,
                   const int     J,
                   const int     P,
                   const int     C,
                   const int     spinProj,
                   const int     refl,
                   const int     index)
	: particleProperties(partName, isospin, G, J, P, C),
	  _spinProj(spinProj),
	  _index   (index),
	  _refl    (refl)
{ }


particle::~particle()
{ }


particle&
particle::operator =(const particle& part)
{
	if (this != &part) {
		particleProperties::operator =(part);
		_spinProj = part._spinProj;
		_lzVecs   = part._lzVecs;
		_index    = part._index;
		_refl     = part._refl;
	}
	return *this;
}


/*void
particle::momenta (std::vector<TVector3>& result) const
{
	parallelLorentzVectorToVector3(_lzVecs, result);
}*/


void
particle::setMomenta(const std::vector<TVector3>& momenta)
{
	const size_t nmbMom = momenta.size();
	_lzVecs.resize(nmbMom);
	// !! EVENT PARALLEL LOOP
	boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
	#pragma omp parallel for
	for(size_t i = 0; i < nmbMom; ++i) {
		const TVector3& mom = momenta[i];
		_lzVecs[i] = TLorentzVector(mom, sqrt(mom.Mag2() + mass2()));
	}
	boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
	uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
	cout << "EPL: particle::setMomenta timediff = " << timeDiff << endl;
}


const std::vector<TLorentzVector>&
particle::transform(const std::vector<TLorentzRotation>& lorentzTransforms)
{
	const size_t nmbLzVec = numParallelEvents();
	if(lorentzTransforms.size() != nmbLzVec) {
		printErr << "size of per-event-data vectors does not match. aborting." << std::endl;
		throw;
	}
	// !! EVENT PARALLEL LOOP
	boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
	#pragma omp parallel for
	for(size_t i = 0; i < nmbLzVec; ++i)
		_lzVecs[i].Transform(lorentzTransforms[i]);
	boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
	uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
	cout << "EPL: particle::transform (rotation) timediff = " << timeDiff << endl;

	return _lzVecs;
}


const std::vector<TLorentzVector>&
particle::transform(const std::vector<TVector3>& boosts)
{
	const size_t nmbLzVec = numParallelEvents();
	if(boosts.size() != nmbLzVec) {
		printErr << "size of per-event-data vectors does not match. aborting." << std::endl;
		throw;
	}
	// !! EVENT PARALLEL LOOP
	boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
	#pragma omp parallel for
	for(size_t i = 0; i < nmbLzVec; ++i)
		_lzVecs[i].Boost(boosts[i]);
	boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
	uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
	std::cout << "EPL: particle::transform (boost) timediff = " << timeDiff << std::endl;

	return _lzVecs;
}


void
particle::scaleLzVecs(double scaleX,
                      double scaleY,
                      double scaleZ,
                      double scaleE)
{
	// !! EVENT PARALLEL LOOP
	boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
	#pragma omp parallel for
	for(unsigned int i = 0; i < numParallelEvents(); ++i) {
		TLorentzVector& v = _lzVecs[i];
		_lzVecs[i].SetXYZT(v.X() * scaleX, v.Y() * scaleY, v.Z() * scaleZ, v.E() * scaleE);
	}
	boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
	uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
	cout << "EPL: particle::scaleLzVecs timediff = " << timeDiff << endl;
}


particle*
particle::doClone() const
{
	particle* particleClone(new particle(*this));
	if (_debug)
		printDebug << "cloned " << *this << "; " << this << " -> " << particleClone << endl;
	return particleClone;
}


bool
particle::isEqualTo(const particleProperties& partProp) const
{
	const particle* part = dynamic_cast<const particle*>(&partProp);
	if (not part)
		return false;
	if (not particleProperties::isEqualTo(partProp))
		return false;
	return (    (spinProj    () == part->spinProj    ())
	        and (index       () == part->index       ())
	        and (reflectivity() == part->reflectivity()));
}


void
particle::setProperties(const particleProperties& prop)
{
	if (this != &prop)
		particleProperties::operator =(prop);
}


string
particle::qnSummary() const
{
	ostringstream out;
	out << name() << "[" << spinQn(isospin()) << parityQn(G())
	    << "(" << spinQn(J()) << parityQn(P()) << parityQn(C()) << ")"
	    << spinQn(spinProj()) << parityQn(reflectivity()) << "]";
	return out.str();
}


ostream&
particle::print(ostream& out) const
{
	particleProperties::print(out);
	out << ", "
	    << "spin proj. = "     << spinQn(_spinProj) << ", "
	    << "reflectivity = "   << _refl             << ", "
	//    << "Lorentz-vectors = [" << _lzVecs            << "] GeV, "
	    << "index = "          << _index;
	return out;
}


string
particle::label() const
{
	ostringstream out;
	out << name() << "[" << spinQn(isospin()) << parityQn(G())
	    << "(" << spinQn(J()) << parityQn(P()) << parityQn(C()) << ")]";
	return out.str();
}
