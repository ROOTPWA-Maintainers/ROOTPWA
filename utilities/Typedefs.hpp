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
//      Defines the types used for vectors, rotations etc.
//
//
// Author List:
//      Armin Gensler          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef RPWA_TYPEDEFS_HPP
#define RPWA_TYPEDEFS_HPP

#include <complex>

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRotation.h"
#include "TLorentzRotation.h"

#include "Complex.hpp"
#include "Vector3.hpp"
#include "LorentzVector.hpp"
#include "Rotation.hpp"
#include "LorentzRotation.hpp"

namespace rpwa {

	typedef double Scalar;

	typedef RpwaComplex<Scalar> Complex;

	typedef RpwaVector3<Scalar> Vector3;

	typedef RpwaLorentzVector<Scalar> LorentzVector;

	typedef RpwaRotation<Scalar> Rotation;

	typedef RpwaLorentzRotation<Scalar> LorentzRotation;

	/*
	typedef double Scalar;

	typedef std::complex<Scalar> Complex;
  
	typedef TVector3 Vector3;
	
	typedef TLorentzVector LorentzVector;
	
	typedef TRotation Rotation;
	
	typedef TLorentzRotation LorentzRotation;
	*/

}  // namespace rpwa


#endif  // RPWA_TYPEDEFS_HPP
