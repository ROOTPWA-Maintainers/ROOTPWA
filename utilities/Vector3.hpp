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
//      An alternative vector implementation similar to TVector3 of ROOT.
//      Usable with CUDA.
//
//
// Author List:
//      Armin Gensler          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef RPWA_VECTOR3_HPP
#define RPWA_VECTOR3_HPP

#include <cmath>
#include <iostream>

#include "TVector3.h"

#include "Rotation.hpp"
#include "cudaUtils.hpp"

namespace rpwa {

	template<typename T>
	class RpwaVector3 {
	
	public:

		typedef T Scalar;
	  
		HOST_DEVICE RpwaVector3(): _x(0), _y(0), _z(0) {}
		HOST_DEVICE RpwaVector3(T x, T y, T z): _x(x), _y(y), _z(z) {}
		HOST        RpwaVector3(const TVector3 & v): _x(v.X()), _y(v.Y()), _z(v.Z()) {}
		
		HOST TVector3 ToROOT() const { return TVector3(_x, _y, _z); }
		
		HOST_DEVICE T X() const { return _x; }
		HOST_DEVICE T Y() const { return _y; }
		HOST_DEVICE T Z() const { return _z; }
		HOST_DEVICE T x() const { return X(); }
		HOST_DEVICE T y() const { return Y(); }
		HOST_DEVICE T z() const { return Z(); }
		HOST_DEVICE T Px() const { return X(); }
		HOST_DEVICE T Py() const { return Y(); }
		HOST_DEVICE T Pz() const { return Z(); }
		
		HOST_DEVICE void SetX(T x) { _x = x; }
		HOST_DEVICE void SetY(T y) { _y = y; }
		HOST_DEVICE void SetZ(T z) { _z = z; }
		HOST_DEVICE void SetXYZ(T x, T y, T z) { _x = x; _y = y; _z = z; }
		
		HOST_DEVICE T Mag2() const { return _x * _x + _y * _y + _z * _z; }
		HOST_DEVICE T Mag() const { return sqrt(Mag2()); }
		
		HOST_DEVICE T Phi() const {
			//return the  azimuth angle. returns phi from -pi to pi
			return (_x == 0 && _y == 0) ? 0 : atan2(_y, _x);
		}
		
		HOST_DEVICE T Theta() const {
			T perp = sqrt(_x * _x + _y * _y);
			return (_x == 0 && _y == 0 && _z == 0) ? 0 : atan2(perp, _z);
		}
		
		HOST_DEVICE bool operator == (const RpwaVector3<T> & other) const {
			return _x == other._x && _y == other._y  && _z == other._z; 
		}
		HOST_DEVICE bool operator != (const RpwaVector3<T> & other) const {
			return ! ((*this) == other); 
		}
		
		HOST_DEVICE RpwaVector3<T>& operator += (const RpwaVector3<T> & other) {
			_x += other._x;
			_y += other._y;
			_z += other._z;
			return *this;
		}
		
		HOST_DEVICE RpwaVector3<T>& operator -= (const RpwaVector3<T> & other) {
			_x -= other._x; 
			_y -= other._y; 
			_z -= other._z; 
			return *this; 
		}
		
		HOST_DEVICE RpwaVector3<T>& operator *= (T f) {
			_x *= f;
			_y *= f; 
			_z *= f; 
			return *this;
		}
		
		HOST_DEVICE RpwaVector3<T>& operator *= (const RpwaRotation<T> & rot) {
			*this = RpwaVector3<T>(
				rot.XX() * _x + rot.XY() * _y + rot.XZ() * _z,
				rot.YX() * _x + rot.YY() * _y + rot.YZ() * _z,
				rot.ZX() * _x + rot.ZY() * _y + rot.ZZ() * _z);
			return *this;
		}

		HOST_DEVICE RpwaVector3<T> operator - () const { return RpwaVector3<T>(- _x, - _y, - _z); }
		
		HOST_DEVICE T Dot(const RpwaVector3<T> & other) const {
			return _x * other._x + _y * other._y + _z * other._z;
		}

		HOST_DEVICE RpwaVector3<T> Cross(const RpwaVector3<T> & other) const {
			return RpwaVector3(
					_y * other._z - _z * other._y,
					_z * other._x - _x * other._z,
					_x * other._y - _y * other._x);
		}

	private:
		T _x, _y, _z;

	};
	
	template <typename T>
	HOST_DEVICE inline RpwaVector3<T> operator + (const RpwaVector3<T>& a, const RpwaVector3<T>& b) {
		return RpwaVector3<T>(a.x() + b.x(), a.y() + b.y(), a.z() + b.z());
	}
	
	template <typename T>
	HOST_DEVICE inline RpwaVector3<T> operator - (const RpwaVector3<T>& a, const RpwaVector3<T>& b) {
		return RpwaVector3<T>(a.x() - b.x(), a.y() - b.y(), a.z() - b.z());
	}
	template <typename T>
	HOST_DEVICE inline RpwaVector3<T> operator * (T f, const RpwaVector3<T>& a) {
		return RpwaVector3<T>(f * a.x(), f * a.y(), f * a.z());
	}
	
	template <typename T>
	HOST_DEVICE inline RpwaVector3<T> operator * (const RpwaVector3<T>& a, T f) {
		return RpwaVector3<T>(f * a.x(), f * a.y(), f * a.z());
	}
	
	template <typename T>
	HOST_DEVICE inline RpwaVector3<T> operator * (const RpwaVector3<T>& a, const RpwaVector3<T>& b) {
		return a.Dot(b);
	}	

	template <typename T>
	HOST inline std::ostream& operator << (std::ostream& stream, const RpwaVector3<T>& vec) {
		return stream << vec.ToROOT();
	}

}  // namespace rpwa


#endif  // VECTOR3_HPP
