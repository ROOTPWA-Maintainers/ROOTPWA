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
//      An alternative vector implementation similar to TRotation of ROOT.
//      Usable with CUDA.
//
//
// Author List:
//      Armin Gensler          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef RPWA_ROTATION_HPP
#define RPWA_ROTATION_HPP

#include "TRotation.h"

#include "cudaUtils.hpp"

namespace rpwa {

	template<typename T>
	class RpwaRotation {
	
	public:

		typedef T Scalar;
		
		HOST_DEVICE RpwaRotation():
			xx(1), xy(0), xz(0),
			yx(0), yy(1), yz(0),
			zx(0), zy(0), zz(1)
		{}
		
		HOST_DEVICE RpwaRotation(T xx, T xy, T xz,
			 T yx, T yy, T yz,
			 T zx, T zy, T zz):
			xx(xx), xy(xy), xz(xz),
			yx(yx), yy(yy), yz(yz),
			zx(zx), zy(zy), zz(zz)
		{}
		
		HOST RpwaRotation(const TRotation & rot):
			xx(rot.XX()), xy(rot.XY()), xz(rot.XZ()),
			yx(rot.YX()), yy(rot.YY()), yz(rot.YZ()),
			zx(rot.ZX()), zy(rot.ZY()), zz(rot.ZZ())
		{}
		
		HOST TRotation ToROOT() const {
			return TRotation(
				xx, xy, xz,
				yx, yy, yz,
				zx, zy, zz);
		}
		
		HOST_DEVICE T XX() const { return xx; }
		HOST_DEVICE T XY() const { return xy; }
		HOST_DEVICE T XZ() const { return xz; }
		HOST_DEVICE T YX() const { return yx; }
		HOST_DEVICE T YY() const { return yy; }
		HOST_DEVICE T YZ() const { return yz; }
		HOST_DEVICE T ZX() const { return zx; }
		HOST_DEVICE T ZY() const { return zy; }
		HOST_DEVICE T ZZ() const { return zz; }
		
		HOST_DEVICE void RotateX(T a) {
			T c = cos(a);
			T s = sin(a);
			T x = yx, y = yy, z = yz;
			yx = c*x - s*zx;
			yy = c*y - s*zy;
			yz = c*z - s*zz;
			zx = s*x + c*zx;
			zy = s*y + c*zy;
			zz = s*z + c*zz;  
		}
		
		HOST_DEVICE void RotateY(T a) {
			T c = cos(a);
			T s = sin(a);
			T x = zx, y = zy, z = zz;
			zx = c*x - s*xx;
			zy = c*y - s*xy;
			zz = c*z - s*xz;
			xx = s*x + c*xx;
			xy = s*y + c*xy;
			xz = s*z + c*xz;
		}
		
		HOST_DEVICE void RotateZ(T a) {
			T c = cos(a);
			T s = sin(a);
			T x = xx, y = xy, z = xz;
			xx = c*x - s*yx;
			xy = c*y - s*yy;
			xz = c*z - s*yz;
			yx = s*x + c*yx;
			yy = s*y + c*yy;
			yz = s*z + c*yz;
		}
		
		HOST_DEVICE RpwaRotation<T> operator* (const RpwaRotation<T> & b) const {
			return RpwaRotation<T>(
					xx*b.xx + xy*b.yx + xz*b.zx,
					xx*b.xy + xy*b.yy + xz*b.zy,
					xx*b.xz + xy*b.yz + xz*b.zz,
					yx*b.xx + yy*b.yx + yz*b.zx,
					yx*b.xy + yy*b.yy + yz*b.zy,
					yx*b.xz + yy*b.yz + yz*b.zz,
					zx*b.xx + zy*b.yx + zz*b.zx,
					zx*b.xy + zy*b.yy + zz*b.zy,
					zx*b.xz + zy*b.yz + zz*b.zz);
		}
		
		HOST_DEVICE void operator *= (const RpwaRotation<T> & m) {
			*this = (*this) * m;
		}
		HOST_DEVICE void Transform(const RpwaRotation<T> & m) {
			*this = m * (*this);
		}
		
	private:
		T xx, xy, xz,
		  yx, yy, yz,
		  zx, zy, zz;
	};

}  // namespace rpwa


#endif  // ROTATION_HPP
