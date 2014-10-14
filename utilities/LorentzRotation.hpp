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
//      An alternative vector implementation similar to TLorentzRotation of ROOT.
//      Usable with CUDA.
//
//
// Author List:
//      Armin Gensler          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef RPWA_LORENTZ_ROTATION_HPP
#define RPWA_LORENTZ_ROTATION_HPP

#include "TLorentzRotation.h"

#include "cudaUtils.hpp"
#include "Vector3.hpp"
#include "Rotation.hpp"

namespace rpwa {

	template<typename T>
	class RpwaLorentzRotation {
	
	public:

		typedef T Scalar;
		
		HOST_DEVICE void SetBoost(T bx, T by, T bz) {
			T bp2 = bx*bx + by*by + bz*bz;
			T gamma = 1 / sqrt(1 - bp2);
			T bgamma = gamma * gamma / (1 + gamma);
			xx = 1 + bgamma * bx * bx;
			yy = 1 + bgamma * by * by;
			zz = 1 + bgamma * bz * bz;
			xy = yx = bgamma * bx * by;
			xz = zx = bgamma * bx * bz;
			yz = zy = bgamma * by * bz;
			xt = tx = gamma * bx;
			yt = ty = gamma * by;
			zt = tz = gamma * bz;
			tt = gamma;
		}
		
		HOST_DEVICE RpwaLorentzRotation():
			xx(1), xy(0), xz(0), xt(0),
			yx(0), yy(1), yz(0), yt(0),
			zx(0), zy(0), zz(1), zt(0),
			tx(0), ty(0), tz(0), tt(1)
		{}
		
		HOST_DEVICE RpwaLorentzRotation(T xx, T xy, T xz, T xt,
			        T yx, T yy, T yz, T yt,
			        T zx, T zy, T zz, T zt,
			        T tx, T ty, T tz, T tt):
			xx(xx), xy(xy), xz(xz), xt(xt),
			yx(yx), yy(yy), yz(yz), yt(yt),
			zx(zx), zy(zy), zz(zz), zt(zt),
			tx(tx), ty(ty), tz(tz), tt(tt)
		{}
		
		HOST_DEVICE RpwaLorentzRotation(const RpwaRotation<T> & rot):
			xx(rot.XX()), xy(rot.XY()), xz(rot.XZ()), xt(0),
			yx(rot.YX()), yy(rot.YY()), yz(rot.YZ()), yt(0),
			zx(rot.ZX()), zy(rot.ZY()), zz(rot.ZZ()), zt(0),
			tx(0), ty(0), tz(0), tt(1)
		{}
		
		HOST_DEVICE RpwaLorentzRotation(const RpwaVector3<T> & boost) {
			SetBoost(boost.X(), boost.Y(), boost.Z());
		}
		
		
		HOST RpwaLorentzRotation(const TLorentzRotation & rot):
			xx(rot.XX()), xy(rot.XY()), xz(rot.XZ()), xt(rot.XT()),
			yx(rot.YX()), yy(rot.YY()), yz(rot.YZ()), yt(rot.YT()),
			zx(rot.ZX()), zy(rot.ZY()), zz(rot.ZZ()), zt(rot.ZT()),
			tx(rot.TX()), ty(rot.TY()), tz(rot.TZ()), tt(rot.TT())
		{}
		
		HOST RpwaLorentzRotation ToROOT() const {
			return RpwaLorentzRotation(
				xx, xy, xz, xt,
				yx, yy, yz, yt,
				zx, zy, zz, zt,
				tx, ty, tz, tt);
		}
		
		HOST_DEVICE T XX() const { return xx; }
		HOST_DEVICE T XY() const { return xy; }
		HOST_DEVICE T XZ() const { return xz; }
		HOST_DEVICE T XT() const { return xt; }
		HOST_DEVICE T YX() const { return yx; }
		HOST_DEVICE T YY() const { return yy; }
		HOST_DEVICE T YZ() const { return yz; }
		HOST_DEVICE T YT() const { return yt; }
		HOST_DEVICE T ZX() const { return zx; }
		HOST_DEVICE T ZY() const { return zy; }
		HOST_DEVICE T ZZ() const { return zz; }
		HOST_DEVICE T ZT() const { return zt; }
		HOST_DEVICE T TX() const { return tx; }
		HOST_DEVICE T TY() const { return ty; }
		HOST_DEVICE T TZ() const { return tz; }
		HOST_DEVICE T TT() const { return tt; }
		
		HOST_DEVICE RpwaLorentzRotation<T> operator * (const RpwaLorentzRotation<T> & b) const {
			return RpwaLorentzRotation<T>(
				xx*b.xx + xy*b.yx + xz*b.zx + xt*b.tx,
				xx*b.xy + xy*b.yy + xz*b.zy + xt*b.ty,
				xx*b.xz + xy*b.yz + xz*b.zz + xt*b.tz,
				xx*b.xt + xy*b.yt + xz*b.zt + xt*b.tt,
				yx*b.xx + yy*b.yx + yz*b.zx + yt*b.tx,
				yx*b.xy + yy*b.yy + yz*b.zy + yt*b.ty,
				yx*b.xz + yy*b.yz + yz*b.zz + yt*b.tz,
				yx*b.xt + yy*b.yt + yz*b.zt + yt*b.tt,
				zx*b.xx + zy*b.yx + zz*b.zx + zt*b.tx,
				zx*b.xy + zy*b.yy + zz*b.zy + zt*b.ty,
				zx*b.xz + zy*b.yz + zz*b.zz + zt*b.tz,
				zx*b.xt + zy*b.yt + zz*b.zt + zt*b.tt,
				tx*b.xx + ty*b.yx + tz*b.zx + tt*b.tx,
				tx*b.xy + ty*b.yy + tz*b.zy + tt*b.ty,
				tx*b.xz + ty*b.yz + tz*b.zz + tt*b.tz,
				tx*b.xt + ty*b.yt + tz*b.zt + tt*b.tt);
		}
		
		HOST_DEVICE void Transform(const RpwaLorentzRotation<T> & rot) {
			*this = rot * (*this);
		}
		HOST_DEVICE void operator *= (const RpwaLorentzRotation<T> & rot) {
			*this = (*this) * rot;  
		}
		
		HOST_DEVICE void Boost(const RpwaVector3<T> & b) {
			Transform(RpwaLorentzRotation<T>(b));
		}
		
		HOST_DEVICE RpwaLorentzRotation<T> Inverse() const {
			return RpwaLorentzRotation<T>(
				 xx,  yx,  zx, -tx,
				 xy,  yy,  zy, -ty,
				 xz,  yz,  zz, -tz,
				-xt, -yt, -zt,  tt);
		}
		
		HOST_DEVICE void Invert() {
			*this = Inverse();
		}
		
	private:
		T xx, xy, xz, xt,
		  yx, yy, yz, yt,
		  zx, zy, zz, zt,
		  tx, ty, tz, tt;
	};	

}  // namespace rpwa


#endif  // LORENTZ_ROTATION_HPP
