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
//    but WITHOUTy ANY WARRANTY; without even the implied warranty of
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
//      An alternative vector implementation similar to TLorentzVector of ROOT.
//      Usable with CUDA.
//
//
// Author List:
//      Armin Gensler          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef RPWA_LORENTZ_VECTOR_HPP
#define RPWA_LORENTZ_VECTOR_HPP

#include <cmath>

#include "TLorentzVector.h"

#include "cudaUtils.hpp"
#include "LorentzRotation.hpp"
#include "Vector3.hpp"

namespace rpwa {

	template<typename Ty> // Ty would shadow method named T
	class RpwaLorentzVector {
	
	public:

		typedef Ty Scalar;
	  
		HOST_DEVICE RpwaLorentzVector(): _xyz(0, 0, 0), _t(0) {}
		HOST_DEVICE RpwaLorentzVector(Ty x, Ty y, Ty z, Ty t): _xyz(x, y, z), _t(t) {}
		HOST_DEVICE RpwaLorentzVector(const RpwaVector3<Ty>& xyz, Ty t): _xyz(xyz), _t(t) {}
		HOST        RpwaLorentzVector(const TLorentzVector & v): _xyz(v.Vect()), _t(v.T()) {}
		
		HOST TLorentzVector ToROOT() const { return TLorentzVector(_xyz.ToROOT(), _t); }
		HOST TLorentzVector operator()() const { return TLorentzVector(_xyz.ToROOT(), _t); }
		
		HOST_DEVICE Ty X() const { return _xyz.X(); }
		HOST_DEVICE Ty Y() const { return _xyz.Y(); }
		HOST_DEVICE Ty Z() const { return _xyz.Z(); }
		HOST_DEVICE Ty T() const { return _t; }
		HOST_DEVICE Ty Px() const { return X(); }
		HOST_DEVICE Ty Py() const { return Y(); }
		HOST_DEVICE Ty Pz() const { return Z(); }
		HOST_DEVICE Ty E() const { return T(); }
		
		HOST_DEVICE void SetX(Ty x) { _xyz.SetX(x); }
		HOST_DEVICE void SetY(Ty y) { _xyz.SetY(y); }
		HOST_DEVICE void SetZ(Ty z) { _xyz.SetZ(z); }
		HOST_DEVICE void SetT(Ty t) { _t = t; }
		HOST_DEVICE void SetPx(Ty x) { SetX(x); }
		HOST_DEVICE void SetPy(Ty y) { SetY(y); }
		HOST_DEVICE void SetPz(Ty z) { SetZ(z); }
		HOST_DEVICE void SetE(Ty t) { SetT(t); }
		
		HOST_DEVICE void SetXYZT(Ty x, Ty y, Ty z, Ty t) { _xyz.SetXYZ(x, y, z); _t = t; }
		HOST_DEVICE void SetPxPyPyE(Ty x, Ty y, Ty z, Ty t) { SetXYZT(x, y, z, t); }
		HOST_DEVICE void SetXYZM(Ty x, Ty y, Ty z, Ty m) {
			if(m >= 0) 
				SetXYZT(x, y, z, sqrt(x*x+y*y+z*z+m*m));
			else 
				SetXYZT(x, y, z, sqrt(max(x*x+y*y+z*z-m*m,0)) );
		}
		
		//HOST_DEVICE void SetVectMag(const Vector3<Ty>& vect, Ty mag) { SetXYZT(vect.X(), vect.Y(), vect.Z(), mag); }
		//HOST_DEVICE void SetVectM(const Vector3<Ty>& vect, Ty m) { SetVectMag(vect, m); }
		
		HOST_DEVICE RpwaVector3<Ty> Vect() const { return _xyz; }
		HOST_DEVICE void SetVect(const RpwaVector3<Ty>& p) { _xyz = p; }
		HOST_DEVICE Ty P() const { return _xyz.Mag(); }	
		
		HOST_DEVICE bool operator == (const RpwaLorentzVector<Ty> & that) const {
			return _xyz == that._xyz && _t == that._t; 
		}
		HOST_DEVICE bool operator != (const RpwaLorentzVector<Ty> & that) const {
			return ! ((*this) == that); 
		}
		
		HOST_DEVICE RpwaLorentzVector<Ty>& operator += (const RpwaLorentzVector<Ty> & that) {
			_xyz += that._xyz; 
			_t += that._t; 
			return *this; 
		}
		
		HOST_DEVICE RpwaLorentzVector<Ty>& operator -= (const RpwaLorentzVector<Ty> & that) {
			_xyz -= that._xyz; 
			_t -= that._t; 
			return *this; 
		}
		
		HOST_DEVICE RpwaLorentzVector<Ty>& operator *= (Ty f) {
			_xyz *= f;
			_t *= f; 
			return *this;
		}
		
		HOST_DEVICE RpwaLorentzVector<Ty> operator + (const RpwaLorentzVector<Ty> & that) const {
			return RpwaLorentzVector<Ty>(_xyz + that._xyz, _t + that._t);
		}
		
		HOST_DEVICE RpwaLorentzVector<Ty> operator - (const RpwaLorentzVector<Ty> & that) const {
			return RpwaLorentzVector<Ty>(_xyz - that._xyz, _t - that._t);
		}
		
		HOST_DEVICE RpwaLorentzVector<Ty> operator * (Ty f) const {
			return RpwaLorentzVector<Ty>(_xyz * f, _t * f);
		}
		
		HOST_DEVICE Ty Dot(const RpwaLorentzVector<Ty> & that) const {
			return T() * that.T() - X() * that.X() - Y() * that.Y() - Y() * that.Y();
		}
		
		HOST_DEVICE Ty operator * (const RpwaLorentzVector<Ty> & that) const {
			return Dot(that);
		}
		
		HOST_DEVICE Ty Mag2() const { 
			return _t * _t - _xyz.Mag2(); 
		}
		
		HOST_DEVICE Ty Mag() const {
			Ty mm = Mag2();
			return mm < 0.0 ? -sqrt(-mm) : sqrt(mm);
		}
		
		HOST_DEVICE Ty Phi() const { return _xyz.Phi(); }
		HOST_DEVICE Ty Theta() const { return _xyz.Theta(); }
		
		HOST_DEVICE Ty M2() const { return Mag2(); }
		HOST_DEVICE Ty M() const { return Mag(); }
		
		HOST_DEVICE void Transform(const RpwaLorentzRotation<Ty> & rot) {
			*this = RpwaLorentzVector(
				rot.XX()*X() + rot.XY()*Y() + rot.XZ()*Z() + rot.XT()*T(),
				rot.YX()*X() + rot.YY()*Y() + rot.YZ()*Z() + rot.YT()*T(),
				rot.ZX()*X() + rot.ZY()*Y() + rot.ZZ()*Z() + rot.ZT()*T(),
				rot.TX()*X() + rot.TY()*Y() + rot.TZ()*Z() + rot.TT()*T());
		}
		
		HOST_DEVICE void operator *= (const RpwaLorentzRotation<Ty> & rot) { Transform(rot); }
		
		//HOST_DEVICE void Transform(const Rotation<Ty> & rot);
		//HOST_DEVICE void operator *= (const Rotation<Ty> & rot) { Transform(rot); }
		
		HOST_DEVICE RpwaVector3<Ty> BoostVector() const {
			return RpwaVector3<Ty>(X()/T(), Y()/T(), Z()/T());
		}
		
		HOST_DEVICE void Boost(const RpwaVector3<Ty> & b) {
			Ty bx = b.X();
			Ty by = b.Y();
			Ty bz = b.Z();
			Ty b2 = bx*bx + by*by + bz*bz;
			Ty gamma = 1 / sqrt(1 - b2);
			Ty bp = bx*X() + by*Y() + bz*Z();
			Ty gamma2 = b2 > 0 ? (gamma - 1)/b2 : 0;
			SetX(X() + gamma2*bp*bx + gamma*bx*T());
			SetY(Y() + gamma2*bp*by + gamma*by*T());
			SetZ(Z() + gamma2*bp*bz + gamma*bz*T());
			SetT(gamma*(T() + bp));
		}
		
	private:
		RpwaVector3<Ty> _xyz;
		Ty _t;
	};
	
	template <typename Ty>
	HOST_DEVICE inline RpwaLorentzVector<Ty> operator * (Ty f, const RpwaLorentzVector<Ty>& b) {
		return b * f;
	}	
	
	template <typename T>
	HOST_DEVICE inline RpwaLorentzVector<T> operator * (const RpwaLorentzRotation<T> rot, const RpwaLorentzVector<T>& v) {
		RpwaLorentzVector<T> result = v;
		result.Transform(rot);
		return result;
	}

	template <typename T>
	HOST inline std::ostream& operator << (std::ostream& stream, const RpwaLorentzVector<T>& vec) {
		return stream << vec.ToROOT();
	}

}  // namespace rpwa


#endif  // LORENTZ_VECTOR_HPP
