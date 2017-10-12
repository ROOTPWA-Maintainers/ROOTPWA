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
//      class hierarchy for mass-dependent part of the amplitude
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef MASSDEPENDENCE_H
#define MASSDEPENDENCE_H


#include <iostream>
#include <complex>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/matrix.hpp>


namespace libconfig {
	class Setting;
}
namespace ublas = boost::numeric::ublas;


namespace rpwa {

	class isobarDecayVertex;
	typedef boost::shared_ptr<isobarDecayVertex> isobarDecayVertexPtr;


	//////////////////////////////////////////////////////////////////////////////
	/// Brief base class for mass dependences
	class massDependence {

	public:

		massDependence()          { }
		virtual ~massDependence() { }

		virtual std::complex<double> amp(const isobarDecayVertex& v) = 0;

		virtual std::complex<double> operator ()(const isobarDecayVertex& v) { return amp(v); }

		virtual std::string name() const = 0;  ///< returns label used in graph visualization, reporting, and key file

		virtual std::string parentLabelForWaveName(const isobarDecayVertex& v) const;  ///< returns label for parent of decay used in wave name

		bool operator ==(const massDependence& rhsMassDep) const { return this->isEqualTo(rhsMassDep); }
		bool operator !=(const massDependence& rhsMassDep) const { return not (*this == rhsMassDep);   }

		virtual std::ostream& print(std::ostream& out) const;

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	protected:

		virtual bool isEqualTo(const massDependence& massDep) const
		{ return typeid(massDep) == typeid(*this); }  ///< returns whether massDep is of same type as this

		static bool _debug;  ///< if set to true, debug messages are printed

	};


	typedef boost::shared_ptr<massDependence> massDependencePtr;


	/// create a mass dependence object as specified by 'massDepType'
	// only marks for mass dependencies not requiring additional arguments,
	// if the mass dependence cannot be created return a NULL pointer
	massDependencePtr createMassDependence(const std::string& massDepType, const libconfig::Setting* setting = nullptr);


	inline
	std::ostream&
	operator <<(std::ostream&         out,
	            const massDependence& massDep)
	{
		return massDep.print(out);
	}


	//////////////////////////////////////////////////////////////////////////////
	/// Brief intermediate class for mass dependencies
	// The idea is to have an intermediate class that provides a couple of
	// static functions that otherwise would have to be added to each class
	// individually, adding a lot of duplicate code. Mass dependencies
	// should inherit from this class and provide a
	// 'static constexpr const char*' member 'cName' indicating the name of
	// the mass dependence used in keyfiles and so on.
	template<class T>
	class massDependenceImpl : public massDependence {

	public:

		massDependenceImpl() : massDependence() { }
		virtual ~massDependenceImpl()           { }

		virtual std::string name() const { return Name(); }  ///< returns label used in graph visualization, reporting, and key file

		static std::string Name() { return T::cName; } ///< returns the name used to trigger the creation of a mass dependence

		template<typename... Args>
		static boost::shared_ptr<T> Create(Args... args) {
			return boost::make_shared<T>(args...);
		}
	};


	//////////////////////////////////////////////////////////////////////////////
	/// Brief trivial flat mass dependence
	class flatMassDependence : public massDependenceImpl<flatMassDependence> {

	public:

		flatMassDependence() : massDependenceImpl<flatMassDependence>() { }
		virtual ~flatMassDependence()                                   { }

		virtual std::complex<double> amp(const isobarDecayVertex&);

		virtual std::string parentLabelForWaveName(const isobarDecayVertex& v) const;  ///< returns label for parent of decay used in wave name

		static constexpr const char* cName = "flat";

	};


	typedef boost::shared_ptr<flatMassDependence> flatMassDependencePtr;


	//////////////////////////////////////////////////////////////////////////////
	/// Brief trivial flat mass dependence over a range
	class binnedMassDependence : public massDependenceImpl<binnedMassDependence> {

	public:

		binnedMassDependence(const double mMin, const double mMax) : massDependenceImpl<binnedMassDependence>(), _mMin(mMin), _mMax(mMax) { }
		virtual ~binnedMassDependence()                                                                                                   { }

		virtual std::complex<double> amp(const isobarDecayVertex&);

		virtual std::string parentLabelForWaveName(const isobarDecayVertex& v) const;  ///< returns label for parent of decay used in wave name

		double getMassMin() const { return _mMin; }
		double getMassMax() const { return _mMax; }

		using massDependenceImpl<binnedMassDependence>::Create;
		static boost::shared_ptr<binnedMassDependence> Create(const libconfig::Setting* settings);

		static constexpr const char* cName = "binned";

	private:
		double _mMin;  ///< Lower limit of the isobar mass bin
		double _mMax;  ///< Upper limit of the isobar mass bin

	};


	typedef boost::shared_ptr<binnedMassDependence> binnedMassDependencePtr;


	//////////////////////////////////////////////////////////////////////////////
	/// Brief relativistic Breit-Wigner with mass-dependent width and Blatt-Weisskopf barrier factors
	class relativisticBreitWigner : public massDependenceImpl<relativisticBreitWigner> {

	public:

		relativisticBreitWigner() : massDependenceImpl<relativisticBreitWigner>() { }
		virtual ~relativisticBreitWigner()           { }

		virtual std::complex<double> amp(const isobarDecayVertex& v);

		virtual std::string parentLabelForWaveName(const isobarDecayVertex& v) const;  ///< returns label for parent of decay used in wave name

		static constexpr const char* cName = "relativisticBreitWigner";

	};


	typedef boost::shared_ptr<relativisticBreitWigner> relativisticBreitWignerPtr;


	//////////////////////////////////////////////////////////////////////////////
	/// Brief relativistic constant-width s-wave Breit-Wigner
	class constWidthBreitWigner : public massDependenceImpl<constWidthBreitWigner> {

	public:

		constWidthBreitWigner() : massDependenceImpl<constWidthBreitWigner>() { }
		virtual ~constWidthBreitWigner()                                      { }

		virtual std::complex<double> amp(const isobarDecayVertex& v);

		static constexpr const char* cName = "constWidthBreitWigner";

	};


	typedef boost::shared_ptr<constWidthBreitWigner> constWidthBreitWignerPtr;


	//////////////////////////////////////////////////////////////////////////////
	/// Brief Breit-Wigner for rho(770) -> pi pi
	/// Breit-Wigner function for L = 1 with the Blatt-Weisskopf barrier
	/// factor replaced by (2 * q^2) / (q^2 + q0^2) so that
	/// Gamma = Gamma0 * m0 / m * (q / q0) * (2 * q^2) / (q^2 + q0^2)
	/// [D. Bisello et al, Phys. Rev. D39 (1989) 701], appendix
	/// http://dx.doi.org/10.1103/PhysRevD.39.701
	class rhoBreitWigner : public massDependenceImpl<rhoBreitWigner> {

	public:

		rhoBreitWigner() : massDependenceImpl<rhoBreitWigner>() { }
		virtual ~rhoBreitWigner()                               { }

		virtual std::complex<double> amp(const isobarDecayVertex& v);

		static constexpr const char* cName = "rhoBreitWigner";

	};


	typedef boost::shared_ptr<rhoBreitWigner> rhoBreitWignerPtr;


	//////////////////////////////////////////////////////////////////////////////
	/// Brief Breit-Wigner for f_0(980) -> pi pi
	/// this is used in piPiSWaveAuMorganPenningtonVes for subtraction of f_0(980)
	/// "Probably this isn't correct S-wave BW form!"
	class f0980BreitWigner : public massDependenceImpl<f0980BreitWigner> {

	public:

		f0980BreitWigner() : massDependenceImpl<f0980BreitWigner>() { }
		virtual ~f0980BreitWigner()                                 { }

		virtual std::complex<double> amp(const isobarDecayVertex& v);

		static constexpr const char* cName = "f0980BreitWigner";

	};


	typedef boost::shared_ptr<f0980BreitWigner> f0980BreitWignerPtr;


	//////////////////////////////////////////////////////////////////////////////
	/// Brief Flatte for f_0(980) -> pi pi
	/// [M. Ablikim et al, Phys. Let. B607, 243] BES II
	class f0980FlatteBesII : public massDependenceImpl<f0980FlatteBesII> {

	public:

		f0980FlatteBesII();
		virtual ~f0980FlatteBesII() { }

		virtual std::complex<double> amp(const isobarDecayVertex& v);

		static constexpr const char* cName = "f0980FlatteBesII";

	private:
		double _piChargedMass;
		double _kaonChargedMass;

	};


	typedef boost::shared_ptr<f0980FlatteBesII> f0980FlatteBesIIPtr;


	//////////////////////////////////////////////////////////////////////////////
	/// Brief Base class for the Au-Morgan-Pennington parameterization of
	///       pi pi s-wave
	/// This class is required to at the same time have a common base class for
	/// all three related mass dependencies and still be able to have access to
	/// the static functions.
	/// See class 'piPiSWaveAuMorganPenningtonM' for references.
	template<class T>
	class piPiSWaveAuMorganPenningtonImpl : public massDependenceImpl<T> {

	public:

		piPiSWaveAuMorganPenningtonImpl();
		virtual ~piPiSWaveAuMorganPenningtonImpl() { }

		virtual std::complex<double> amp(const isobarDecayVertex& v);

	protected:

		ublas::matrix<std::complex<double> >               _T;
		std::vector<ublas::matrix<std::complex<double> > > _a;
		std::vector<ublas::matrix<std::complex<double> > > _c;
		ublas::matrix<double>                              _sP;
		int                                                _vesSheet;

		double _piChargedMass;
		double _piNeutralMass;
		double _kaonChargedMass;
		double _kaonNeutralMass;
		double _kaonMeanMass;

	};


	//////////////////////////////////////////////////////////////////////////////
	/// Brief Au-Morgan-Pennington parameterization of pi pi s-wave
	/// [K.L. Au et al, Phys. Rev. D35, 1633] M solution.
	/// we have introduced a small modification by setting the
	/// off-diagonal elements of the M-matrix to zero.
	class piPiSWaveAuMorganPenningtonM : public piPiSWaveAuMorganPenningtonImpl<piPiSWaveAuMorganPenningtonM> {

	public:

		piPiSWaveAuMorganPenningtonM() : piPiSWaveAuMorganPenningtonImpl<piPiSWaveAuMorganPenningtonM>() { }
		virtual ~piPiSWaveAuMorganPenningtonM()                                                          { }

		static constexpr const char* cName = "piPiSWaveAuMorganPenningtonM";

	};


	typedef boost::shared_ptr<piPiSWaveAuMorganPenningtonM> piPiSWaveAuMorganPenningtonMPtr;


	//////////////////////////////////////////////////////////////////////////////
	/// Brief old VES pi pi s-wave parameterization
	/// [K.L. Au et al, Phys. Rev. D35, 1633] M solution.
	/// brute force subtraction of the f0(980)
	class piPiSWaveAuMorganPenningtonVes : public piPiSWaveAuMorganPenningtonImpl<piPiSWaveAuMorganPenningtonVes> {

	public:

		piPiSWaveAuMorganPenningtonVes();
		virtual ~piPiSWaveAuMorganPenningtonVes() { }

		virtual std::complex<double> amp(const isobarDecayVertex& v);

		static constexpr const char* cName = "piPiSWaveAuMorganPenningtonVes";

	};


	typedef boost::shared_ptr<piPiSWaveAuMorganPenningtonVes> piPiSWaveAuMorganPenningtonVesPtr;


	//////////////////////////////////////////////////////////////////////////////
	/// Brief Kachaev's version of the AMP pi pi s-wave parameterization
	///
	/// from the original fortran code:
	/// [K.L. Au et al, Phys. Rev. D35, 1633] M solution.
	/// 04-Mar-2003 See eps_k1.for for description.
	/// here matrix M = K^{-1} is parametrized with one pole.
	/// misprint in the article (other than in K1--K3 solutions)
	/// was corrected.
	///
	/// 14-Mar-2003 nice amplitude for pi-pi S-wave without f0(975).
	/// it is smooth and nicely tends to zero after approx 1.5 GeV.
	/// f0(975) pole excluded; coupling to KK zeroed; set C411=C422=0.
	/// the largest effect from C411, zeroing of C422 looks insignificant.
	class piPiSWaveAuMorganPenningtonKachaev : public piPiSWaveAuMorganPenningtonImpl<piPiSWaveAuMorganPenningtonKachaev> {

	public:

		piPiSWaveAuMorganPenningtonKachaev();
		virtual ~piPiSWaveAuMorganPenningtonKachaev() { }

		static constexpr const char* cName = "piPiSWaveAuMorganPenningtonKachaev";

	};


	typedef boost::shared_ptr<piPiSWaveAuMorganPenningtonKachaev> piPiSWaveAuMorganPenningtonKachaevPtr;


	//////////////////////////////////////////////////////////////////////////////
	/// combined amplitude for rho(1450)/rho(1700)
	/// [A. Donnachie et al, Z. Phys. C 33 (1987) 407] http://dx.doi.org/10.1007/BF01552547, sec. 4
	class rhoPrimeMassDep : public massDependenceImpl<rhoPrimeMassDep> {

	public:

		rhoPrimeMassDep() : massDependenceImpl<rhoPrimeMassDep>() { }
		virtual ~rhoPrimeMassDep()                                { }

		virtual std::complex<double> amp(const isobarDecayVertex& v);

		static constexpr const char* cName = "rhoPrime";

	};


	typedef boost::shared_ptr<rhoPrimeMassDep> rhoPrimeMassDepPtr;


}  // namespace rpwa


#endif  // MASSDEPENDENCE_H
