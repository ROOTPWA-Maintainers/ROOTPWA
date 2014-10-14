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

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/matrix.hpp>


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

		virtual std::vector<Complex> amp(const isobarDecayVertex& v) = 0;

		virtual std::vector<Complex> operator ()(const isobarDecayVertex& v) { return amp(v); }

		virtual std::string name() const { return "massDependence"; }  ///< returns label used in graph visualization, reporting, and key file

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


	inline
	std::ostream&
	operator <<(std::ostream&         out,
	            const massDependence& massDep)
	{
		return massDep.print(out);
	}


	//////////////////////////////////////////////////////////////////////////////
	/// Brief trivial flat mass dependence
	class flatMassDependence : public massDependence {

	public:

		flatMassDependence() : massDependence() { }
		virtual ~flatMassDependence()           { }

		virtual std::vector<Complex> amp(const isobarDecayVertex&);

		virtual std::string name() const { return "flatMassDependence"; }  ///< returns label used in graph visualization, reporting, and key file

	};


	typedef boost::shared_ptr<flatMassDependence> flatMassDependencePtr;


	inline
	flatMassDependencePtr
	createFlatMassDependence()
	{
		flatMassDependencePtr massDep(new flatMassDependence());
		return massDep;
	}


	//////////////////////////////////////////////////////////////////////////////
	/// Brief trivial flat mass dependence over a range
	class flatRangeMassDependence : public massDependence {

	public:

		flatRangeMassDependence() : massDependence() { }
		virtual ~flatRangeMassDependence()           { }

		virtual std::vector<Complex> amp(const isobarDecayVertex&);

		virtual std::string name() const { return "flatRangeMassDependence"; }  ///< returns label used in graph visualization, reporting, and key file

	};


	typedef boost::shared_ptr<flatRangeMassDependence> flatRangeMassDependencePtr;


	inline
	flatRangeMassDependencePtr
	createFlatRangeMassDependence()
	{
		flatRangeMassDependencePtr massDep(new flatRangeMassDependence());
		return massDep;
	}


	//////////////////////////////////////////////////////////////////////////////
	/// Brief relativistic Breit-Wigner with mass-dependent width and Blatt-Weisskopf barrier factors
	class relativisticBreitWigner : public massDependence {

	public:

		relativisticBreitWigner() : massDependence() { }
		virtual ~relativisticBreitWigner()           { }

		virtual std::vector<Complex> amp(const isobarDecayVertex& v);

		virtual std::string name() const { return "relativisticBreitWigner"; }  ///< returns label used in graph visualization, reporting, and key file

	};


	typedef boost::shared_ptr<relativisticBreitWigner> relativisticBreitWignerPtr;


	inline
	relativisticBreitWignerPtr
	createRelativisticBreitWigner()
	{
		relativisticBreitWignerPtr massDep(new relativisticBreitWigner());
		return massDep;
	}


	//////////////////////////////////////////////////////////////////////////////
	/// Brief relativistic constant-width s-wave Breit-Wigner
	class constWidthBreitWigner : public massDependence {

	public:

		constWidthBreitWigner() : massDependence() { }
		virtual ~constWidthBreitWigner()           { }

		virtual std::vector<Complex> amp(const isobarDecayVertex& v);

		virtual std::string name() const { return "constWidthBreitWigner"; }  ///< returns label used in graph visualization, reporting, and key file

	};


	typedef boost::shared_ptr<constWidthBreitWigner> constWidthBreitWignerPtr;


	inline
	constWidthBreitWignerPtr
	createConstWidthBreitWigner()
	{
		constWidthBreitWignerPtr massDep(new constWidthBreitWigner());
		return massDep;
	}


	//////////////////////////////////////////////////////////////////////////////
	/// Brief Breit-Wigner for rho(770) -> pi pi
	/// Breit-Wigner function for L = 1 with the Blatt-Weisskopf barrier
	/// factor replaced by (2 * q^2) / (q^2 + q0^2) so that
	/// Gamma = Gamma0 * m0 / m * (q / q0) * (2 * q^2) / (q^2 + q0^2)
	/// [D. Bisello et al, Phys. Rev. D39 (1989) 701], appendix
	/// http://dx.doi.org/10.1103/PhysRevD.39.701
	class rhoBreitWigner : public massDependence {

	public:

		rhoBreitWigner() : massDependence() { }
		virtual ~rhoBreitWigner()           { }

		virtual std::vector<Complex> amp(const isobarDecayVertex& v);

		virtual std::string name() const { return "rhoBreitWigner"; }  ///< returns label used in graph visualization, reporting, and key file

	};


	typedef boost::shared_ptr<rhoBreitWigner> rhoBreitWignerPtr;


	inline
	rhoBreitWignerPtr
	createRhoBreitWigner()
	{
		rhoBreitWignerPtr massDep(new rhoBreitWigner());
		return massDep;
	}


	//////////////////////////////////////////////////////////////////////////////
	/// Brief Breit-Wigner for f_0(980) -> pi pi
	/// this is used in piPiSWaveAuMorganPenningtonVes for subtraction of f_0(980)
	/// "Probably this isn't correct S-wave BW form!"
	class f0980BreitWigner : public massDependence {

	public:

		f0980BreitWigner() : massDependence() { }
		virtual ~f0980BreitWigner()           { }

		virtual std::vector<Complex> amp(const isobarDecayVertex& v);

		virtual std::string name() const { return "f0980BreitWigner"; }  ///< returns label used in graph visualization, reporting, and key file

	};


	typedef boost::shared_ptr<f0980BreitWigner> f0980BreitWignerPtr;


	inline
	f0980BreitWignerPtr
	createF0980BreitWigner()
	{
		f0980BreitWignerPtr massDep(new f0980BreitWigner());
		return massDep;
	}


	//////////////////////////////////////////////////////////////////////////////
	/// Brief Au-Morgan-Pennington parameterization of pi pi s-wave
	/// [K.L. Au et al, Phys. Rev. D35, 1633] M solution.
	/// we have introduced a small modification by setting the
	/// off-diagonal elements of the M-matrix to zero.
	class piPiSWaveAuMorganPenningtonM : public massDependence {

	public:

		piPiSWaveAuMorganPenningtonM();
		virtual ~piPiSWaveAuMorganPenningtonM() { }

		virtual std::vector<Complex> amp(const isobarDecayVertex& v);

		virtual std::string name() const { return "piPiSWaveAuMorganPenningtonM"; }  ///< returns label used in graph visualization, reporting, and key file

	protected:

		ublas::matrix<Complex>               _T;
		std::vector<ublas::matrix<Complex> > _a;
		std::vector<ublas::matrix<Complex> > _c;
		ublas::matrix<double>                              _sP;
		int                                                _vesSheet;

		double _piChargedMass;
		double _piNeutralMass;
		double _kaonChargedMass;
		double _kaonNeutralMass;
		double _kaonMeanMass;

	};


	typedef boost::shared_ptr<piPiSWaveAuMorganPenningtonM> piPiSWaveAuMorganPenningtonMPtr;


	inline
	piPiSWaveAuMorganPenningtonMPtr
	createPiPiSWaveAuMorganPenningtonM()
	{
		piPiSWaveAuMorganPenningtonMPtr massDep(new piPiSWaveAuMorganPenningtonM());
		return massDep;
	}


	//////////////////////////////////////////////////////////////////////////////
	/// Brief old VES pi pi s-wave parameterization
	/// [K.L. Au et al, Phys. Rev. D35, 1633] M solution.
	/// brute force subtraction of the f0(980)
	class piPiSWaveAuMorganPenningtonVes : public piPiSWaveAuMorganPenningtonM {

	public:

		piPiSWaveAuMorganPenningtonVes();
		virtual ~piPiSWaveAuMorganPenningtonVes() { }

		virtual std::vector<Complex> amp(const isobarDecayVertex& v);

		virtual std::string name() const { return "piPiSWaveAuMorganPenningtonVes"; }  ///< returns label used in graph visualization, reporting, and key file

	};


	typedef boost::shared_ptr<piPiSWaveAuMorganPenningtonVes> piPiSWaveAuMorganPenningtonVesPtr;


	inline
	piPiSWaveAuMorganPenningtonVesPtr
	createPiPiSWaveAuMorganPenningtonVes()
	{
		piPiSWaveAuMorganPenningtonVesPtr massDep(new piPiSWaveAuMorganPenningtonVes());
		return massDep;
	}


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
	class piPiSWaveAuMorganPenningtonKachaev : public piPiSWaveAuMorganPenningtonM {

	public:

		piPiSWaveAuMorganPenningtonKachaev();
		virtual ~piPiSWaveAuMorganPenningtonKachaev() { }

		virtual std::string name() const { return "piPiSWaveAuMorganPenningtonKachaev"; }  ///< returns label used in graph visualization, reporting, and key file

	};


	typedef boost::shared_ptr<piPiSWaveAuMorganPenningtonKachaev> piPiSWaveAuMorganPenningtonKachaevPtr;


	inline
	piPiSWaveAuMorganPenningtonKachaevPtr
	createPiPiSWaveAuMorganPenningtonKachaev()
	{
		piPiSWaveAuMorganPenningtonKachaevPtr massDep(new piPiSWaveAuMorganPenningtonKachaev());
		return massDep;
	}


	//////////////////////////////////////////////////////////////////////////////
	/// combined amplitude for rho(1450)/rho(1700)
	/// [A. Donnachie et al, Z. Phys. C 33 (1987) 407] http://dx.doi.org/10.1007/BF01552547, sec. 4
	class rhoPrimeMassDep : public massDependence {

	public:

		rhoPrimeMassDep() : massDependence() { }
		virtual ~rhoPrimeMassDep()           { }

		virtual std::vector<Complex> amp(const isobarDecayVertex& v);

		virtual std::string name() const { return "rhoPrimeMassDep"; }  ///< returns label used in graph visualization, reporting, and key file

	};


	typedef boost::shared_ptr<rhoPrimeMassDep> rhoPrimeMassDepPtr;


	inline
	rhoPrimeMassDepPtr
	createRhoPrimeMassDep()
	{
		rhoPrimeMassDepPtr massDep(new rhoPrimeMassDep());
		return massDep;
	}


}  // namespace rpwa


#endif  // MASSDEPENDENCE_H
