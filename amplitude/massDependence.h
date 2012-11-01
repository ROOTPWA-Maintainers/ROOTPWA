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

		virtual std::complex<double> amp(const isobarDecayVertex& v) = 0;

		virtual std::complex<double> operator ()(const isobarDecayVertex& v) { return amp(v); }

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

		virtual std::complex<double> amp(const isobarDecayVertex&);

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
	/// Brief relativistic Breit-Wigner with mass-dependent width and Blatt-Weisskopf barrier factors
	class relativisticBreitWigner : public massDependence {

	public:

		relativisticBreitWigner() : massDependence() { }
		virtual ~relativisticBreitWigner()           { }

		virtual std::complex<double> amp(const isobarDecayVertex& v);

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
	/// Brief Au-Morgan-Pennington parameterization of pi pi s-wave
	/// source: K.L. Au et al, Phys.Rev. D35, P 1633. M solution.
	/// we have introduced a small modification by setting the
	/// off-diagonal elements of the M-matrix to zero.
	class piPiSWaveAuMorganPenningtonM : public massDependence {

	public:

		piPiSWaveAuMorganPenningtonM();
		virtual ~piPiSWaveAuMorganPenningtonM() { }

		virtual std::complex<double> amp(const isobarDecayVertex& v);

		virtual std::string name() const { return "piPiSWaveAuMorganPenningtonM"; }  ///< returns label used in graph visualization, reporting, and key file

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
	/// source: K.L. Au et al, Phys.Rev. D35, P 1633. M solution.
	/// brute force subtraction of the f0(980)
	class piPiSWaveAuMorganPenningtonVes : public piPiSWaveAuMorganPenningtonM {

	public:

		piPiSWaveAuMorganPenningtonVes();
		virtual ~piPiSWaveAuMorganPenningtonVes() { }

		virtual std::complex<double> amp(const isobarDecayVertex& v);

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
	/// source: K.L. Au et al, Phys.Rev. D35, P 1633. M solution.
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
	/// see DOI: 10.1007/BF01552547
	class rhoPrimeMassDep : public massDependence {

	public:

		rhoPrimeMassDep() : massDependence() { }
		virtual ~rhoPrimeMassDep()           { }

		virtual std::complex<double> amp(const isobarDecayVertex& v);

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
