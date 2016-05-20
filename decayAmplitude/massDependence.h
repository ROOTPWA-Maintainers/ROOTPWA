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

#include <TFormula.h>


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

		virtual std::string name() const { return "flat"; }  ///< returns label used in graph visualization, reporting, and key file

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
	class binnedMassDependence : public massDependence {

	public:

		binnedMassDependence(const double mMin, const double mMax) : massDependence(), _mMin(mMin), _mMax(mMax) { }
		virtual ~binnedMassDependence()                                                                         { }

		virtual std::complex<double> amp(const isobarDecayVertex&);

		virtual std::string name() const { return "binned"; }  ///< returns label used in graph visualization, reporting, and key file

		double getMassMin() const { return _mMin; }
		double getMassMax() const { return _mMax; }

	private:

		double _mMin;  ///< Lower limit of the isobar mass bin
		double _mMax;  ///< Upper limit of the isobar mass bin

	};


	typedef boost::shared_ptr<binnedMassDependence> binnedMassDependencePtr;


	inline
	binnedMassDependencePtr
	createBinnedMassDependence(const double mMin, const double mMax)
	{
		binnedMassDependencePtr massDep(new binnedMassDependence(mMin, mMax));
		return massDep;
	}


	//////////////////////////////////////////////////////////////////////////////
	/// Brief polynomial mass dependence with flexible coefficients
	class polynomialMassDependence : public massDependence {

	public:

		polynomialMassDependence(const std::vector<std::complex<double> >& coefficients, const double mMin, const double mMax)
			: massDependence(),
			  _mMin(mMin),
			  _mMax(mMax),
			  _coefficients(coefficients) { }
		virtual ~polynomialMassDependence()   { }

		virtual std::complex<double> amp(const isobarDecayVertex&);

		virtual std::string name() const { return "polynomial"; }  ///< returns label used in graph visualization, reporting, and key file

		const std::vector<std::complex<double> >& getCoefficients() const { return _coefficients; }
		double                                    getMassMin()      const { return _mMin;         }
		double                                    getMassMax()      const { return _mMax;         }

	private:

		double                             _mMin;          ///< Mass point that gets scaled to -1. // mMin = -1, mMax = 1 -> no scaling
		double                             _mMax;          ///< Mass point that gets scaled to  1.

		std::vector<std::complex<double> > _coefficients;  ///< Coefficients of the polynomial

	};


	typedef boost::shared_ptr<polynomialMassDependence> polynomialMassDependencePtr;


	inline
	polynomialMassDependencePtr
	createPolynomialMassDependencePtr(const std::vector<std::complex<double> >& coefficients, const double mMin, const double mMax)
	{
		polynomialMassDependencePtr massDep(new polynomialMassDependence(coefficients, mMin, mMax));
		return massDep;
	}


	//////////////////////////////////////////////////////////////////////////////
	/// Brief complex exponential mass dependence with certain degree
	class complexExponentialMassDependence : public massDependence {

	public:

		complexExponentialMassDependence(const int degree, const double mMin, const double mMax) : massDependence(), _degree(degree), _mMin(mMin), _mMax(mMax) { }
		virtual ~complexExponentialMassDependence()                                                                                                            { }

		virtual std::complex<double> amp(const isobarDecayVertex&);

		virtual std::string name() const { return "complexExponential"; }  ///< returns label used in graph visualization, reporting, and key file

		int    getDegree()  const { return _degree; }
		double getMassMin() const { return _mMin;   }
		double getMassMax() const { return _mMax;   }

	private:

		int    _degree;  ///< Degree of the expansion
		double _mMin;    ///< Mass point that gets scaled to 0.
		double _mMax;    ///< Mass point that gets scaled to 2.

	};


	typedef boost::shared_ptr<complexExponentialMassDependence> complexExponentialMassDependencePtr;


	inline
	complexExponentialMassDependencePtr
	createComplexExponentialMassDependencePtr(const int degree, const double mMin, const double mMax)
	{
		complexExponentialMassDependencePtr massDep(new complexExponentialMassDependence(degree, mMin, mMax));
		return massDep;
	}


	//////////////////////////////////////////////////////////////////////////////
	/// Brief complex TF1 mass dependence with arbitrary function
	class arbitraryFunctionMassDependence : public massDependence {

	public:

		arbitraryFunctionMassDependence(const std::string& name, const std::string& realPart, const std::string& imagPart)
			: massDependence(),
			  _name(name),
			  _realFunctionString(realPart),
			  _imagFunctionString(imagPart),
			  _realPart((name+"_real").c_str(), realPart.c_str()),
			  _imagPart((name+"_imag").c_str(), imagPart.c_str()) { }
		virtual ~arbitraryFunctionMassDependence()                    { }

		virtual std::complex<double> amp(const isobarDecayVertex&);

		virtual std::string name() const { return "arbitraryFunction"; }  ///< returns label used in graph visualization, reporting, and key file

		const std::string& getName()               const { return _name;               }
		const std::string& getRealFunctionString() const { return _realFunctionString; }
		const std::string& getImagFunctionString() const { return _imagFunctionString; }

	private:

		std::string _name;                ///< Function name
		std::string _realFunctionString;  ///< Function string of the real part
		std::string _imagFunctionString;  ///< Function string of the imag part
		TFormula    _realPart;            ///< Function describing the real part
		TFormula    _imagPart;            ///< Function describing the imag part

	};


	typedef boost::shared_ptr<arbitraryFunctionMassDependence> arbitraryFunctionMassDependencePtr;


	inline
	arbitraryFunctionMassDependencePtr
	createArbitraryFunctionMassDependencePtr(const std::string& name, const std::string& realPart, const std::string& imagPart)
	{
		arbitraryFunctionMassDependencePtr massDep(new arbitraryFunctionMassDependence(name, realPart, imagPart));
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
	/// Brief relativistic constant-width s-wave Breit-Wigner
	class constWidthBreitWigner : public massDependence {

	public:

		constWidthBreitWigner() : massDependence() { }
		virtual ~constWidthBreitWigner()           { }

		virtual std::complex<double> amp(const isobarDecayVertex& v);

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

		virtual std::complex<double> amp(const isobarDecayVertex& v);

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

		virtual std::complex<double> amp(const isobarDecayVertex& v);

		virtual std::string name() const { return "f_0(980)"; }  ///< returns label used in graph visualization, reporting, and key file

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
	/// Brief Flatte for f_0(980) -> pi pi
	/// [M. Ablikim et al, Phys. Let. B607, 243] BES II
	class f0980Flatte : public massDependence {

	public:

		f0980Flatte();
		virtual ~f0980Flatte()           { }

		virtual std::complex<double> amp(const isobarDecayVertex& v);

		virtual std::string name() const { return "f_0(980)Flatte"; }  ///< returns label used in graph visualization, reporting, and key file

	private:
		double _piChargedMass;
		double _kaonChargedMass;

	};


	typedef boost::shared_ptr<f0980Flatte> f0980FlattePtr;


	inline
	f0980FlattePtr
	createF0980Flatte()
	{
		f0980FlattePtr massDep(new f0980Flatte());
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
	/// [K.L. Au et al, Phys. Rev. D35, 1633] M solution.
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

		virtual std::complex<double> amp(const isobarDecayVertex& v);

		virtual std::string name() const { return "rhoPrime"; }  ///< returns label used in graph visualization, reporting, and key file

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
