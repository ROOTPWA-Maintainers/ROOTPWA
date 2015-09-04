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
//      functor class hierarchy for mass-dependent part of the amplitude
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <boost/numeric/ublas/io.hpp>

#include "physUtils.hpp"
#include "isobarDecayVertex.h"
#include "particleDataTable.h"
#include "phaseSpaceIntegral.h"
#include "massDependence.h"


using namespace std;
using namespace boost::numeric::ublas;
using namespace rpwa;


////////////////////////////////////////////////////////////////////////////////
bool massDependence::_debug = false;


ostream&
massDependence::print(ostream& out) const
{
	out << name();
	return out;
}


////////////////////////////////////////////////////////////////////////////////
complex<double>
flatMassDependence::amp(const isobarDecayVertex&)
{
	if (_debug)
		printDebug << name() << " = 1" << endl;
	return 1;
}


////////////////////////////////////////////////////////////////////////////////
complex<double>
flatRangeMassDependence::amp(const isobarDecayVertex& v)
{
	complex<double> amp = 0.;

	const particlePtr& parent = v.parent();
	if (fabs(parent->lzVec().M() - parent->mass()) < parent->width()/2.)
		amp = 1.;

	if (_debug)
		printDebug << name() << " M = " << parent->lzVec().M()
		                     << ", M0 = " << parent->mass()
		                     << ", G0 = " << parent->width()
		                     << ", amp = " << amp << endl;

	return amp;
}

////////////////////////////////////////////////////////////////////////////////
complex<double>
steplikeMassDependence::amp(const isobarDecayVertex& v)
{
	const particlePtr& parent = v.parent();
	double mIsob = parent->lzVec().M();
	if (mIsob>_mMin and mIsob<_mMax){
		return complex<double>(1.,0.);
	};
	return complex<double>(0.,0.);
};

////////////////////////////////////////////////////////////////////////////////
complex<double>
relativisticBreitWigner::amp(const isobarDecayVertex& v)
{
	const particlePtr& parent = v.parent();
	const particlePtr& daughter1 = v.daughter1();
	const particlePtr& daughter2 = v.daughter2();

	complex<double> bw;
	if(daughter1->isStable() and daughter2->isStable()) {

		// get Breit-Wigner parameters
		const double       M      = parent->lzVec().M();         // parent mass
		const double       m1     = daughter1->lzVec().M();  // daughter 1 mass
		const double       m2     = daughter2->lzVec().M();  // daughter 2 mass
		const double       q      = breakupMomentum(M,  m1, m2);
		const double       M0     = parent->mass();              // resonance peak position
		const double       q02    = breakupMomentumSquared(M0, m1, m2, true);
		// !NOTE! the following is incorrect but this is how it was done in PWA2000
		const double       q0     = sqrt(fabs(q02));
		const double       Gamma0 = parent->width();             // resonance peak width
		const unsigned int L      = v.L();

		bw = breitWigner(M, M0, Gamma0, L, q, q0);
		if (_debug)
			printDebug << name() << "(m = " << maxPrecision(M) << " GeV/c^2, m_0 = " << maxPrecision(M0)
			           << " GeV/c^2, Gamma_0 = " << maxPrecision(Gamma0) << " GeV/c^2, L = " << spinQn(L)
			           << ", q = " << maxPrecision(q) << " GeV/c, q0 = "
			           << maxPrecision(q0) << " GeV/c) = " << maxPrecisionDouble(bw) << endl;
	} else {

		bw = (*phaseSpaceIntegral::instance())(v);

	}

	return bw;
}


////////////////////////////////////////////////////////////////////////////////
complex<double>
constWidthBreitWigner::amp(const isobarDecayVertex& v)
{
	const particlePtr& parent = v.parent();

	// get Breit-Wigner parameters
	const double M      = parent->lzVec().M();  // parent mass
	const double M0     = parent->mass();       // resonance peak position
	const double Gamma0 = parent->width();      // resonance peak width

	// A / (B - iA) = (A / (B^2 + A^2)) * (B + iA)
	const double          A  = M0 * Gamma0;
	const double          B  = M0 * M0 - M * M;
	const complex<double> bw = (A / (B * B + A * A)) * complex<double>(B, A);
	// const complex<double> bw = (M0 * Gamma0) / (M0 * M0 - M * M - imag * M0 * Gamma0);
	if (_debug)
		printDebug << name() << "(m = " << maxPrecision(M) << " GeV/c^2, m_0 = " << maxPrecision(M0)
		           << " GeV/c^2, Gamma_0 = " << maxPrecision(Gamma0) << " GeV/c^2) = "
		           << maxPrecisionDouble(bw) << endl;
	return bw;
}


////////////////////////////////////////////////////////////////////////////////
complex<double>
rhoBreitWigner::amp(const isobarDecayVertex& v)
{
	const particlePtr& parent = v.parent();

	// get Breit-Wigner parameters
	const double M      = parent->lzVec().M();         // parent mass
	const double m1     = v.daughter1()->lzVec().M();  // daughter 1 mass
	const double m2     = v.daughter2()->lzVec().M();  // daughter 2 mass
	const double q2     = breakupMomentumSquared(M,  m1, m2);
	const double q      = sqrt(q2);
	const double M0     = parent->mass();              // resonance peak position
	const double q02    = breakupMomentumSquared(M0, m1, m2);
	const double q0     = sqrt(q02);
	const double Gamma0 = parent->width();             // resonance peak width

	const double F      = 2 * q2 / (q02 + q2);
	const double Gamma  = Gamma0 * (M0 / M) * (q / q0) * F;
	// in the original publication the width reads
	// Gamma = Gamma0 * (q / q0) * F

	// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
	const double          A  = M0 * Gamma0 * sqrt(F);
	// in the original publication A reads
	// A = sqrt(M0 * Gamma0 * (m / q0) * F)
	const double          B  = M0 * M0 - M * M;
	const double          C  = M0 * Gamma;
	const complex<double> bw = (A / (B * B + C * C)) * std::complex<double>(B, C);
	// return (M0 * Gamma0 * sqrt(F)) / (M0 * M0 - M * M - imag * M0 * Gamma);
	if (_debug)
		printDebug << name() << "(m = " << maxPrecision(M) << " GeV/c^2, m_0 = " << maxPrecision(M0)
		           << " GeV/c^2, Gamma_0 = " << maxPrecision(Gamma0) << " GeV/c^2, "
		           << "q = " << maxPrecision(q) << " GeV/c, q0 = " << maxPrecision(q0) << " GeV/c) "
		           << "= " << maxPrecisionDouble(bw) << endl;
	return bw;
}


////////////////////////////////////////////////////////////////////////////////
complex<double>
f0980BreitWigner::amp(const isobarDecayVertex& v)
{
	const particlePtr& parent = v.parent();

	// get Breit-Wigner parameters
	const double M      = parent->lzVec().M();         // parent mass
	const double m1     = v.daughter1()->lzVec().M();  // daughter 1 mass
	const double m2     = v.daughter2()->lzVec().M();  // daughter 2 mass
	const double q      = breakupMomentum(M,  m1, m2);
	const double M0     = parent->mass();              // resonance peak position
	const double q0     = breakupMomentum(M0, m1, m2);
	const double Gamma0 = parent->width();             // resonance peak width

	const double Gamma  = Gamma0 * (q / q0);

	// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
	const double          C  = M0 * Gamma;
	const double          A  = C * M / q;
	const double          B  = M0 * M0 - M * M;
	const complex<double> bw = (A / (B * B + C * C)) * std::complex<double>(B, C);
	// return ((M0 * Gamma0 * M / q) / (M0 * M0 - M * M - imag * M0 * Gamma);
	if (_debug)
		printDebug << name() << "(m = " << maxPrecision(M) << " GeV/c^2, m_0 = " << maxPrecision(M0)
		           << " GeV/c^2, Gamma_0 = " << maxPrecision(Gamma0) << " GeV/c^2, "
		           << "q = " << maxPrecision(q) << " GeV/c, q0 = " << maxPrecision(q0) << " GeV/c) "
		           << "= " << maxPrecisionDouble(bw) << endl;
	return bw;
}


////////////////////////////////////////////////////////////////////////////////
f0980Flatte::f0980Flatte()
	: massDependence()
{
	particleDataTable& pdt = particleDataTable::instance();
	const string partList[] = {"pi+", "pi0", "K+", "K0"};
	for (unsigned int i = 0; i < sizeof(partList) / sizeof(partList[0]); ++i)
		if (not pdt.isInTable(partList[i])) {
			printErr << "cannot find particle " << partList[i] << " in particle data table. "
			         << "Aborting..." << endl;
			throw;
		}
	_piChargedMass   = pdt.entry("pi+")->mass();
	_kaonChargedMass = pdt.entry("K+" )->mass();
}


complex<double>
f0980Flatte::amp(const isobarDecayVertex& v)
{
	const complex<double> imag(0, 1);

	const particlePtr& parent = v.parent();

	// get masses of states involved
	const double M      = parent->lzVec().M();         // parent mass

	// Flatte parameters
	const double M0   = 0.965;
	const double g1   = 0.165;
	const double g2g1 = 4.21;

	// break-up momenta
	const std::complex<double> kPi = breakupMomentumComplex(M, _piChargedMass,   _piChargedMass);
	const std::complex<double> kK  = breakupMomentumComplex(M, _kaonChargedMass, _kaonChargedMass);

	// phase space factors
	const std::complex<double> rhoPi = 2. * kPi / M;
	const std::complex<double> rhoK  = 2. * kK  / M;

	std::complex<double> denom(M0*M0 - M*M, 0.);

	denom -= imag * g1 * (rhoPi + g2g1*rhoK);

	const std::complex<double> amp(1. / denom);

	if (_debug)
		printDebug << name() << "(m = " << maxPrecision(M) << " GeV/c^2, m_0 = " << maxPrecision(M0)
		           << " GeV/c^2, g1 = " << maxPrecision(g1) << " GeV/c^2, g2/g1 = " << maxPrecision(g2g1)
		           << ") = " << maxPrecisionDouble(amp) << endl;
	return amp;
}


////////////////////////////////////////////////////////////////////////////////
piPiSWaveAuMorganPenningtonM::piPiSWaveAuMorganPenningtonM()
	: massDependence(),
	  _T       (2, 2),
	  _a       (2, matrix<complex<double> >(2, 2)),
	  _c       (5, matrix<complex<double> >(2, 2)),
	  _sP      (1, 2),
	  _vesSheet(0)
{
	const double f[2] = {0.1968, -0.0154};  // AMP Table 1, M solution: f_1^1 and f_2^1

	_a[0](0, 0) =  0.1131;  // AMP Table 1, M solution: f_2^2
	_a[0](0, 1) =  0.0150;  // AMP Table 1, M solution: f_1^3
	_a[0](1, 0) =  0.0150;  // AMP Table 1, M solution: f_1^3
	_a[0](1, 1) = -0.3216;  // AMP Table 1, M solution: f_2^3
	_a[1](0, 0) = f[0] * f[0];
	_a[1](0, 1) = f[0] * f[1];
	_a[1](1, 0) = f[1] * f[0];
	_a[1](1, 1) = f[1] * f[1];

	_c[0](0, 0) =  0.0337;                // AMP Table 1, M solution: c_11^0
	_c[1](0, 0) = -0.3185;                // AMP Table 1, M solution: c_11^1
	_c[2](0, 0) = -0.0942;                // AMP Table 1, M solution: c_11^2
	_c[3](0, 0) = -0.5927;                // AMP Table 1, M solution: c_11^3
	_c[4](0, 0) =  0.1957;                // AMP Table 1, M solution: c_11^4
	_c[0](0, 1) = _c[0](1, 0) = -0.2826;  // AMP Table 1, M solution: c_12^0
	_c[1](0, 1) = _c[1](1, 0) =  0.0918;  // AMP Table 1, M solution: c_12^1
	_c[2](0, 1) = _c[2](1, 0) =  0.1669;  // AMP Table 1, M solution: c_12^2
	_c[3](0, 1) = _c[3](1, 0) = -0.2082;  // AMP Table 1, M solution: c_12^3
	_c[4](0, 1) = _c[4](1, 0) = -0.1386;  // AMP Table 1, M solution: c_12^4
	_c[0](1, 1) =  0.3010;                // AMP Table 1, M solution: c_22^0
	_c[1](1, 1) = -0.5140;                // AMP Table 1, M solution: c_22^1
	_c[2](1, 1) =  0.1176;                // AMP Table 1, M solution: c_22^2
	_c[3](1, 1) =  0.5204;                // AMP Table 1, M solution: c_22^3
	_c[4](1, 1) = -0.3977;                // AMP Table 1, M solution: c_22^4

	_sP(0, 0) = -0.0074;  // AMP Table 1, M solution: s_0
	_sP(0, 1) =  0.9828;  // AMP Table 1, M solution: s_1

	particleDataTable& pdt = particleDataTable::instance();
	const string partList[] = {"pi+", "pi0", "K+", "K0"};
	for (unsigned int i = 0; i < sizeof(partList) / sizeof(partList[0]); ++i)
		if (not pdt.isInTable(partList[i])) {
			printErr << "cannot find particle " << partList[i] << " in particle data table. "
			         << "Aborting..." << endl;
			throw;
		}
	_piChargedMass   = pdt.entry("pi+")->mass();
	_piNeutralMass   = pdt.entry("pi0")->mass();
	_kaonChargedMass = pdt.entry("K+" )->mass();
	_kaonNeutralMass = pdt.entry("K0" )->mass();
	_kaonMeanMass    = (_kaonChargedMass + _kaonNeutralMass) / 2;
}


complex<double>
piPiSWaveAuMorganPenningtonM::amp(const isobarDecayVertex& v)
{
	const complex<double> imag(0, 1);

	double mass = v.parent()->lzVec().M();
	double s    = mass * mass;
	if (fabs(s - _sP(0, 1)) < 1e-6) {
		mass += 1e-6;
		s     = mass * mass;
	}

	const complex<double> qPiPi   = breakupMomentumComplex(mass, _piChargedMass,   _piChargedMass  );
	const complex<double> qPi0Pi0 = breakupMomentumComplex(mass, _piNeutralMass,   _piNeutralMass  );
	const complex<double> qKK     = breakupMomentumComplex(mass, _kaonChargedMass, _kaonChargedMass);
	const complex<double> qK0K0   = breakupMomentumComplex(mass, _kaonNeutralMass, _kaonNeutralMass);
	complex<double>       qKmKm   = breakupMomentumComplex(mass, _kaonMeanMass,    _kaonMeanMass   );

	matrix<complex<double> > rho(2, 2);
	if (_vesSheet) {
		if (qKmKm.imag() > 0)
			qKmKm *= -1;
		rho(0, 0) = (2. * qPiPi) / mass;
		rho(1, 1) = (2. * qKmKm) / mass;
	} else {
		rho(0, 0) = ((2. * qPiPi) / mass + (2. * qPi0Pi0) / mass) / 2.;
		rho(1, 1) = ((2. * qKK)   / mass + (2. * qK0K0)   / mass) / 2.;
	}
	rho(0, 1) = rho(1, 0) = 0;

	const double scale = (s / (4 * _kaonMeanMass * _kaonMeanMass)) - 1;

	matrix<complex<double> > M(zero_matrix<complex<double> >(2, 2));
	for (unsigned int i = 0; i < _sP.size2(); ++i) {
		const complex<double> fa = 1. / (s - _sP(0, i));
		M += fa * _a[i];
	}
	for (unsigned int i = 0; i < _c.size(); ++i) {
		const complex<double> sc = pow(scale, (int)i);
		M += sc *_c[i];
	}

	// modification: off-diagonal terms set to 0
	M(0, 1) = 0;
	M(1, 0) = 0;

	invertMatrix<complex<double> >(M - imag * rho, _T);
	const complex<double> amp = _T(0, 0);
	if (_debug)
		printDebug << name() << "(m = " << maxPrecision(mass) << " GeV) = "
		           << maxPrecisionDouble(amp) << endl;

	return amp;
}


////////////////////////////////////////////////////////////////////////////////
piPiSWaveAuMorganPenningtonVes::piPiSWaveAuMorganPenningtonVes()
	: piPiSWaveAuMorganPenningtonM()
{
	_vesSheet = 1;
}


complex<double>
piPiSWaveAuMorganPenningtonVes::amp(const isobarDecayVertex& v)
{
	const double M = v.parent()->lzVec().M();

	const double          f0Mass  = 0.9837;  // [GeV]
	const double          f0Width = 0.0376;  // [GeV]
	const complex<double> coupling(-0.3743, 0.3197);

	const complex<double> ampM = piPiSWaveAuMorganPenningtonM::amp(v);

	complex<double> bw = 0;
	if (M > 2 * _piChargedMass) {
		const double q     = breakupMomentum(M,      _piChargedMass, _piChargedMass);
		const double q0    = breakupMomentum(f0Mass, _piChargedMass, _piChargedMass);
		const double Gamma = f0Width * (q / q0);
		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double C     = f0Mass * Gamma;
		const double A     = C * (M / q);
		const double B     = f0Mass * f0Mass - M * M;
		bw = (A / (B * B + C * C)) * complex<double>(B, C);
	}

	const complex<double> amp = ampM - coupling * bw;
	if (_debug)
		printDebug << name() << "(m = " << maxPrecision(M) << " GeV) = "
		           << maxPrecisionDouble(amp) << endl;

	return amp;
}


////////////////////////////////////////////////////////////////////////////////
piPiSWaveAuMorganPenningtonKachaev::piPiSWaveAuMorganPenningtonKachaev()
	: piPiSWaveAuMorganPenningtonM()
{
	// change parameters according to Kachaev's prescription
	_c[4](0, 0) = 0; // was 0.1957;
	_c[4](1, 1) = 0; // was -0.3977;

	_a[0](0, 1) = 0; // was 0.0150
	_a[0](1, 0) = 0; // was 0.0150

	// _a[1] are the f's from the AMP paper
	_a[1](0, 0) = 0;
	_a[1](0, 1) = 0;
	_a[1](1, 0) = 0;
	_a[1](1, 1) = 0;
}


////////////////////////////////////////////////////////////////////////////////
complex<double>
rhoPrimeMassDep::amp(const isobarDecayVertex& v)
{

	// get Breit-Wigner parameters
	const double M   = v.parent()->lzVec().M();                 // parent mass

	// rho' parameters
	const double M01     = 1.465;  // rho(1450) mass [GeV/c^]
	const double Gamma01 = 0.235;  // rho(1450) width [GeV/c^]
	const double M02     = 1.700;  // rho(1700) mass [GeV/c^]
	const double Gamma02 = 0.220;  // rho(1700) width [GeV/c^]

	// const complex<double> bw1 = breitWigner(M, M01, Gamma01, L, q, q0);
	// const complex<double> bw2 = breitWigner(M, M02, Gamma02, L, q, q0);

	// A / (B - iA) = (A / (B^2 + A^2)) * (B + iA)
	const double          A1  = M01 * Gamma01;
	const double          B1  = M01 * M01 - M * M;
	const complex<double> bw1 = (A1 / (B1 * B1 + A1 * A1)) * complex<double>(B1, A1);
	// const complex<double> bw = (M0 * Gamma0) / (M0 * M0 - M * M - imag * M0 * Gamma0);

	// A / (B - iA) = (A / (B^2 + A^2)) * (B + iA)
	const double          A2  = M02 * Gamma02;
	const double          B2  = M02 * M02 - M * M;
	const complex<double> bw2 = (A2 / (B2 * B2 + A2 * A2)) * complex<double>(B2, A2);
	// const complex<double> bw = (M0 * Gamma0) / (M0 * M0 - M * M - imag * M0 * Gamma0);

	const complex<double> bw = (4 * bw1 - 3 * bw2) / 7;

	if (_debug)
		printDebug << name() << "(m = " << maxPrecision(M) << " GeV/c^2, "
		           << "GeV/c) = " << maxPrecisionDouble(bw) << endl;

	return bw;

}
