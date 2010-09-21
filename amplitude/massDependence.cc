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
//      functor class hierarchy for mass-dependent part of the amplitude
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "utilities.h"
#include "mathUtils.hpp"
#include "isobarDecayVertex.h"
#include "particleDataTable.h"
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
		printInfo << name() << " = 1" << endl;
	return 1;
}


ostream&
flatMassDependence::print(ostream& out) const
{
	out << name();
	return out;
}


////////////////////////////////////////////////////////////////////////////////
complex<double>
relativisticBreitWigner::amp(const isobarDecayVertex& v)
{
	const particlePtr& parent = v.parent();

	const double M      = parent->lzVec().M();          // parent mass
	const double m1     = v.daughter1()->lzVec().M();   // daughter 1 mass
	const double m2     = v.daughter2()->lzVec().M();   // daughter 2 mass
	const double M0     = parent->mass();               // resonance peak position
	const double Gamma0 = parent->width();              // resonance peak width
	const double q      = breakupMomentum(M,  m1, m2);  // breakup momentum
	//const double q0     = breakupMomentum(M0, m1, m2);  // breakup momentum at peak position
	const unsigned int L      = v.L();

	// this is how it is done in PWA2000
	const double M02    = M0 * M0;
	const double m12    = m1 * m1;
	const double m22    = m2 * m2;
	// const double m12    = v.daughter1()->mass() * v.daughter1()->mass();
	// const double m22    = v.daughter2()->mass() * v.daughter2()->mass();
	const double lambda = M02 * M02 + m12 * m12 + m22 * m22 - 2 * (M02 * m12 + m12 * m22 + m22 * M02);
	const double q0     = sqrt(fabs(lambda / (4 * M02)));  //!!! the fabs is probably wrong

	const complex<double> bw = breitWigner(M, M0, Gamma0, L, q, q0);
	if (_debug)
		printInfo << name() << "(m = " << M << " GeV, m_0 = " << M0 << "GeV, "
		          << "Gamma_0 = " << Gamma0 << "GeV, L = " << 0.5 * L << ", q = " << q << "GeV, "
		          << q0 << "GeV) = " << maxPrecisionDouble(bw) << endl;
	return bw;
}


ostream&
relativisticBreitWigner::print(ostream& out) const
{
	out << name();
	return out;
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
	const double f[2] = {0.1968, -0.0154};

	_a[0](0, 0) =  0.1131;
	_a[0](0, 1) =  0.0150;
	_a[0](1, 0) =  0.0150;
	_a[0](1, 1) = -0.3216;
	_a[1](0, 0) = f[0] * f[0];
	_a[1](0, 1) = f[0] * f[1];
	_a[1](1, 0) = f[1] * f[0];
	_a[1](1, 1) = f[1] * f[1];
	
	_c[0](0, 0) =  0.0337;
	_c[1](0, 0) = -0.3185;
	_c[2](0, 0) = -0.0942;
	_c[3](0, 0) = -0.5927;
	_c[4](0, 0) =  0.1957; 
	_c[0](0, 1) = _c[0](1, 0) = -0.2826;
	_c[1](0, 1) = _c[1](1, 0) =  0.0918;
	_c[2](0, 1) = _c[2](1, 0) =  0.1669;
	_c[3](0, 1) = _c[3](1, 0) = -0.2082;
	_c[4](0, 1) = _c[4](1, 0) = -0.1386;
	_c[0](1, 1) =  0.3010;
	_c[1](1, 1) = -0.5140;
	_c[2](1, 1) =  0.1176;
	_c[3](1, 1) =  0.5204;
	_c[4](1, 1) = -0.3977; 

	_sP(0, 0) = -0.0074;
	_sP(0, 1) =  0.9828;
  
	particleDataTable& pdt = particleDataTable::instance();
	const string partList[] = {"pi", "pi0", "K", "K0"};
	for (unsigned int i = 0; i < sizeof(partList) / sizeof(partList[0]); ++i)
		if (not pdt.isInTable(partList[i])) {
			printErr << "cannot find particle " << partList[i] << " in particle data table. "
			         << "was the table initiatlized properly?" << endl;
			throw;
		}
	_piChargedMass   = pdt.entry("pi" )->mass();
	_piNeutralMass   = pdt.entry("pi0")->mass();
	_kaonChargedMass = pdt.entry("K"  )->mass();
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

	const complex<double> qPiPi   = q(mass, _piChargedMass,   _piChargedMass);
	const complex<double> qPi0Pi0 = q(mass, _piNeutralMass,   _piNeutralMass);
	const complex<double> qKK     = q(mass, _kaonChargedMass, _kaonChargedMass);
	const complex<double> qK0K0   = q(mass, _kaonNeutralMass, _kaonNeutralMass);
	complex<double>       qKmKm   = q(mass, _kaonMeanMass,    _kaonMeanMass);

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
	for (unsigned int p = 0; p < _sP.size2(); ++p) {
		const complex<double> fa = 1. / (s - _sP(0, p));
		M += fa * _a[p];
	}
	for (unsigned int n = 0; n < _c.size(); ++n) {
	  const complex<double> sc = pow(scale, (int) n);
		M += sc *_c[n];
	}
	
	// modification: off-diagonal terms set to 0
	M(0, 1) = 0;
	M(1, 0) = 0;

	invertMatrix<complex<double> >(M - imag * rho, _T);
	const complex<double> amp = _T(0, 0);
	if (_debug)
		printInfo << name() << "(m = " << mass << "GeV) = " << maxPrecisionDouble(amp) << endl;

	return amp;
}


ostream&
piPiSWaveAuMorganPenningtonM::print(ostream& out) const
{
	out << name();
	return out;
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
	double mass = v.parent()->lzVec().M();

	const double          f0Mass  = 0.9837;  // [GeV]
	const double          f0Width = 0.0376;  // [GeV]
	const complex<double> coupling(-0.3743, 0.3197);

	const complex<double> ampM = piPiSWaveAuMorganPenningtonM::amp(v);

	complex<double> bw;
	if (mass > 2 * _piChargedMass) {
		const double p     = q(mass,   _piChargedMass, _piChargedMass).real();
		const double p0    = q(f0Mass, _piChargedMass, _piChargedMass).real();
		const double Gamma = f0Width * (p / p0);
		const double A     = f0Mass * f0Mass - mass * mass;
		const double B     = f0Mass * Gamma;
		const double C     = B * (mass / p);
		const double denom = C / (A * A + B * B);
		bw = denom * complex<double>(A, B);
	}

	const complex<double> amp = ampM - coupling * bw;
	if (_debug)
		printInfo << name() << "(m = " << mass << "GeV) = " << maxPrecisionDouble(amp) << endl;

	return amp;
}


ostream&
piPiSWaveAuMorganPenningtonVes::print(ostream& out) const
{
	out << name();
	return out;
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


ostream&
piPiSWaveAuMorganPenningtonKachaev::print(ostream& out) const
{
	out << name();
	return out;
}
