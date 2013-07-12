#include <iostream>
#include <sstream>
#include <complex>

#include "TApplication.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TF1.h"

#include "physUtils.hpp"
#include "particleDataTable.h"
#include "waveDescription.h"


using namespace std;
using namespace rpwa;


////////////////////////////////////////////////////////////////////////////////
// computes relativistic Breit-Wigner amplitude with mass-dependent
// width for 2-body decay using various forms
// !NOTE! L is units of hbar/2
inline
complex<double>
myBreitWigner(const double       M,       // mass
              const double       M0,      // peak position
              const double       Gamma0,  // total width
              const int          L,       // relative orbital angular momentum
              const double       q,       // 2-body breakup momentum
              const double       q0,      // 2-body breakup momentum at peak position
              const unsigned int algo = 0)
{
	if (q0 == 0)
		return 0;
	const double dyn   = (M0 / M) * (q / q0)
		                   * (barrierFactorSquared(L, q) / barrierFactorSquared(L, q0));
	const double Gamma = Gamma0 * dyn;
	switch (algo) {

	case 0: {
		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double A = M0 * Gamma0;
		const double B = M0 * M0 - M * M;
		const double C = M0 * Gamma;
		return (A / (B * B + C * C)) * complex<double>(B, C);
		// return (M0 * Gamma0) / (M0 * M0 - M * M - imag * M0 * Gamma);
	}

	// >>> lies on unitarity circle
	case 1: {
		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double A = M0 * Gamma;
		const double B = M0 * M0 - M * M;
		const double C = M0 * Gamma;
		return (A / (B * B + C * C)) * complex<double>(B, C);
		// return (M0 * Gamma) / (M0 * M0 - M * M - imag * M0 * Gamma);
	}

	case 2: {
		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double A = M * Gamma;
		const double B = M0 * M0 - M * M;
		const double C = M0 * Gamma;
		return (A / (B * B + C * C)) * complex<double>(B, C);
		// return (M * Gamma) / (M0 * M0 - M * M - imag * M0 * Gamma);
	}

	case 3: {
		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double A = sqrt(M * Gamma) * sqrt(M0 * Gamma0);
		const double B = M0 * M0 - M * M;
		const double C = M0 * Gamma;
		return (A / (B * B + C * C)) * complex<double>(B, C);
		// return (sqrt(M * Gamma)) / (M0 * M0 - M * M - imag * M0 * Gamma);
	}

	case 4: {
		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double A = M0 * Gamma0;
		const double B = M0 * M0 - M * M;
		const double C = M0 * Gamma;
		return (A / (B * B + C * C)) * (barrierFactor(L, q) / barrierFactor(L, q0))
			* complex<double>(B, C);
		// return (M0 * Gamma0 * F_L(q) / F_L(q0)) / (M0 * M0 - M * M - imag * M0 * Gamma);
	}

	case 5: {
		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double A = M0 * Gamma0;
		const double B = M0 * M0 - M * M;
		const double C = M0 * Gamma;
		return (A / (B * B + C * C)) * (barrierFactorSquared(L, q) / barrierFactorSquared(L, q0))
			* complex<double>(B, C);
		// return (M0 * Gamma0 * F_L(q) / F_L(q0)) / (M0 * M0 - M * M - imag * M0 * Gamma);
	}

	// >>> lies on unitarity circle
	case 6: {
		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double A = M0 * Gamma0;
		const double B = M0 * M0 - M * M;
		const double C = M0 * Gamma0;
		return (A / (B * B + C * C)) * complex<double>(B, C);
		// return (M0 * Gamma0) / (M0 * M0 - M * M - imag * M0 * Gamma0);
	}

	case 7: {
		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double A = M0 * Gamma0;
		const double B = M0 * M0 - M * M;
		const double C = M0 * Gamma0;
		return (A / (B * B + C * C)) * dyn * complex<double>(B, C);
		// return (M0 * Gamma0 * dyn) / (M0 * M0 - M * M - imag * M0 * Gamma0);
	}
	}

	return 0;
}


////////////////////////////////////////////////////////////////////////////////
// dynamic amplitude for two-body decays with (quasi-)stable daughters
class dynamicAmplitude {

public:

	dynamicAmplitude(const double       mass          = 0,
	                 const double       width         = 0,
	                 const unsigned int L             = 0,
	                 const double       daughterMass0 = 0,
	                 const double       daughterMass1 = 0,
	                 const unsigned int algo          = 0)
		: _mass         (mass),
		  _width        (width),
		  _L            (L),
		  _daughterMass0(daughterMass0),
		  _daughterMass1(daughterMass1),
		  _algo         (algo)
	{ }

	virtual ~dynamicAmplitude() { }

	virtual complex<double> amp(const double m) const
	{
		// const double q  = breakupMomentum(m,     _daughterMass0, _daughterMass1);
		const double q  = sqrt(fabs(breakupMomentumSquared(m, _daughterMass0, _daughterMass1, true)));
		const double q0 = breakupMomentum(_mass, _daughterMass0, _daughterMass1);
		return myBreitWigner(m, _mass, _width, _L, q, q0, _algo);
	}

	double       _mass;
	double       _width;
	unsigned int _L;
	double       _daughterMass0;
	double       _daughterMass1;
	unsigned int _algo;

};  // dynamicAmplitude


// intensity functor
class dynAmpIntensity : public dynamicAmplitude {

public:

	dynAmpIntensity(const double       mass          = 0,
	                const double       width         = 0,
	                const unsigned int L             = 0,
	                const double       daughterMass0 = 0,
	                const double       daughterMass1 = 0,
	                const unsigned int algo          = 0)
		: dynamicAmplitude(mass, width, L, daughterMass0, daughterMass1, algo)
	{ }

	virtual ~dynAmpIntensity() { }

	double operator() (double* x, double*) const { return norm(amp(x[0])); }

};  // dynAmpIntensity


////////////////////////////////////////////////////////////////////////////////
// dynamic amplitude for three-body decay in isobar model where
// integral over isobar phase space is performed numerically
class threeBodyDynAmpInt : public dynamicAmplitude {

public:

	threeBodyDynAmpInt(const dynamicAmplitude& isobarAmp,
	                   const double            mass         = 0,
	                   const double            width        = 0,
	                   const unsigned int      L            = 0,
	                   const double            bachelorMass = 0,
	                   const double            isobarMass   = 0,
	                   const unsigned int      algo         = 1)
		: dynamicAmplitude(mass, width, L, bachelorMass, isobarMass, 1),
		  _isobarAmp(isobarAmp),
		  _algo     (algo)
	{ }
	virtual ~threeBodyDynAmpInt() { }

	virtual double dyn(const double m) const
	{
		switch (_algo) {
		case 1:
			{
				// const double q  = breakupMomentum(m,     _daughterMass0, _daughterMass1);
				const double q  = sqrt(fabs(breakupMomentumSquared(m, _daughterMass0, _daughterMass1, true)));
				const double q0 = breakupMomentum(_mass, _daughterMass0, _daughterMass1);
				if (q0 == 0)
					return 0;
				return (_mass / m) * (q / q0) * (barrierFactorSquared(_L, q) / barrierFactorSquared(_L, q0));
			}
			break;
		case 2:
			return (_mass / m) * (psIntegral(m) / psIntegral(_mass));
			break;
		default:
			return 1;
		}
	}

	virtual complex<double> amp(const double m) const
	{
		const double Gamma = _width * dyn(m);
		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double A = _mass * Gamma;
		const double B = _mass * _mass - m * m;
		const double C = _mass * Gamma;
		return (A / (B * B + C * C)) * complex<double>(B, C);
		// return (M0 * Gamma) / (M0 * M0 - M * M - imag * M0 * Gamma);
	}

	virtual double phaseSpace(double* x,
	                  double* p) const
	{
		const double twoBodyMass   = x[0];
		const double threeBodyMass = p[0];
		// breakup momentum between isobar and bachelor particle
		const double q   = breakupMomentum(threeBodyMass, _daughterMass0, twoBodyMass);
		// breakup momentum between daughter particles of isobar
		const double q12 = breakupMomentum(twoBodyMass, _isobarAmp._daughterMass0, _isobarAmp._daughterMass1);
		const complex<double> bwAmp = _isobarAmp.amp(twoBodyMass);
		return q * barrierFactorSquared(_L, q) * q12 * norm(bwAmp);
	}

	virtual double psIntegral(const double threeBodyMass) const
	{
		TF1 f("f", this, &threeBodyDynAmpInt::phaseSpace, 0, 1, 1, "threeBodyDynAmpInt", "phaseSpace");
		f.SetParameter(0, threeBodyMass);
		return f.Integral(_isobarAmp._daughterMass0 + _isobarAmp._daughterMass1, threeBodyMass - _daughterMass0);
	}

	const dynamicAmplitude& _isobarAmp;
	unsigned int            _algo;

};  // threeBodyDynAmpInt


////////////////////////////////////////////////////////////////////////////////
// dynamic amplitude for three-body decay in isobar model where
// integral over isobar phase space is performed by MC simulation of decay amplitude
class threeBodyDynAmpMc : public threeBodyDynAmpInt {

public:

	threeBodyDynAmpMc(const dynamicAmplitude& isobarAmp,
	                  const double            mass         = 0,
	                  const double            width        = 0,
	                  const unsigned int      L            = 0,
	                  const double            bachelorMass = 0,
	                  const double            isobarMass   = 0,
	                  const unsigned int      algo         = 1)
		: threeBodyDynAmpInt(mass, width, L, bachelorMass, isobarMass, algo)
	{ }
	virtual ~threeBodyDynAmpMc() { }


};  // threeBodyDynAmpMc


void
drawHists(const unsigned int  nmbGraphs,
          const unsigned int* colors,
          TGraph**            intensities,
          TGraph**            argands,
          TGraph**            phases,
          TGraph**            dyn = 0)
{
	if ((not intensities) or (not argands) or (not phases)) {
		printErr << "null pointer: intensities = " << intensities << ", argands = " << argands << ", "
		         << "phases = " << phases << endl;
		throw;
	}
	TCanvas* canv = 0;
	if (dyn) {
		canv = new TCanvas("peak_shape", "Peak Shape", 1980, 470);
		canv->Divide(4, 1);
	} else {
		canv = new TCanvas("peak_shape", "Peak Shape", 1800, 600);
		canv->Divide(3, 1);
	}

	// draw intensity
	canv->cd(1);
	intensities[0]->SetTitle("Intensity");
	intensities[0]->Draw("APC");
	intensities[0]->SetLineColor  (colors[0]);
	intensities[0]->SetMarkerColor(colors[0]);
	for (unsigned int i = 1; i < nmbGraphs; ++i) {
		intensities[i]->SetTitle("Intensity");
		intensities[i]->Draw("PC SAME");
		intensities[i]->SetLineColor  (colors[i]);
		intensities[i]->SetMarkerColor(colors[i]);
	}

	// draw Argand plot
	canv->cd(2);
	argands[0]->SetTitle("Argand Plot");
	argands[0]->Draw("APC");
	argands[0]->SetLineColor  (colors[0]);
	argands[0]->SetMarkerColor(colors[0]);
	for (unsigned int i = 1; i < nmbGraphs; ++i) {
		argands[i]->SetTitle("Argand Plot");
		argands[i]->Draw("PC SAME");
		argands[i]->SetLineColor  (colors[i]);
		argands[i]->SetMarkerColor(colors[i]);
	}
	TEllipse* circ = new TEllipse(0, 0.5, 0.5);
	circ->SetLineStyle(2);
	circ->SetLineColor(16);
	circ->SetFillStyle(0);
	circ->Draw();

	// draw phase
	canv->cd(3);
	phases[0]->SetTitle("Phase");
	phases[0]->Draw("APC");
	phases[0]->SetLineColor  (colors[0]);
	phases[0]->SetMarkerColor(colors[0]);
	for (unsigned int i = 1; i < nmbGraphs; ++i) {
		phases[i]->SetTitle("Phase");
		phases[i]->Draw("PC SAME");
		phases[i]->SetLineColor  (colors[i]);
		phases[i]->SetMarkerColor(colors[i]);
	}

	// draw dynamic factor
	if (dyn) {
		canv->cd(4);
		dyn[0]->SetTitle("Dynamic Factor");
		dyn[0]->Draw("APC");
		dyn[0]->SetLineColor  (colors[0]);
		dyn[0]->SetMarkerColor(colors[0]);
		for (unsigned int i = 1; i < nmbGraphs; ++i) {
			dyn[i]->SetTitle("Dynamic Factor");
			dyn[i]->Draw("PC SAME");
			dyn[i]->SetLineColor  (colors[i]);
			dyn[i]->SetMarkerColor(colors[i]);
		}
	}

	canv->ForceUpdate();
	canv->Flush();
}


int
main()
{
	TApplication app("", 0, 0);

	// test two-body dynamic amplitude
	if (0) {
		// rho(770)
		// const double bwMass          = 0.7665;
		// const double bwWidth         = 0.1502;
		// const double bwL             = 2;
		// const double bwDaughterMass0 = 0.13957018;
		// const double bwDaughterMass1 = bwDaughterMass0;
		// const double massMin         = bwDaughterMass0 + bwDaughterMass1;

		// a2(1320) -> rho [D] pi
		const double bwMass          = 1.3183;
		const double bwWidth         = 0.107;
		const double bwL             = 4;
		const double bwDaughterMass0 = 0.7665;
		const double bwDaughterMass1 = 0.13957018;
		const double massMin         = 3 * bwDaughterMass1;

		const double       massMax   = 3;
		const unsigned int nmbPoints = 1000;

		// const unsigned int algos[]  = {1};
		// const unsigned int nmbAlgos = sizeof(algos) / sizeof(algos[0]);
		// const unsigned int nmbGraphs = nmbAlgos;

		const unsigned int nmbGraphs         = 5;
		const double       daughterMassDelta = 0.05;

		const unsigned int colors[] = {1, 2, 3, 4, 6, 28};
		TGraph* intensities[nmbGraphs];
		TGraph* argands    [nmbGraphs];
		TGraph* phases     [nmbGraphs];

		for (unsigned int iGraph = 0; iGraph < nmbGraphs; ++iGraph) {

			// dynamicAmplitude dynAmp(bwMass, bwWidth, bwL, bwDaughterMass0, bwDaughterMass1, algos[iGraph]);
			dynamicAmplitude dynAmp(bwMass, bwWidth, bwL, bwDaughterMass0 + iGraph * daughterMassDelta,
			                        bwDaughterMass1, 1);

			intensities[iGraph] = new TGraph(nmbPoints);
			argands    [iGraph] = new TGraph(nmbPoints);
			phases     [iGraph] = new TGraph(nmbPoints);

			// calculate amplitude and fill graphs
			const double massStep = (massMax - massMin) / (double)nmbPoints;
			for (unsigned int i = 0; i < nmbPoints; ++i) {

				const double          mass  = massMin + i * massStep;
				const complex<double> amp   = dynAmp.amp(mass);
				const double          phase = arg(amp);
				// phase w.r.t. to center of unitarity circle is identical to
				// phase; at least for amplitudes that lie on the unitarity circle
				// double                phase2 = arg((amp - complex<double>(0, 0.5))
				//                                   * exp(complex<double>(0, piHalf))) / 2;
				// if (phase2 < 0)
				// 	phase2 += pi;
				intensities[iGraph]->SetPoint(i, mass,       norm(amp));
				argands    [iGraph]->SetPoint(i, amp.real(), amp.imag());
				phases     [iGraph]->SetPoint(i, mass,       phase);
			}

		}

		drawHists(nmbGraphs, colors, intensities, argands, phases);

		// test integrability
		if (0) {
			dynAmpIntensity    dynAmpInt(bwMass, bwWidth, bwL, bwDaughterMass0, bwDaughterMass1, 1);
			TF1*               intFunc = new TF1("intFunc", &dynAmpInt, massMin, massMax, 0, "dynAmpIntensity");
			const double       intLimits[]  = {massMax, 10, 100, 1000, 10000};
			const unsigned int nmbIntLimits = sizeof(intLimits) / sizeof(intLimits[0]);
			for (unsigned int i = 0; i < nmbIntLimits; ++i)
				printInfo << "integral in range [" << intFunc->GetXmin() << ", " << intLimits[i] << "] = "
				          << intFunc->Integral(intFunc->GetXmin(), intLimits[i]) << endl;
		}
	}


	// test three-body dynamic amplitude
	if (0) {
		// rho(770) -> pi+ pi-
		const double isobarMass          = 0.7665;
		const double isobarWidth         = 0.1502;
		const double isobarL             = 2;
		const double isobarDaughterMass0 = 0.13957018;
		const double isobarDaughterMass1 = isobarDaughterMass0;
		const dynamicAmplitude isobarDynAmp(isobarMass, isobarWidth, isobarL, isobarDaughterMass0,
		                                    isobarDaughterMass1, 1);

		// a1(1260) -> rho [S] pi
		// const double bwMass       = 1.267;
		// const double bwWidth      = 0.365;
		// const double bwL          = 0;
		// const double bachelorMass = isobarDaughterMass0;
		// const double massMin      = 3 * bachelorMass;

		// a2(1320) -> rho [D] pi
		const double bwMass       = 1.3183;
		const double bwWidth      = 0.107;
		const double bwL          = 4;
		const double bachelorMass = isobarDaughterMass0;
		const double massMin      = 3 * bachelorMass;

		const double       massMax   = 3;
		const unsigned int nmbPoints = 1000;

		const unsigned int algos[]  = {1, 2};
		const unsigned int nmbAlgos = sizeof(algos) / sizeof(algos[0]);
		const unsigned int nmbGraphs = nmbAlgos;

		// const unsigned int nmbGraphs         = 1;
		// const double       daughterMassDelta = 0.05;

		const unsigned int colors[] = {1, 2, 3, 4, 6, 28};
		TGraph* intensities[nmbGraphs];
		TGraph* argands    [nmbGraphs];
		TGraph* phases     [nmbGraphs];
		TGraph* dyn        [nmbGraphs];

		for (unsigned int iGraph = 0; iGraph < nmbGraphs; ++iGraph) {

			threeBodyDynAmpInt dynAmp(isobarDynAmp, bwMass, bwWidth, bwL, bachelorMass, isobarMass, algos[iGraph]);

			intensities[iGraph] = new TGraph(nmbPoints);
			argands    [iGraph] = new TGraph(nmbPoints);
			phases     [iGraph] = new TGraph(nmbPoints);
			dyn        [iGraph] = new TGraph(nmbPoints);

			// calculate amplitude and fill graphs
			const double massStep = (massMax - massMin) / (double)nmbPoints;
			for (unsigned int i = 0; i < nmbPoints; ++i) {

				const double          mass      = massMin + i * massStep;
				const complex<double> amp       = dynAmp.amp(mass);
				const double          phase     = arg(amp);
				const double          dynFactor = dynAmp.dyn(mass);
				intensities[iGraph]->SetPoint(i, mass,       norm(amp));
				argands    [iGraph]->SetPoint(i, amp.real(), amp.imag());
				phases     [iGraph]->SetPoint(i, mass,       phase);
				dyn        [iGraph]->SetPoint(i, mass,       dynFactor);
			}

		}

		drawHists(nmbGraphs, colors, intensities, argands, phases, dyn);
	}

	// gApplication->SetReturnFromRun(kFALSE);
	// gSystem->Run();

	if (1) {
    // initialize particle data table
    particleDataTable::readFile("./particleDataTable.txt");
    stringstream keyFileContent;
    keyFileContent << "productionVertex : {type = \"diffractiveDissVertex\";"
                   << "beam : {name = \"pi-\";}; target : {name = \"p+\";};};" << endl
                   << "decayVertex : {" << endl
                   << "XQuantumNumbers : {isospin =  2; G = +1; J = 2; P = -1; C = -1; M = 0; refl = +1;};" << endl
                   << "XDecay : {fsParticles = ({name = \"pi+\";}, {name = \"pi-\";}); L = 2; S = 0;};" << endl
                   << "};" << endl;

		waveDescription    waveDesc;
		isobarAmplitudePtr amplitude;
		if (   not waveDesc.parseKeyFileContent(keyFileContent.str())
		    or not waveDesc.constructAmplitude(amplitude)) {
			printErr << "problems constructing decay topology from key file content string:" << endl;
			waveDesc.printKeyFileContent(cout, keyFileContent.str());
			cout << "aborting." << endl;
			exit(1);
		}
		printInfo << *amplitude;
	}

	return 0;
}
