#include <iostream>
#include <complex>

#include "TGraph.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TF1.h"


using namespace std;


const double pi     = 2 * asin((double)1);
const double piHalf = pi / 2;
const double twoPi  = 2 * pi;
const double fourPi = 4 * pi;


// computes squared breakup momentum of 2-body decay
inline
double
breakupMomentumSquared(const double M,   // mass of mother particle
                       const double m1,  // mass of daughter particle 1
                       const double m2,  // mass of daughter particle 2
                       const bool   allowSubThr = false)  // if set sub-threshold decays with negative return values are allowed
{
	const double mSum  = m1 + m2;
	if (not allowSubThr and (M < mSum)) {
		cerr << "mother mass " << M << " GeV/c^2 is smaller than sum of daughter masses "
		     << m1 << " + " << m2 << " GeV/c^2. this should never happen. aborting." << endl;
		throw;
	}
	const double mDiff = m1 - m2;
	return (M - mSum) * (M + mSum) * (M - mDiff) * (M + mDiff) / (4 * M * M);
}


// computes breakup momentum of 2-body decay
inline
double
breakupMomentum(const double M,   // mass of mother particle
                const double m1,  // mass of daughter particle 1
                const double m2)  // mass of daughter particle 2
{
	return sqrt(fabs(breakupMomentumSquared(M, m1, m2, true)));
	// return sqrt(breakupMomentumSquared(M, m1, m2, false));
}


// computes square of Blatt-Weisskopf barrier factor for 2-body decay
// !NOTE! L is units of hbar/2
inline
double
barrierFactorSquared(const int    L,               // relative orbital angular momentum
                     const double breakupMom,      // breakup momentum of 2-body decay [GeV/c]
                     const bool   debug = false,
                     const double Pr    = 0.1973)  // momentum scale 0.1973 GeV/c corresponds to 1 fm interaction radius
{
	const double z   = (breakupMom * breakupMom) / (Pr * Pr);
	double       bf2 = 0;
	switch (L) {
	case 0:  // L = 0
		bf2 = 1;
		break;
	case 2:  // L = 1
		bf2 = (2 * z) / (z + 1);
		break;
	case 4:  // L = 2
		bf2 = (13 * z * z) / (z * (z + 3) + 9);
		break;
	case 6:  // L = 3
		bf2 = (277 * z * z * z) / (z * (z * (z + 6) + 45) + 225);
		break;
	case 8:  // L = 4
		{
			const double z2 = z * z;
			bf2 = (12746 * z2 * z2) / (z * (z * (z * (z + 10) + 135) + 1575) + 11025);
		}
		break;
	case 10:  // L = 5
		{
			const double z2 = z * z;
			bf2 = (998881 * z2 * z2 * z)
				/ (z * (z * (z * (z * (z + 15) + 315) + 6300) + 99225) + 893025);
		}
		break;
	case 12:  // L = 6
		{
			const double z3 = z * z * z;
			bf2 = (118394977 * z3 * z3)
				/ (z * (z * (z * (z * (z * (z + 21) + 630) + 18900) + 496125) + 9823275) + 108056025);
		}
		break;
	case 14:  // L = 7
		{
			const double z3 = z * z * z;
			bf2 = (19727003738LL * z3 * z3 * z)
				/ (z * (z * (z * (z * (z * (z * (z + 28) + 1134) + 47250) + 1819125) + 58939650)
				        + 1404728325L) + 18261468225LL);
		}
		break;
	default:
		cout << "calculation of Blatt-Weisskopf barrier factor is not (yet) implemented for L = "
		     << L << ". returning 0." << endl;
		return 0;
	}
	if (debug)
		cout << "squared Blatt-Weisskopf barrier factor(L = " << L << ", "
		     << "q = " << breakupMom << " GeV/c; P_r = " << Pr << " GeV/c) = "
		     << bf2 << endl;
	return bf2;
}


// computes Blatt-Weisskopf barrier factor for 2-body decay
// !NOTE! L is units of hbar/2
inline
double
barrierFactor(const int    L,               // relative orbital angular momentum
              const double breakupMom,      // breakup momentum of 2-body decay [GeV/c]
              const bool   debug = false,
              const double Pr    = 0.1973)  // momentum scale 0.1973 GeV/c corresponds to 1 fm interaction radius
{
	const double bf = sqrt(barrierFactorSquared(L, breakupMom, false, Pr));
	if (debug)
		cout << "Blatt-Weisskopf barrier factor(L = " << L << ", "
		     << "q = " << breakupMom << " GeV/c; P_r = " << Pr << " GeV/c) = "
		     << bf << endl;
	return bf;
}


// computes relativistic Breit-Wigner amplitude with mass-dependent width for 2-body decay
// !NOTE! L is units of hbar/2
inline
complex<double>
breitWigner(const double       M,       // mass
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
		const double q  = breakupMomentum(m,     _daughterMass0, _daughterMass1);
		const double q0 = breakupMomentum(_mass, _daughterMass0, _daughterMass1);
		return breitWigner(m, _mass, _width, _L, q, q0, _algo);
	}

	double       _mass;
	double       _width;
	unsigned int _L;
	double       _daughterMass0;
	double       _daughterMass1;
	unsigned int _algo;

};


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

};


// dynamic amplitude for three-body decay in isobar model
class threeBodyDynAmp : public dynamicAmplitude {

public:

	threeBodyDynAmp(const dynamicAmplitude& isobarAmp,
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
	virtual ~threeBodyDynAmp() { }

	double dyn(const double m) const
	{
		switch (_algo) {
		case 1:
			{
				const double q  = breakupMomentum(m,     _daughterMass0, _daughterMass1);
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

	double phaseSpace(double* x,
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

	double psIntegral(const double threeBodyMass) const
	{
		TF1 f("f", this, &threeBodyDynAmp::phaseSpace, 0, 1, 1, "threeBodyDynAmp", "phaseSpace");
		f.SetParameter(0, threeBodyMass);
		return f.Integral(_isobarAmp._daughterMass0 + _isobarAmp._daughterMass1, threeBodyMass - _daughterMass0);
	}

	const dynamicAmplitude& _isobarAmp;
	unsigned int            _algo;

};


void
testPeakShapeInstableIsobar()
{

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

		// draw amplitude
		TCanvas* canv = new TCanvas("peak_shape", "Peak Shape", 1800, 600);
		canv->Divide(3, 1);
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

		canv->ForceUpdate();
		canv->Flush();

		// test integrability
		if (0) {
			dynAmpIntensity    dynAmpInt(bwMass, bwWidth, bwL, bwDaughterMass0, bwDaughterMass1, 1);
			TF1*               intFunc = new TF1("intFunc", &dynAmpInt, massMin, massMax, 0, "dynAmpIntensity");
			const double       intLimits[]  = {massMax, 10, 100, 1000, 10000};
			const unsigned int nmbIntLimits = sizeof(intLimits) / sizeof(intLimits[0]);
			for (unsigned int i = 0; i < nmbIntLimits; ++i)
				cout << "!!! integral in range [" << intFunc->GetXmin() << ", " << intLimits[i] << "] = "
				     << intFunc->Integral(intFunc->GetXmin(), intLimits[i]) << endl;
		}
	}


	// test three-body dynamic amplitude
	if (1) {
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
		// const double bwMass       = 1.3183;
		const double bwMass       = 0.8;
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

			threeBodyDynAmp dynAmp(isobarDynAmp, bwMass, bwWidth, bwL, bachelorMass, isobarMass, algos[iGraph]);

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

		// draw amplitude
		TCanvas* canv = new TCanvas("peak_shape", "Peak Shape", 2000, 500);
		canv->Divide(4, 1);
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

		canv->ForceUpdate();
		canv->Flush();

	}

}
