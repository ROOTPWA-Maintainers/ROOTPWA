
#include<TH2D.h>
#include<TCanvas.h>
#include<TApplication.h>
#include<physUtils.hpp>
#include<complex>
#include<TSystem.h>


#include <iostream>
#include <sstream>
#include <complex>

#include <boost/progress.hpp>

#include "TApplication.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TH2.h"

#include "physUtils.hpp"
#include "particleDataTable.h"
#include "waveDescription.h"

////////////////////////////////////////////////////////////////////////////////
// computes relativistic Breit-Wigner amplitude with mass-dependent
// width for 2-body decay using various forms
// !NOTE! L is units of hbar/2
inline
std::complex<double>
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
		                   * (rpwa::barrierFactorSquared(L, q) / rpwa::barrierFactorSquared(L, q0));
	const double Gamma = Gamma0 * dyn;
	switch (algo) {

	case 0: {
		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double A = M0 * Gamma0;
		const double B = M0 * M0 - M * M;
		const double C = M0 * Gamma;
		return (A / (B * B + C * C)) * std::complex<double>(B, C);
		// return (M0 * Gamma0) / (M0 * M0 - M * M - imag * M0 * Gamma);
	}

	// >>> lies on unitarity circle
	case 1: {
		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double A = M0 * Gamma;
		const double B = M0 * M0 - M * M;
		const double C = M0 * Gamma;
		return (A / (B * B + C * C)) * std::complex<double>(B, C);
		// return (M0 * Gamma) / (M0 * M0 - M * M - imag * M0 * Gamma);
	}

	case 2: {
		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double A = M * Gamma;
		const double B = M0 * M0 - M * M;
		const double C = M0 * Gamma;
		return (A / (B * B + C * C)) * std::complex<double>(B, C);
		// return (M * Gamma) / (M0 * M0 - M * M - imag * M0 * Gamma);
	}

	case 3: {
		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double A = sqrt(M * Gamma) * sqrt(M0 * Gamma0);
		const double B = M0 * M0 - M * M;
		const double C = M0 * Gamma;
		return (A / (B * B + C * C)) * std::complex<double>(B, C);
		// return (sqrt(M * Gamma)) / (M0 * M0 - M * M - imag * M0 * Gamma);
	}

	case 4: {
		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double A = M0 * Gamma0;
		const double B = M0 * M0 - M * M;
		const double C = M0 * Gamma;
		return (A / (B * B + C * C)) * (rpwa::barrierFactor(L, q) / rpwa::barrierFactor(L, q0))
			* std::complex<double>(B, C);
		// return (M0 * Gamma0 * F_L(q) / F_L(q0)) / (M0 * M0 - M * M - imag * M0 * Gamma);
	}

	case 5: {
		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double A = M0 * Gamma0;
		const double B = M0 * M0 - M * M;
		const double C = M0 * Gamma;
		return (A / (B * B + C * C)) * (rpwa::barrierFactorSquared(L, q) / rpwa::barrierFactorSquared(L, q0))
			* std::complex<double>(B, C);
		// return (M0 * Gamma0 * F_L(q) / F_L(q0)) / (M0 * M0 - M * M - imag * M0 * Gamma);
	}

	// >>> lies on unitarity circle
	case 6: {
		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double A = M0 * Gamma0;
		const double B = M0 * M0 - M * M;
		const double C = M0 * Gamma0;
		return (A / (B * B + C * C)) * std::complex<double>(B, C);
		// return (M0 * Gamma0) / (M0 * M0 - M * M - imag * M0 * Gamma0);
	}

	case 7: {
		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double A = M0 * Gamma0;
		const double B = M0 * M0 - M * M;
		const double C = M0 * Gamma0;
		return (A / (B * B + C * C)) * dyn * std::complex<double>(B, C);
		// return (M0 * Gamma0 * dyn) / (M0 * M0 - M * M - imag * M0 * Gamma0);
	}
	}

	return 0;
}

#include<TFile.h>
#include<iostream>
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

	virtual std::complex<double> amp(const double m) const
	{
		// const double q  = breakupMomentum(m,     _daughterMass0, _daughterMass1);
		const double q  = sqrt(fabs(rpwa::breakupMomentumSquared(m, _daughterMass0, _daughterMass1, true)));
		const double q0 = rpwa::breakupMomentum(_mass, _daughterMass0, _daughterMass1);
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

	double operator() (double* x, double*) const { return std::norm(amp(x[0])); }

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
		const double psint = psIntegral(m);
		const double psint0 = psIntegral(_mass);
		const double dyn = psint / psint0;
/*		std::cout<<"psint="<<psint<<" psint0="<<psint0<<" dyn="<<dyn<<std::endl;
		std::cout<<"integrating from "<<_isobarAmp._daughterMass0 + _isobarAmp._daughterMass1<<" to "<<m - _daughterMass0<<std::endl;
*/
		return (_mass / m) * dyn;
	}

	virtual std::complex<double> amp(const double m) const
	{
		const double Gamma = _width * dyn(m);
		// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
		const double A = _mass * _width/*<-Gamma0*/;
		const double B = _mass * _mass - m * m;
		const double C = _mass * Gamma;
		return (A / (B * B + C * C)) * std::complex<double>(B, C);
		// return (M0 * Gamma0) / (M0 * M0 - M * M - imag * M0 * Gamma);
	}

	virtual double phaseSpace(double* x,
	                  double* p) const
	{
		const double twoBodyMass   = x[0];
		const double threeBodyMass = p[0];
		// breakup momentum between isobar and bachelor particle
		const double q   = rpwa::breakupMomentum(threeBodyMass, _daughterMass0, twoBodyMass);
		// breakup momentum between daughter particles of isobar
		const double q12 = rpwa::breakupMomentum(twoBodyMass, _isobarAmp._daughterMass0, _isobarAmp._daughterMass1);
		const std::complex<double> bwAmp = _isobarAmp.amp(twoBodyMass);
		return q * rpwa::barrierFactorSquared(_L, q) * q12 * std::norm(bwAmp);
	}

	virtual double psIntegral(const double threeBodyMass) const
	{
		TF1 f("f", this, &threeBodyDynAmpInt::phaseSpace, 0, 1, 1, "threeBodyDynAmpInt", "phaseSpace");
		f.SetParameter(0, threeBodyMass);
		const double result = f.Integral(_isobarAmp._daughterMass0 + _isobarAmp._daughterMass1, threeBodyMass - _daughterMass0);
/*if(threeBodyMass != _mass) {
std::cout<<"integrating from "<<_isobarAmp._daughterMass0 + _isobarAmp._daughterMass1<<" to "<<threeBodyMass - _daughterMass0<<std::endl;
std::cout<<"result = "<<result<<std::endl;
}*/
		return result;
	}

	const dynamicAmplitude& _isobarAmp;
	unsigned int            _algo;

};  // threeBodyDynAmpInt


int main() {

	TApplication app("", 0, 0);

	const unsigned int N_HISTS = 1;

	TFile* outFile = TFile::Open("amplitudeTests.root", "RECREATE");

	TCanvas* c1 = new TCanvas("c1");
//	c1->Divide(3, 3);

	TH2D* h1[N_HISTS];
	for(unsigned int i = 0; i < N_HISTS; ++i) {
		std::stringstream sstr;
		sstr<<"a1(1260)_1_"<<i;
		h1[i] = new TH2D(sstr.str().c_str(), sstr.str().c_str(), 300, 0, 2.5, 300, 0, 2.5);
		h1[i]->GetXaxis()->SetTitle("isobar mass (GeV)");
		h1[i]->GetYaxis()->SetTitle("m1+m2 (GeV)");
	}

	TCanvas* c2 = new TCanvas("c2");
//	c2->Divide(3, 3);

	TH2D* h2[N_HISTS];
	for(unsigned int i = 0; i < N_HISTS; ++i) {
		std::stringstream sstr;
		sstr<<"a1(1260)_2_"<<i;
		h2[i] = new TH2D(sstr.str().c_str(), sstr.str().c_str(), 300, 0, 2.5, 300, 0, 2.5);
		h2[i]->GetXaxis()->SetTitle("isobar mass (GeV)");
		h2[i]->GetYaxis()->SetTitle("m1+m2 (GeV)");
	}

	for(int i = 0; i < h1[0]->GetNbinsX(); ++i) {
		for(int j = 0; j < h1[0]->GetNbinsY(); ++j) {

			static bool first = true;
			static bool doneSomething = false;

			for(unsigned int canI = 0; canI < N_HISTS; ++canI) {

				const double x = h1[canI]->GetXaxis()->GetBinCenter(i);
				const double y = h1[canI]->GetYaxis()->GetBinCenter(j);

				if(y < 0.556) continue;
				if(x < 0.556) continue;
				if(y >= x) continue;

//				std::cout<<"x="<<x<<" y="<<y<<std::endl;

//				const double massPart = y / (2 * N_HISTS);

				// get Breit-Wigner parameters
				const double       M      = x;         // parent mass
				const double       m1     = 0.13957018;  // daughter 1 mass
				const double       m2     = y - m1;  // daughter 2 mass
				const double       q      = rpwa::breakupMomentum(M,  m1, m2);
				const double       M0     = 1.23;              // resonance peak position
				const double       q02    = rpwa::breakupMomentumSquared(M0, m1, m2, true);
				// !NOTE! the following is incorrect but this is how it was done in PWA2000
				const double       q0     = sqrt(fabs(q02));
				const double       Gamma0 = 0.425;             // resonance peak width
				const unsigned int L      = 0;

				const std::complex<double> bw = rpwa::breitWigner<std::complex<double> >(M, M0, Gamma0, L, q, q0);

//				std::cout<<"x="<<x<<" y="<<y<<" z="<<bw<<std::endl;

				h1[canI]->Fill(x, y, std::norm(bw));
				if(first) {
					std::stringstream sstr;
					sstr.precision(3);
					sstr<<"a1(1260)_1_"<<m1/y<<"*M_m2="<<m2/y<<"*M";
//					std::cout<<"setting histogram title to "<<sstr.str()<<std::endl;
					h1[canI]->SetTitle(sstr.str().c_str());
					h1[canI]->SetName(sstr.str().c_str());
					doneSomething = true;
				}

			}

			if(first and doneSomething) {
				first = false;
			}

		}
	}

	const int nBins = h1[0]->GetNbinsX() * h1[0]->GetNbinsY();
	boost::progress_display* progressIndicator = new boost::progress_display(nBins, std::cout, "");

	const double isobarMass          = 0.7665;
	const double isobarWidth         = 0.1502;
	const double isobarL             = 2;
	const double isobarDaughterMass0 = 0.13957018;
	const double isobarDaughterMass1 = isobarDaughterMass0;

	const dynamicAmplitude isobarDynAmp(isobarMass, isobarWidth, isobarL, isobarDaughterMass0,
			                                    isobarDaughterMass1, 0);

	for(int i = 0; i < h1[0]->GetNbinsX(); ++i) {
		for(int j = 0; j < h1[0]->GetNbinsY(); ++j) {

			static bool first = true;
			static bool doneSomething = false;

			++(*progressIndicator);

			for(unsigned int canI = 0; canI < N_HISTS; ++canI) {

				const double x = h1[canI]->GetXaxis()->GetBinCenter(i);
				const double y = h1[canI]->GetYaxis()->GetBinCenter(j);

				if(y < 0.556) continue;
				if(x < 0.556) continue;
				if(y >= x) continue;

//				std::cout<<"x="<<x<<" y="<<y<<std::endl;

//				const double massPart = y / (2 * N_HISTS);

				// get Breit-Wigner parameters
				const double       M      = x;         // parent mass
				const double       m1     = 0.13957018;  // daughter 1 mass
				const double       m2     = y - m1;  // daughter 2 mass
				const double       M0     = 1.35;              // resonance peak position
				const double       Gamma0 = 0.35;             // resonance peak width
				const unsigned int L      = 0;

				threeBodyDynAmpInt dynAmp(isobarDynAmp, M0, Gamma0, L, m1, m2, 2);
				const std::complex<double> bw = dynAmp.amp(M);

//				std::cout<<"x="<<x<<" y="<<y<<" z="<<bw<<std::endl;

				h2[canI]->Fill(x, y, std::norm(bw));
				if(first) {
					std::stringstream sstr;
					sstr.precision(3);
					sstr<<"a1(1260)_2_m1="<<m1/y<<"*M_m2="<<m2/y<<"*M";
//					std::cout<<"setting histogram title to "<<sstr.str()<<std::endl;
					h2[canI]->SetTitle(sstr.str().c_str());
					h2[canI]->SetName(sstr.str().c_str());
					doneSomething = true;
				}

			}

			if(first and doneSomething) {
				first = false;
			}

		}
	}

	outFile->cd();

	for(unsigned int i = 0; i < N_HISTS; ++i) {
		c1->cd();
//		c1->cd(i+1);
		h1[i]->SetMaximum(3);
		h1[i]->Draw("colz");
		h1[i]->SetStats(false);
		c2->cd();
//		c2->cd(i+1);
		h2[i]->SetMaximum(3);
		h2[i]->Draw("colz");
		h2[i]->SetStats(false);
	}
	c1->Draw();
	c2->Draw();
	c1->Write();
	c2->Write();
	outFile->Write();
	outFile->Close();
	gApplication->SetReturnFromRun(kFALSE);
	gSystem->Run();
	return 0;

}
