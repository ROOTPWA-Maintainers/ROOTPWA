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
//      tests nBodyPhaseSpaceGen
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <iostream>
#include <string>
#include <algorithm>

#include "TStopwatch.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TGenPhaseSpace.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include "TColor.h"

#include "reportingUtils.hpp"
#include "physUtils.hpp"
#include "conversionUtils.hpp"
#include "nBodyPhaseSpaceGen.h"


using namespace std;
using namespace rpwa;


const double gChargedPionMass = 0.13957018;  // charged pion rest mass [GeV/c^2] [PDG08]


// constucts Lorentz vector of mother particle
TLorentzVector constructMother(TRandom3&    random,
                               const double mass,            // [GeV/c^2]
                               const double maxPt  = 0.150,  // [GeV/c]
                               const double maxEta = 1)
{
	double pt      = maxPt * random.Rndm();                 // n-body transverse momentum in lab frame [0, 150] MeV/c
	double eta     = 2 * maxEta * random.Rndm() - maxEta;   // n-body pseudorapidity in lab frame [-1, 1]
	double phi     = 4 * asin(1.) * random.Rndm();          // n-body azimuth in lab frame [0, 2 pi]
	double tanhEta = tanh(eta);
	TLorentzVector mother;
	mother.SetPx(pt * cos(phi));
	mother.SetPy(pt * sin(phi));
	mother.SetPz(tanhEta * sqrt((mass * mass + pt * pt) / (1 - tanhEta * tanhEta)));
	mother.SetVectM(mother.Vect(), mass);
	return mother;
}


// signature for TF1
double
dalitzKinematicBorder(double* var,
                      double* par)
{
	// variables
	const double mass_2 = var[0];  // 2-body mass squared on x-axis

	// parameters
	const double  M   = par[0];                       // 3-body mass
	const double* m   = &par[1];                      // array with the 3 daughter masses
	const bool    min = (par[4] < 0) ? true : false;  // switches between curves for minimum and maximum mass squared on y-axis

	return dalitzKinematicBorder(mass_2, M, m, min);
}


void setResidualPalette(const unsigned int nmbSteps = 99)
{
	const unsigned int nmbCol           = 3;
	Double_t           rChannel[nmbCol] = {0, 1,   1};
	Double_t           gChannel[nmbCol] = {0, 1,   0};
	Double_t           bChannel[nmbCol] = {1, 1,   0};
	Double_t           zPos    [nmbCol] = {0, 0.5, 1};
	TColor::CreateGradientColorTable(nmbCol, zPos, rChannel, gChannel, bChannel, nmbSteps);
}


// fit function for flat Dalitz plot
double
dalitzFitFunc(double* var,
              double* par)
{
	const double m01_2 = var[0];
	const double m12_2 = var[1];

	// parameters
	const double  M   = par[1];   // 3-body mass
	const double* m   = &par[2];  // array with the 3 daughter masses

	if (   (m12_2 > dalitzKinematicBorder(m01_2, M, m, true ))
	    && (m12_2 < dalitzKinematicBorder(m01_2, M, m, false)))
		return par[0];
	else {
		TF2::RejectPoint();
		return 0;
	}
}


// generalization of TF2::Integral
double
nDimIntegral(TF1*          f,
             int           n,                         // number of variables
             const double* min,                       // array of lower integration borders
             const double* max,                       // array of upper integration borders
             const double* funcPar,                   // function parameters
             int           maxNmbOfPoints = 10000,    // maximum allowed number of function evaluations
             double        epsilon        = 0.00001)  // relative accuracy
{
	if (n == 1)
		return f->Integral(min[0], max[0], funcPar);
	f->SetParameters(funcPar);
	double resRelErr;      // estimate of the relative accuracy of the result
	int    nmbOfFuncEval;  // number of function evaluations performed
	int    errCode;        // 0 = normal exit,  1 = maxNmbOfPoints too small, 3 =  n < 2 or n > 15
	double result = f->IntegralMultiple(n, min, max, -1, maxNmbOfPoints, epsilon, resRelErr, nmbOfFuncEval, errCode);
	return result;
}


double
nBodyPhaseSpaceElement(double* var,  // effective masses of the intermediate (i + 2)-body systems
                       double* par)  // # of daughters, n daughter masses, M_n
{
	const int     n = static_cast<int>(par[0]);  // number of decay daughters > 2
	const double* m = &par[1];                   // masses of the daughter particles
	double        M[n];                          // (i + 1)-body invariant mass
	M[0] = m[0];                                 // 1-body mass is mass of first daughter
	for (int i = 1; i < (n - 1); ++i)            // intermediate (i + 1)-body effective masses
		M[i] = var[i - 1];
	M[n - 1] = par[n + 1];                       // n-body mass

	// reject points outside of kinematical range
	for (int i = 1; i < n; ++i)
		if (M[i] < M[i - 1] + m[i]) {
			TF1::RejectPoint();
			return 0;
		}

	double momProd = 1;          // product of breakup momenta
	for (int i = 1; i < n; ++i)  // loop over 2- to n-bodies
		momProd *= breakupMomentum(M[i], M[i - 1], m[i]);
	return momProd;  // the term that needs to be intergrated; normalization in integration routine
}


TF1* dLips;  // n-body Lorentz-invariant phase space element


double
nBodyPhaseSpace(double* var,  // n-body mass
                double* par)  // # of daughters, n daughter masses
{
	const double  M_n = var[0];  // n-body mass

	const int     n   = static_cast<int>(par[0]);  // number of decay daughters > 2
	const double* m   = &par[1];                   // masses of the daughter particles
	const double  A   = par[n + 1];                // scale factor

	// reject points outside of kinematical range
	double mSum = 0;
	for (int i = 0; i < n; ++i)
		mSum += m[i];
	if (M_n < mSum)
		return 0;

	// normalization
	const double f_n  = 1 / (2 * pow(twoPi, 2 * n - 3));  // S. U. Chung's normalization
	//const double f_n  = pi * pow(twoPi, n - 2);           // genbod's normalization
	//const double f_n  = 1;
	const double norm = f_n / M_n;

	if (n == 2)
		return A * norm * breakupMomentum(M_n, m[0], m[1]);

	// prepare parameter array for phase space element
	double lipsPar[n + 2];
	lipsPar[0] = n;
	for (int i = 1; i <= n; ++i)
		lipsPar[i] = m[i - 1];
	lipsPar[n + 1] = M_n;

	// set integration limits
	double min[n - 2];  // lower integration limits for effective (i + 2)-body masses
	min[0] = m[0] + m[1];
	for (int i = 1; i < n - 2; ++i)
		min[i] = min[i - 1] + m[i + 1];
	double max[n - 2];  // upper integration limits for effective (i + 2)-body masses
	mSum = 0;
	for (int i = n - 3; i >= 0; --i) {
		mSum  += m[i + 2];
		max[i] = M_n - mSum;
	}

	// perform integration
	return A * norm * nDimIntegral(dLips, n - 2, min, max, lipsPar, 500000, 0.00001);
}


void testNBodyPhaseSpaceGen(const unsigned int nmbEvents   = 1000000,
                            const string&      outFileName = "./testNBodyPhaseSpaceGen.root")
{
	TStopwatch timer;
	timer.Start();

	// set parameters
	const unsigned int n            = 3;                                                       // number of daughters
	double             m[n]         = {gChargedPionMass, gChargedPionMass, gChargedPionMass};  // daughter masses
	const double       massRange[2] = {m[0] + m[1] + m[2], 2};                                 // [GeV/c^2]
	const double       mass         = (massRange[0] + massRange[1]) / 2;
	const unsigned int seed         = 1234567890;
	const bool         hitMissMode  = false;

	// create output file
	printInfo << "creating output file '" << outFileName << "' ... " << endl;
	TFile* outFile = TFile::Open(outFileName.c_str(), "RECREATE");
	if (!outFile) {
		printErr << "could not create output file. aborting..." << endl;
		return;
	}
	cout << "    ... success" << endl;

	// setup nBodyPhaseSpaceGen generators
	const unsigned int nmbGen = 5;
	nBodyPhaseSpaceGen* gen[5];
	for (unsigned int i = 0; i < nmbGen; ++i)
		gen[i] = new nBodyPhaseSpaceGen();
	gen[0]->setWeightType(nBodyPhaseSpaceGen::S_U_CHUNG);
	gen[1]->setWeightType(nBodyPhaseSpaceGen::NUPHAZ);
	gen[2]->setWeightType(nBodyPhaseSpaceGen::GENBOD);
	gen[3]->setWeightType(nBodyPhaseSpaceGen::FLAT);
	gen[4]->setWeightType(nBodyPhaseSpaceGen::S_U_CHUNG);
	for (unsigned int i = 0; i < 4; ++i)
		gen[i]->setKinematicsType(nBodyPhaseSpaceGen::RAUBOLD_LYNCH);
	gen[4]->setKinematicsType(nBodyPhaseSpaceGen::BLOCK);
	for (unsigned int i = 0; i < nmbGen; ++i) {
		gen[i]->setDecay(n, m);
		gen[i]->setSeed(seed);
		gen[i]->setMaxWeight(1.01 * gen[i]->estimateMaxWeight(mass));
	}
	// setup TGenPhaseSpace generator
	TGenPhaseSpace TGen;
	gRandom->SetSeed(seed);  // TGenPhaseSpace uses gRandom
	double maxTGenWeight = 0;
	for (unsigned int i = 0; i < 10000; ++i) {
		TLorentzVector mother(0, 0, 0, mass);
		TGen.SetDecay(mother, (int)n, m);
		const double weight = TGen.Generate();
		maxTGenWeight = (weight > maxTGenWeight) ? weight : maxTGenWeight;
	}
	maxTGenWeight *= 1.01;
	double   maxTGenWeightObserved = 0;
	TRandom3 random(seed);

	// book histograms
	TH2F* hDalitz [nmbGen + 1];
	TH1F* hMass   [nmbGen + 1];
	for (unsigned int i = 0; i < nmbGen; ++i) {
		hDalitz[i] = new TH2F(("hDalitz" + toString(i)).c_str(),
		                      ("Dalitz Plot Residuum nBodyPhaseSpaceGen " + toString(i) + ";m_{01}^{2} [(GeV/c^{2})^{2}];m_{12}^{2} [(GeV/c^{2})^{2}]").c_str(),
		                      100, 0, 1.5, 100, 0, 1.5);
		hMass[i]   = new TH1F(("hMass" + toString(i)).c_str(),
		                      ("Mother Mass nBodyPhaseSpaceGen " + toString(i) + ";m [GeV/c^{2}]").c_str(),
		                      1000, 0, 2.5);
	}
	hDalitz [nmbGen] = new TH2F("hDalitzTGen", "Dalitz Plot Residuum TGen;m_{01}^{2} [(GeV/c^{2})^{2}];m_{12}^{2} [(GeV/c^{2})^{2}]", 100, 0, 1.5, 100, 0, 1.5);
	hMass   [nmbGen] = new TH1F("hMassTGen",   "Mother Mass TGen;m [GeV/c^{2}]", 1000, 0, 2.5);
	TH1F* hMassRaw   = new TH1F("hMassRaw",    "Mother Mass Raw;m [GeV/c^{2}]",  1000, 0, 2.5);

	if (0) {
		printInfo << "testing generators for constant 3-body mass = " << mass << " GeV/c^2" << endl
		          << "    generating " << nmbEvents << " events ... " << endl;
		unsigned int countAtt[nmbGen + 1];
		unsigned int countAcc[nmbGen + 1];
		bool         done = false;
		for (unsigned int i = 0; i <= nmbGen; ++i)
			countAtt[i] = countAcc[i] = 0;
		while (!done) {
			TLorentzVector mother = constructMother(random, mass);
			for (unsigned int i = 0; i < nmbGen; ++i)
				if (countAcc[i] < nmbEvents){
					++countAtt[i];
					if (gen[i]->generateDecayAccepted(mother)) {
						++countAcc[i];
						const TLorentzVector m01 = gen[i]->daughter(0) + gen[i]->daughter(1);
						const TLorentzVector m12 = gen[i]->daughter(1) + gen[i]->daughter(2);
						hDalitz[i]->Fill(m01.M2(), m12.M2());
					}
				}
			TGen.SetDecay(mother, (int)n, m);
			const double TGenWeight = TGen.Generate();
			maxTGenWeightObserved   = (TGenWeight > maxTGenWeightObserved) ? TGenWeight : maxTGenWeightObserved;
			if (countAcc[nmbGen] < nmbEvents) {
				++countAtt[nmbGen];
				if ((TGenWeight / maxTGenWeight) > random.Rndm()) {
					++countAcc[nmbGen];
					const TLorentzVector m01 = *TGen.GetDecay(0) + *TGen.GetDecay(1);
					const TLorentzVector m12 = *TGen.GetDecay(1) + *TGen.GetDecay(2);
					hDalitz[nmbGen]->Fill(m01.M2(), m12.M2());
				}
			}
			done = true;
			for (unsigned int i = 0; i <= nmbGen; ++i)
				if (countAcc[i] < nmbEvents)
					done = false;
		}
		cout << "    ... done." << endl << endl;
		for (unsigned int i = 0; i < nmbGen; ++i)
			printInfo << "gen[" << i << "]: " << countAtt[i] << " attempts, " << ((double)nmbEvents / countAtt[i]) * 100 << " % efficiency" << endl
			          << *gen[i] << endl;
		printInfo << "TGenPhaseSpace: " << countAtt[nmbGen] << " attempts, " << ((double)nmbEvents / countAtt[nmbGen]) * 100 << " % efficiency" << endl
		          << "TGenPhaseSpace maximum weight = " << maxTGenWeight
		          << " vs. maximum weight observed = " << maxTGenWeightObserved << endl << endl;

		// plot residuum of Dalitz plots after subtracting expected flat component
		setResidualPalette(99);
		for (unsigned int i = 0; i <= nmbGen; ++i) {
			new TCanvas();
			const double xMin = (m[0] + m[1]) * (m[0] + m[1]);
			const double xMax = (mass - m[2]) * (mass - m[2]);
			const double yMin = (m[1] + m[2]) * (m[1] + m[2]);
			const double yMax = (mass - m[0]) * (mass - m[0]);
			TF2* dalitzFit = new TF2("dalitzFit", dalitzFitFunc, xMin, xMax, yMin, yMax, 5);
			dalitzFit->SetParameter(0, 1);
			dalitzFit->FixParameter(1, mass);
			dalitzFit->FixParameter(2, m[0]);
			dalitzFit->FixParameter(3, m[1]);
			dalitzFit->FixParameter(4, m[2]);
			const double dalitzInt = dalitzFit->Integral(xMin, xMax, yMin, yMax, 1e-6)
				/ (hDalitz[i]->GetXaxis()->GetBinWidth(1) * hDalitz[i]->GetYaxis()->GetBinWidth(1));
			dalitzFit->SetParameter(0, nmbEvents / dalitzInt);
			hDalitz[i]->Sumw2();
			hDalitz[i]->Add(dalitzFit, -1, "I");  // ROOT does not support I flag for 2D functions; increased residual at borders
			// center z-range
			const double maxVal = max(fabs(hDalitz[i]->GetMinimum()), fabs(hDalitz[i]->GetMaximum()));
			hDalitz[i]->SetMinimum(-maxVal);
			hDalitz[i]->SetMaximum( maxVal);
			hDalitz[i]->SetContour(99);
			hDalitz[i]->Draw("COLZ");
			// draw border
			TF1* dalitzBorder = new TF1("dalitzBorder", dalitzKinematicBorder, xMin, xMax, 5);
			dalitzBorder->SetParameter(0, mass);
			dalitzBorder->SetParameter(1, m[0]);
			dalitzBorder->SetParameter(2, m[1]);
			dalitzBorder->SetParameter(3, m[2]);
			dalitzBorder->SetParameter(4, -1);
			dalitzBorder->SetNpx(1000);
			dalitzBorder->SetLineWidth(1);
			dalitzBorder->SetLineColor(3);
			dalitzBorder->DrawCopy("SAME");
			dalitzBorder->SetParameter(4, +1);
			dalitzBorder->DrawCopy("SAME");
			gPad->Update();
			gPad->RedrawAxis();
		}
	}

	if (1) {
		printInfo << "testing generators for 3-body in mass range [" << massRange[0] << ", " << massRange[1] << "] GeV/c^2" << endl
		          << "    generating " << nmbEvents << " events ... " << endl;
		// recalculate maximum weights
		gen[0]->setMaxWeight(1.01 * gen[0]->estimateMaxWeight(massRange[1]));
		gen[1]->setMaxWeight(1.01 * gen[1]->estimateMaxWeight(massRange[1]));
		gen[2]->setMaxWeight(1.01 * gen[2]->estimateMaxWeight(massRange[0] + 0.0001));
		gen[3]->setMaxWeight(1.01 * gen[3]->estimateMaxWeight(massRange[1]));
		gen[4]->setMaxWeight(1.01 * gen[4]->estimateMaxWeight(massRange[1]));
		maxTGenWeight = 0;
		for (unsigned int i = 0; i < 10000; ++i) {
			TLorentzVector mother(0, 0, 0, massRange[0] + 0.0001);
			TGen.SetDecay(mother, (int)n, m);
			const double weight = TGen.Generate();
			maxTGenWeight = (weight > maxTGenWeight) ? weight : maxTGenWeight;
		}
		maxTGenWeight        *= 1.01;
		maxTGenWeightObserved = 0;

		unsigned int countAtt[nmbGen + 1];
		unsigned int countAcc[nmbGen + 1];
		bool         done = false;
		for (unsigned int i = 0; i <= nmbGen; ++i)
			countAtt[i] = countAcc[i] = 0;
		while (!done) {
			TLorentzVector mother = constructMother(random, massRange[0] + random.Rndm() * (massRange[1] - massRange[0]));
			hMassRaw->Fill(mother.M());
			for (unsigned int i = 0; i < nmbGen; ++i)
				if (countAcc[i] < nmbEvents) {
					++countAtt[i];
					if (hitMissMode) {
						if (gen[i]->generateDecayAccepted(mother)) {
							++countAcc[i];
							hMass[i]->Fill(mother.M());
						}
					} else {
						hMass[i]->Fill(mother.M(), gen[i]->generateDecay(mother));
						if (gen[i]->eventAccepted())
							++countAcc[i];
					}
				}
			TGen.SetDecay(mother, (int)n, m);
			const double TGenWeight = TGen.Generate();
			maxTGenWeightObserved   = (TGenWeight > maxTGenWeightObserved) ? TGenWeight : maxTGenWeightObserved;
			if (countAcc[nmbGen] < nmbEvents) {
				++countAtt[nmbGen];
				if (hitMissMode) {
					if ((TGenWeight / maxTGenWeight) > random.Rndm())
						hMass[nmbGen]->Fill(mother.M());
				} else
					hMass[nmbGen]->Fill(mother.M(), TGenWeight);
				if ((TGenWeight / maxTGenWeight) > random.Rndm())
					++countAcc[nmbGen];
			}
			done = true;
			for (unsigned int i = 0; i <= nmbGen; ++i)
				if (countAcc[i] < nmbEvents)
					done = false;
		}
		cout << "    ... done." << endl << endl;
		for (unsigned int i = 0; i < nmbGen; ++i)
			printInfo << "gen[" << i << "]: " << countAtt[i] << " attempts, " << ((double)nmbEvents / countAtt[i]) * 100 << " % efficiency" << endl
			          << *gen[i] << endl;
		printInfo << "TGenPhaseSpace: " << countAtt[nmbGen] << " attempts, " << ((double)nmbEvents / countAtt[nmbGen]) * 100 << " % efficiency" << endl
		          << "maximum weight = " << maxTGenWeight
		          << " vs. maximum weight observed = " << maxTGenWeightObserved << endl << endl;

		// draw histograms
		dLips     = new TF1("nBodyPhaseSpaceElement", nBodyPhaseSpaceElement, 0,            100,          n + 2);
		TF1* Lips = new TF1("nBodyPhaseSpace",        nBodyPhaseSpace,        massRange[0], massRange[1], n + 2);
		Lips->SetParameter(0, n);
		for (unsigned int i = 1; i <= n; ++i)
			Lips->SetParameter(i, m[i - 1]);
		Lips->SetLineColor(2);
		Lips->SetLineWidth(2);
		for (unsigned int i = 0; i <= nmbGen; ++i) {
			Lips->SetParameter(n + 1, 1);
			double A;
			if (hitMissMode)
				// scale function to last bin
				A = hMass[i]->GetBinContent(hMass[i]->FindBin(massRange[1]) - 1) / Lips->Eval(massRange[1]);
			else {
				// scale histograms by number of events per bin
				for (int j = 1; j < hMassRaw->GetNbinsX(); ++j) {
					const double nmbEvBin = hMassRaw->GetBinContent(j);
					if (nmbEvBin > 0)
						hMass[i]->SetBinContent(j, hMass[i]->GetBinContent(j) / nmbEvBin);
				}
				A = 1;
			}
			Lips->SetParameter(n + 1, A);
			new TCanvas();
			hMass[i]->Draw();
			Lips->DrawCopy("LSAME");
			if (hitMissMode)
				hMass[i]->Draw("SAME E");
		}
		new TCanvas;
		hMassRaw->Draw();
	}

	if (0) {
		const unsigned int n2            = 4;                                     // number of daughters
		double             m2[n2]        = {gChargedPionMass, gChargedPionMass,
		                                    gChargedPionMass, gChargedPionMass};  // daughter masses
		const double       massRange2[2] = {m2[0] + m2[1] + m2[2] + m2[3], 2};    // [GeV/c^2]
		printInfo << "comparing run time of generators for 4-body in mass range [" << massRange2[0] << ", " << massRange2[1] << "] GeV/c^2" << endl
		          << "    generating " << nmbEvents << " events using TGenPhaseSpace ... " << endl;

		// setup nBodyPhaseSpaceGen generator
		nBodyPhaseSpaceGen gen2;
		gen2.setWeightType(nBodyPhaseSpaceGen::GENBOD);
		//gen2.setWeightType(nBodyPhaseSpaceGen::S_U_CHUNG);
		gen2.setKinematicsType(nBodyPhaseSpaceGen::BLOCK);
		gen2.setDecay(n2, m2);
		gen2.setSeed(seed);
		gen2.setMaxWeight(1.01 * gen2.estimateMaxWeight(massRange2[0] + 0.0001));
		//gen2.setMaxWeight(1.01 * gen2.estimateMaxWeight(massRange2[1]));
		// recalculate maximum weight for TGenPhaseSpace generator
		maxTGenWeight = 0;
		for (unsigned int i = 0; i < 10000; ++i) {
			TLorentzVector mother(0, 0, 0, massRange2[0] + 0.0001);
			TGen.SetDecay(mother, (int)n2, m2);
			const double weight = TGen.Generate();
			maxTGenWeight = (weight > maxTGenWeight) ? weight : maxTGenWeight;
		}
		maxTGenWeight        *= 1.01;
		maxTGenWeightObserved = 0;

		TStopwatch timer2;
		timer2.Start();
		unsigned int countAtt = 0;
		unsigned int countAcc = 0;
		random.SetSeed(seed);
		while (countAcc < nmbEvents) {
			TLorentzVector mother = constructMother(random, massRange2[0] + random.Rndm() * (massRange2[1] - massRange2[0]));
			TGen.SetDecay(mother, (int)n2, m2);
			const double TGenWeight = TGen.Generate();
			maxTGenWeightObserved   = (TGenWeight > maxTGenWeightObserved) ? TGenWeight : maxTGenWeightObserved;
			++countAtt;
			if ((TGenWeight / maxTGenWeight) > random.Rndm())
				++countAcc;
		}
		timer2.Stop();
		const double timeTGen = timer2.RealTime();
		cout << "    ... done." << endl
		     << "    " << countAtt << " attempts, " << ((double)nmbEvents / countAtt) * 100 << " % efficiency" << endl
		     << "    consumed time: ";
		timer2.Print();
		timer2.Reset();
		cout << "    TGenPhaseSpace maximum weight = " << maxTGenWeight
		     << " vs. maximum weight observed = " << maxTGenWeightObserved << endl;

		cout << endl << "    generating " << nmbEvents << " events using nBodyPhaseSpaceGen ... " << endl;
		timer2.Start();
		countAtt = 0;
		countAcc = 0;
		random.SetSeed(seed);
		while (countAcc < nmbEvents) {
			TLorentzVector mother = constructMother(random, massRange2[0] + random.Rndm() * (massRange2[1] - massRange2[0]));
			++countAtt;
			if (gen2.generateDecayAccepted(mother))
				++countAcc;
		}
		timer2.Stop();
		const double timeNBody = timer2.RealTime();
		cout << "    ... done." << endl
		     << "    " <<  countAtt << " attempts, " << ((double)nmbEvents / countAtt) * 100 << " % efficiency" << endl
		     << "    consumed time: ";
		timer2.Print();
		timer2.Reset();
		printInfo << gen2 << endl;
		printInfo << "nBodyPhaseSpaceGen speed gain w.r.t. TGenPhaseSpace = " << 100 * timeTGen / timeNBody << " %" << endl << endl;

		printInfo << "checking mass dependence of nBodyPhaseSpaceGen 4-body phase space in mass range [" << massRange2[0] << ", " << massRange2[1] << "] GeV/c^2" << endl;
		gen2.setWeightType(nBodyPhaseSpaceGen::S_U_CHUNG);
		gen2.setKinematicsType(nBodyPhaseSpaceGen::BLOCK);
		gen2.setDecay(n2, m2);
		gen2.setSeed(seed);
		gen2.setMaxWeight(1.01 * gen2.estimateMaxWeight(massRange2[1]));
		TH1F* hMass2 = new TH1F("hMass4Body", "4-Body Mass nBodyPhaseSpaceGen;m [GeV/c^{2}]", 1000, 0, 2.5);
		countAtt = 0;
		countAcc = 0;
		random.SetSeed(seed);
		while (countAcc < nmbEvents) {
			TLorentzVector mother = constructMother(random, massRange2[0] + random.Rndm() * (massRange2[1] - massRange2[0]));
			++countAtt;
			if (gen2.generateDecayAccepted(mother)) {
				++countAcc;
				hMass2->Fill(mother.M());
			}
		}
		cout << "    ... done." << endl
		     << "    " <<  countAtt << " attempts, " << ((double)nmbEvents / countAtt) * 100 << " % efficiency" << endl
		     << gen2 << endl;
		dLips     = new TF1("nBodyPhaseSpaceElement", nBodyPhaseSpaceElement, 0,             100,           n2 + 2);
		TF1* Lips = new TF1("nBodyPhaseSpace",        nBodyPhaseSpace,        massRange2[0], massRange2[1], n2 + 2);
		Lips->SetParameter(0, n2);
		for (unsigned int i = 1; i <= n2; ++i)
			Lips->SetParameter(i, m2[i - 1]);
		Lips->SetLineColor(2);
		Lips->SetLineWidth(2);
		Lips->SetParameter(n2 + 1, 1);
		const double A = hMass2->GetBinContent(hMass2->FindBin(massRange2[1]) - 1) / Lips->Eval(massRange2[1]);
		Lips->SetParameter(n2 + 1, A);
		new TCanvas();
		hMass2->Draw();
		Lips->DrawCopy("LSAME");
		hMass2->Draw("SAME E");
	}

	timer.Stop();
	printInfo << endl << "this job consumed: ";
	timer.Print();
}
