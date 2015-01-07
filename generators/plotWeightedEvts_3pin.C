/*
 * IMPORTANT:
 *
 * The naming scheme of histogram files is critical in this program.
 * More complex diagrams, that are automatically generated, rely on the specific scheme defined below.
 * The program requests names for the Monte Carlo histograms to contain "MC_"
 * or end with "MC". (the reason for _ is that the string doesnt no accidentally contain another "MC")
 * The Data histograms must have the same filename as the corresponding MC histogram but "Data" replaced for "MC".
 *
 * Also for each massbin create a directory and create the same filenames for each variable
 */

#include <iostream>
#include <string>
#include <vector>

#include <boost/progress.hpp>

#include <TClonesArray.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TKey.h>
#include <TLorentzRotation.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TTree.h>

#ifndef __CINT__
#include <particleDataTable.h>
#endif // __CINT__
#include <reportingUtils.hpp>

using namespace std;

const int    HISTLIMITS_COSTHETA_BINS = 100;
const double HISTLIMITS_COSTHETA_MAX  =   1.;
const double HISTLIMITS_COSTHETA_MIN  =  -1.;
const int    HISTLIMITS_3PI_MASS_BINS = 52;
const double HISTLIMITS_3PI_MASS_MAX  =  2.54;
const double HISTLIMITS_3PI_MASS_MIN  =  0.46;
const int    HISTLIMITS_2PI_MASS_BINS = 250;
const double HISTLIMITS_2PI_MASS_MAX  =   2.5;
const double HISTLIMITS_2PI_MASS_MIN  =   0.;
const int    HISTLIMITS_2PI_MASS2_BINS = 100;
const double HISTLIMITS_2PI_MASS2_MAX  =   5.;
const double HISTLIMITS_2PI_MASS2_MIN  =   0.;
const int    HISTLIMITS_PHI_BINS = 100;
const double HISTLIMITS_PHI_MAX  =  TMath::Pi();
const double HISTLIMITS_PHI_MIN  = -TMath::Pi();
const int    HISTLIMITS_TPRIME_BINS = 400;
const double HISTLIMITS_TPRIME_MAX  =   4.;
const double HISTLIMITS_TPRIME_MIN  =   0.;

struct GJHistBunch {
	// base histograms
	std::vector<TH1D*> isobar_mass;
	std::vector<TH1D*> costheta_GJF;
	std::vector<TH1D*> phi_GJF;
	std::vector<TH2D*> costheta_GJF_tprime;
	std::vector<TH2D*> phi_costheta_GJF;

	// resolved for positive and negative reflectivity, and flat wave
	std::vector<THStack*> costheta_GJF_Stack;
	std::vector<TH1D*> costheta_GJF_PosRef;
	std::vector<TH1D*> costheta_GJF_NegRef;
	std::vector<TH1D*> costheta_GJF_Flat;
};

struct HelicityHistBunch {
	// base histograms
	std::vector<TH1D*> costheta_HF;
	std::vector<TH1D*> phi_HF;
	std::vector<TH2D*> phi_costheta_HF;
};

struct HelicityAngles {
	double cosTheta;
	double phi;
};

int getTotalCharge(std::pair<int, int> p) {
	return p.first + p.second;
}

GJHistBunch GJHistBunchFactory(const std::string& name_prefix, const bool twoMc) {
	GJHistBunch temp;
	TH1D* hMIsobarData = new TH1D(("hMIsobarData_" + name_prefix).c_str(), (name_prefix + " Isobar Mass (Data)").c_str(), HISTLIMITS_2PI_MASS_BINS, HISTLIMITS_2PI_MASS_MIN, HISTLIMITS_2PI_MASS_MAX);
	hMIsobarData->SetXTitle("isobar mass [GeV]");
	hMIsobarData->SetYTitle("# of events");
	temp.isobar_mass.push_back(hMIsobarData);
	if (twoMc) {
		TH1D* hMIsobarMcPsp = new TH1D(("hMIsobarMcPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Mass (McPsp)").c_str(), HISTLIMITS_2PI_MASS_BINS, HISTLIMITS_2PI_MASS_MIN, HISTLIMITS_2PI_MASS_MAX);
		hMIsobarMcPsp->SetXTitle("isobar mass [GeV]");
		hMIsobarMcPsp->SetYTitle("# of events");
		temp.isobar_mass.push_back(hMIsobarMcPsp);
		TH1D* hMIsobarMcAcc = new TH1D(("hMIsobarMcAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Mass (McAcc)").c_str(), HISTLIMITS_2PI_MASS_BINS, HISTLIMITS_2PI_MASS_MIN, HISTLIMITS_2PI_MASS_MAX);
		hMIsobarMcAcc->SetXTitle("isobar mass [GeV]");
		hMIsobarMcAcc->SetYTitle("# of events");
		temp.isobar_mass.push_back(hMIsobarMcAcc);
	} else {
		TH1D* hMIsobarMc = new TH1D(("hMIsobarMc_" + name_prefix).c_str(), (name_prefix + " Isobar Mass (Mc)").c_str(), HISTLIMITS_2PI_MASS_BINS, HISTLIMITS_2PI_MASS_MIN, HISTLIMITS_2PI_MASS_MAX);
		hMIsobarMc->SetXTitle("isobar mass [GeV]");
		hMIsobarMc->SetYTitle("# of events");
		temp.isobar_mass.push_back(hMIsobarMc);
	}

	TH1D* hGJData = new TH1D(("hGJData_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (Data)").c_str(), HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
	hGJData->SetXTitle("isobar cos(#theta_{GJ})");
	hGJData->SetYTitle("# of events");
	temp.costheta_GJF.push_back(hGJData);
	if (twoMc) {
		TH1D* hGJMcPsp = new TH1D(("hGJMcPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (McPsp)").c_str(), HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
		hGJMcPsp->SetXTitle("isobar cos(#theta_{GJ})");
		hGJMcPsp->SetYTitle("# of events");
		temp.costheta_GJF.push_back(hGJMcPsp);
		TH1D* hGJMcAcc = new TH1D(("hGJMcAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (McAcc)").c_str(), HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
		hGJMcAcc->SetXTitle("isobar cos(#theta_{GJ})");
		hGJMcAcc->SetYTitle("# of events");
		temp.costheta_GJF.push_back(hGJMcAcc);
	} else {
		TH1D* hGJMc = new TH1D(("hGJMc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (Mc)").c_str(), HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
		hGJMc->SetXTitle("isobar cos(#theta_{GJ})");
		hGJMc->SetYTitle("# of events");
		temp.costheta_GJF.push_back(hGJMc);
	}

	TH2D* hGJtData = new TH2D(("hGJtData_" + name_prefix).c_str(), (name_prefix + " Isobar Cos GJ Theta vs t' (Data)").c_str(), HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX, HISTLIMITS_TPRIME_BINS, HISTLIMITS_TPRIME_MIN, HISTLIMITS_TPRIME_MAX);
	hGJtData->SetXTitle("isobar cos(#theta_{GJ})");
	hGJtData->SetYTitle("t' [GeV]");
	hGJtData->SetZTitle("# events");
	hGJtData->SetOption("COLZ");
	temp.costheta_GJF_tprime.push_back(hGJtData);
	if (twoMc) {
		TH2D* hGJtMcPsp = new TH2D(("hGJtMcPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Cos GJ Theta vs t' (McPsp)").c_str(), HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX, HISTLIMITS_TPRIME_BINS, HISTLIMITS_TPRIME_MIN, HISTLIMITS_TPRIME_MAX);
		hGJtMcPsp->SetXTitle("isobar cos(#theta_{GJ})");
		hGJtMcPsp->SetYTitle("t' [GeV]");
		hGJtMcPsp->SetZTitle("# events");
		hGJtMcPsp->SetOption("COLZ");
		temp.costheta_GJF_tprime.push_back(hGJtMcPsp);
		TH2D* hGJtMcAcc = new TH2D(("hGJtMcAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos GJ Theta vs t' (McAcc)").c_str(), HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX, HISTLIMITS_TPRIME_BINS, HISTLIMITS_TPRIME_MIN, HISTLIMITS_TPRIME_MAX);
		hGJtMcAcc->SetXTitle("isobar cos(#theta_{GJ})");
		hGJtMcAcc->SetYTitle("t' [GeV]");
		hGJtMcAcc->SetZTitle("# events");
		hGJtMcAcc->SetOption("COLZ");
		temp.costheta_GJF_tprime.push_back(hGJtMcAcc);
	} else {
		TH2D* hGJtMc = new TH2D(("hGJtMc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos GJ Theta vs t' (Mc)").c_str(), HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX, HISTLIMITS_TPRIME_BINS, HISTLIMITS_TPRIME_MIN, HISTLIMITS_TPRIME_MAX);
		hGJtMc->SetXTitle("isobar cos(#theta_{GJ})");
		hGJtMc->SetYTitle("t' [GeV]");
		hGJtMc->SetZTitle("# events");
		hGJtMc->SetOption("COLZ");
		temp.costheta_GJF_tprime.push_back(hGJtMc);
	}

	TH1D* hTYData = new TH1D(("hTYData_" + name_prefix).c_str(), (name_prefix + " Isobar Treiman-Yang Phi (Data)").c_str(), HISTLIMITS_PHI_BINS, HISTLIMITS_PHI_MIN, HISTLIMITS_PHI_MAX);
	hTYData->SetXTitle("isobar #phi_{TY} [rad]");
	hTYData->SetYTitle("# of events");
	temp.phi_GJF.push_back(hTYData);
	if (twoMc) {
		TH1D* hTYMcPsp = new TH1D(("hTYMcPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Treiman-Yang Phi (McPsp)").c_str(), HISTLIMITS_PHI_BINS, HISTLIMITS_PHI_MIN, HISTLIMITS_PHI_MAX);
		hTYMcPsp->SetXTitle("isobar #phi_{TY} [rad]");
		hTYMcPsp->SetYTitle("# of events");
		temp.phi_GJF.push_back(hTYMcPsp);
		TH1D* hTYMcAcc = new TH1D(("hTYMcAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Treiman-Yang Phi (McAcc)").c_str(), HISTLIMITS_PHI_BINS, HISTLIMITS_PHI_MIN, HISTLIMITS_PHI_MAX);
		hTYMcAcc->SetXTitle("isobar #phi_{TY} [rad]");
		hTYMcAcc->SetYTitle("# of events");
		temp.phi_GJF.push_back(hTYMcAcc);
	} else {
		TH1D* hTYMc = new TH1D(("hTYMc_" + name_prefix).c_str(), (name_prefix + " Isobar Treiman-Yang Phi (Mc)").c_str(), HISTLIMITS_PHI_BINS, HISTLIMITS_PHI_MIN, HISTLIMITS_PHI_MAX);
		hTYMc->SetXTitle("isobar #phi_{TY} [rad]");
		hTYMc->SetYTitle("# of events");
		temp.phi_GJF.push_back(hTYMc);
	}

	TH2D* hTYVsGJData = new TH2D(("hTYVsGJData_" + name_prefix).c_str(), (name_prefix + " Isobar Treiman-Yang Phi vs Cos GJ Theta (Data)").c_str(), HISTLIMITS_PHI_BINS, HISTLIMITS_PHI_MIN, HISTLIMITS_PHI_MAX, HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
	hTYVsGJData->SetXTitle("isobar #phi_{TY} [rad]");
	hTYVsGJData->SetYTitle("isobar cos(#theta_{GJ})");
	hTYVsGJData->SetZTitle("# of events");
	hTYVsGJData->SetOption("COLZ");
	temp.phi_costheta_GJF.push_back(hTYVsGJData);
	if (twoMc) {
		TH2D* hTYVsGJMcPsp = new TH2D(("hTYVsGJMcPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Treiman-Yang Phi vs Cos GJ Theta (McPsp)").c_str(), HISTLIMITS_PHI_BINS, HISTLIMITS_PHI_MIN, HISTLIMITS_PHI_MAX, HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
		hTYVsGJMcPsp->SetXTitle("isobar #phi_{TY} [rad]");
		hTYVsGJMcPsp->SetYTitle("isobar cos(#theta_{GJ})");
		hTYVsGJMcPsp->SetZTitle("# of events");
		hTYVsGJMcPsp->SetOption("COLZ");
		temp.phi_costheta_GJF.push_back(hTYVsGJMcPsp);
		TH2D* hTYVsGJMcAcc = new TH2D(("hTYVsGJMcAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Treiman-Yang Phi vs Cos GJ Theta (McAcc)").c_str(), HISTLIMITS_PHI_BINS, HISTLIMITS_PHI_MIN, HISTLIMITS_PHI_MAX, HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
		hTYVsGJMcAcc->SetXTitle("isobar #phi_{TY} [rad]");
		hTYVsGJMcAcc->SetYTitle("isobar cos(#theta_{GJ})");
		hTYVsGJMcAcc->SetZTitle("# of events");
		hTYVsGJMcAcc->SetOption("COLZ");
		temp.phi_costheta_GJF.push_back(hTYVsGJMcAcc);
	} else {
		TH2D* hTYVsGJMc = new TH2D(("hTYVsGJMc_" + name_prefix).c_str(), (name_prefix + " Isobar Treiman-Yang Phi vs Cos GJ Theta (Mc)").c_str(), HISTLIMITS_PHI_BINS, HISTLIMITS_PHI_MIN, HISTLIMITS_PHI_MAX, HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
		hTYVsGJMc->SetXTitle("isobar #phi_{TY} [rad]");
		hTYVsGJMc->SetYTitle("isobar cos(#theta_{GJ})");
		hTYVsGJMc->SetZTitle("# of events");
		hTYVsGJMc->SetOption("COLZ");
		temp.phi_costheta_GJF.push_back(hTYVsGJMc);
	}

	if (twoMc) {
		THStack* hGJF_Stack_McPsp = new THStack(("hGJF_Stack_McPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (McPsp)").c_str());

		TH1D* hGJF_Flat_McPsp = new TH1D(("hGJF_Flat_McPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (McPsp)").c_str(), HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
		hGJF_Flat_McPsp->SetXTitle("isobar cos(#theta_{GJ})");
		hGJF_Flat_McPsp->SetYTitle("# of events");
		temp.costheta_GJF_Flat.push_back(hGJF_Flat_McPsp);
		hGJF_Stack_McPsp->Add(hGJF_Flat_McPsp);

		TH1D* hGJF_Neg_McPsp = new TH1D(("hGJF_Neg_McPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (McPsp)").c_str(), HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
		hGJF_Neg_McPsp->SetXTitle("isobar cos(#theta_{GJ})");
		hGJF_Neg_McPsp->SetYTitle("# of events");
		temp.costheta_GJF_NegRef.push_back(hGJF_Neg_McPsp);
		hGJF_Stack_McPsp->Add(hGJF_Neg_McPsp);

		TH1D* hGJF_Pos_McPsp = new TH1D(("hGJF_Pos_McPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (McPsp)").c_str(), HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
		hGJF_Pos_McPsp->SetXTitle("isobar cos(#theta_{GJ})");
		hGJF_Pos_McPsp->SetYTitle("# of events");
		temp.costheta_GJF_PosRef.push_back(hGJF_Pos_McPsp);
		hGJF_Stack_McPsp->Add(hGJF_Pos_McPsp);

		temp.costheta_GJF_Stack.push_back(hGJF_Stack_McPsp);

		THStack* hGJF_Stack_McAcc = new THStack(("hGJF_Stack_McAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (McAcc)").c_str());

		TH1D* hGJF_Flat_McAcc = new TH1D(("hGJF_Flat_McAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (McAcc)").c_str(), HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
		hGJF_Flat_McAcc->SetXTitle("isobar cos(#theta_{GJ})");
		hGJF_Flat_McAcc->SetYTitle("# of events");
		temp.costheta_GJF_Flat.push_back(hGJF_Flat_McAcc);
		hGJF_Stack_McAcc->Add(hGJF_Flat_McAcc);

		TH1D* hGJF_Neg_McAcc = new TH1D(("hGJF_Neg_McAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (McAcc)").c_str(), HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
		hGJF_Neg_McAcc->SetXTitle("isobar cos(#theta_{GJ})");
		hGJF_Neg_McAcc->SetYTitle("# of events");
		temp.costheta_GJF_NegRef.push_back(hGJF_Neg_McAcc);
		hGJF_Stack_McAcc->Add(hGJF_Neg_McAcc);

		TH1D* hGJF_Pos_McAcc = new TH1D(("hGJF_Pos_McAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (McAcc)").c_str(), HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
		hGJF_Pos_McAcc->SetXTitle("isobar cos(#theta_{GJ})");
		hGJF_Pos_McAcc->SetYTitle("# of events");
		temp.costheta_GJF_PosRef.push_back(hGJF_Pos_McAcc);
		hGJF_Stack_McAcc->Add(hGJF_Pos_McAcc);

		temp.costheta_GJF_Stack.push_back(hGJF_Stack_McAcc);
	} else {
		THStack* hGJF_Stack_Mc = new THStack(("hGJF_Stack_Mc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (Mc)").c_str());

		TH1D* hGJF_Flat_Mc = new TH1D(("hGJF_Flat_Mc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (Mc)").c_str(), HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
		hGJF_Flat_Mc->SetXTitle("isobar cos(#theta_{GJ})");
		hGJF_Flat_Mc->SetYTitle("# of events");
		temp.costheta_GJF_Flat.push_back(hGJF_Flat_Mc);
		hGJF_Stack_Mc->Add(hGJF_Flat_Mc);

		TH1D* hGJF_Neg_Mc = new TH1D(("hGJF_Neg_Mc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (Mc)").c_str(), HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
		hGJF_Neg_Mc->SetXTitle("isobar cos(#theta_{GJ})");
		hGJF_Neg_Mc->SetYTitle("# of events");
		temp.costheta_GJF_NegRef.push_back(hGJF_Neg_Mc);
		hGJF_Stack_Mc->Add(hGJF_Neg_Mc);

		TH1D* hGJF_Pos_Mc = new TH1D(("hGJF_Pos_Mc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (Mc)").c_str(), HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
		hGJF_Pos_Mc->SetXTitle("isobar cos(#theta_{GJ})");
		hGJF_Pos_Mc->SetYTitle("# of events");
		temp.costheta_GJF_PosRef.push_back(hGJF_Pos_Mc);
		hGJF_Stack_Mc->Add(hGJF_Pos_Mc);

		temp.costheta_GJF_Stack.push_back(hGJF_Stack_Mc);
	}

	return temp;
}

HelicityHistBunch HelicityHistBunchFactory(const std::string& name_prefix, const bool twoMc) {
	HelicityHistBunch temp;
	TH1D* hHThetaData = new TH1D(("hHThetaData_" + name_prefix).c_str(), (name_prefix + " Isobar Helicity Cos Theta (Data)").c_str(), HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
	hHThetaData->SetXTitle("cos(#theta_{hel}) of #pi^{0} from isobar");
	hHThetaData->SetYTitle("# of events");
	temp.costheta_HF.push_back(hHThetaData);
	if (twoMc) {
		TH1D* hHThetaMcPsp = new TH1D(("hHThetaMcPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Helicity Cos Theta (McPsp)").c_str(), HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
		hHThetaMcPsp->SetXTitle("cos(#theta_{hel} of #pi^{0} from isobar)");
		hHThetaMcPsp->SetYTitle("# of events");
		temp.costheta_HF.push_back(hHThetaMcPsp);
		TH1D* hHThetaMcAcc = new TH1D(("hHThetaMcAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Helicity Cos Theta (McAcc)").c_str(), HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
		hHThetaMcAcc->SetXTitle("cos(#theta_{hel} of #pi^{0} from isobar)");
		hHThetaMcAcc->SetYTitle("# of events");
		temp.costheta_HF.push_back(hHThetaMcAcc);
	} else {
		TH1D* hHThetaMc = new TH1D(("hHThetaMc_" + name_prefix).c_str(), (name_prefix + " Isobar Helicity Cos Theta (Mc)").c_str(), HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
		hHThetaMc->SetXTitle("cos(#theta_{hel} of #pi^{0} from isobar)");
		hHThetaMc->SetYTitle("# of events");
		temp.costheta_HF.push_back(hHThetaMc);
	}

	TH1D* hHPhiData = new TH1D(("hHPhiData_" + name_prefix).c_str(), (name_prefix + " Isobar Helicity Phi (Data)").c_str(), HISTLIMITS_PHI_BINS, HISTLIMITS_PHI_MIN, HISTLIMITS_PHI_MAX);
	hHPhiData->SetXTitle("#phi_{hel} [rad] of #pi^{0} from isobar");
	hHPhiData->SetYTitle("# of events");
	temp.phi_HF.push_back(hHPhiData);
	if (twoMc) {
		TH1D* hHPhiMcPsp = new TH1D(("hHPhiMcPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Helicity Phi (McPsp)").c_str(), HISTLIMITS_PHI_BINS, HISTLIMITS_PHI_MIN, HISTLIMITS_PHI_MAX);
		hHPhiMcPsp->SetXTitle("#phi_{hel} [rad] of #pi^{0} from isobar");
		hHPhiMcPsp->SetYTitle("# of events");
		temp.phi_HF.push_back(hHPhiMcPsp);
		TH1D* hHPhiMcAcc = new TH1D(("hHPhiMcAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Helicity Phi (McAcc)").c_str(), HISTLIMITS_PHI_BINS, HISTLIMITS_PHI_MIN, HISTLIMITS_PHI_MAX);
		hHPhiMcAcc->SetXTitle("#phi_{hel} [rad] of #pi^{0} from isobar");
		hHPhiMcAcc->SetYTitle("# of events");
		temp.phi_HF.push_back(hHPhiMcAcc);
	} else {
		TH1D* hHPhiMc = new TH1D(("hHPhiMc_" + name_prefix).c_str(), (name_prefix + " Isobar Helicity Phi (Mc)").c_str(), HISTLIMITS_PHI_BINS, HISTLIMITS_PHI_MIN, HISTLIMITS_PHI_MAX);
		hHPhiMc->SetXTitle("#phi_{hel} [rad] of #pi^{0} from isobar");
		hHPhiMc->SetYTitle("# of events");
		temp.phi_HF.push_back(hHPhiMc);
	}

	TH2D* hHPhiVsHThetaData = new TH2D(("hHPhiVsHThetaData_" + name_prefix).c_str(), (name_prefix + " Isobar Helicity Phi vs Cos Theta (Data)").c_str(), HISTLIMITS_PHI_BINS, HISTLIMITS_PHI_MIN, HISTLIMITS_PHI_MAX, HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
	hHPhiVsHThetaData->SetXTitle("isobar #phi_{hel} [rad]");
	hHPhiVsHThetaData->SetYTitle("cos(#theta_{hel} of #pi^{0} from isobar)");
	hHPhiVsHThetaData->SetZTitle("# of events");
	hHPhiVsHThetaData->SetOption("COLZ");
	temp.phi_costheta_HF.push_back(hHPhiVsHThetaData);
	if (twoMc) {
		TH2D* hHPhiVsHThetaMcPsp = new TH2D(("hHPhiVsHThetaMcPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Helicity Phi vs Cos Theta (McPsp)").c_str(), HISTLIMITS_PHI_BINS, HISTLIMITS_PHI_MIN, HISTLIMITS_PHI_MAX, HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
		hHPhiVsHThetaMcPsp->SetXTitle("isobar #phi_{hel} [rad]");
		hHPhiVsHThetaMcPsp->SetYTitle("cos(#theta_{hel} of #pi^{0} from isobar)");
		hHPhiVsHThetaMcPsp->SetZTitle("# of events");
		hHPhiVsHThetaMcPsp->SetOption("COLZ");
		temp.phi_costheta_HF.push_back(hHPhiVsHThetaMcPsp);
		TH2D* hHPhiVsHThetaMcAcc = new TH2D(("hHPhiVsHThetaMcAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Helicity Phi vs Cos Theta (McAcc)").c_str(), HISTLIMITS_PHI_BINS, HISTLIMITS_PHI_MIN, HISTLIMITS_PHI_MAX, HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
		hHPhiVsHThetaMcAcc->SetXTitle("isobar #phi_{hel} [rad]");
		hHPhiVsHThetaMcAcc->SetYTitle("cos(#theta_{hel} of #pi^{0} from isobar)");
		hHPhiVsHThetaMcAcc->SetZTitle("# of events");
		hHPhiVsHThetaMcAcc->SetOption("COLZ");
		temp.phi_costheta_HF.push_back(hHPhiVsHThetaMcAcc);
	} else {
		TH2D* hHPhiVsHThetaMc = new TH2D(("hHPhiVsHThetaMc_" + name_prefix).c_str(), (name_prefix + " Isobar Helicity Phi vs Cos Theta (Mc)").c_str(), HISTLIMITS_PHI_BINS, HISTLIMITS_PHI_MIN, HISTLIMITS_PHI_MAX, HISTLIMITS_COSTHETA_BINS, HISTLIMITS_COSTHETA_MIN, HISTLIMITS_COSTHETA_MAX);
		hHPhiVsHThetaMc->SetXTitle("isobar #phi_{hel} [rad]");
		hHPhiVsHThetaMc->SetYTitle("cos(#theta_{hel} of #pi^{0} from isobar)");
		hHPhiVsHThetaMc->SetZTitle("# of events");
		hHPhiVsHThetaMc->SetOption("COLZ");
		temp.phi_costheta_HF.push_back(hHPhiVsHThetaMc);
	}

	return temp;
}

void fillWeightedHelicityAnglePlots(const HelicityAngles &ha, double weight, unsigned int tree_index, HelicityHistBunch &hhb) {
	hhb.costheta_HF[tree_index]->Fill(ha.cosTheta, weight);
	hhb.phi_HF[tree_index]->Fill(ha.phi, weight);
	hhb.phi_costheta_HF[tree_index]->Fill(ha.phi, ha.cosTheta, weight);
}

void fillWeightedGJAnglePlots(const TLorentzVector &isobar, double weight, double weightPosRef, double weightNegRef, double weightFlat, double tprime, unsigned int tree_index, GJHistBunch &hBunch) {
	hBunch.costheta_GJF[tree_index]->Fill(isobar.CosTheta(), weight);
	if (tree_index != 0) {
		hBunch.costheta_GJF_PosRef[tree_index-1]->Fill(isobar.CosTheta(), weightPosRef);
		hBunch.costheta_GJF_NegRef[tree_index-1]->Fill(isobar.CosTheta(), weightNegRef);
		hBunch.costheta_GJF_Flat[tree_index-1]->Fill(isobar.CosTheta(), weightFlat);
	}
	hBunch.costheta_GJF_tprime[tree_index]->Fill(isobar.CosTheta(), tprime, weight);
	hBunch.phi_GJF[tree_index]->Fill(isobar.Phi(), weight);
	hBunch.isobar_mass[tree_index]->Fill(isobar.M(), weight);
	hBunch.phi_costheta_GJF[tree_index]->Fill(isobar.Phi(), isobar.CosTheta(), weight);
}

HelicityAngles calculateHelicityAngles(const TLorentzVector& isobar, TLorentzVector particle, TLorentzVector *beam = NULL, bool first = true) {
	HelicityAngles temp;

	TVector3 zaxis_gjf;
	if(beam != NULL) {
		zaxis_gjf = beam->Vect().Unit();
	}
	else {
		zaxis_gjf = TVector3(0,0,1);
	}
	// create helicity frame coordinate system
	TVector3 zaxis = isobar.Vect().Unit();
	TVector3 yaxis = zaxis_gjf.Cross(zaxis).Unit();
	TVector3 xaxis = yaxis.Cross(zaxis);

	// boost NParticleState into isobar rest frame
	const TVector3 boost_vec = isobar.BoostVector();
	particle.Boost(-boost_vec);
	// theta
	double theta = zaxis.Angle(particle.Vect());
	temp.cosTheta = TMath::Cos(theta);
	double fX = xaxis.Dot(particle.Vect());
	double fY = yaxis.Dot(particle.Vect());
	if(fX == 0.0 && fY == 0.0)
		temp.phi = 0.0;
	else
		temp.phi = TMath::ATan2(fY,fX);

	return temp;
}

void makeDifferencePlots(TDirectory* dir) {
	// create list of histograms to compare data and mc, and to calculate
	// the acceptance
	TList plotsDifference;
	TList plotsAcceptance;
	TIter histiter(dir->GetListOfKeys());
	TKey* key;
	while ((key = dynamic_cast<TKey*>(histiter()))) {
		const std::string s(key->GetName());
		if ((s.length() >= 2 && s.substr(s.length()-2, 2) == "Mc") ||
		    (s.length() >= 5 && s.substr(s.length()-5, 5) == "McPsp") ||
		    (s.length() >= 5 && s.substr(s.length()-5, 5) == "McAcc") ||
		    s.find("Mc_") != std::string::npos ||
		    s.find("McPsp_") != std::string::npos ||
		    s.find("McAcc_") != std::string::npos) {
			plotsDifference.Add(key);
		}
		if ((s.length() >= 5 && s.substr(s.length()-5, 5) == "McAcc") ||
		    s.find("McAcc_") != std::string::npos) {
			plotsAcceptance.Add(key);
		}
	}

	// process histograms for acceptance
	// process those first as some histograms are rescaled for creating the
	// difference plots
	histiter = TIter(&plotsAcceptance);
	while ((key = dynamic_cast<TKey*>(histiter()))) {
		// generate difference histograms
		const std::string nameMcAcc(key->GetName());
		const size_t pos = nameMcAcc.find("McAcc");

		std::string nameMcPsp(nameMcAcc);
		nameMcPsp.erase(pos, 5);
		nameMcPsp.insert(pos, "McPsp");

		std::string nameAcceptance(nameMcAcc);
		nameAcceptance.erase(pos, 5);
		nameAcceptance.insert(pos, "Acceptance");

		TH1* histMcAcc;
		dir->GetObject(nameMcAcc.c_str(), histMcAcc);
		TH1* histMcPsp;
		dir->GetObject(nameMcPsp.c_str(), histMcPsp);

		if (histMcAcc && histMcPsp) {
			dir->cd();

			TH1* histAcceptance = dynamic_cast<TH1*>(histMcAcc->Clone(nameAcceptance.c_str()));
			assert(histAcceptance != NULL);
			histAcceptance->SetTitle("");
			histAcceptance->Divide(histMcPsp);
			if (dynamic_cast<TH2*>(histAcceptance) != NULL)
				histAcceptance->SetZTitle("acceptance (acc. MC / PS MC)");
			else
				histAcceptance->SetYTitle("acceptance (acc. MC / PS MC)");
			histAcceptance->Write(NULL, TObject::kOverwrite);
		}
	}

	// process histograms for difference
	histiter = TIter(&plotsDifference);
	while ((key = dynamic_cast<TKey*>(histiter()))) {
		// generate difference histograms
		std::string hnamemc(key->GetName());
		int pos = hnamemc.find("Mc");
		// create new string with MC exchanged for Diff
		std::string hnamediff(hnamemc);
		hnamediff.erase(pos, 2);
		hnamediff.insert(pos, "Diff");
		// create new string with MC exchanged for RelDiff
		std::string hnamereldiff(hnamemc);
		hnamereldiff.erase(pos, 2);
		hnamereldiff.insert(pos, "RelDiff");
		// create new string with MC exchanged for Data
		std::string hnamedata(hnamemc);
		if (hnamemc.substr(pos, 5) == "McPsp" || hnamemc.substr(pos, 5) == "McAcc") {
			hnamedata.erase(pos, 5);
		} else {
			hnamedata.erase(pos, 2);
		}
		hnamedata.insert(pos, "Data");

		TH1* mchist;
		dir->GetObject(hnamemc.c_str(), mchist);
		TH1* datahist;
		dir->GetObject(hnamedata.c_str(), datahist);

		if (datahist && mchist) {
			dir->cd();

			TH1* diffhist = dynamic_cast<TH1*>(mchist->Clone(hnamediff.c_str()));
			assert(diffhist != NULL);
			diffhist->SetTitle("");
			diffhist->Add(datahist, -1.);
			if (dynamic_cast<TH2*>(diffhist) != NULL)
				diffhist->SetZTitle("# of events difference(MC-Data)");
			else
				diffhist->SetYTitle("# of events difference(MC-Data)");
			diffhist->SetMaximum();
			diffhist->SetMinimum();
			double max = std::max(std::abs(diffhist->GetMaximum()), std::abs(diffhist->GetMinimum()));
			if (max == 0.)
				max = 1. / 1.1;
			diffhist->SetMaximum( 1.1 * max);
			diffhist->SetMinimum(-1.1 * max);
			diffhist->Write(NULL, TObject::kOverwrite);

			TH1* reldiffhist = dynamic_cast<TH1*>(diffhist->Clone(hnamereldiff.c_str()));
			assert(reldiffhist != NULL);
			reldiffhist->Divide(datahist);
			if (dynamic_cast<TH2*>(reldiffhist) != NULL)
				reldiffhist->SetZTitle("relative difference((MC-Data)/Data)");
			else
				reldiffhist->SetYTitle("relative difference((MC-Data)/Data)");
			reldiffhist->SetMaximum();
			reldiffhist->SetMinimum();
			max = std::max(std::abs(reldiffhist->GetMaximum()), std::abs(reldiffhist->GetMinimum()));
			if (max == 0.)
				max = 1. / 1.1;
			reldiffhist->SetMaximum( 1.1 * max);
			reldiffhist->SetMinimum(-1.1 * max);
			reldiffhist->Write(NULL, TObject::kOverwrite);
		}
	}
}

void
createWeightedPlots(const std::string& dataFileName,
                    const std::string& dataTreeName,
                    const std::string& dataProdKinPartNamesName,
                    const std::string& dataProdKinMomentaName,
                    const std::string& dataDecayKinPartNamesName,
                    const std::string& dataDecayKinMomentaName,
                    const std::string& mcPspFileName,
                    const std::string& mcPspTreeName,
                    const std::string& mcPspProdKinPartNamesName,
                    const std::string& mcPspProdKinMomentaName,
                    const std::string& mcPspDecayKinPartNamesName,
                    const std::string& mcPspDecayKinMomentaName,
                    const std::string& mcPspWeightFileName,
                    const std::string& mcPspWeightTreeName,
                    const std::string& mcAccFileName,
                    const std::string& mcAccTreeName,
                    const std::string& mcAccProdKinPartNamesName,
                    const std::string& mcAccProdKinMomentaName,
                    const std::string& mcAccDecayKinPartNamesName,
                    const std::string& mcAccDecayKinMomentaName,
                    const std::string& mcAccWeightFileName,
                    const std::string& mcAccWeightTreeName,
                    const std::string& massBin,
                    const std::string& outFileName,
                    const std::string& pdgFileName,
                    const long int     treeCacheSize)
{
	// keep track of the processing time
	TStopwatch timer;
	timer.Start();

	// initialize particle data table
	rpwa::particleDataTable::readFile(pdgFileName);

	// open data file
	TFile* dataFile = TFile::Open(dataFileName.c_str());

	// tree containg the real data events
	TTree* dataTree;
	dataFile->GetObject(dataTreeName.c_str(), dataTree);

	// names of particles in tree
	TClonesArray* dataProdKinPartNames(NULL);
	dataFile->GetObject(dataProdKinPartNamesName.c_str(), dataProdKinPartNames);
	assert(dataProdKinPartNames->GetEntries() == 1);
	for (Int_t i=0; i<dataProdKinPartNames->GetEntries(); i++) {
		if (!rpwa::particleDataTable::isInTable(((TObjString*)(dataProdKinPartNames->At(i)))->String().Data())) {
			printErr << "Unknown particle '" << ((TObjString*)(dataProdKinPartNames->At(i)))->String() << "' found in input tree '" << dataFileName << ":" << dataTreeName << "'." << std::endl;
			return;
		}
	}

	TClonesArray* dataDecayKinPartNames(NULL);
	dataFile->GetObject(dataDecayKinPartNamesName.c_str(), dataDecayKinPartNames);
	assert(dataDecayKinPartNames->GetEntries() == 3);
	for (Int_t i=0; i<dataDecayKinPartNames->GetEntries(); i++) {
		if (!rpwa::particleDataTable::isInTable(((TObjString*)(dataDecayKinPartNames->At(i)))->String().Data())) {
			printErr << "Unknown particle '" << ((TObjString*)(dataDecayKinPartNames->At(i)))->String() << "' found in input tree '" << dataFileName << ":" << dataTreeName << "'." << std::endl;
			return;
		}
	}

	// open weighted MC file
	TFile* mcPspFile = TFile::Open(mcPspFileName.c_str());

	// tree containing the phase space events
	TTree* mcPspTree;
	mcPspFile->GetObject(mcPspTreeName.c_str(), mcPspTree);

	// add friend containing weights
	mcPspTree->AddFriend(mcPspWeightTreeName.c_str(), mcPspWeightFileName.c_str());

	// names of particles in tree
	TClonesArray* mcPspProdKinPartNames(NULL);
	mcPspFile->GetObject(mcPspProdKinPartNamesName.c_str(), mcPspProdKinPartNames);
	assert(mcPspProdKinPartNames->GetEntries() == 1);
	for (Int_t i=0; i<mcPspProdKinPartNames->GetEntries(); i++) {
		if (!rpwa::particleDataTable::isInTable(((TObjString*)(mcPspProdKinPartNames->At(i)))->String().Data())) {
			printErr << "Unknown particle '" << ((TObjString*)(mcPspProdKinPartNames->At(i)))->String() << "' found in input tree '" << mcPspFileName << ":" << mcPspTreeName << "'." << std::endl;
			return;
		}
	}

	TClonesArray* mcPspDecayKinPartNames(NULL);
	mcPspFile->GetObject(mcPspDecayKinPartNamesName.c_str(), mcPspDecayKinPartNames);
	assert(mcPspDecayKinPartNames->GetEntries() == 3);
	for (Int_t i=0; i<mcPspDecayKinPartNames->GetEntries(); i++) {
		if (!rpwa::particleDataTable::isInTable(((TObjString*)(mcPspDecayKinPartNames->At(i)))->String().Data())) {
			printErr << "Unknown particle '" << ((TObjString*)(mcPspDecayKinPartNames->At(i)))->String() << "' found in input tree '" << mcPspFileName << ":" << mcPspTreeName << "'." << std::endl;
			return;
		}
	}

	// open weighted MC file
	TFile* mcAccFile(NULL);
	if (mcAccFileName != "") {
		mcAccFile = TFile::Open(mcAccFileName.c_str());
	}

	// tree containing the phase space events
	TTree* mcAccTree(NULL);
	if (mcAccFile != NULL) {
		mcAccFile->GetObject(mcAccTreeName.c_str(), mcAccTree);
	}

	// add friend containing weights
	if (mcAccFile != NULL) {
		mcAccTree->AddFriend(mcAccWeightTreeName.c_str(), mcAccWeightFileName.c_str());
	}

	// names of particles in tree
	TClonesArray* mcAccProdKinPartNames(NULL);
	if (mcAccFile != NULL) {
		mcAccFile->GetObject(mcAccProdKinPartNamesName.c_str(), mcAccProdKinPartNames);
		assert(mcAccProdKinPartNames->GetEntries() == 1);
		for (Int_t i=0; i<mcAccProdKinPartNames->GetEntries(); i++) {
			if (!rpwa::particleDataTable::isInTable(((TObjString*)(mcAccProdKinPartNames->At(i)))->String().Data())) {
				printErr << "Unknown particle '" << ((TObjString*)(mcAccProdKinPartNames->At(i)))->String() << "' found in input tree '" << mcAccFileName << ":" << mcAccTreeName << "'." << std::endl;
				return;
			}
		}
	}

	TClonesArray* mcAccDecayKinPartNames(NULL);
	if (mcAccFile != NULL) {
		mcAccFile->GetObject(mcAccDecayKinPartNamesName.c_str(), mcAccDecayKinPartNames);
		assert(mcAccDecayKinPartNames->GetEntries() == 3);
		for (Int_t i=0; i<mcAccDecayKinPartNames->GetEntries(); i++) {
			if (!rpwa::particleDataTable::isInTable(((TObjString*)(mcAccDecayKinPartNames->At(i)))->String().Data())) {
				printErr << "Unknown particle '" << ((TObjString*)(mcAccDecayKinPartNames->At(i)))->String() << "' found in input tree '" << mcAccFileName << ":" << mcAccTreeName << "'." << std::endl;
				return;
			}
		}
	}

	double massval = 0.0;
	unsigned int datatreeentries = 0;

	unsigned int pointpos = massBin.find(".");
	if(pointpos == 0 || pointpos == massBin.size())
		std::cout<<"Warning: Bad massbin name!"<<std::endl;
	std::string masshigh = massBin.substr(pointpos+1);
	massval = atof(masshigh.c_str());
	massval /=1000;
	datatreeentries = dataTree->GetEntries();


	gROOT->SetStyle("Plain");
	TFile* outfile = TFile::Open(outFileName.c_str(), "UPDATE");
	outfile->cd();

	TDirectory* outdir = gDirectory->GetDirectory(massBin.c_str());
	if (outdir == NULL) {
		outdir = gDirectory->mkdir(massBin.c_str());
	}
	outdir->cd();

	// --------------- global diagrams
	std::vector<TH1D*> hM;
	TH1D* hMData = new TH1D("hResMassData", "Mass (Data)", HISTLIMITS_3PI_MASS_BINS, HISTLIMITS_3PI_MASS_MIN, HISTLIMITS_3PI_MASS_MAX);
	hM.push_back(hMData);
	if (mcAccFile != NULL) {
		TH1D* hMMcPsp = new TH1D("hResMassMcPsp", "Mass (McPsp)", HISTLIMITS_3PI_MASS_BINS, HISTLIMITS_3PI_MASS_MIN, HISTLIMITS_3PI_MASS_MAX);
		hM.push_back(hMMcPsp);
		TH1D* hMMcAcc = new TH1D("hResMassMcAcc", "Mass (McAcc)", HISTLIMITS_3PI_MASS_BINS, HISTLIMITS_3PI_MASS_MIN, HISTLIMITS_3PI_MASS_MAX);
		hM.push_back(hMMcAcc);
	} else {
		TH1D* hMMc = new TH1D("hResMassMc", "Mass (Mc)", HISTLIMITS_3PI_MASS_BINS, HISTLIMITS_3PI_MASS_MIN, HISTLIMITS_3PI_MASS_MAX);
		hM.push_back(hMMc);
	}

	// Dalitz plots
	std::vector<TH2D*> dalitz_neutral;
	TH2D* hDalitzData = new TH2D("hDalitzData", "Dalitz Plot #pi^{0}#pi^{0} vs. #pi^{-}#pi^{0} (Data)", HISTLIMITS_2PI_MASS2_BINS, HISTLIMITS_2PI_MASS2_MIN, HISTLIMITS_2PI_MASS2_MAX, HISTLIMITS_2PI_MASS2_BINS, HISTLIMITS_2PI_MASS2_MIN, HISTLIMITS_2PI_MASS2_MAX);
	hDalitzData->SetXTitle("mass^{2}(#pi^{0}#pi^{0}) [GeV^{2}/c^{4}]");
	hDalitzData->SetYTitle("mass^{2}(#pi^{-}#pi^{0}) [GeV^{2}/c^{4}]");
	hDalitzData->SetZTitle("# of events");
	hDalitzData->SetOption("COLZ");
	hDalitzData->SetStats(0);
	dalitz_neutral.push_back(hDalitzData);
	if (mcAccFile != NULL) {
		TH2D* hDalitzMcPsp = new TH2D("hDalitzMcPsp", "Dalitz Plot #pi^{0}#pi^{0} vs. #pi^{-}#pi^{0} (McPsp)", HISTLIMITS_2PI_MASS2_BINS, HISTLIMITS_2PI_MASS2_MIN, HISTLIMITS_2PI_MASS2_MAX, HISTLIMITS_2PI_MASS2_BINS, HISTLIMITS_2PI_MASS2_MIN, HISTLIMITS_2PI_MASS2_MAX);
		hDalitzMcPsp->SetXTitle("mass^{2}(#pi^{0}#pi^{0}) [GeV^{2}/c^{4}]");
		hDalitzMcPsp->SetYTitle("mass^{2}(#pi^{-}#pi^{0}) [GeV^{2}/c^{4}]");
		hDalitzMcPsp->SetZTitle("# of events");
		hDalitzMcPsp->SetOption("COLZ");
		hDalitzMcPsp->SetStats(0);
		dalitz_neutral.push_back(hDalitzMcPsp);
		TH2D* hDalitzMcAcc = new TH2D("hDalitzMcAcc", "Dalitz Plot #pi^{0}#pi^{0} vs. #pi^{-}#pi^{0} (McAcc)", HISTLIMITS_2PI_MASS2_BINS, HISTLIMITS_2PI_MASS2_MIN, HISTLIMITS_2PI_MASS2_MAX, HISTLIMITS_2PI_MASS2_BINS, HISTLIMITS_2PI_MASS2_MIN, HISTLIMITS_2PI_MASS2_MAX);
		hDalitzMcAcc->SetXTitle("mass^{2}(#pi^{0}#pi^{0}) [GeV^{2}/c^{4}]");
		hDalitzMcAcc->SetYTitle("mass^{2}(#pi^{-}#pi^{0}) [GeV^{2}/c^{4}]");
		hDalitzMcAcc->SetZTitle("# of events");
		hDalitzMcAcc->SetOption("COLZ");
		hDalitzMcAcc->SetStats(0);
		dalitz_neutral.push_back(hDalitzMcAcc);
	} else {
		TH2D* hDalitzMc = new TH2D("hDalitzMc", "Dalitz Plot #pi^{0}#pi^{0} vs. #pi^{-}#pi^{0} (Mc)", HISTLIMITS_2PI_MASS2_BINS, HISTLIMITS_2PI_MASS2_MIN, HISTLIMITS_2PI_MASS2_MAX, HISTLIMITS_2PI_MASS2_BINS, HISTLIMITS_2PI_MASS2_MIN, HISTLIMITS_2PI_MASS2_MAX);
		hDalitzMc->SetXTitle("mass^{2}(#pi^{0}#pi^{0}) [GeV^{2}/c^{4}]");
		hDalitzMc->SetYTitle("mass^{2}(#pi^{-}#pi^{0}) [GeV^{2}/c^{4}]");
		hDalitzMc->SetZTitle("# of events");
		hDalitzMc->SetOption("COLZ");
		hDalitzMc->SetStats(0);
		dalitz_neutral.push_back(hDalitzMc);
	}

	// --------------- generate histogram bunches
	// neutral isobar
	// GJ Histogram Bunch
	GJHistBunch GJHB_neutral_isobar = GJHistBunchFactory("Neutral", mcAccFile!=NULL);
	// Helicity Histogram Bunch
	HelicityHistBunch HHB_neutral_isobar = HelicityHistBunchFactory("Neutral", mcAccFile!=NULL);

	// charged isobar
	// GJ Histogram Bunch MC
	GJHistBunch GJHB_charged_isobar = GJHistBunchFactory("Charged", mcAccFile!=NULL);
	// Helicity Histogram Bunch
	HelicityHistBunch HHB_charged_isobar = HelicityHistBunchFactory("Charged", mcAccFile!=NULL);

	// charged isobar with rho mass cut
	// GJ Histogram Bunch MC
	GJHistBunch GJHB_rho = GJHistBunchFactory("ChargedRho", mcAccFile!=NULL);
	// Helicity Histogram Bunch
	HelicityHistBunch HHB_rho = HelicityHistBunchFactory("ChargedRho", mcAccFile!=NULL);

	//Loop both over data and mc tree
	// itree = 0: data tree
	// itree = 1: mc tree
	// itree = 1: (accepted) mc tree
	for (unsigned int itree = 0; itree < 3; ++itree) {
		TTree* tree(NULL);
		if      (itree == 0) { tree = dataTree; }
		else if (itree == 1) { tree = mcPspTree; }
		else if (itree == 2) { tree = mcAccTree; }
		if (tree == NULL)
			continue;

		if (itree == 0) {
			printInfo << "Step 1: creating plots for data events" << std::endl;
		} else if (itree == 1) {
			printInfo << "Step 2: creating plots for phase-space events" << std::endl;
		} else {
			printInfo << "Step 3: creating plots for accepted phase-space events" << std::endl;
		}

		// in case of MC tree connect the weights
		double weight = 1;
		double weightPosRef = 1.;
		double weightNegRef = 1.;
		double weightFlat = 1.;
		double impweight = 1;
		if (itree != 0) {
			tree->SetBranchAddress("weight", &weight);
			tree->SetBranchAddress("weightPosRef", &weightPosRef);
			tree->SetBranchAddress("weightNegRef", &weightNegRef);
			tree->SetBranchAddress("weightFlat", &weightFlat);
			tree->SetBranchAddress("impweight", &impweight);
		}

		TClonesArray* prodKinPartNames  = NULL;
		TClonesArray* decayKinPartNames = NULL;
		std::string   prodKinMomentaName;
		std::string   decayKinMomentaName;
		TClonesArray* prodKinMomenta    = NULL;
		TClonesArray* decayKinMomenta   = NULL;

		if (itree == 0) {
			prodKinPartNames    = dataProdKinPartNames;
			decayKinPartNames   = dataDecayKinPartNames;
			prodKinMomentaName  = dataProdKinMomentaName;
			decayKinMomentaName = dataDecayKinMomentaName;
		} else if (itree == 1) {
			prodKinPartNames    = mcPspProdKinPartNames;
			decayKinPartNames   = mcPspDecayKinPartNames;
			prodKinMomentaName  = mcPspProdKinMomentaName;
			decayKinMomentaName = mcPspDecayKinMomentaName;
		} else {
			prodKinPartNames    = mcAccProdKinPartNames;
			decayKinPartNames   = mcAccDecayKinPartNames;
			prodKinMomentaName  = mcAccProdKinMomentaName;
			decayKinMomentaName = mcAccDecayKinMomentaName;
		}

		tree->SetBranchAddress(prodKinMomentaName.c_str(),  &prodKinMomenta);
		tree->SetBranchAddress(decayKinMomentaName.c_str(), &decayKinMomenta);

		tree->SetCacheSize(treeCacheSize);
		if (itree != 0) {
			tree->AddBranchToCache("weight", true);
			tree->AddBranchToCache("weightPosRef", true);
			tree->AddBranchToCache("weightNegRef", true);
			tree->AddBranchToCache("weightFlat", true);
			tree->AddBranchToCache("impweight", true);
		}
		tree->AddBranchToCache(prodKinMomentaName.c_str(),  true);
		tree->AddBranchToCache(decayKinMomentaName.c_str(), true);
		tree->StopCacheLearningPhase();

		TLorentzVector lvBeam;
		int qBeam;

		TLorentzVector lvPart[3];
		int qPart[3];

		TLorentzVector lvOut;

		// loop over tree entries
		const long int nmbEvents = tree->GetEntries();
		double maxweight = 0;
		double avweight = 0;
		boost::progress_display progressIndicator(nmbEvents, cout, "");
		for (long int i = 0; i < nmbEvents; ++i) {
			++progressIndicator;
			tree->GetEntry(i);

			assert(prodKinMomenta->GetEntries() == 1);
			for (Int_t j=0; j<prodKinMomenta->GetEntries(); j++) {
				const rpwa::particleProperties* pp = rpwa::particleDataTable::entry(((TObjString*)(prodKinPartNames->At(j)))->String().Data());
				if (pp == NULL) {
					printErr << "Unknown particle '" << ((TObjString*)(prodKinPartNames->At(j)))->String() << "' found in input tree." << std::endl;
					return;
				}
				lvBeam.SetVectM(*((TVector3*)(prodKinMomenta->At(j))), pp->mass());
				qBeam = pp->charge();
			}

			lvOut.SetXYZM(0., 0., 0., 0.);
			assert(decayKinMomenta->GetEntries() == 3);
			for (Int_t j=0; j<decayKinMomenta->GetEntries(); j++) {
				const rpwa::particleProperties* pp = rpwa::particleDataTable::entry(((TObjString*)(decayKinPartNames->At(j)))->String().Data());
				if (pp == NULL) {
					printErr << "Unknown particle '" << ((TObjString*)(decayKinPartNames->At(j)))->String() << "' found in input tree." << std::endl;
					return;
				}

				lvPart[j].SetVectM(*((TVector3*)(decayKinMomenta->At(j))), pp->mass());
				qPart[j] = pp->charge();

				lvOut += lvPart[j];
			}
			// in case its data tree (itree=0) put weights to 1
			if (itree == 0) {
				weight = 1;
				impweight = 1;
			}
			if (weight == 0)
				continue;

			if (impweight != 0)
				weight /= impweight;

			if (weight > maxweight)
				maxweight = weight;
			avweight += weight;

			double tprime;
			{
				TVector3 dir=lvBeam.Vect();
				double const mpi=lvBeam.M();
				double k=sqrt(lvOut.E()*lvOut.E()-mpi*mpi)/dir.Mag();
				dir*=k;
				TLorentzVector beam;
				beam.SetVectM(dir,mpi);
				tprime = -(beam-lvOut).M2();
			}

			unsigned int npart = 3;

			// ----dalitz plots
			// make all combinations
			std::vector<std::pair<std::pair<int, int>, double> > comb;
			for(unsigned int i=0; i < npart; i++) {
				for(unsigned int j=0; j < i; j++) {
					TLorentzVector lv = lvPart[i] + lvPart[j];
					double mass2 = lv.M2();
					std::pair<int, int> charge(qPart[i], qPart[j]);
					std::pair<std::pair<int, int>, double> temp;
					temp = make_pair(charge, mass2);
					comb.push_back(temp);
				}
			}

			// now fill corresponding dalitz plots
			// this part is pi-pi0pi0 specific

			// put pi0 pi0 combination to front of vector to make combination picks easier later on
			for(unsigned int i = 0; i < comb.size(); i++) {
				if(getTotalCharge(comb[i].first) == 0) {
					if(i == 0) break;
					else {
						std::pair<std::pair<int, int>, double> temp = comb[0];
						comb[0] = comb[i];
						comb[i] = temp;
						break;
					}
				}
			}
			// now actually fill the histograms
			dalitz_neutral[itree]->Fill(comb[0].second, comb[1].second, weight);
			//dalitz_charged[itree]->Fill(comb[1].second, comb[2].second, weight);

			// transform into GJ
			{
				TLorentzVector tempX=lvOut;
				// rotate event into scattering plane
				// get normal vector
				TVector3 y(0,1,0);
				TVector3 N=lvBeam.Vect().Cross(tempX.Vect());
				TVector3 rot=N.Cross(y);
				TRotation t;
				double a=N.Angle(y);
				t.Rotate(a,rot);

				TLorentzRotation T(t);
				TLorentzRotation L1(T);
				tempX*=T;
				lvBeam.Transform(T);

				// boost to X rest frame
				TVector3 boost=-tempX.BoostVector();
				TLorentzRotation b;
				b.Boost(boost);
				tempX*=b;
				lvBeam.Transform(b);

				// put beam along z-axis
				TVector3 beamdir=lvBeam.Vect();
				a=beamdir.Angle(TVector3(0,0,1));
				TRotation t2;
				t2.Rotate(a,TVector3(0,1,0));
				T=TLorentzRotation(t2);
				lvBeam.Transform(T);

				// transform pions
				lvOut.SetXYZM(0., 0., 0., 0.);
				for(unsigned int i=0; i<npart; ++i){
					lvPart[i].Transform(L1);
					lvPart[i].Transform(b);
					lvPart[i].Transform(T);

					lvOut += lvPart[i];
				}
			}

			hM[itree]->Fill(lvOut.M(), weight);

			// loop over all states that contain 2 final state particles
			// and plot angles
			unsigned int nstates = 3;
			const unsigned int mapstates[3][2] = { {0, 1}, {0, 2}, {1, 2} };
			//cout<<"number of substates: "<<event.nStates()<<endl;
			for (unsigned int is = 0; is < nstates; ++is) {
				const TLorentzVector lvIsobar = lvPart[mapstates[is][0]] + lvPart[mapstates[is][1]];

				if ((qPart[mapstates[is][0]] + qPart[mapstates[is][1]]) == 0) {
					// this is a neutral isobar state with 2 final state particles
					fillWeightedGJAnglePlots(lvIsobar, weight, weightPosRef, weightNegRef, weightFlat, tprime, itree, GJHB_neutral_isobar);
					fillWeightedHelicityAnglePlots(calculateHelicityAngles(lvIsobar, lvPart[mapstates[is][0]], NULL,  true), weight, itree, HHB_neutral_isobar);
					fillWeightedHelicityAnglePlots(calculateHelicityAngles(lvIsobar, lvPart[mapstates[is][1]], NULL, false), weight, itree, HHB_neutral_isobar);
				}
				else if ((qPart[mapstates[is][0]] + qPart[mapstates[is][1]]) == -1) {
					int neutralidx = 0;
					if (qPart[mapstates[is][1]] == 0) neutralidx = 1;

					// this is a negativly charged isobar state with 2 final state particles
					fillWeightedGJAnglePlots(lvIsobar, weight, weightPosRef, weightNegRef, weightFlat, tprime, itree, GJHB_charged_isobar);
					fillWeightedHelicityAnglePlots(calculateHelicityAngles(lvIsobar, lvPart[mapstates[is][neutralidx]]), weight, itree, HHB_charged_isobar);

					// fill the angles only if the isobar mass is close the the rho mass
					const rpwa::particleProperties* pp = rpwa::particleDataTable::entry("rho(770)-");
					if (pp == NULL) {
						printErr << "Unknown particle 'rho(770)-' which is required for rho mass cut on charged isobar." << std::endl;
						return;
					}

					if (std::abs(lvIsobar.M() - pp->mass()) < pp->width()/2.) {
						fillWeightedGJAnglePlots(lvIsobar, weight, weightPosRef, weightNegRef, weightFlat, tprime, itree, GJHB_rho);
						fillWeightedHelicityAnglePlots(calculateHelicityAngles(lvIsobar, lvPart[mapstates[is][neutralidx]]), weight, itree, HHB_rho);
					}
				}
			}

		}// end loop over events
		if (itree != 0) {
			avweight /= (double) nmbEvents;
			cout << "Maxweight=" << maxweight << endl;
			cout << "Average weight=" << avweight << endl;
		}
		std::cout << std::endl;
	}// end loop over trees

	GJHB_neutral_isobar.costheta_GJF_Stack[0]->Write(NULL, TObject::kOverwrite);
	GJHB_charged_isobar.costheta_GJF_Stack[0]->Write(NULL, TObject::kOverwrite);
	if (mcAccFile != NULL) {
		GJHB_neutral_isobar.costheta_GJF_Stack[1]->Write(NULL, TObject::kOverwrite);
		GJHB_charged_isobar.costheta_GJF_Stack[1]->Write(NULL, TObject::kOverwrite);
	}

	outfile->Write(NULL, TObject::kOverwrite);
	makeDifferencePlots(outdir);

	TList* Hlist = gDirectory->GetList();
	Hlist->Remove(dataTree);
	Hlist->Remove(mcPspTree);
	Hlist->Remove(mcAccTree);
	//Hlist->Remove("hWeights");
	int nobj = Hlist->GetEntries();
	std::cout << "Found " << nobj << " Objects in HList" << std::endl;

	outfile->Close();

	gROOT->cd();

	// print information on the processing time
	timer.Stop();
	printInfo << "this job consumed: ";
	timer.Print();
}

void
plotWeightedEvts_3pin(const std::string& dataFileName,
                      const std::string& mcFileName,
                      const std::string& mcWeightFileName,
                      const std::string& massBin,
                      const std::string& outFileName,
                      const std::string& pdgFileName      = std::string(getenv("ROOTPWA")) + "/amplitude/particleDataTable.txt")
{
	// set some default values for tree names and so on
	const std::string dataTreeName              = "rootPwaEvtTree";
	const std::string dataProdKinPartNamesName  = "prodKinParticles";
	const std::string dataProdKinMomentaName    = "prodKinMomenta";
	const std::string dataDecayKinPartNamesName = "decayKinParticles";
	const std::string dataDecayKinMomentaName   = "decayKinMomenta";
	const std::string mcTreeName                = "rootPwaEvtTree";
	const std::string mcProdKinPartNamesName    = "prodKinParticles";
	const std::string mcProdKinMomentaName      = "prodKinMomenta";
	const std::string mcDecayKinPartNamesName   = "decayKinParticles";
	const std::string mcDecayKinMomentaName     = "decayKinMomenta";
	const std::string mcWeightTreeName          = "rootPwaWeightTree";

	// set 25 MByte ROOT tree read cache
	const long int    treeCacheSize             = 25000000;

	createWeightedPlots(dataFileName, dataTreeName, dataProdKinPartNamesName, dataProdKinMomentaName, dataDecayKinPartNamesName, dataDecayKinMomentaName,
	                    mcFileName,   mcTreeName,   mcProdKinPartNamesName,   mcProdKinMomentaName,   mcDecayKinPartNamesName,   mcDecayKinMomentaName,
	                    mcWeightFileName, mcWeightTreeName,
	                    "",           "",           "",                       "",                     "",                        "",
	                    "",               "",
	                    massBin, outFileName,
	                    pdgFileName, treeCacheSize);
}

void
plotWeightedEvts_3pin(const std::string& dataFileName,
                      const std::string& mcPspFileName,
                      const std::string& mcPspWeightFileName,
                      const std::string& mcAccFileName,
                      const std::string& mcAccWeightFileName,
                      const std::string& massBin,
                      const std::string& outFileName,
                      const std::string& pdgFileName         = std::string(getenv("ROOTPWA")) + "/amplitude/particleDataTable.txt")
{
	// set some default values for tree names and so on
	const std::string dataTreeName               = "rootPwaEvtTree";
	const std::string dataProdKinPartNamesName   = "prodKinParticles";
	const std::string dataProdKinMomentaName     = "prodKinMomenta";
	const std::string dataDecayKinPartNamesName  = "decayKinParticles";
	const std::string dataDecayKinMomentaName    = "decayKinMomenta";
	const std::string mcPspTreeName              = "rootPwaEvtTree";
	const std::string mcPspProdKinPartNamesName  = "prodKinParticles";
	const std::string mcPspProdKinMomentaName    = "prodKinMomenta";
	const std::string mcPspDecayKinPartNamesName = "decayKinParticles";
	const std::string mcPspDecayKinMomentaName   = "decayKinMomenta";
	const std::string mcPspWeightTreeName        = "rootPwaWeightTree";
	const std::string mcAccTreeName              = "rootPwaEvtTree";
	const std::string mcAccProdKinPartNamesName  = "prodKinParticles";
	const std::string mcAccProdKinMomentaName    = "prodKinMomenta";
	const std::string mcAccDecayKinPartNamesName = "decayKinParticles";
	const std::string mcAccDecayKinMomentaName   = "decayKinMomenta";
	const std::string mcAccWeightTreeName        = "rootPwaWeightTree";

	// set 25 MByte ROOT tree read cache
	const long int    treeCacheSize             = 25000000;

	createWeightedPlots(dataFileName,  dataTreeName,  dataProdKinPartNamesName,  dataProdKinMomentaName,  dataDecayKinPartNamesName,  dataDecayKinMomentaName,
	                    mcPspFileName, mcPspTreeName, mcPspProdKinPartNamesName, mcPspProdKinMomentaName, mcPspDecayKinPartNamesName, mcPspDecayKinMomentaName,
	                    mcPspWeightFileName, mcPspWeightTreeName,
	                    mcAccFileName, mcAccTreeName, mcAccProdKinPartNamesName, mcAccProdKinMomentaName, mcAccDecayKinPartNamesName, mcAccDecayKinMomentaName,
	                    mcAccWeightFileName, mcAccWeightTreeName,
	                    massBin, outFileName,
	                    pdgFileName, treeCacheSize);
}
