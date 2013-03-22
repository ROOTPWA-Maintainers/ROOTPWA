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
#include <TLorentzVector.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TTree.h>

#include <NParticleEvent.h>
#include <reportingUtils.hpp>

using namespace std;

int nbninsm = 144;
int nbinsang = 80;



struct GJHistBunch {
	// base histograms
	std::vector<TH1D*> isobar_mass;
	std::vector<TH1D*> costheta_GJF;
	std::vector<TH1D*> phi_GJF;
	std::vector<TH1D*> costheta_GJF_MC_raw;
	std::vector<TH2D*> costheta_GJF_tprime;

	// resolved for positive and negative reflectivity, and flat wave
	std::vector<THStack*> costheta_GJF_Stack;
	std::vector<TH1D*> costheta_GJF_PosRef;
	std::vector<TH1D*> costheta_GJF_NegRef;
	std::vector<TH1D*> costheta_GJF_Flat;
	
	// difference histograms
	TH1D* isobar_mass_diff;
	TH1D* costheta_GJF_diff;
	TH1D* phi_GJF_diff;
	TH2D* costheta_GJF_tprime_diff;
};

struct HelicityHistBunch {
	// base histograms
	std::vector<TH1D*> costheta_HF;
	std::vector<TH1D*> phi_HF;
	
	// difference histograms
	TH1D* costheta_HF_diff;
	TH1D* phi_HF_diff;
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
	TH1D* hMIsobarData = new TH1D(("hMIsobarData_" + name_prefix).c_str(), (name_prefix + " Isobar Mass (Data)").c_str(), nbninsm,
				      0.0, 1.5);
	hMIsobarData->SetXTitle("isobar mass [GeV]");
	hMIsobarData->SetYTitle("# of events");
	temp.isobar_mass.push_back(hMIsobarData);
	if (twoMc) {
		TH1D* hMIsobarMcPsp = new TH1D(("hMIsobarMcPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Mass (McPsp)").c_str(), nbninsm, 0.0, 1.5);
		hMIsobarMcPsp->SetXTitle("isobar mass [GeV]");
		hMIsobarMcPsp->SetYTitle("# of events");
		temp.isobar_mass.push_back(hMIsobarMcPsp);
		TH1D* hMIsobarMcAcc = new TH1D(("hMIsobarMcAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Mass (McAcc)").c_str(), nbninsm, 0.0, 1.5);
		hMIsobarMcAcc->SetXTitle("isobar mass [GeV]");
		hMIsobarMcAcc->SetYTitle("# of events");
		temp.isobar_mass.push_back(hMIsobarMcAcc);
	} else {
		TH1D* hMIsobarMc = new TH1D(("hMIsobarMc_" + name_prefix).c_str(), (name_prefix + " Isobar Mass (Mc)").c_str(), nbninsm, 0.0, 1.5);
		hMIsobarMc->SetXTitle("isobar mass [GeV]");
		hMIsobarMc->SetYTitle("# of events");
		temp.isobar_mass.push_back(hMIsobarMc);
	}
	
	TH1D* hGJData = new TH1D(("hGJData_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (Data)").c_str(), nbinsang, -1, 1);
	hGJData->SetXTitle("isobar cos(#theta_{GJ})");
	hGJData->SetYTitle("# of events");
	temp.costheta_GJF.push_back(hGJData);
	if (twoMc) {
		TH1D* hGJMcPsp = new TH1D(("hGJMcPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (McPsp)").c_str(), nbinsang, -1, 1);
		hGJMcPsp->SetXTitle("isobar cos(#theta_{GJ})");
		hGJMcPsp->SetYTitle("# of events");
		temp.costheta_GJF.push_back(hGJMcPsp);
		TH1D* hGJMcAcc = new TH1D(("hGJMcAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (McAcc)").c_str(), nbinsang, -1, 1);
		hGJMcAcc->SetXTitle("isobar cos(#theta_{GJ})");
		hGJMcAcc->SetYTitle("# of events");
		temp.costheta_GJF.push_back(hGJMcAcc);
	} else {
		TH1D* hGJMc = new TH1D(("hGJMc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (Mc)").c_str(), nbinsang, -1, 1);
		hGJMc->SetXTitle("isobar cos(#theta_{GJ})");
		hGJMc->SetYTitle("# of events");
		temp.costheta_GJF.push_back(hGJMc);
	}

	if (twoMc) {
		TH1D* hGJMcPsp_raw = new TH1D(("hGJMcPsp_raw" + name_prefix).c_str(), "Cos Gottfried-Jackson Theta (unweighted McPsp)", nbinsang, -1, 1);
		hGJMcPsp_raw->SetXTitle("isobar cos(#theta_{GJ})");
		hGJMcPsp_raw->SetYTitle("# of events");
		temp.costheta_GJF_MC_raw.push_back(hGJMcPsp_raw);
		TH1D* hGJMcAcc_raw = new TH1D(("hGJMcAcc_raw" + name_prefix).c_str(), "Cos Gottfried-Jackson Theta (unweighted McAcc)", nbinsang, -1, 1);
		hGJMcAcc_raw->SetXTitle("isobar cos(#theta_{GJ})");
		hGJMcAcc_raw->SetYTitle("# of events");
		temp.costheta_GJF_MC_raw.push_back(hGJMcAcc_raw);
	} else {
		TH1D* hGJMc_raw = new TH1D(("hGJMc_raw" + name_prefix).c_str(),
					   "Cos Gottfried-Jackson Theta (unweighted Mc)", nbinsang, -1, 1);
		hGJMc_raw->SetXTitle("isobar cos(#theta_{GJ})");
		hGJMc_raw->SetYTitle("# of events");
		temp.costheta_GJF_MC_raw.push_back(hGJMc_raw);
	}
	
	TH2D* hGJtData = new TH2D(("hGJtData_" + name_prefix).c_str(), (name_prefix + " Isobar Cos GJ Theta vs t' (Data)").c_str(), nbinsang, -1, 1, 40, 0., 2.);
	hGJtData->SetXTitle("isobar cos(#theta_{GJ})");
	hGJtData->SetYTitle("t' [GeV]");
	hGJtData->SetOption("COLZ");
	temp.costheta_GJF_tprime.push_back(hGJtData);
	if (twoMc) {
		TH2D* hGJtMcPsp = new TH2D(("hGJtMcPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Cos GJ Theta vs t' (McPsp)").c_str(), nbinsang, -1, 1, 40, 0., 2.);
		hGJtMcPsp->SetXTitle("isobar cos(#theta_{GJ})");
		hGJtMcPsp->SetYTitle("t' [GeV]");
		hGJtMcPsp->SetOption("COLZ");
		temp.costheta_GJF_tprime.push_back(hGJtMcPsp);
		TH2D* hGJtMcAcc = new TH2D(("hGJtMcAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos GJ Theta vs t' (McAcc)").c_str(), nbinsang, -1, 1, 40, 0., 2.);
		hGJtMcAcc->SetXTitle("isobar cos(#theta_{GJ})");
		hGJtMcAcc->SetYTitle("t' [GeV]");
		hGJtMcAcc->SetOption("COLZ");
		temp.costheta_GJF_tprime.push_back(hGJtMcAcc);
	} else {
		TH2D* hGJtMc = new TH2D(("hGJtMc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos GJ Theta vs t' (Mc)").c_str(), nbinsang, -1, 1, 40, 0., 2.);
		hGJtMc->SetXTitle("isobar cos(#theta_{GJ})");
		hGJtMc->SetYTitle("t' [GeV]");
		hGJtMc->SetOption("COLZ");
		temp.costheta_GJF_tprime.push_back(hGJtMc);
	}
	
	TH1D* hTYData = new TH1D(("hTYData_" + name_prefix).c_str(), (name_prefix + " Isobar Treiman-Yang Phi (Data)").c_str(), nbinsang, -TMath::Pi(), TMath::Pi());
	hTYData->SetXTitle("isobar #phi_{TY} [rad]");
	hTYData->SetYTitle("# of events");
	temp.phi_GJF.push_back(hTYData);
	if (twoMc) {
		TH1D* hTYMcPsp = new TH1D(("hTYMcPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Treiman-Yang Phi (McPsp)").c_str(), nbinsang, -TMath::Pi(), TMath::Pi());
		hTYMcPsp->SetXTitle("isobar #phi_{TY} [rad]");
		hTYMcPsp->SetYTitle("# of events");
		temp.phi_GJF.push_back(hTYMcPsp);
		TH1D* hTYMcAcc = new TH1D(("hTYMcAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Treiman-Yang Phi (McAcc)").c_str(), nbinsang, -TMath::Pi(), TMath::Pi());
		hTYMcAcc->SetXTitle("isobar #phi_{TY} [rad]");
		hTYMcAcc->SetYTitle("# of events");
		temp.phi_GJF.push_back(hTYMcAcc);
	} else {
		TH1D* hTYMc = new TH1D(("hTYMc_" + name_prefix).c_str(), (name_prefix + " Isobar Treiman-Yang Phi (Mc)").c_str(), nbinsang, -TMath::Pi(), TMath::Pi());
		hTYMc->SetXTitle("isobar #phi_{TY} [rad]");
		hTYMc->SetYTitle("# of events");
		temp.phi_GJF.push_back(hTYMc);
	}

	if (twoMc) {
		THStack* hGJF_Stack_McPsp = new THStack(("hGJF_Stack_McPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (McPsp)").c_str());
	
		TH1D* hGJF_Flat_McPsp = new TH1D(("hGJF_Flat_McPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (McPsp)").c_str(), nbinsang, -1, 1);
		hGJF_Flat_McPsp->SetXTitle("isobar cos(#theta_{GJ})");
		hGJF_Flat_McPsp->SetYTitle("# of events");
		temp.costheta_GJF_Flat.push_back(hGJF_Flat_McPsp);
		hGJF_Stack_McPsp->Add(hGJF_Flat_McPsp);
		
		TH1D* hGJF_Neg_McPsp = new TH1D(("hGJF_Neg_McPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (McPsp)").c_str(), nbinsang, -1, 1);
		hGJF_Neg_McPsp->SetXTitle("isobar cos(#theta_{GJ})");
		hGJF_Neg_McPsp->SetYTitle("# of events");
		temp.costheta_GJF_NegRef.push_back(hGJF_Neg_McPsp);
		hGJF_Stack_McPsp->Add(hGJF_Neg_McPsp);
		
		TH1D* hGJF_Pos_McPsp = new TH1D(("hGJF_Pos_McPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (McPsp)").c_str(), nbinsang, -1, 1);
		hGJF_Pos_McPsp->SetXTitle("isobar cos(#theta_{GJ})");
		hGJF_Pos_McPsp->SetYTitle("# of events");
		temp.costheta_GJF_PosRef.push_back(hGJF_Pos_McPsp);
		hGJF_Stack_McPsp->Add(hGJF_Pos_McPsp);
	
		temp.costheta_GJF_Stack.push_back(hGJF_Stack_McPsp);

		THStack* hGJF_Stack_McAcc = new THStack(("hGJF_Stack_McAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (McAcc)").c_str());
	
		TH1D* hGJF_Flat_McAcc = new TH1D(("hGJF_Flat_McAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (McAcc)").c_str(), nbinsang, -1, 1);
		hGJF_Flat_McAcc->SetXTitle("isobar cos(#theta_{GJ})");
		hGJF_Flat_McAcc->SetYTitle("# of events");
		temp.costheta_GJF_Flat.push_back(hGJF_Flat_McAcc);
		hGJF_Stack_McAcc->Add(hGJF_Flat_McAcc);
		
		TH1D* hGJF_Neg_McAcc = new TH1D(("hGJF_Neg_McAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (McAcc)").c_str(), nbinsang, -1, 1);
		hGJF_Neg_McAcc->SetXTitle("isobar cos(#theta_{GJ})");
		hGJF_Neg_McAcc->SetYTitle("# of events");
		temp.costheta_GJF_NegRef.push_back(hGJF_Neg_McAcc);
		hGJF_Stack_McAcc->Add(hGJF_Neg_McAcc);
		
		TH1D* hGJF_Pos_McAcc = new TH1D(("hGJF_Pos_McAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (McAcc)").c_str(), nbinsang, -1, 1);
		hGJF_Pos_McAcc->SetXTitle("isobar cos(#theta_{GJ})");
		hGJF_Pos_McAcc->SetYTitle("# of events");
		temp.costheta_GJF_PosRef.push_back(hGJF_Pos_McAcc);
		hGJF_Stack_McAcc->Add(hGJF_Pos_McAcc);
	
		temp.costheta_GJF_Stack.push_back(hGJF_Stack_McAcc);
	} else {
		THStack* hGJF_Stack_Mc = new THStack(("hGJF_Stack_Mc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (Mc)").c_str());
	
		TH1D* hGJF_Flat_Mc = new TH1D(("hGJF_Flat_Mc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (Mc)").c_str(), nbinsang, -1, 1);
		hGJF_Flat_Mc->SetXTitle("isobar cos(#theta_{GJ})");
		hGJF_Flat_Mc->SetYTitle("# of events");
		temp.costheta_GJF_Flat.push_back(hGJF_Flat_Mc);
		hGJF_Stack_Mc->Add(hGJF_Flat_Mc);
		
		TH1D* hGJF_Neg_Mc = new TH1D(("hGJF_Neg_Mc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (Mc)").c_str(), nbinsang, -1, 1);
		hGJF_Neg_Mc->SetXTitle("isobar cos(#theta_{GJ})");
		hGJF_Neg_Mc->SetYTitle("# of events");
		temp.costheta_GJF_NegRef.push_back(hGJF_Neg_Mc);
		hGJF_Stack_Mc->Add(hGJF_Neg_Mc);
		
		TH1D* hGJF_Pos_Mc = new TH1D(("hGJF_Pos_Mc_" + name_prefix).c_str(), (name_prefix + " Isobar Cos Gottfried-Jackson Theta (Mc)").c_str(), nbinsang, -1, 1);
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
	TH1D* hHThetaData = new TH1D(("hHThetaData_" + name_prefix).c_str(), (name_prefix + " Isobar Helicity Cos Theta (Data)").c_str(), nbinsang, -1,
				     1);
	hHThetaData->SetXTitle("cos(#theta_{hel}) of #pi^{0} from isobar");
	hHThetaData->SetYTitle("# of events");
	temp.costheta_HF.push_back(hHThetaData);
	if (twoMc) {
		TH1D* hHThetaMcPsp = new TH1D(("hHThetaMcPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Helicity Cos Theta (McPsp)").c_str(), nbinsang, -1, 1);
		hHThetaMcPsp->SetXTitle("cos(#theta_{hel} of #pi^{0} from isobar)");
		hHThetaMcPsp->SetYTitle("# of events");
		temp.costheta_HF.push_back(hHThetaMcPsp);
		TH1D* hHThetaMcAcc = new TH1D(("hHThetaMcAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Helicity Cos Theta (McAcc)").c_str(), nbinsang, -1, 1);
		hHThetaMcAcc->SetXTitle("cos(#theta_{hel} of #pi^{0} from isobar)");
		hHThetaMcAcc->SetYTitle("# of events");
		temp.costheta_HF.push_back(hHThetaMcAcc);
	} else {
		TH1D* hHThetaMc = new TH1D(("hHThetaMc_" + name_prefix).c_str(), (name_prefix + " Isobar Helicity Cos Theta (Mc)").c_str(), nbinsang, -1, 1);
		hHThetaMc->SetXTitle("cos(#theta_{hel} of #pi^{0} from isobar)");
		hHThetaMc->SetYTitle("# of events");
		temp.costheta_HF.push_back(hHThetaMc);
	}
	
	TH1D* hHPhiData = new TH1D(("hHPhiData_" + name_prefix).c_str(), (name_prefix + " Isobar Helicity Phi (Data)").c_str(), nbinsang, -TMath::Pi(), TMath::Pi());
	hHPhiData->SetXTitle("#phi_{hel} [rad] of #pi^{0} from isobar");
	hHPhiData->SetYTitle("# of events");
	temp.phi_HF.push_back(hHPhiData);
	if (twoMc) {
		TH1D* hHPhiMcPsp = new TH1D(("hHPhiMcPsp_" + name_prefix).c_str(), (name_prefix + " Isobar Helicity Phi (McPsp)").c_str(), nbinsang, -TMath::Pi(), TMath::Pi());
		hHPhiMcPsp->SetXTitle("#phi_{hel} [rad] of #pi^{0} from isobar");
		hHPhiMcPsp->SetYTitle("# of events");
		temp.phi_HF.push_back(hHPhiMcPsp);
		TH1D* hHPhiMcAcc = new TH1D(("hHPhiMcAcc_" + name_prefix).c_str(), (name_prefix + " Isobar Helicity Phi (McAcc)").c_str(), nbinsang, -TMath::Pi(), TMath::Pi());
		hHPhiMcAcc->SetXTitle("#phi_{hel} [rad] of #pi^{0} from isobar");
		hHPhiMcAcc->SetYTitle("# of events");
		temp.phi_HF.push_back(hHPhiMcAcc);
	} else {
		TH1D* hHPhiMc = new TH1D(("hHPhiMc_" + name_prefix).c_str(), (name_prefix + " Isobar Helicity Phi (Mc)").c_str(), nbinsang, -TMath::Pi(), TMath::Pi());
		hHPhiMc->SetXTitle("#phi_{hel} [rad] of #pi^{0} from isobar");
		hHPhiMc->SetYTitle("# of events");
		temp.phi_HF.push_back(hHPhiMc);
	}
	
	return temp;
}

void fillWeightedHelicityAnglePlots(const HelicityAngles &ha, double weight, unsigned int tree_index, HelicityHistBunch &hhb) {
	hhb.costheta_HF[tree_index]->Fill(ha.cosTheta, weight);
	hhb.phi_HF[tree_index]->Fill(ha.phi, weight);
}

void fillWeightedGJAnglePlots(const TLorentzVector &isobar, double weight, double weightPosRef, double weightNegRef, double weightFlat, double tprime, unsigned int tree_index, GJHistBunch &hBunch) {
	hBunch.costheta_GJF[tree_index]->Fill(isobar.CosTheta(), weight);
	if (tree_index != 0) {
		hBunch.costheta_GJF_MC_raw[tree_index-1]->Fill(isobar.CosTheta());
		hBunch.costheta_GJF_PosRef[tree_index-1]->Fill(isobar.CosTheta(), weightPosRef);
		hBunch.costheta_GJF_NegRef[tree_index-1]->Fill(isobar.CosTheta(), weightNegRef);
		hBunch.costheta_GJF_Flat[tree_index-1]->Fill(isobar.CosTheta(), weightFlat);
	}
	hBunch.costheta_GJF_tprime[tree_index]->Fill(isobar.CosTheta(), tprime, weight);
	hBunch.phi_GJF[tree_index]->Fill(isobar.Phi(), weight);
	hBunch.isobar_mass[tree_index]->Fill(isobar.M(), weight);
}

HelicityAngles calculateHelicityAngles(const NParticleState &isobar, TLorentzVector *beam = NULL, bool first = true) {
	HelicityAngles temp;
	
	TVector3 zaxis_gjf;
	if(beam != NULL) {
		zaxis_gjf = beam->Vect().Unit();
	}
	else {
		zaxis_gjf = TVector3(0,0,1);
	}
	// create helicity frame coordinate system
	TVector3 zaxis = isobar.p().Vect().Unit();
	TVector3 yaxis = zaxis_gjf.Cross(zaxis).Unit();
	TVector3 xaxis = yaxis.Cross(zaxis);
	
	// boost NParticleState into isobar rest frame
	TLorentzVector particle;
	if(isobar.getParticle(0)->q() == 0 && first)
		particle = isobar.getParticle(0)->p();
	else
		particle = isobar.getParticle(1)->p();
	const TVector3 boost_vec = isobar.p().BoostVector();
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
	// make list of MC Histograms
	TList mclist;
	TList *histlist = dir->GetListOfKeys();
	TIter histiter(histlist);
	TObject *obj;
	while ((obj = histiter())) {
		const std::string s(obj->GetName());
		if ((s.length() >= 2 && s.substr(s.length()-2, 2) == "Mc") ||
		    (s.length() >= 5 && s.substr(s.length()-5, 5) == "McPsp") ||
		    (s.length() >= 5 && s.substr(s.length()-5, 5) == "McAcc") ||
		    s.find("Mc_") != std::string::npos ||
		    s.find("McPsp_") != std::string::npos ||
		    s.find("McAcc_") != std::string::npos) {
			mclist.Add(obj);
		}
	}
	histiter = TIter(&mclist);
	TH1D *diffhist, *reldiffhist, *mchist, *datahist;
	double scale;
	while ((mchist = (TH1D*) histiter())) {
		// generate difference histograms
		std::string hnamemc(mchist->GetName());
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
		
		dir->GetObject(hnamediff.c_str(), diffhist);
		dir->GetObject(hnamereldiff.c_str(), reldiffhist);
		dir->GetObject(hnamedata.c_str(), datahist);
		dir->GetObject(hnamemc.c_str(), mchist);
		if (!diffhist) {
			if (datahist) {
				dir->cd();
				scale = datahist->Integral();
				scale = scale / (mchist->Integral());
				mchist->Scale(scale);
				
				diffhist = new TH1D(*mchist);
				diffhist->SetName(hnamediff.c_str());
				diffhist->SetTitle("");
				diffhist->Add(datahist, -1.);
				diffhist->SetYTitle("# of events difference(MC-Data)");
				diffhist->Write();
				
				reldiffhist = new TH1D(*diffhist);
				reldiffhist->SetName(hnamereldiff.c_str());
				reldiffhist->Divide(datahist);
				reldiffhist->SetYTitle("relative difference((MC-Data)/Data)");
				reldiffhist->Write();
			}
		}
	}
	
	histiter = TIter(&mclist);
	TH2D *diffhist2d, *reldiffhist2d, *mchist2d, *datahist2d;
	scale = 1.0;
	while ((mchist2d = (TH2D*) histiter())) {
		// generate difference histograms
		std::string hnamemc(mchist2d->GetName());
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
		
		dir->GetObject(hnamediff.c_str(), diffhist2d);
		dir->GetObject(hnamereldiff.c_str(), reldiffhist2d);
		dir->GetObject(hnamedata.c_str(), datahist2d);
		dir->GetObject(hnamemc.c_str(), mchist2d);
		if (!diffhist2d) {
			if (datahist2d) {
				dir->cd();
				scale = datahist2d->Integral();
				scale = scale / (mchist2d->Integral());
				mchist2d->Scale(scale);
				
				diffhist2d = new TH2D(*mchist2d);
				diffhist2d->SetName(hnamediff.c_str());
				diffhist2d->SetTitle("diff(MC-Data)");
				diffhist2d->Add(datahist2d, -1.);
				double max = diffhist2d->GetMaximum();
				if(max < TMath::Abs(diffhist2d->GetMinimum()))
					max = TMath::Abs(diffhist2d->GetMinimum());
				diffhist2d->SetMaximum(max);
				diffhist2d->SetMinimum(-max);
				diffhist2d->Write();
				
				reldiffhist2d = new TH2D(*diffhist2d);
				reldiffhist2d->SetName(hnamereldiff.c_str());
				reldiffhist2d->SetTitle("rel. diff((MC-Data)/Data)");
				reldiffhist2d->Divide(datahist2d);
				reldiffhist2d->SetMaximum();
				reldiffhist2d->SetMinimum();
				max = reldiffhist2d->GetMaximum();
				if (max < TMath::Abs(reldiffhist2d->GetMinimum()))
					max = TMath::Abs(reldiffhist2d->GetMinimum());
				reldiffhist2d->SetMaximum(max);
				reldiffhist2d->SetMinimum(-max);
				reldiffhist2d->Write();
			}
		}
	}
}

TH2D* createDalitzHistogram(const std::string& name, const std::string& title, double mass, unsigned int treeentries) {
	double rangelow = 0.0;
	double rangehigh = pow(mass-0.1, 2);
	unsigned int nbins = 0;
	nbins = 0.6*TMath::Sqrt(treeentries)*(rangehigh-rangelow);
	std::cout<<"nbins for mass "<<mass<<": "<<nbins<<std::endl;
	// lower bound is 20 bins
	if(nbins < 20) nbins = 20;
	// upper bound is 130 bins
	if(nbins > 130) nbins = 130;
	TH2D* phist = new TH2D(name.c_str(), title.c_str(), nbins, rangelow, rangehigh, nbins, rangelow, rangehigh);
	return phist;
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
		    const std::string& mcAccFileName,
		    const std::string& mcAccTreeName,
		    const std::string& mcAccProdKinPartNamesName,
		    const std::string& mcAccProdKinMomentaName,
		    const std::string& mcAccDecayKinPartNamesName,
		    const std::string& mcAccDecayKinMomentaName,
		    const std::string& massBin,
		    const std::string& outFileName,
		    const long int     treeCacheSize)
{
	// keep track of the processing time
	TStopwatch timer;
	timer.Start();

	// open data file
	TFile* dataFile = TFile::Open(dataFileName.c_str());

	// tree containg the real data events
	TTree* dataTree;
	dataFile->GetObject(dataTreeName.c_str(), dataTree);

	// names of particles in tree
	TClonesArray* dataProdKinPartNames(NULL);
	dataFile->GetObject(dataProdKinPartNamesName.c_str(), dataProdKinPartNames);
	assert(dataProdKinPartNames->GetEntries() == 1);
	for (Int_t i=0; i<dataProdKinPartNames->GetEntries(); i++)
		assert(((TObjString*)(dataProdKinPartNames->At(i)))->String().EqualTo("pi-"));

	TClonesArray* dataDecayKinPartNames(NULL);
	dataFile->GetObject(dataDecayKinPartNamesName.c_str(), dataDecayKinPartNames);
	assert(dataDecayKinPartNames->GetEntries() == 3);
	for (Int_t i=0; i<dataDecayKinPartNames->GetEntries(); i++)
		assert(((TObjString*)(dataDecayKinPartNames->At(i)))->String().EqualTo("pi-") || ((TObjString*)(dataDecayKinPartNames->At(i)))->String().EqualTo("pi0"));

	// open weighted MC file
	TFile* mcPspFile = TFile::Open(mcPspFileName.c_str());

	// tree containing the phase space events and weights
	TTree* mcPspTree;
	mcPspFile->GetObject(mcPspTreeName.c_str(), mcPspTree);

	// names of particles in tree
	TClonesArray* mcPspProdKinPartNames(NULL);
	mcPspFile->GetObject(mcPspProdKinPartNamesName.c_str(), mcPspProdKinPartNames);
	assert(mcPspProdKinPartNames->GetEntries() == 1);
	for (Int_t i=0; i<mcPspProdKinPartNames->GetEntries(); i++)
		assert(((TObjString*)(mcPspProdKinPartNames->At(i)))->String().EqualTo("pi-"));

	TClonesArray* mcPspDecayKinPartNames(NULL);
	mcPspFile->GetObject(mcPspDecayKinPartNamesName.c_str(), mcPspDecayKinPartNames);
	assert(mcPspDecayKinPartNames->GetEntries() == 3);
	for (Int_t i=0; i<mcPspDecayKinPartNames->GetEntries(); i++)
		assert(((TObjString*)(mcPspDecayKinPartNames->At(i)))->String().EqualTo("pi-") || ((TObjString*)(mcPspDecayKinPartNames->At(i)))->String().EqualTo("pi0"));

	// open weighted MC file
	TFile* mcAccFile(NULL);
	if (mcAccFileName != "") {
		mcAccFile = TFile::Open(mcAccFileName.c_str());
	}

	// tree containing the phase space events and weights
	TTree* mcAccTree(NULL);
	if (mcAccFile != NULL) {
		mcAccFile->GetObject(mcAccTreeName.c_str(), mcAccTree);
	}

	// names of particles in tree
	TClonesArray* mcAccProdKinPartNames(NULL);
	if (mcAccFile != NULL) {
		mcAccFile->GetObject(mcAccProdKinPartNamesName.c_str(), mcAccProdKinPartNames);
		assert(mcAccProdKinPartNames->GetEntries() == 1);
		for (Int_t i=0; i<mcAccProdKinPartNames->GetEntries(); i++)
			assert(((TObjString*)(mcAccProdKinPartNames->At(i)))->String().EqualTo("pi-"));
	}

	TClonesArray* mcAccDecayKinPartNames(NULL);
	if (mcAccFile != NULL) {
		mcAccFile->GetObject(mcAccDecayKinPartNamesName.c_str(), mcAccDecayKinPartNames);
		assert(mcAccDecayKinPartNames->GetEntries() == 3);
		for (Int_t i=0; i<mcAccDecayKinPartNames->GetEntries(); i++)
			assert(((TObjString*)(mcAccDecayKinPartNames->At(i)))->String().EqualTo("pi-") || ((TObjString*)(mcAccDecayKinPartNames->At(i)))->String().EqualTo("pi0"));
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
	TH1D* hMData = new TH1D("hResMassData", "Mass (Data)", 60, 0.5, 2.5);
	hM.push_back(hMData);
	if (mcAccFile != NULL) {
		TH1D* hMMcPsp = new TH1D("hResMassMcPsp", "Mass (McPsp)", 60, 0.5, 2.5);
		hM.push_back(hMMcPsp);
		TH1D* hMMcAcc = new TH1D("hResMassMcAcc", "Mass (McAcc)", 60, 0.5, 2.5);
		hM.push_back(hMMcAcc);
	} else {
		TH1D* hMMc = new TH1D("hResMassMc", "Mass (Mc)", 60, 0.5, 2.5);
		hM.push_back(hMMc);
	}
	
	vector<TH1D*> hThetaLab;
	TH1D* hThetaLabData = new TH1D("hThetaLabData", "Cos Theta Lab (Data)", nbninsm, 0.997, 1);
	hThetaLab.push_back(hThetaLabData);
	if (mcAccFile != NULL) {
		TH1D* hThetaLabMcPsp = new TH1D("hThetaLabMcPsp", "Cos Theta Lab (McPsp)", nbninsm, 0.997, 1);
		hThetaLab.push_back(hThetaLabMcPsp);
		TH1D* hThetaLabMcAcc = new TH1D("hThetaLabMcAcc", "Cos Theta Lab (McAcc)", nbninsm, 0.997, 1);
		hThetaLab.push_back(hThetaLabMcAcc);
	} else {
		TH1D* hThetaLabMc = new TH1D("hThetaLabMc", "Cos Theta Lab (Mc)", nbninsm, 0.997, 1);
		hThetaLab.push_back(hThetaLabMc);
	}
	
	// Dalitz plots
	std::vector<TH2D*> dalitz_neutral;
	TH2D* hDalitzData = createDalitzHistogram("hDalitzData",
						  "Dalitz Plot #pi^{0}#pi^{0} vs. #pi^{-}#pi^{0} (Data)", massval, datatreeentries);
	hDalitzData->SetXTitle("mass^{2}(#pi^{0}#pi^{0}) [GeV^{2}/c^{4}]");
	hDalitzData->SetYTitle("mass^{2}(#pi^{-}#pi^{0}) [GeV^{2}/c^{4}]");
	hDalitzData->SetOption("COLZ");
	hDalitzData->SetStats(0);
	dalitz_neutral.push_back(hDalitzData);
	if (mcAccFile != NULL) {
		TH2D* hDalitzMcPsp = createDalitzHistogram("hDalitzMcPsp",
							   "Dalitz Plot #pi^{0}#pi^{0} vs. #pi^{-}#pi^{0} (McPsp)", massval, datatreeentries);
		hDalitzMcPsp->SetXTitle("mass^{2}(#pi^{0}#pi^{0}) [GeV^{2}/c^{4}]");
		hDalitzMcPsp->SetYTitle("mass^{2}(#pi^{-}#pi^{0}) [GeV^{2}/c^{4}]");
		hDalitzMcPsp->SetOption("COLZ");
		hDalitzMcPsp->SetStats(0);
		dalitz_neutral.push_back(hDalitzMcPsp);
		TH2D* hDalitzMcAcc = createDalitzHistogram("hDalitzMcAcc",
							   "Dalitz Plot #pi^{0}#pi^{0} vs. #pi^{-}#pi^{0} (McAcc)", massval, datatreeentries);
		hDalitzMcAcc->SetXTitle("mass^{2}(#pi^{0}#pi^{0}) [GeV^{2}/c^{4}]");
		hDalitzMcAcc->SetYTitle("mass^{2}(#pi^{-}#pi^{0}) [GeV^{2}/c^{4}]");
		hDalitzMcAcc->SetOption("COLZ");
		hDalitzMcAcc->SetStats(0);
		dalitz_neutral.push_back(hDalitzMcAcc);
	} else {
		TH2D* hDalitzMc = createDalitzHistogram("hDalitzMc",
							"Dalitz Plot #pi^{0}#pi^{0} vs. #pi^{-}#pi^{0} (Mc)", massval, datatreeentries);
		hDalitzMc->SetXTitle("mass^{2}(#pi^{0}#pi^{0}) [GeV^{2}/c^{4}]");
		hDalitzMc->SetYTitle("mass^{2}(#pi^{-}#pi^{0}) [GeV^{2}/c^{4}]");
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
	
	double avweight = 1;
	
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
		double maxweight = 0;
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

		TLorentzVector* beam = new TLorentzVector;
		int qbeam;
		TClonesArray* p = new TClonesArray("TLorentzVector");
		std::vector<int>* q = new std::vector<int>;

		TVector3 vertex;

		NParticleEvent event(p, q, beam, &qbeam, &vertex);

		// loop over tree entries
		const long int nmbEvents = tree->GetEntries();
		boost::progress_display progressIndicator(nmbEvents, cout, "");
		for (long int i = 0; i < nmbEvents; ++i) {
			++progressIndicator;
			tree->GetEntry(i);
			p->Delete();
			q->clear();

			assert(prodKinMomenta->GetEntries() == 1);
			for (Int_t j=0; j<prodKinMomenta->GetEntries(); j++) {
				if (((TObjString*)(prodKinPartNames->At(j)))->String().EqualTo("pi-")) {
					beam->SetVectM(*((TVector3*)(prodKinMomenta->At(j))), 0.13957018);
					qbeam = -1;
				} else
					assert(false);
			}
			assert(decayKinMomenta->GetEntries() == 3);
			for (Int_t j=0; j<decayKinMomenta->GetEntries(); j++) {
				new((*p)[j]) TLorentzVector;
				if (((TObjString*)(decayKinPartNames->At(j)))->String().EqualTo("pi-")) {
					((TLorentzVector*)(p->At(j)))->SetVectM(*((TVector3*)(decayKinMomenta->At(j))), 0.13957018);
					q->push_back(-1);
				} else if (((TObjString*)(decayKinPartNames->At(j)))->String().EqualTo("pi0")) {
					((TLorentzVector*)(p->At(j)))->SetVectM(*((TVector3*)(decayKinMomenta->At(j))), 0.1349766);
					q->push_back(0);
				} else
					assert(false);
			}
			// in case its data tree (itree=0) put weights to 1
			if (itree == 0) {
				weight = 1;
				impweight = 1;
			}
			if (weight == 0)
				continue;

			// this builds all subsystems
			event.refresh();
			
			if (impweight != 0)
				weight /= impweight;
			
			if (weight > maxweight)
				maxweight = weight;
			if (itree != 0)
				avweight += weight;

			double tprime = event.tprime();
			//cout << tprime << endl;
			unsigned int npart = event.nParticles();
			// loop over pions
			for (unsigned int ipart = 0; ipart < npart; ++ipart) {
				TLorentzVector p = event.getParticle(ipart).p();
				hThetaLab[itree]->Fill(p.CosTheta(), weight);
			}
			
			// ----dalitz plots
			// make all combinations
			std::vector<std::pair<std::pair<int, int>, double> > comb;
			for(unsigned int i=0; i < npart; i++) {
				for(unsigned int j=0; j < i; j++) {
					TLorentzVector lv = event.getParticle(i).p() + event.getParticle(j).p();
					double mass2 = lv.M2();
					std::pair<int, int> charge(event.getParticle(i).q(), event.getParticle(j).q());
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
			event.toGJ();
			// build again to get particles in GJF
			event.build();
			
			// loop over all states that contain n-1 final state particles
			// and plot angles
			unsigned int nstates = event.nStates();
			//cout<<"number of substates: "<<event.nStates()<<endl;
			for (unsigned int is = 0; is < nstates; ++is) {
				
				const NParticleState& state = event.getState(is);
				if (state.n() == npart) {
					hM[itree]->Fill(state.p().M());
				}
				if (state.n() == npart - 1 && state.q() == 0) {
					// this is a neutral isobar state with n-1 (here 2) final state particles
					fillWeightedGJAnglePlots(state.p(), weight, weightPosRef, weightNegRef, weightFlat, tprime, itree, GJHB_neutral_isobar);
					fillWeightedHelicityAnglePlots(calculateHelicityAngles(state, NULL,  true), weight, itree, HHB_neutral_isobar);
					fillWeightedHelicityAnglePlots(calculateHelicityAngles(state, NULL, false), weight, itree, HHB_neutral_isobar);
				}
				else if (state.n() == npart - 1 && state.q() == -1) {
					// this is a negativly charged isobar state with n-1 (here 2) final state particles
					fillWeightedGJAnglePlots(state.p(), weight, weightPosRef, weightNegRef, weightFlat, tprime, itree, GJHB_charged_isobar);
					fillWeightedHelicityAnglePlots(calculateHelicityAngles(state), weight, itree, HHB_charged_isobar);
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
	
	GJHB_neutral_isobar.costheta_GJF_Stack[0]->Write();
	GJHB_charged_isobar.costheta_GJF_Stack[0]->Write();
	if (mcAccFile != NULL) {
		GJHB_neutral_isobar.costheta_GJF_Stack[1]->Write();
		GJHB_charged_isobar.costheta_GJF_Stack[1]->Write();
	}
	
	outfile->Write();
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
                      const std::string& massBin,
                      const std::string& outFileName)
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

	// set 25 MByte ROOT tree read cache
	const long int    treeCacheSize             = 25000000;

	createWeightedPlots(dataFileName, dataTreeName, dataProdKinPartNamesName, dataProdKinMomentaName, dataDecayKinPartNamesName, dataDecayKinMomentaName,
	                    mcFileName,   mcTreeName,   mcProdKinPartNamesName,   mcProdKinMomentaName,   mcDecayKinPartNamesName,   mcDecayKinMomentaName,
			    "",           "",           "",                       "",                     "",                        "",
	                    massBin, outFileName,
			    treeCacheSize);
}

void
plotWeightedEvts_3pin(const std::string& dataFileName,
                      const std::string& mcPspFileName,
                      const std::string& mcAccFileName,
                      const std::string& massBin,
                      const std::string& outFileName)
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
	const std::string mcAccTreeName              = "rootPwaEvtTree";
	const std::string mcAccProdKinPartNamesName  = "prodKinParticles";
	const std::string mcAccProdKinMomentaName    = "prodKinMomenta";
	const std::string mcAccDecayKinPartNamesName = "decayKinParticles";
	const std::string mcAccDecayKinMomentaName   = "decayKinMomenta";

	// set 25 MByte ROOT tree read cache
	const long int    treeCacheSize             = 25000000;

	createWeightedPlots(dataFileName,  dataTreeName,  dataProdKinPartNamesName,  dataProdKinMomentaName,  dataDecayKinPartNamesName,  dataDecayKinMomentaName,
	                    mcPspFileName, mcPspTreeName, mcPspProdKinPartNamesName, mcPspProdKinMomentaName, mcPspDecayKinPartNamesName, mcPspDecayKinMomentaName,
	                    mcAccFileName, mcAccTreeName, mcAccProdKinPartNamesName, mcAccProdKinMomentaName, mcAccDecayKinPartNamesName, mcAccDecayKinMomentaName,
	                    massBin, outFileName,
			    treeCacheSize);
}
