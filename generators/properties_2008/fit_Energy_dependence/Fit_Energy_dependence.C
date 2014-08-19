// reading BeamPar_70507 in order to retrieve the depencence of the energy on the vertex position
// the histograms will be fitted and the fitted values will be filled into 2D histograms
// depending on the vertex position
//
// to do: for each x,y bin set one has produce the energy distribution fit it and
// write the mean value and sigma to a TH2(x,y)
//
// do > root -l fit_Energy_dependence.C+ // to run it
// by P. Jasinski Promme@web.de jasinski@kph.uni-mainz.de
// BeamPar_70507.root was produced by Hauke Wöhrmann Hauke.Woehrmann@physik.uni-muenchen.de
// big thanx for that!


#include <iostream>
#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TTree.h"
#include <string>

using namespace std;

const string filename("BeamPar_70507.root");

void Fit_Energy_dependence(){
	// set some style options
	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	// try to read the file containing the energy distributions
	TFile* beampar = new TFile("BeamPar_70507.root");
	if (beampar->IsZombie()){
		cout << " Error: could not open file! " << endl;
		return;
	}
	// retrieve the histograms
	TH2* hEVsDX = (TH2*) gDirectory->Get("hEVsDX");
	TH2* hEVsDY = (TH2*) gDirectory->Get("hEVsDY");
	TH2* hEVsX = (TH2*) gDirectory->Get("hEVsX");
	TH2* hEVsY = (TH2*) gDirectory->Get("hEVsY");
	TH2* hdxdyB = (TH2*) gDirectory->Get("hdxdyB");
	TH2* hxdxB = (TH2*) gDirectory->Get("hxdxB");
	TH2* hxdyB = (TH2*) gDirectory->Get("hxdyB");
	TH2* hydxB = (TH2*) gDirectory->Get("hydxB");
	TH2* hydyB = (TH2*) gDirectory->Get("hydyB");
	TH2* hxyB = (TH2*) gDirectory->Get("hxyB");
	TTree* BpT = (TTree*) gDirectory->Get("BpT");

	if (!hEVsDX || !hEVsDY || !hEVsX || !hEVsY || !hdxdyB ||
			!hxdxB || !hxdyB || !hydxB || !hydyB || !hxyB || !BpT ){
		cout << " Error: At least one histogram could not be retrieved " << endl;
		return;
	}
	TCanvas* canvas = new TCanvas ("canvas", "loaded histograms", 800, 800);
	canvas->Divide(3,3);
	canvas->cd(1);
	hEVsDX->Draw("COLZ");
	canvas->cd(2);
	hEVsDY->Draw("COLZ");
	canvas->cd(3);
	hEVsX->Draw("COLZ");
	canvas->cd(4);
	hEVsY->Draw("COLZ");
	canvas->cd(5);
	hdxdyB->Draw("COLZ");
	canvas->cd(6);
	hxdxB->Draw("COLZ");
	canvas->cd(7);
	hxdyB->Draw("COLZ");
	canvas->cd(8);
	hydxB->Draw("COLZ");
	canvas->cd(9);
	hydyB->Draw("COLZ");
	canvas->Update();
	canvas->Print("Loaded_histograms.pdf");

	TFile* outputfile = new TFile("energy_distributions.root", "Recreate");

	// this part needs to get implemented...
	// for (int x = 0; x < 21; x++){
	// 		for (int y = 0; y < 21; y++){
	//			BpT->Draw("EBeamApprox", " .. < x && x < .. & .. < y && y < .. ");
	// 			// get, fit, save
	//			see void UserEvent_basics_root_scripts::Analyze_beam_properties_output()
	// 			to get a clue how to do it ...

	outputfile->Write();
	outputfile->Close();
}
