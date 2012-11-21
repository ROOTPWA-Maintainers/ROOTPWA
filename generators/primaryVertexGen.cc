/*
 * author: Prometeusz (Promme) jasinski
 *         jasinski@kph.uni-mainz.de Promme@web.de
 *
 * This file contains some methods simulating the Primary Vertex distribution
 * in the target cell as well as the incoming beam properties
 *
 * or use the methods here to adapt it in your event generator
 *
 * (2010.03.03)
 * - implementing 2008 h- beam simulation
 *
 * (2010.06.16)
 * - moved to rootpwa svn repository and changed to class structure
 * - changed to a class
 *
 */

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include "TColor.h"
#include "TROOT.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom.h"

#include "randomNumberGenerator.h"
#include "reportingUtils.hpp"
#include "primaryVertexGen.h"


using namespace std;
using namespace rpwa;

primaryVertexGen::primaryVertexGen(string histfilename,
                                   double beam_part_mass,
                                   double mean_beam_energy,
                                   double mean_beam_energy_spread)
	: _histogramfile(NULL),
	  _hist_angles_vert_mean(NULL),
	  _hist_angles_horiz_mean(NULL),
	  _hist_angles_vert_sigma(NULL),
	  _hist_angles_horiz_sigma(NULL),
	  _hist_vertex_distr_xy(NULL),
	  _hist_vertex_distr_z(NULL),
	  _beam_part_mass(beam_part_mass),
	  _beam_energy_mean(mean_beam_energy),
	  _beam_energy_sigma(mean_beam_energy_spread)
{
	_histograms_loaded = loadHistograms(histfilename, true);
}


primaryVertexGen::~primaryVertexGen() {
	if(_histogramfile) {
		_histogramfile->Close();
	}
	delete _histogramfile;
	delete _hist_angles_vert_mean;
	delete _hist_angles_horiz_mean;
	delete _hist_angles_vert_sigma;
	delete _hist_angles_horiz_sigma;
	delete _hist_vertex_distr_xy;
	delete _hist_vertex_distr_z;
}


bool primaryVertexGen::loadHistograms(string filename, bool plot) {
	if(plot) {
		gROOT->SetStyle("Plain");
		gStyle->SetPalette(1);
		gesPalette();
	}
	TFile* histogramfile = new TFile(filename.c_str());
	if(histogramfile->IsZombie()) {
		cout << " Error: Could not read the given file containing histograms! " << endl;
		return false;
	}
	_hist_angles_vert_mean   = (TH2*)gDirectory->Get("hist_angles_vert_mean");
	_hist_angles_horiz_mean  = (TH2*)gDirectory->Get("hist_angles_horiz_mean");
	_hist_angles_vert_sigma  = (TH2*)gDirectory->Get("hist_angles_vert_sigma");
	_hist_angles_horiz_sigma = (TH2*)gDirectory->Get("hist_angles_horiz_sigma");
	_hist_vertex_distr_xy    = (TH2*)gDirectory->Get("hist_vertex_distr_xy");
	_hist_vertex_distr_z     = (TH1*)gDirectory->Get("hist_vertex_distr_z");
	// check if it worked
	if(!_hist_angles_vert_mean ||
	   !_hist_angles_horiz_mean ||
	   !_hist_angles_vert_sigma ||
	   !_hist_angles_horiz_sigma ||
	   !_hist_vertex_distr_xy ||
	   !_hist_vertex_distr_z)
	{
		cout << " Error: histograms not found!" << endl;
		return false;
	}
	// draw the output
	if(plot) {
		TCanvas *analyze_beam_properties_output = new TCanvas("analyze_beam_properties_output", "analyze_beam_properties_output method output", 600, 900);
		analyze_beam_properties_output->Divide(2,3);
		analyze_beam_properties_output->cd(1);
		_hist_angles_vert_mean->Draw("COLZ");
		gPad->Update();
		analyze_beam_properties_output->cd(2);
		_hist_angles_vert_sigma->Draw("COLZ");
		gPad->Update();
		analyze_beam_properties_output->cd(3);
		_hist_angles_horiz_mean->Draw("COLZ");
		gPad->Update();
		analyze_beam_properties_output->cd(4);
		_hist_angles_horiz_sigma->Draw("COLZ");
		gPad->Update();
		analyze_beam_properties_output->cd(5);
		gPad->SetLogz();
		_hist_vertex_distr_xy->Draw("COLZ");
		gPad->Update();
		analyze_beam_properties_output->cd(6);
		_hist_vertex_distr_z->Draw("");
		gPad->Update();
		analyze_beam_properties_output->Print("Fitted_beam_property_distributions.pdf");
	}
	return true;
}


bool primaryVertexGen::check() {
	return _histograms_loaded;
}


TVector3& primaryVertexGen::getVertex(const float cutR,
                                      const float cutZ_low,
                                      const float cutZ_high)
{
	if(!_histograms_loaded) {
		printErr << "cannot produce vertex before histograms are loaded correctly." << endl;
		throw;
	}
	double x;
	double y;
	double z;
	do {
		_hist_vertex_distr_xy->GetRandom2(x, y);
		z = _hist_vertex_distr_z->GetRandom();
	} while((sqrt(x*x + y*y) > cutR) || (z < cutZ_low) || (z > cutZ_high));
	TVector3* result = new TVector3(x,y,z);
	return *result;
}


TVector3& primaryVertexGen::getBeamDir(const TVector3 vertex) {
	TVector3* result = new TVector3(0.,0.,1.);
	if(!_histograms_loaded) {
		printErr << "cannot produce beam before histograms are loaded correctly." << endl;
		throw;
	}
	// check whether we are in the valid ranges of the histograms
	double sigma_horiz = _hist_angles_horiz_sigma->GetBinContent(
	          _hist_angles_horiz_sigma->FindBin(vertex.X(), vertex.Y()));
	double sigma_vert = _hist_angles_vert_sigma->GetBinContent(
	          _hist_angles_vert_sigma->FindBin(vertex.X(), vertex.Y()));
	if((sigma_horiz == 0.) || (sigma_vert == 0.)) {
		printErr << "vertex position out of range." << endl;
		throw;
	}
	// if check passed retrieve interpolated values
	double mean_horiz = _hist_angles_horiz_mean->Interpolate(vertex.X(), vertex.Y());
	double mean_vert  = _hist_angles_vert_mean->Interpolate(vertex.X(), vertex.Y());
	sigma_horiz = _hist_angles_horiz_sigma->Interpolate(vertex.X(), vertex.Y());
	sigma_vert  = _hist_angles_vert_sigma->Interpolate(vertex.X(), vertex.Y());
	// now we can retrieve the direction randomly
	TRandom3* random = randomNumberGenerator::instance()->getGenerator();
	double angle_horiz = random->Gaus(mean_horiz, sigma_horiz);
	double angle_vert  = random->Gaus(mean_vert , sigma_vert);
	result->RotateX(-angle_vert);
	result->RotateY(angle_horiz);
	return *result;
}

TLorentzVector& primaryVertexGen::getBeamPart(const TVector3 beam_dir) {
	TRandom3* random = randomNumberGenerator::instance()->getGenerator();
	double energy = random->Gaus(_beam_energy_mean, _beam_energy_sigma);
	TVector3 _beam_dir(beam_dir);
	// m² = E² - p² -> |p| = sqrt(E²-m²)
	_beam_dir.SetMag(sqrt(energy*energy-_beam_part_mass*_beam_part_mass));
	TLorentzVector* result = new TLorentzVector(_beam_dir.X(), _beam_dir.Y(), _beam_dir.Z(), energy);
	return *result;
}

void primaryVertexGen::gesPalette(int i) {
	TPad foo; // never remove this line :-)))
	if(i == 0) {
		const Int_t NRGBs = 5;
		const Int_t NCont = 255;
		Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
		Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
		Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
		Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
		TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
		gStyle->SetNumberContours(NCont);
	}

	if(i == 1) {
		const UInt_t Number = 3;
		Double_t Red[Number]    = { 1.00, 0.00, 0.00 };
		Double_t Green[Number]  = { 0.00, 1.00, 0.00 };
		Double_t Blue[Number]   = { 1.00, 0.00, 1.00 };
		Double_t Length[Number] = { 0.00, 0.50, 1.00 };
		Int_t nb=255;
		TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
		gStyle->SetNumberContours(nb);
	}
}

void Beam_simulation() {
	primaryVertexGen primaryVertexGen("properties_2008/primary_vertex_properties.root",
	                                  0.13957018,
	                                  191.29,
	                                  1.94476);
	if(!primaryVertexGen.check()) {
		return;
	}
	// some cross check histograms
	//const int n_steps_x = 21;
	//const int n_steps_y = 21;
	//const float cell_size_x = 4.;
	//const float cell_size_y = 4.;

	TH2F *hist_vertex_distr_xy_sim   = new TH2F("hist_vertex_distr_xy_sim", "simulated vertex X Y distribution", 200, -2., 2., 200, -2., 2.);
	TH1F *hist_vertex_distr_z_sim    = new TH1F("hist_vertex_distr_z_sim","simulated vertex Z distribution",1000, -100, 0);
	TH1F *hist_angle_distr_horiz     = new TH1F("hist_angle_distr_horiz","horizontal angle distribution",1000, -0.001, 0.001);
	TH1F *hist_angle_distr_vert      = new TH1F("hist_angle_distr_vert ","vertical angle distribution",1000, -0.001, 0.001);
	TH1F *hist_energy_distr_sim      = new TH1F("hist_energy_distr_sim", "simulated energy distribution",1000, 170, 210);
	TH2F* hist_beamtrack_horiz_plane = new TH2F("hist_beamtrack_horiz_plane", "beamtrack in horizontal plane", 1000, 0, 6000, 100, -5, 5);
	TH2F* hist_beamtrack_vert_plane  = new TH2F("hist_beamtrack_vert_plane", "beamtrack in vertical plane", 1000, 0, 6000, 100, -5, 5);
	TH1F* hist_horiz_beamdivergence_at_CEDAR = new TH1F("hist_horiz_beamdivergence_at_CEDAR", "horizontal beam divergence at CEDAR region", 1000, -400e-6, 400e-6);
	TH1F* hist_vert_beamdivergence_at_CEDAR  = new TH1F("hist_vert_beamdivergence_at_CEDAR", "vertical beam divergence at CEDAR region", 1000, -400e-6, 400e-6);

	int counter = 0;
	TCanvas* canvas = new TCanvas("canvas", "simulated distributions", 900, 600);
	canvas->Divide(3,2);
	TCanvas* canvas2 = new TCanvas("canvas2", "simulated tracks", 600, 600);
	canvas2->Divide(1,2);
	TCanvas* canvas3 = new TCanvas("canvas3", "simulated distributions at CEDAR region", 600, 400);
	canvas3->Divide(2,1);
	for(int i = 0; i < 1000000; i++) {
		TVector3 vertex = primaryVertexGen.getVertex();
		TVector3 beam_dir = primaryVertexGen.getBeamDir(vertex);
		if(beam_dir.Mag() == 0) {
			//cout << " skipping " << endl;
			continue;
		}
		hist_vertex_distr_xy_sim->Fill(vertex.X(), vertex.Y());
		hist_vertex_distr_z_sim->Fill(vertex.Z());
		TLorentzVector beam_part = primaryVertexGen.getBeamPart(beam_dir);
		//cout << beam_part.M() << endl;
		double azi = beam_part.Px()/beam_part.Pz();
		double dip = beam_part.Py()/beam_part.Pz();
		hist_angle_distr_horiz->Fill(azi);
		hist_angle_distr_vert->Fill(dip);
		hist_energy_distr_sim->Fill(beam_part.E());

		// check the beam spot 30m downstream by extrapolation from 0 to 6000 cm
		for(int i = 0; i < 1000; i++) {
			// get the distance between pos z of the vertex to the point
			// to be extrapolated
			float pos_z = i*6;
			float extrapol_dist_z =  pos_z - vertex.Z();
			// now we can calculate the factor for the pointing vector for extrapolation
			float extrapol_fac = extrapol_dist_z/beam_dir.Z();
			// this factor has to be applied in all 3 dimensions
			hist_beamtrack_horiz_plane->Fill(pos_z, vertex.X()+beam_dir.X()*extrapol_fac);
			hist_beamtrack_vert_plane ->Fill(pos_z, vertex.Y()+beam_dir.Y()*extrapol_fac);
		}
		// have a look on the properties at the CEDAR region
		// transport the vertex position 30m downstream
		float pos_z = 3000;
		float extrapol_dist_z =  pos_z - vertex.Z();
		float extrapol_fac = extrapol_dist_z/beam_dir.Z();
		vertex.SetXYZ(vertex.X()+beam_dir.X()*extrapol_fac,
		              vertex.Y()+beam_dir.Y()*extrapol_fac,
		              pos_z);

		// transport the particle track back to the CEDAR postion upstream
		// by using the transportation matrix output given by Lau 26.01.2009
		/*
		 * 1COMPASS HADRON OPTICS 200 GEV V.2007                                                      TRANSPORT RUN26/01/09
		0POSITION TYPE      STRENGTH *         H O R I Z O N T A L      *           V E R T I C A L        *         D I S P E R S I O N
		METERS          T*M,T/M*M *     R11     R12     R21     R22  *     R33     R34     R43     R44  *     R16     R26     R36     R46
		                 T/M**2*M *    MM/MM   MM/MR   MR/MM   MR/MR *    MM/MM   MM/MR   MR/MM   MR/MR *    MM/PC   MR/PC   MM/PC   MR/PC
		************************************************************************************************************************************
		30.000  3                 *    1.000  30.000   0.000   1.000 *    1.000  30.000   0.000   1.000 *    0.000   0.000   0.000   0.000
		69.396  3                 *    0.670  47.160  -0.021   0.000 *    0.924  73.374  -0.014   0.000 *    0.000   0.000   0.772   0.022
		75.713  3  CED2           *    0.536  47.160  -0.021   0.000 *    0.838  73.374  -0.014   0.000 *    0.000   0.000   0.908   0.022
		75.913  3                 *    0.531  47.160  -0.021   0.000 *    0.835  73.374  -0.014   0.000 *    0.000   0.000   0.912   0.022
		82.230  3  CED1           *    0.398  47.160  -0.021   0.000 *    0.749  73.374  -0.014   0.000 *    0.000   0.000   1.048   0.022
		                                         ^very parallel beam               ^very parallel beam :)
		 */

		// compass is measuring in cm!
		// x is horizontal plane
		// y is vertical plane
		// the azimuth lies in x
		// the dip in y
		// at 82.230 m upstream the CEDAR in Lau's table
		float x_displacement_upstream = vertex.X()*(0.398) + azi*(47.160*100); // mm/mrad -> cm/rad = * 1000/10 = * 100
		float y_displacement_upstream = vertex.Y()*(0.749) + dip*(73.374*100);
		// at 69.396 m downstream the CEDAR in Lau's table
		float x_displacement_downstream = vertex.X()*(0.670) + azi*(47.160*100);
		float y_displacement_downstream = vertex.Y()*(0.924) + dip*(73.374*100);
		// and calculate the beam divergence
		float divergence_horiz = (x_displacement_downstream-x_displacement_upstream) / ((82.230-69.396) * 100);
		float divergence_vert  = (y_displacement_downstream-y_displacement_upstream) / ((82.230-69.396) * 100);
		hist_horiz_beamdivergence_at_CEDAR->Fill(divergence_horiz);
		hist_vert_beamdivergence_at_CEDAR->Fill(divergence_vert);

		++counter;
		if((counter % 10000) == 0) {
			canvas->cd(1);
			gPad->Clear();
			gPad->SetLogz();
			hist_vertex_distr_xy_sim->Draw("COLZ");
			gPad->Update();
			canvas->cd(2);
			gPad->Clear();
			hist_vertex_distr_z_sim->Draw("");
			gPad->Update();
			canvas->cd(3);
			gPad->Clear();
			hist_angle_distr_horiz->Draw("");
			gPad->Update();
			canvas->cd(4);
			gPad->Clear();
			hist_angle_distr_vert->Draw("");
			gPad->Update();
			canvas->cd(5);
			gPad->Clear();
			hist_energy_distr_sim->Draw("");
			gPad->Update();
			canvas2->cd(1);
			gPad->Clear();
			hist_beamtrack_horiz_plane->Draw("COLZ");
			gPad->Update();
			canvas2->cd(2);
			gPad->Clear();
			hist_beamtrack_vert_plane->Draw("COLZ");
			gPad->Update();
			canvas3->cd(1);
			gPad->Clear();
			hist_horiz_beamdivergence_at_CEDAR->Draw();
			gPad->Update();
			canvas3->cd(2);
			gPad->Clear();
			hist_vert_beamdivergence_at_CEDAR->Draw();
			gPad->Update();
		}
	}
	canvas2->Print("Simulated_tracks.pdf");
	canvas->Print("Simulated_distributions.pdf");
	canvas3->Print("Simulated_distributions_at_CEDAR.pdf");
}

