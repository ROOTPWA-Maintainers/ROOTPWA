#include "TCanvas.h"
#include "TTree.h"
#include <string>
#include <vector>
#include <utils>
#include <iostream>
#include <map>

using namespace std;

void plot_all_phase_WA03(){
	TTree* tree = loadFitResult("/home/Promme/disk2/analysis/Kp_Kppipi/PWA2/FIT_SET_WA03/*.result.root");
	TCanvas* canvas_result = new TCanvas();

	// phi 2-S0+(K**pi) vs. 1+S0+(rhoK)
	plot4(tree, 13, 4, 1.4, 2.2,"fitResult_v2", canvas_result);
	string filename(canvas_result->GetName());
	canvas_result->Print(("phase_"+filename+".pdf").c_str());

	delete canvas_result;
	canvas_result = new TCanvas();

	// 2-S0+(fK) vs. 2-S0+(K**pi)
	plot4(tree, 12, 13, 1.4, 2.2,"fitResult_v2", canvas_result);
	filename = canvas_result->GetName();
	canvas_result->Print(("phase_"+filename+".pdf").c_str());

	delete canvas_result;
	canvas_result = new TCanvas();

	// 2-P0+(K*pi) vs. 2-S0+(K**pi)
	plot4(tree, 14, 13, 1.4, 2.2,"fitResult_v2", canvas_result);
	filename = canvas_result->GetName();
	canvas_result->Print(("phase_"+filename+".pdf").c_str());

	// 2+D1+(K*pi) vs. 1+S0+(K*pi)
	plot4(tree, 16, 10, 1., 2.,"fitResult_v2", canvas_result);
	filename = canvas_result->GetName();
	canvas_result->Print(("phase_"+filename+".pdf").c_str());

	delete canvas_result;
	canvas_result = new TCanvas();

	// 2+D1+(rhoK) vs. 1+S0+(K*pi)
	plot4(tree, 17, 10, 1., 2.,"fitResult_v2", canvas_result);
	filename = canvas_result->GetName();
	canvas_result->Print(("phase_"+filename+".pdf").c_str());

	delete canvas_result;
	canvas_result = new TCanvas();

	// 1+S1+(K*pi) vs. 1+S0+(K*pi)
	plot4(tree, 13, 10, 1., 2.,"fitResult_v2", canvas_result);
	filename = canvas_result->GetName();
	canvas_result->Print(("phase_"+filename+".pdf").c_str());

	delete canvas_result;
	canvas_result = new TCanvas();

	// 1+S1+(rhoK) vs. 1+S0+(K*pi)
	plot4(tree, 11, 10, 1., 2.,"fitResult_v2", canvas_result);
	filename = canvas_result->GetName();
	canvas_result->Print(("phase_"+filename+".pdf").c_str());

	delete canvas_result;
	canvas_result = new TCanvas();

	// 0-P0+(K*pi) vs. 1+S0+(K*pi)
	plot4(tree, 0, 10, 1., 2.,"fitResult_v2", canvas_result);
	filename = canvas_result->GetName();
	canvas_result->Print(("phase_"+filename+".pdf").c_str());

	delete canvas_result;
	canvas_result = new TCanvas();

	// 0-P0+(rhoK) vs. 1+S0+(K*pi)
	plot4(tree, 17, 10, 1., 2.,"fitResult_v2", canvas_result);
	filename = canvas_result->GetName();
	canvas_result->Print(("phase_"+filename+".pdf").c_str());

	delete canvas_result;
	canvas_result = new TCanvas();

	// 0-S0+(eK) vs. 1+S0+(K*pi)
	plot4(tree, 2, 10, 1., 2.,"fitResult_v2", canvas_result);
	filename = canvas_result->GetName();
	canvas_result->Print(("phase_"+filename+".pdf").c_str());

	delete canvas_result;
	canvas_result = new TCanvas();

	// 1+S0+(K*pi) vs. 1+S0+(K*pi)
	plot4(tree, 10, 10, 1., 2.,"fitResult_v2", canvas_result);
	filename = canvas_result->GetName();
	canvas_result->Print(("phase_"+filename+".pdf").c_str());

	delete canvas_result;
	canvas_result = new TCanvas();

	// 1+S0+(rhoK) vs. 1+S0+(K*pi)
	plot4(tree, 6, 10, 1., 2.,"fitResult_v2", canvas_result);
	filename = canvas_result->GetName();
	canvas_result->Print(("phase_"+filename+".pdf").c_str());

	delete canvas_result;
	canvas_result = new TCanvas();

	// 1+P0+(kpi) vs. 1+S0+(K*pi)
	plot4(tree, 3, 10, 1., 2.,"fitResult_v2", canvas_result);
	filename = canvas_result->GetName();
	canvas_result->Print(("phase_"+filename+".pdf").c_str());

	delete canvas_result;
	canvas_result = new TCanvas();

	// 1+P0+(epi) vs. 1+S0+(K*pi)
	plot4(tree, 8, 10, 1., 2.,"fitResult_v2", canvas_result);
	filename = canvas_result->GetName();
	canvas_result->Print(("phase_"+filename+".pdf").c_str());

	delete canvas_result;
	canvas_result = new TCanvas();

	// 1+D0+(K*pi) vs. 1+S0+(K*pi)
	plot4(tree, 5, 10, 1., 2.,"fitResult_v2", canvas_result);
	filename = canvas_result->GetName();
	canvas_result->Print(("phase_"+filename+".pdf").c_str());
}
