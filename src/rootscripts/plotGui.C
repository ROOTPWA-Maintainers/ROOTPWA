
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>

#include <TApplication.h>
#include <TGClient.h>
#include <TGButton.h>
#include <TGListBox.h>
#include <TList.h>
#include <TFile.h>
#include <TTree.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TGraphErrors.h>

#include <fitResult.h>
#include <reportingUtils.hpp>


using namespace rpwa;
using namespace std;

vector<string>
readWaveList(const string& waveListFileName)
{
	printInfo << "reading amplitude names and thresholds from wave list file "
	          << "'" << waveListFileName << "'." << endl;
	ifstream waveListFile(waveListFileName.c_str());
	if (not waveListFile) {
		printErr << "cannot open file '" << waveListFileName << "'. aborting." << endl;
		throw;
	}
	vector<string> waveNames;
	unsigned int         countWave = 0;
	unsigned int         lineNmb   = 0;
	string               line;
	while (getline(waveListFile, line)) {
		if (line[0] == '#')  // comments start with #
			continue;
		stringstream lineStream;
		lineStream.str(line);
		string waveName;
		if (lineStream >> waveName) {
			double threshold;
			// !!! it would be safer to make the threshold value in the wave list file mandatory
			if (not (lineStream >> threshold))
				threshold = 0;
			waveNames.push_back(waveName);
			++countWave;
		} else
			printWarn << "cannot parse line '" << line << "' in wave list file "
			          << "'" << waveListFileName << "'" << endl;
		++lineNmb;
	}
	waveListFile.close();
	return waveNames;
}

class MyMainFrame : public TGMainFrame {

	private:
		TGListBox           *_listBox;
		TGListBox           *_listBox2;
		TGCheckButton       *_checkMulti;
		TList               *_selected;
		std::vector<std::string> _waveNames;
		std::vector<const rpwa::fitResult*> _fitResults;
		const double _binWidth;

	public:
		MyMainFrame(const TGWindow *p,
		            UInt_t w,
		            UInt_t h,
		            TTree* tree,
		            const double& binWidth,
		            const double& intensityThreshold,
		            const string& waveListFileName);
		virtual ~MyMainFrame();
		void DoExit();
		void DoSelect();
		void HandleButtons();
		void PrintSelected();

		ClassDef(MyMainFrame, 0)
};

void MyMainFrame::DoSelect()
{
	Printf("Slot DoSelect()");
}

void MyMainFrame::DoExit()
{
	Printf("Slot DoExit()");
	gApplication->Terminate(0);
}

MyMainFrame::MyMainFrame(const TGWindow *p,
                         UInt_t w,
                         UInt_t h,
                         TTree* tree,
                         const double& binWidth,
                         const double& intensityThreshold,
                         const string& waveListFileName) :
	TGMainFrame(p, w, h),
	_binWidth(binWidth)
{

	{
		fitResult* res = 0;//new fitResult();
		string branchname = "fitResult_v2";
		if(tree->FindBranch(branchname.c_str()) == 0){
			cerr << "Invalid branch ." << branchname << endl;
			throw;
		}
		tree->SetBranchAddress(branchname.c_str(),&res);
		long entries = tree->GetEntries();
		if(entries < 1) {
			cerr << "Tree has no entries" << endl;
			throw;
		}
		cout << "Mass bins: " << entries << endl;
		_fitResults.resize(entries, 0);
		for(long i = 0; i < entries; ++i) {
			cout << "loading result " << i << " of " << entries << "... " << std::flush;
			tree->GetEntry(i);
			_fitResults[i] = new rpwa::fitResult(*res);
			cout << "done!" << endl;
		}
	}

	vector<string> whitelistedWaves;
	if(waveListFileName != "") {
		cout << "reading wave whitelist '" << waveListFileName << "'." << endl;
		whitelistedWaves = readWaveList(waveListFileName);
		cout << "found " << whitelistedWaves.size() << " white-listed waves." << endl;
	}

	// get list of waves
	_waveNames = _fitResults[0]->waveNames();
	if(intensityThreshold > 0.) {
		cout << "found threshold (" << intensityThreshold << "), starting with " << _waveNames.size() << " waves." << endl;
		vector<string> wavesToBeRemoved;
		for(unsigned int i = 0; i < _waveNames.size(); ++i) {
			const string& waveName = _waveNames[i];
			if(whitelistedWaves.size() > 0 and std::find(whitelistedWaves.begin(), whitelistedWaves.end(), waveName) == whitelistedWaves.end()) {
				wavesToBeRemoved.push_back(waveName);
				continue;
			}
			bool aboveThreshold = false;
			cout << "checking wave " << waveName << " (" << i+1 << "/" << _waveNames.size() << ")" << endl;
			for(unsigned int bin_i = 0; bin_i < _fitResults.size(); ++bin_i) {
				const double intensity = _fitResults[bin_i]->intensity(waveName.c_str());
				cout << "intensity[" << bin_i+1 << "/" << _fitResults.size() << "] = " << intensity << endl;
				if(intensity >= intensityThreshold) {
					aboveThreshold = true;
					cout <<"accepted!"<<endl;
					break;
				}
			}
			if(not aboveThreshold) {
				wavesToBeRemoved.push_back(waveName);
			}
		}
		for(unsigned int i = 0; i < wavesToBeRemoved.size(); ++i) {
			_waveNames.erase(std::find(_waveNames.begin(), _waveNames.end(), wavesToBeRemoved[i]));
		}
		cout << "after thresholding, " << _waveNames.size() << " waves are left." << endl;
	}


	// Create main frame
	_listBox = new TGListBox(this, 89);
	_listBox2 = new TGListBox(this, 88);
	_selected = new TList;

	for (unsigned int i = 0; i < _waveNames.size(); ++i) {
		cout << _waveNames[i] << endl;
		_listBox->AddEntry(_waveNames[i].c_str(), i);
		_listBox2->AddEntry(_waveNames[i].c_str(), i);
	}
	_listBox->Resize(400,250);
	_listBox2->Resize(400,250);

	AddFrame(_listBox, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 5, 5, 5, 5));
	AddFrame(_listBox2, new TGLayoutHints(kLHintsTop | kLHintsRight| kLHintsExpandX, 5, 5, 5, 5));

	_checkMulti = new TGCheckButton(this, "&Mutliple selection", 10);
	AddFrame(_checkMulti, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));
	_checkMulti->Connect("Clicked()", "MyMainFrame", this, "HandleButtons()");
	// Create a horizontal frame containing button(s)
	TGHorizontalFrame *hframe = new TGHorizontalFrame(this, 400, 20, kFixedWidth);
	TGTextButton *show = new TGTextButton(hframe, "&Show");
	show->SetToolTipText("Click here to print the selection you made");
	show->Connect("Pressed()", "MyMainFrame", this, "PrintSelected()");
	hframe->AddFrame(show, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));
	TGTextButton *exit = new TGTextButton(hframe, "&Exit ");
	exit->Connect("Pressed()", "MyMainFrame", this, "DoExit()");
	hframe->AddFrame(exit, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));
	AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 2, 2, 5, 1));

	// Set a name to the main frame
	SetWindowName("List Box");
	MapSubwindows();

	// Initialize the layout algorithm via Resize()
	Resize(GetDefaultSize());

	// Map main frame
	MapWindow();
	_listBox->Select(0);
	_listBox2->Select(0);
}

MyMainFrame::~MyMainFrame()
{
	// Clean up main frame...
	Cleanup();
	if (_selected) {
		_selected->Delete();
		delete _selected;
	}
}

void MyMainFrame::HandleButtons()
{
	// Handle check button.
	Int_t id;
	TGButton *btn = (TGButton *) gTQSender;
	id = btn->WidgetId();

	printf("HandleButton: id = %d\n", id);

	//if (id == 10)
	//   fListBox->SetMultipleSelections(fCheckMulti->GetState());
}


void MyMainFrame::PrintSelected()
{
	// Writes selected entries in TList if multiselection.
	string w1=_waveNames[_listBox->GetSelected()];
	string w2=_waveNames[_listBox2->GetSelected()];
	cout << w1 << endl;
	cout << w2 << endl;
	// Produce plots
	unsigned int nFitResults = _fitResults.size();

	Int_t colour=1;

	TGraphErrors* gph = new TGraphErrors(nFitResults);
	stringstream graphName;
	graphName << "PHI"<<w1<<"---"<<"PHI"<<w2;
	cout << "creating graph   " << graphName.str() << endl;
	gph->SetName (graphName.str().c_str());
	gph->SetTitle(graphName.str().c_str());
	gph->SetMarkerStyle(21);
	gph->SetMarkerSize(0.5);
	gph->SetMarkerColor(colour);
	gph->SetLineColor(colour);

	TGraphErrors* gphP1 = (TGraphErrors*)gph->Clone("gph+1");
	TGraphErrors* gphM1 = (TGraphErrors*)gph->Clone("gph-1");
	gphP1->SetMarkerColor(2);
	gphP1->SetLineColor(2);
	gphM1->SetMarkerColor(3);
	gphM1->SetLineColor(3);

	TGraphErrors* gRe = new TGraphErrors(nFitResults);
	graphName.str("");
	graphName.clear();
	graphName << "RE_"<<w1<<"---"<<""<<w2;
	cout << "creating graph   " << graphName.str() << endl;
	gRe->SetName (graphName.str().c_str());
	gRe->SetTitle(graphName.str().c_str());
	gRe->SetMarkerStyle(21);
	gRe->SetMarkerSize(0.5);
	gRe->SetMarkerColor(colour);
	gRe->SetLineColor(colour);

	TGraphErrors* gIm = new TGraphErrors(nFitResults);
	graphName.str("");
	graphName.clear();
	graphName << "IM_"<<w1<<"---"<<""<<w2;
	cout << "creating graph   " << graphName.str() << endl;
	gIm->SetName (graphName.str().c_str());
	gIm->SetTitle(graphName.str().c_str());
	gIm->SetMarkerStyle(21);
	gIm->SetMarkerSize(0.5);
	gIm->SetMarkerColor(colour);
	gIm->SetLineColor(colour);

	TGraphErrors* g1 = new TGraphErrors(nFitResults);
	graphName.str("");
	graphName.clear();
	graphName<<"g"<<w1;
	g1->SetName (graphName.str().c_str());
	g1->SetTitle(graphName.str().c_str());
	g1->SetMarkerStyle(21);
	g1->SetMarkerSize(0.5);
	TGraphErrors* g2 = new TGraphErrors(nFitResults);
	graphName.str("");
	graphName.clear();
	graphName<<"g"<<w2;
	g2->SetName (graphName.str().c_str());
	g2->SetTitle(graphName.str().c_str());
	g2->SetMarkerStyle(21);
	g2->SetMarkerSize(0.5);

	for(unsigned int i = 0; i < nFitResults; ++i){
		const fitResult* result = _fitResults[i];

		if(not result->converged()) {
			continue;
		}

		const double intensity1=result->intensity(w1.c_str());
		if((numeric_limits<double>::has_infinity and intensity1 == numeric_limits<double>::infinity()) or intensity1!=intensity1) {
			continue;
		}

		g1->SetPoint(i, result->massBinCenter()*0.001, intensity1);
		g1->SetPointError(i, _binWidth*0.5, result->intensityErr(w1.c_str()));

		const double intensity2=result->intensity(w2.c_str());
		if((numeric_limits<double>::has_infinity and intensity2 == numeric_limits<double>::infinity()) or intensity2!=intensity2) {
			continue;
		}

		g2->SetPoint(i, result->massBinCenter()*0.001, intensity2);
		g2->SetPointError(i, _binWidth*0.5, result->intensityErr(w2.c_str()));

		double ph=result->phase(w1.c_str(),w2.c_str());
		const double pherr=result->phaseErr(w1.c_str(),w2.c_str());
		// check if we should make a transformation by 2pi
		// this is needed because of cyclical variable phi
		if(i>11){
			double mpre;
			double phpre;
			gph->GetPoint(i-1,mpre,phpre);
			double diff1=fabs(ph-phpre);
			double diff2=fabs(ph+360-phpre);
			double diff3=fabs(ph-360-phpre);
			if(diff2<diff1 && diff2<diff3)ph+=360;
			else if(diff3<diff1 && diff3<diff2)ph-=360;
		}

		gph->SetPoint(i, result->massBinCenter()*0.001, ph);
		gph->SetPointError(i, _binWidth*0.5, pherr);

		// add point +- 360 degree
		gphP1->SetPoint(i, result->massBinCenter()*0.001, ph+360);
		gphP1->SetPointError(i, _binWidth*0.5, pherr);

		gphM1->SetPoint(i, result->massBinCenter()*0.001, ph-360);
		gphM1->SetPointError(i, _binWidth*0.5, pherr);

		unsigned int wi1=result->waveIndex(w1);
		unsigned int wi2=result->waveIndex(w2);
		complex<double> rho=result->spinDensityMatrixElem(wi1,wi2);
		TMatrixT<double> rhoCov=result->spinDensityMatrixElemCov(wi1,wi2);
		gRe->SetPoint(i, result->massBinCenter()*0.001, rho.real());
		gRe->SetPointError(i, _binWidth*0.5, sqrt(rhoCov[0][0]));
		gIm->SetPoint(i, result->massBinCenter()*0.001, rho.imag());
		gIm->SetPointError(i, _binWidth*0.5, sqrt(rhoCov[1][1]));
	}// end loop over bins

	// plot graphs
	TCanvas*c=new TCanvas("c","c",10,10,1200,800);
	c->Divide(2,3);
	c->cd(1);
	gph->Draw("AP");
	gph->GetXaxis()->SetTitle("5#pi mass (GeV/c^2)");
	gph->GetYaxis()->SetTitle("Phase difference");
	gph->GetYaxis()->SetRangeUser(-270,270);

	gphP1->Draw("PSAME");
	gphM1->Draw("PSAME");

	c->cd(3);
	g1->Draw("AP");
	g1->GetXaxis()->SetTitle("5#pi mass (GeV/c^2)");
	g1->GetYaxis()->SetTitle("Intensity");
	c->cd(5);
	g2->Draw("AP");
	g2->GetXaxis()->SetTitle("5#pi mass (GeV/c^2)");
	g2->GetYaxis()->SetTitle("Intensity");

	c->cd(2);
	gRe->Draw("AP");
	gRe->GetXaxis()->SetTitle("5#pi mass (GeV/c^2)");
	gRe->GetYaxis()->SetTitle("Re(#rho_{ij})");
	c->cd(4);
	gIm->Draw("AP");
	gIm->GetXaxis()->SetTitle("5#pi mass (GeV/c^2)");
	gIm->GetYaxis()->SetTitle("Im(#rho_{ij})");

}

void plotGui(const std::string& infilename,
             const double binWidth = 0.03,
             const double intensityThreshold = -1.)
{

	// load fitResult tree
	TFile* infile=TFile::Open(infilename.c_str(),"READ");
	if(infile==NULL){
		cerr << "File " << infilename << " not found!"<< endl;
		return;
	}
	TTree* pwa=(TTree*)infile->FindObjectAny("pwa");
	if(pwa==NULL){
		cerr << "Tree not found!"<< endl;
		return;
	}

	// Popup the GUI...
	new MyMainFrame(gClient->GetRoot(), 20, 20, pwa, binWidth, intensityThreshold, "");
}
