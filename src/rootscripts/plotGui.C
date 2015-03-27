
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <vector>

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


set<string>
readWaveList(const string& waveListFileName)
{
	printInfo << "reading amplitude names and thresholds from wave list file "
	          << "'" << waveListFileName << "'." << endl;
	ifstream waveListFile(waveListFileName.c_str());
	if (not waveListFile) {
		printErr << "cannot open file '" << waveListFileName << "'. Aborting..." << endl;
		throw;
	}
	set<string> waveNames;
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
			waveNames.insert(waveName);
			++countWave;
		} else
			printWarn << "cannot parse line '" << line << "' in wave list file "
			          << "'" << waveListFileName << "'" << endl;
		++lineNmb;
	}
	waveListFile.close();
	return waveNames;
}


class plotGuiMainFrame : public TGMainFrame {

	private:
		TGListBox*          _listBox1;
		TGListBox*          _listBox2;
		TGCheckButton*      _checkDrawNewCanvas;
		std::vector<std::string> _waveNames;
		std::vector<const rpwa::fitResult*> _fitResults;
		const double _binWidth;
		unsigned int _canvasCounter;
		TCanvas* _currentCanvas;

	public:
		plotGuiMainFrame(const TGWindow *p,
		            UInt_t w,
		            UInt_t h,
		            TTree* tree,
		            const double& binWidth,
		            const double& intensityThreshold,
		            const string& waveListFileName);
		virtual ~plotGuiMainFrame();
		void DoExit();
		void DoSelect();
		void ActiveCanvasClosed();
		void HandleButtons();
		void PrintSelected();

		static bool _debug;

		ClassDef(plotGuiMainFrame, 0)

};


bool plotGuiMainFrame::_debug = false;


void plotGuiMainFrame::DoSelect()
{
	if(_debug) {
		printDebug << "Slot DoSelect()" << endl;
	}
}


void plotGuiMainFrame::DoExit()
{
	if(_debug) {
		printDebug << "Slot DoExit()" << endl;
	}
	gApplication->Terminate(0);
}


plotGuiMainFrame::plotGuiMainFrame(const TGWindow *p,
                         UInt_t w,
                         UInt_t h,
                         TTree* tree,
                         const double& binWidth,
                         const double& intensityThreshold,
                         const string& waveListFileName) :
	TGMainFrame(p, w, h),
	_binWidth(binWidth),
	_canvasCounter(0),
	_currentCanvas(0)
{

	{
		fitResult* res = 0;//new fitResult();
		string branchname = "fitResult_v2";
		if(tree->FindBranch(branchname.c_str()) == 0){
			printErr << "invalid branch '" << branchname << "'. Aborting..." << endl;
			throw;
		}
		tree->SetBranchAddress(branchname.c_str(),&res);
		long entries = tree->GetEntries();
		if(entries < 1) {
			printErr << "tree has no entries. Aborting..." << endl;
			throw;
		}
		printInfo << "Mass bins: " << entries << endl;
		_fitResults.resize(entries, 0);
		for(long i = 0; i < entries; ++i) {
			printInfo << "loading result " << i << " of " << entries << "... " << std::flush;
			tree->GetEntry(i);
			_fitResults[i] = new rpwa::fitResult(*res);
			cout << "done!" << endl;
		}
	}

	{
		// get list of waves
		set<string> whitelistedWaves;
		if(waveListFileName != "") {
			printInfo << "reading wave whitelist '" << waveListFileName << "'." << endl;
			whitelistedWaves = readWaveList(waveListFileName);
			printInfo << "found " << whitelistedWaves.size() << " white-listed waves." << endl;
		}
		set<string> waveNameCollector;
		for(unsigned int i = 0; i < _fitResults.size(); ++i) {
			const vector<string>& waveNames = _fitResults[i]->waveNames();
			waveNameCollector.insert(waveNames.begin(), waveNames.end());
		}
		vector<string> tentativeWaveNames(waveNameCollector.size());
		std::copy(waveNameCollector.begin(), waveNameCollector.end(), tentativeWaveNames.begin());
		printInfo << "found " << tentativeWaveNames.size() << " waves in input fit results." << endl;
		if(intensityThreshold > 0.) {
			printInfo << "checking thresholds for " << tentativeWaveNames.size() << " waves." << endl;
			for(unsigned int i = 0; i < tentativeWaveNames.size(); ++i) {
				const string& waveName = tentativeWaveNames[i];
				if(_debug) {
					printDebug << "checking if intensity is above threshold for wave '" << waveName << "'." << endl;
				}
				for(unsigned int bin_i = 0; bin_i < _fitResults.size(); ++bin_i) {
					const double intensity = _fitResults[bin_i]->intensity(waveName.c_str());
					if(intensity >= intensityThreshold) {
						whitelistedWaves.insert(waveName);
						break;
					}
				}
			}
			printInfo << "thresholds checked." << endl;
		}
		if(whitelistedWaves.empty()) {
			if(waveListFileName != "" or intensityThreshold > 0) {
				printErr << "thresholds and/or wave white-list did not leave any waves to plot. Aborting..." << endl;
				throw;
			}
			printInfo << "all waves will be added to the selection list." << endl;
			_waveNames = tentativeWaveNames;
		} else {
			for(set<string>::const_iterator it = whitelistedWaves.begin(); it != whitelistedWaves.end(); ++it) {
				if(std::find(tentativeWaveNames.begin(), tentativeWaveNames.end(), *it) != tentativeWaveNames.end()) {
					_waveNames.push_back(*it);
				} else {
					printWarn << "could not find white listed wave '" << *it << "' in any of the fit result files." << endl;
				}
			}
		}
	}
	printInfo << "ending up with " << _waveNames.size() << " waves to plot." << endl;

	// Create main frame
	_listBox1 = new TGListBox(this, 89);
	_listBox2 = new TGListBox(this, 88);

	for (unsigned int i = 0; i < _waveNames.size(); ++i) {
		if(_debug) {
			printDebug << _waveNames[i] << endl;
		}
		_listBox1->AddEntry(_waveNames[i].c_str(), i);
		_listBox2->AddEntry(_waveNames[i].c_str(), i);
	}
	_listBox1->Resize(400,250);
	_listBox2->Resize(400,250);

	AddFrame(_listBox1, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 5, 5, 5, 5));
	AddFrame(_listBox2, new TGLayoutHints(kLHintsTop | kLHintsRight| kLHintsExpandX, 5, 5, 5, 5));

	_checkDrawNewCanvas = new TGCheckButton(this, "&Draw in new canvas", 11);
	AddFrame(_checkDrawNewCanvas, new TGLayoutHints(kLHintsTop | kLHintsRight, 5, 5, 5, 5));
	// Create a horizontal frame containing button(s)
	TGHorizontalFrame *hframe = new TGHorizontalFrame(this, 400, 20, kFixedWidth);
	TGTextButton *show = new TGTextButton(hframe, "&Show");
	show->SetToolTipText("Click here to print the selection you made");
	show->Connect("Pressed()", "plotGuiMainFrame", this, "PrintSelected()");
	hframe->AddFrame(show, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));
	TGTextButton *exit = new TGTextButton(hframe, "&Exit ");
	exit->Connect("Pressed()", "plotGuiMainFrame", this, "DoExit()");
	hframe->AddFrame(exit, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));
	AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 2, 2, 5, 1));

	// Set a name to the main frame
	SetWindowName("List Box");
	MapSubwindows();

	// Initialize the layout algorithm via Resize()
	Resize(GetDefaultSize());

	// Map main frame
	MapWindow();
	_listBox1->Select(0);
	_listBox2->Select(0);
}


plotGuiMainFrame::~plotGuiMainFrame()
{
	// Clean up main frame...
	Cleanup();
}


void plotGuiMainFrame::HandleButtons() { }


void plotGuiMainFrame::ActiveCanvasClosed() {
	if(_debug) {
		printDebug << "canvas closed" << endl;
	}
	_currentCanvas = 0;
}


void plotGuiMainFrame::PrintSelected()
{
	// Writes selected entries in TList if multiselection.

	string w1 = "";
	string w2 = "";
	try {
		w1 = _waveNames.at(_listBox1->GetSelected());
		w2 = _waveNames.at(_listBox2->GetSelected());
	} catch (std::out_of_range&) {
		printErr << "selection lists seem to be corrupted." << endl;
		cout     << "   _listBox1->GetSelected() = " << _listBox1->GetSelected() << endl;
		cout     << "   _listBox2->GetSelected() = " << _listBox2->GetSelected() << endl;
		cout     << "    _waveNames.size()       = " << _waveNames.size()        << endl;
		return;
	}

	// Produce plots
	unsigned int nFitResults = _fitResults.size();

	Int_t colour=1;

	TGraphErrors* gph = new TGraphErrors(nFitResults);
	stringstream graphName;
	graphName << "PHI"<<w1<<"---"<<"PHI"<<w2;
	if(_debug) {
		printDebug << "creating graph   " << graphName.str() << endl;
	}
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
	if(_debug) {
		printDebug << "creating graph   " << graphName.str() << endl;
	}
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
	if(_debug) {
		printDebug << "creating graph   " << graphName.str() << endl;
	}
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

		int wi1 = result->waveIndex(w1);
		int wi2 = result->waveIndex(w2);
		double intensity1 = 0.;
		double intensity2 = 0.;
		double intensityErr1 = 0.;
		double intensityErr2 = 0.;
		double ph = 0.;
		double phErr = 0.;
		complex<double> rho(0., 0.);
		TMatrixT<double> rhoCov(2, 2);
		if(wi1 >= 0) {
			intensity1 = result->intensity(w1.c_str());
			intensityErr1 = result->intensityErr(w1.c_str());
		}
		if(wi2 >= 0) {
			intensity2 = result->intensity(w2.c_str());
			intensityErr2 = result->intensityErr(w2.c_str());
		}
		if(wi1 >= 0 and wi2 >= 0) {
			ph = result->phase(w1.c_str(), w2.c_str());
			phErr = result->phaseErr(w1.c_str(), w2.c_str());
			rho = result->spinDensityMatrixElem(wi1, wi2);
			rhoCov = result->spinDensityMatrixElemCov(wi1, wi2);
		}

		if((numeric_limits<double>::has_infinity and intensity1 == numeric_limits<double>::infinity()) or intensity1!=intensity1) {
			continue;
		}

		g1->SetPoint(i, result->massBinCenter()*0.001, intensity1);
		g1->SetPointError(i, _binWidth*0.5, intensityErr1);

		if((numeric_limits<double>::has_infinity and intensity2 == numeric_limits<double>::infinity()) or intensity2!=intensity2) {
			continue;
		}

		g2->SetPoint(i, result->massBinCenter()*0.001, intensity2);
		g2->SetPointError(i, _binWidth*0.5, intensityErr2);

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
		gph->SetPointError(i, _binWidth*0.5, phErr);

		// add point +- 360 degree
		gphP1->SetPoint(i, result->massBinCenter()*0.001, ph+360);
		gphP1->SetPointError(i, _binWidth*0.5, phErr);

		gphM1->SetPoint(i, result->massBinCenter()*0.001, ph-360);
		gphM1->SetPointError(i, _binWidth*0.5, phErr);

		gRe->SetPoint(i, result->massBinCenter()*0.001, rho.real());
		gRe->SetPointError(i, _binWidth*0.5, sqrt(rhoCov[0][0]));
		gIm->SetPoint(i, result->massBinCenter()*0.001, rho.imag());
		gIm->SetPointError(i, _binWidth*0.5, sqrt(rhoCov[1][1]));
	}// end loop over bins

	// plot graphs
	if(not _currentCanvas or _checkDrawNewCanvas->IsOn()) {
		if(_currentCanvas) {
			_currentCanvas->Disconnect("Closed()", this, "ActiveCanvasClosed()");
		}
		stringstream sstr;
		sstr << "plotGui_c" << _canvasCounter++;
		_currentCanvas= new TCanvas(sstr.str().c_str(), sstr.str().c_str(), 10, 10, 1200, 800);
		_currentCanvas->Connect("Closed()", "plotGuiMainFrame", this, "ActiveCanvasClosed()");
	} else {
		_currentCanvas->Clear();
	}
	_currentCanvas->Divide(2,3);
	_currentCanvas->cd(1);

	gph->Draw("AP");
	gph->GetXaxis()->SetTitle("5#pi mass (GeV/c^2)");
	gph->GetYaxis()->SetTitle("Phase difference");
	gph->GetYaxis()->SetRangeUser(-270,270);

	gphP1->Draw("PSAME");
	gphM1->Draw("PSAME");

	_currentCanvas->cd(3);
	g1->Draw("AP");
	g1->GetXaxis()->SetTitle("5#pi mass (GeV/c^2)");
	g1->GetYaxis()->SetTitle("Intensity");

	_currentCanvas->cd(5);
	g2->Draw("AP");
	g2->GetXaxis()->SetTitle("5#pi mass (GeV/c^2)");
	g2->GetYaxis()->SetTitle("Intensity");

	_currentCanvas->cd(2);
	gRe->Draw("AP");
	gRe->GetXaxis()->SetTitle("5#pi mass (GeV/c^2)");
	gRe->GetYaxis()->SetTitle("Re(#rho_{ij})");

	_currentCanvas->cd(4);
	gIm->Draw("AP");
	gIm->GetXaxis()->SetTitle("5#pi mass (GeV/c^2)");
	gIm->GetYaxis()->SetTitle("Im(#rho_{ij})");

	_currentCanvas->Draw();
	_currentCanvas->Update();

}

void plotGui(const std::string& inFileName,
             const double& binWidth = 0.03,
             const double& intensityThreshold = -1.,
             const string& whiteListFileName = "")
{

	printInfo << "opening file '" << inFileName << "'." << endl;
	TFile* infile=TFile::Open(inFileName.c_str(),"READ");
	if(infile==NULL){
		printErr << "file '" << inFileName << " not found. Aborting..."<< endl;
		return;
	}
	TTree* pwa=(TTree*)infile->FindObjectAny("pwa");
	if(pwa==NULL){
		printErr << "tree not found in input file '" << inFileName << "'. Aborting..."<< endl;
		return;
	}

	// Popup the GUI...
	new plotGuiMainFrame(gClient->GetRoot(), 20, 20, pwa, binWidth, intensityThreshold, whiteListFileName);
}
