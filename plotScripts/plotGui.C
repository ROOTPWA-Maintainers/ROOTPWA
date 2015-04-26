
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
#include <TMultiGraph.h>
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
		std::vector<std::vector<const rpwa::fitResult*> > _fitResults;
		const double _binWidth;
		unsigned int _canvasCounter;
		TCanvas* _currentCanvas;
		vector<int> _usefulColors;

	public:
		plotGuiMainFrame(const TGWindow *p,
		            UInt_t w,
		            UInt_t h,
		            vector<TTree*> trees,
		            const double& binWidth,
		            const double& intensityThreshold,
		            const string& waveListFileName);
		virtual ~plotGuiMainFrame();
		void DoExit();
		void DoSelect();
		void ActiveCanvasClosed();
		void HandleButtons();
		void ListBoxSelectionChangedByMouse();
		void ListBox1SelectionChangedByKeyboard(TGFrame* frame);
		void ListBox2SelectionChangedByKeyboard(TGFrame* frame);
		void PrintSelected(const int sel1 = -1, const int sel2 = -1);

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
                         vector<TTree*> trees,
                         const double& binWidth,
                         const double& intensityThreshold,
                         const string& waveListFileName) :
	TGMainFrame(p, w, h),
	_fitResults(trees.size(), std::vector<const rpwa::fitResult*>()),
	_binWidth(binWidth),
	_canvasCounter(0),
	_currentCanvas(0)
{

	for(unsigned int treeIndex = 0; treeIndex < trees.size(); ++treeIndex)
	{
		printInfo << "reading tree " << treeIndex << " of " << trees.size() << "..." << endl;
		fitResult* res = 0;//new fitResult();
		string branchname = "fitResult_v2";
		if(trees[treeIndex]->FindBranch(branchname.c_str()) == 0){
			printErr << "invalid branch '" << branchname << "'. Aborting..." << endl;
			throw;
		}
		trees[treeIndex]->SetBranchAddress(branchname.c_str(),&res);
		long entries = trees[treeIndex]->GetEntries();
		if(entries < 1) {
			printErr << "tree has no entries. Aborting..." << endl;
			throw;
		}
		printInfo << "Mass bins: " << entries << endl;
		_fitResults[treeIndex].resize(entries, 0);
		for(long i = 0; i < entries; ++i) {
			printInfo << "loading result " << i << " of " << entries << "... " << std::flush;
			trees[treeIndex]->GetEntry(i);
			_fitResults[treeIndex][i] = new rpwa::fitResult(*res);
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
		for(unsigned int treeIndex = 0; treeIndex < _fitResults.size(); ++treeIndex) {
			for(unsigned int i = 0; i < _fitResults[treeIndex].size(); ++i) {
				const vector<string>& waveNames = _fitResults[treeIndex][i]->waveNames();
				waveNameCollector.insert(waveNames.begin(), waveNames.end());
			}
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
				for(unsigned int treeIndex = 0; treeIndex < _fitResults.size(); ++treeIndex) {
					for(unsigned int bin_i = 0; bin_i < _fitResults[treeIndex].size(); ++bin_i) {
						const double intensity = _fitResults[treeIndex][bin_i]->intensity(waveName.c_str());
						if(intensity >= intensityThreshold) {
							whitelistedWaves.insert(waveName);
							break;
						}
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

	_listBox1->Connect("Selected(Int_t)", "plotGuiMainFrame", this, "ListBoxSelectionChangedByMouse()");
	_listBox2->Connect("Selected(Int_t)", "plotGuiMainFrame", this, "ListBoxSelectionChangedByMouse()");
	((TGLBContainer *)_listBox1->GetContainer())->Connect("CurrentChanged(TGFrame*)", "plotGuiMainFrame", this, "ListBox1SelectionChangedByKeyboard(TGFrame*)");
	((TGLBContainer *)_listBox2->GetContainer())->Connect("CurrentChanged(TGFrame*)", "plotGuiMainFrame", this, "ListBox2SelectionChangedByKeyboard(TGFrame*)");

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

	vector<int> badColors;
	int badColorsArray[12] = {5, 10, 18, 19, 21, 22, 23, 25, 26, 34, 37, 31};
	badColors.assign(badColorsArray, badColorsArray+12);
	for(unsigned int i = 1; i < 41; ++i) {
		if(std::find(badColors.begin(), badColors.end(), i)==badColors.end()) {
			_usefulColors.push_back(i);
		}
	}

}


plotGuiMainFrame::~plotGuiMainFrame()
{
	// Clean up main frame...
	Cleanup();
}


void plotGuiMainFrame::HandleButtons() { }


void plotGuiMainFrame::ListBox1SelectionChangedByKeyboard(TGFrame* frame) {
	if(frame) {
		PrintSelected(_listBox1->FindEntry(((TGTextLBEntry *)frame)->GetText()->Data())->EntryId());
	}
}


void plotGuiMainFrame::ListBox2SelectionChangedByKeyboard(TGFrame* frame) {
	if(frame) {
		PrintSelected(-1, _listBox2->FindEntry(((TGTextLBEntry *)frame)->GetText()->Data())->EntryId());
	}
}


void plotGuiMainFrame::ListBoxSelectionChangedByMouse() {
	PrintSelected();
}


void plotGuiMainFrame::ActiveCanvasClosed() {
	if(_debug) {
		printDebug << "canvas closed" << endl;
	}
	_currentCanvas = 0;
}


void plotGuiMainFrame::PrintSelected(const int sel1, const int sel2)
{
	// Writes selected entries in TList if multiselection.

	string w1 = "";
	string w2 = "";
	{
		int selectionList1 = sel1 < 0 ? _listBox1->GetSelected() : sel1;
		int selectionList2 = sel2 < 0 ? _listBox2->GetSelected() : sel2;
		try {
			w1 = _waveNames.at(selectionList1);
			w2 = _waveNames.at(selectionList2);
		} catch (std::out_of_range&) {
			printErr << "selection lists seem to be corrupted." << endl;
			cout     << "   _listBox1->GetSelected() = " << selectionList1 << endl;
			cout     << "   _listBox2->GetSelected() = " << selectionList2 << endl;
			cout     << "    _waveNames.size()       = " << _waveNames.size()        << endl;
			return;
		}
		if(_debug) {
			printDebug << "PrintSelected called with list box entry ids (" << selectionList1 << ", " << selectionList2 << ")." << endl;
		}
	}

	stringstream sstr;
	sstr << "g" << w1;
	const string intensity1GraphName = sstr.str();
	sstr.str("");
	sstr << "g" << w2;
	const string intensity2GraphName = sstr.str();
	sstr.str("");
	sstr << "PHI"<< w1 << "---" << "PHI" << w2;
	const string phaseGraphName = sstr.str();
	sstr.str("");
	sstr << "RE_" << w1 << "---" << w2;
	const string reGraphName = sstr.str();
	sstr.str("");
	sstr << "IM_" << w1 << "---" << w2;
	const string imGraphName = sstr.str();

	TMultiGraph* intensity1Graph = new TMultiGraph();
	TMultiGraph* intensity2Graph = new TMultiGraph();
	TMultiGraph* phaseGraph = new TMultiGraph();
	TMultiGraph* reGraph = new TMultiGraph();
	TMultiGraph* imGraph = new TMultiGraph();

	intensity1Graph->SetName(intensity1GraphName.c_str());
	intensity1Graph->SetTitle(intensity1GraphName.c_str());
	intensity2Graph->SetName(intensity2GraphName.c_str());
	intensity2Graph->SetTitle(intensity2GraphName.c_str());
	phaseGraph->SetName(phaseGraphName.c_str());
	phaseGraph->SetTitle(phaseGraphName.c_str());
	reGraph->SetName(reGraphName.c_str());
	reGraph->SetTitle(reGraphName.c_str());
	imGraph->SetName(imGraphName.c_str());
	imGraph->SetTitle(imGraphName.c_str());

	for(unsigned int treeIndex = 0; treeIndex < _fitResults.size(); ++treeIndex) {

		const unsigned int& nPoints = _fitResults[treeIndex].size();

		TGraphErrors* gph = new TGraphErrors(nPoints);
		stringstream graphName;
		graphName << phaseGraphName << "_" << treeIndex;
		if(_debug) {
			printDebug << "creating graph   " << graphName.str() << endl;
		}
		gph->SetName (graphName.str().c_str());
		gph->SetTitle(graphName.str().c_str());
		gph->SetMarkerStyle(21);
		gph->SetMarkerSize(0.5);
		gph->SetMarkerColor(_usefulColors[treeIndex % _usefulColors.size()]);
		gph->SetLineColor(_usefulColors[treeIndex % _usefulColors.size()]);

		TGraphErrors* gphP1 = (TGraphErrors*)gph->Clone("gph+1");
		TGraphErrors* gphM1 = (TGraphErrors*)gph->Clone("gph-1");
		if(_fitResults.size() == 1) {
			gphP1->SetMarkerColor(2);
			gphP1->SetLineColor(2);
			gphM1->SetMarkerColor(3);
			gphM1->SetLineColor(3);
		} else {
			gphP1->SetMarkerColor(_usefulColors[treeIndex % _usefulColors.size()]);
			gphP1->SetLineColor(_usefulColors[treeIndex % _usefulColors.size()]);
			gphM1->SetMarkerColor(_usefulColors[treeIndex % _usefulColors.size()]);
			gphM1->SetLineColor(_usefulColors[treeIndex % _usefulColors.size()]);
		}

		TGraphErrors* gRe = new TGraphErrors(nPoints);
		graphName.str("");
		graphName << reGraphName << "_" << treeIndex;
		if(_debug) {
			printDebug << "creating graph   " << graphName.str() << endl;
		}
		gRe->SetName (graphName.str().c_str());
		gRe->SetTitle(graphName.str().c_str());
		gRe->SetMarkerStyle(21);
		gRe->SetMarkerSize(0.5);
		gRe->SetMarkerColor(_usefulColors[treeIndex % _usefulColors.size()]);
		gRe->SetLineColor(_usefulColors[treeIndex % _usefulColors.size()]);

		TGraphErrors* gIm = new TGraphErrors(nPoints);
		graphName.str("");
		graphName << imGraphName << "_" << treeIndex;
		if(_debug) {
			printDebug << "creating graph   " << graphName.str() << endl;
		}
		gIm->SetName (graphName.str().c_str());
		gIm->SetTitle(graphName.str().c_str());
		gIm->SetMarkerStyle(21);
		gIm->SetMarkerSize(0.5);
		gIm->SetMarkerColor(_usefulColors[treeIndex % _usefulColors.size()]);
		gIm->SetLineColor(_usefulColors[treeIndex % _usefulColors.size()]);

		TGraphErrors* g1 = new TGraphErrors(nPoints);
		graphName.str("");
		graphName << intensity1GraphName << "_" << treeIndex;
		g1->SetName (graphName.str().c_str());
		g1->SetTitle(graphName.str().c_str());
		g1->SetMarkerStyle(21);
		g1->SetMarkerSize(0.5);
		g1->SetMarkerColor(_usefulColors[treeIndex % _usefulColors.size()]);
		g1->SetLineColor(_usefulColors[treeIndex % _usefulColors.size()]);
		TGraphErrors* g2 = new TGraphErrors(nPoints);
		graphName.str("");
		graphName << intensity1GraphName << "_" << treeIndex;
		g2->SetName (graphName.str().c_str());
		g2->SetTitle(graphName.str().c_str());
		g2->SetMarkerStyle(21);
		g2->SetMarkerSize(0.5);
		g2->SetMarkerColor(_usefulColors[treeIndex % _usefulColors.size()]);
		g2->SetLineColor(_usefulColors[treeIndex % _usefulColors.size()]);

		for(unsigned int i = 0; i < nPoints; ++i){
			const fitResult* result = _fitResults[treeIndex][i];

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

		intensity1Graph->Add(g1);
		intensity2Graph->Add(g2);
		phaseGraph->Add(gph);
		phaseGraph->Add(gphP1);
		phaseGraph->Add(gphM1);
		reGraph->Add(gRe);
		imGraph->Add(gIm);

	}

	// plot graphs
	if(not _currentCanvas or _checkDrawNewCanvas->IsOn()) {
		if(_currentCanvas) {
			_currentCanvas->Disconnect("Closed()", this, "ActiveCanvasClosed()");
		}
		stringstream sstr;
		sstr << "plotGui_c" << _canvasCounter++;
		_currentCanvas= new TCanvas(sstr.str().c_str(), sstr.str().c_str(), 10, 10, 1600, 900);
		_currentCanvas->Connect("Closed()", "plotGuiMainFrame", this, "ActiveCanvasClosed()");
	} else {
		_currentCanvas->Clear();
	}
	_currentCanvas->Divide(2,2);

	_currentCanvas->cd(3);
	phaseGraph->Draw("AP");
	phaseGraph->GetXaxis()->SetTitle("5#pi mass (GeV/c^2)");
	phaseGraph->GetYaxis()->SetTitle("Phase difference");
	phaseGraph->GetYaxis()->SetRangeUser(-270,270);

	_currentCanvas->cd(1);
	intensity1Graph->Draw("AP");
	intensity1Graph->GetXaxis()->SetTitle("5#pi mass (GeV/c^2)");
	intensity1Graph->GetYaxis()->SetTitle("Intensity");

	_currentCanvas->cd(2);
	intensity2Graph->Draw("AP");
	intensity2Graph->GetXaxis()->SetTitle("5#pi mass (GeV/c^2)");
	intensity2Graph->GetYaxis()->SetTitle("Intensity");

	_currentCanvas->Draw();
	_currentCanvas->Update();

}

void plotGui(const vector<std::string>& inFileNames,
             const double& binWidth = 0.03,
             const double& intensityThreshold = -1.,
             const string& whiteListFileName = "")
{

	printInfo << "opening file '" << inFileNames << "'." << endl;
	vector<TTree*> trees(inFileNames.size(), 0);
	for(unsigned int i = 0; i < inFileNames.size(); ++i) {
		TFile* infile=TFile::Open(inFileNames[i].c_str(),"READ");
		if(not infile){
			printErr << "file '" << inFileNames << " not found. Aborting..."<< endl;
			return;
		}
		TTree* pwa=(TTree*)infile->FindObjectAny("pwa");
		if(not pwa){
			printErr << "tree not found in input file '" << inFileNames << "'. Aborting..."<< endl;
			return;
		}
		trees[i] = pwa;
	}

	// Popup the GUI...
	new plotGuiMainFrame(gClient->GetRoot(), 20, 20, trees, binWidth, intensityThreshold, whiteListFileName);
}

void plotGui(const std::string& inFileName,
             const double& binWidth = 0.03,
             const double& intensityThreshold = -1.,
             const string& whiteListFileName = "")
{
	const std::vector<string> fileNames(1, inFileName);
	plotGui(fileNames, binWidth, intensityThreshold, whiteListFileName);
}
