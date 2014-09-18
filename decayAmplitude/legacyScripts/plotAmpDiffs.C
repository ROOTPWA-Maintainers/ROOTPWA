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
//
// Description:
//      plots differences between amplitude files
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <string>
#include <set>

#include <boost/progress.hpp>

#include "TFile.h"
#include "TChain.h"
#include "TTreeFormula.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"

#include "reportingUtils.hpp"


using namespace std;
using namespace boost;
using namespace rpwa;


bool
plotAmpDiffs(const string&  inFileNamePattern,
             const string&  outFileName  = "ampDiffPlots.root",
             const long int maxNmbEvents = -1,
             const string&  inTreeName   = "ampDiffTree")
{
	// open input file(s)
	printInfo << "opening input file(s) '" << inFileNamePattern << "'" << endl;
	TChain inTree(inTreeName.c_str());
	if (inTree.Add(inFileNamePattern.c_str()) < 1) {
		printWarn << "no events in input file(s) '" << inFileNamePattern << "'" << endl;
		return false;
	}
	inTree.GetListOfFiles()->ls();

	// open output file
	TFile* outFile = 0;
	if (outFileName != "") {
		printInfo << "writing difference plots to '" << outFileName << "'" << endl;
		outFile = TFile::Open(outFileName.c_str(), "RECREATE");
	}

	// connect tree leafs
	TObjString* ampName = 0;
	UInt_t      eventNmb;
	UInt_t      massBinMin, massBinMax;  // [MeV/c^2]
	Double_t    valReal[2], valImag[2];
	Double_t    absDiffReal, absDiffImag;
	Double_t    relDiffReal, relDiffImag;
	inTree.SetBranchAddress("ampName",     &ampName);
	inTree.SetBranchAddress("eventNmb",    &eventNmb);
	inTree.SetBranchAddress("massBinMin",  &massBinMin);
	inTree.SetBranchAddress("massBinMax",  &massBinMax);
	inTree.SetBranchAddress("valReal",     valReal);
	inTree.SetBranchAddress("valImag",     valImag);
	inTree.SetBranchAddress("absDiffReal", &absDiffReal);
	inTree.SetBranchAddress("absDiffImag", &absDiffImag);
	inTree.SetBranchAddress("relDiffReal", &relDiffReal);
	inTree.SetBranchAddress("relDiffImag", &relDiffImag);

	// book histograms
	const long int nmbEventsTree = inTree.GetEntries();
	const long int nmbEvents     = ((maxNmbEvents > 0) ? min(maxNmbEvents, nmbEventsTree)
	                                : nmbEventsTree);
	// TH1F* hAmps1Re   = new TH1F("amps1Real",   "amps1Real;Event Number;#Rgothic[Amp 1]", nmbEvents, -0.5, nmbEvents - 0.5);
	// TH1F* hAmps1Im   = new TH1F("amps1Imag",   "amps1Imag;Event Number;#Jgothic[Amp 1]", nmbEvents, -0.5, nmbEvents - 0.5);
	// TH1F* hAmps2Re   = new TH1F("amps2Real",   "amps2Real;Event Number;#Rgothic[Amp 2]", nmbEvents, -0.5, nmbEvents - 0.5);
	// TH1F* hAmps2Im   = new TH1F("amps2Imag",   "amps2Imag;Event Number;#Jgothic[Amp 2]", nmbEvents, -0.5, nmbEvents - 0.5);
	TH1F* hAbsDiffRe = new TH1F("absDiffReal", "absDiffReal;#Rgothic[Amp 1 - Amp 2];Count", 1000000, -1e-9, 1e-9);
	TH1F* hAbsDiffIm = new TH1F("absDiffImag", "absDiffImag;#Jgothic[Amp 1 - Amp 2];Count", 1000000, -1e-8, 1e-8);
	TH1F* hRelDiffRe = new TH1F("relDiffReal", "relDiffReal;(#Rgothic[Ampl 1] - #Rgothic[Ampl 2]) / #Rgothic[Ampl 1];Count", 1000000, -1e-4, 1e-4);
	TH1F* hRelDiffIm = new TH1F("relDiffImag", "relDiffImag;(#Jgothic[Ampl 1] - #Jgothic[Ampl 2]) / #Jgothic[Ampl 1];Count", 1000000, -1e-4, 1e-4);
	// TH2F* hCorrRe    = new TH2F("corrReal",    "corrReal;#Rgothic[Amp 1];#Rgothic[Amp 2]", 5000, -5, 5, 5000, -5, 5);
	// TH2F* hCorrIm    = new TH2F("corrImag",    "corrImag;#Jgothic[Amp 1];#Jgothic[Amp 2]", 5000, -5, 5, 5000, -5, 5);

	// loop over tree
	set<string>      ampNames;
	double           maxAbsDiff      = 0;
	long int         maxAbsDiffIndex = -1;
	double           maxRelDiff      = 0;
	long int         maxRelDiffIndex = -1;
	progress_display progressIndicator(nmbEvents, cout, "");
	for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
		++progressIndicator;

		if (inTree.LoadTree(eventIndex) < 0)
			break;
		inTree.GetEntry(eventIndex);

		ampNames.insert(ampName->GetString().Data());

		// hAmps1Re->SetBinContent(eventIndex + 1, valReal[0]);
		// hAmps1Im->SetBinContent(eventIndex + 1, valImag[0]);
		// hAmps2Re->SetBinContent(eventIndex + 1, valReal[1]);
		// hAmps2Im->SetBinContent(eventIndex + 1, valImag[1]);
		hAbsDiffRe->Fill(absDiffReal);
		hAbsDiffIm->Fill(absDiffImag);
		hRelDiffRe->Fill(relDiffReal);
		hRelDiffIm->Fill(relDiffImag);
		// hCorrRe->Fill(valReal[0], valReal[1]);
		// hCorrIm->Fill(valImag[0], valImag[1]);
		// compute maximum deviations
		const double absMax = max(fabs(absDiffReal), fabs(absDiffImag));
		const double relMax = max(fabs(relDiffReal), fabs(relDiffImag));
		if (absMax > maxAbsDiff) {
			maxAbsDiff      = absMax;
			maxAbsDiffIndex = eventIndex;
		}
		if (relMax > maxRelDiff) {
			maxRelDiff      = relMax;
			maxRelDiffIndex = eventIndex;
		}
	}
	printInfo << "maximum observed deviations: absolute = " << maxPrecision(maxAbsDiff) << " "
	          << "(event " << maxAbsDiffIndex << "), relative = " << maxPrecision(maxRelDiff) << " "
	          << "(event " << maxRelDiffIndex << ")" << endl;

	// print events with max deviation
	inTree.GetEntry(maxAbsDiffIndex);
	printInfo << "event[" << maxAbsDiffIndex << "]:"                 << endl
	          << "    ampName ....... " << ampName->GetString()      << endl
	          << "    eventNmb ...... " << eventNmb                  << endl
	          << "    massBinMin .... " << massBinMin                << endl
	          << "    massBinMax .... " << massBinMax                << endl
	          << "    valReal[0] .... " << maxPrecision(valReal[0] ) << endl
	          << "    valReal[1] .... " << maxPrecision(valReal[1] ) << endl
	          << "    valImag[0] .... " << maxPrecision(valImag[0] ) << endl
	          << "    valImag[1] .... " << maxPrecision(valImag[1] ) << endl
	          << "    absDiffReal ... " << maxPrecision(absDiffReal) << endl
	          << "    absDiffImag ... " << maxPrecision(absDiffImag) << endl
	          << "    relDiffReal ... " << maxPrecision(relDiffReal) << endl
	          << "    relDiffImag ... " << maxPrecision(relDiffImag) << endl;
	inTree.GetEntry(maxRelDiffIndex);
	printInfo << "event[" << maxRelDiffIndex << "]:"                 << endl
	          << "    ampName ....... " << ampName->GetString()      << endl
	          << "    eventNmb ...... " << eventNmb                  << endl
	          << "    massBinMin .... " << massBinMin                << endl
	          << "    massBinMax .... " << massBinMax                << endl
	          << "    valReal[0] .... " << maxPrecision(valReal[0] ) << endl
	          << "    valReal[1] .... " << maxPrecision(valReal[1] ) << endl
	          << "    valImag[0] .... " << maxPrecision(valImag[0] ) << endl
	          << "    valImag[1] .... " << maxPrecision(valImag[1] ) << endl
	          << "    absDiffReal ... " << maxPrecision(absDiffReal) << endl
	          << "    absDiffImag ... " << maxPrecision(absDiffImag) << endl
	          << "    relDiffReal ... " << maxPrecision(relDiffReal) << endl
	          << "    relDiffImag ... " << maxPrecision(relDiffImag) << endl;

	const string leafsToDraw[] = {"absDiffReal",
	                              "absDiffImag",
	                              "relDiffReal",
	                              "relDiffImag"};
	const unsigned int nmbLeafs = sizeof(leafsToDraw) / sizeof(leafsToDraw[0]);
	TH2F* hLeafs[nmbLeafs];
	for (unsigned int i = 0; i < nmbLeafs; ++i) {
		printInfo << "filling histogram for leaf '" << leafsToDraw[i] << "'" << endl;
		const double min = inTree.GetMinimum(leafsToDraw[i].c_str());
		const double max = inTree.GetMaximum(leafsToDraw[i].c_str());
		hLeafs[i] = new TH2F((leafsToDraw[i] + "Wave").c_str(),
		                     (leafsToDraw[i] + " vs. Wave;Wave;" + leafsToDraw[i]).c_str(),
		                     1, 0, 1, 100000, min - 0.1 * fabs(min), max + 0.1 * fabs(max));
		hLeafs[i]->SetBit(TH1::kCanRebin);
		TTreeFormula*    leafVal = new TTreeFormula("", leafsToDraw[i].c_str(), &inTree);
		progress_display progressIndicator(nmbEvents, cout, "");
		for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
			++progressIndicator;
			if (inTree.LoadTree(eventIndex) < 0)
				break;
			inTree.GetEntry(eventIndex);
			leafVal->UpdateFormulaLeaves();
			hLeafs[i]->Fill(ampName->GetString().Data(), leafVal->EvalInstance(), 1.);
		}
		hLeafs[i]->LabelsDeflate("X");
		hLeafs[i]->SetLabelSize(0.02, "X");
		hLeafs[i]->SetDrawOption("COLZ");
	}

	if (outFile) {
		outFile->Write();
		outFile->Close();
		delete outFile;
		printInfo << "wrote difference plots to '" << outFileName << "'" << endl;
	}
	return true;
}
