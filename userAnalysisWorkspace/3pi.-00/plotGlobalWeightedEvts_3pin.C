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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <sstream>

#include <TCanvas.h>
#include <TClass.h>
#include <TColor.h>
#include <TExec.h>
#include <TFile.h>
#include <TFrame.h>
#include <TGaxis.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>

const unsigned int NUMBER_CONTOURS = 31;

class MassBin {
public:
	MassBin(TDirectory* d, const double mc, const double ml, const double mu) : dir (d), massBinCenter(mc), massBinLower(ml), massBinUpper(mu) {}

	TDirectory* GetDir() const { return dir; }
	double GetMassBinCenter() const { return massBinCenter; }
	double GetMassBinLower() const { return massBinLower; }
	double GetMassBinUpper() const { return massBinUpper; }

	bool operator<(const MassBin& rhs) const { return massBinCenter < rhs.massBinCenter; }

private:
	TDirectory* dir;
	double massBinCenter;
	double massBinLower;
	double massBinUpper;
};

class BookyDefinition {
public:
	BookyDefinition(const std::string& b) : bookyName(b) {
		canvas = new TCanvas(bookyName.c_str(), bookyName.c_str());
		canvas->Print((bookyName + ".ps[").c_str());
	}

	BookyDefinition(const std::string& b, const unsigned int w, const unsigned int h) : bookyName(b) {
		canvas = new TCanvas(bookyName.c_str(), bookyName.c_str(), w, h);
		canvas->Print((bookyName + ".ps[").c_str());
	}

	virtual ~BookyDefinition() {
		canvas->Print((bookyName + ".ps]").c_str());
		delete canvas;
	}

	virtual void Process(const MassBin& massBin, const bool twoMc) = 0;
	virtual std::vector<TCanvas*> Finalize() = 0;

	TCanvas* GetCanvas() const { return canvas; }
	const std::string& GetName() const { return bookyName; }

private:
	std::string bookyName;
	TCanvas* canvas;
};

class Default2DComparison : public BookyDefinition {
public:
	Default2DComparison(const std::string& b, const std::string& t) : BookyDefinition(b), templateName(t) {
	}

	~Default2DComparison() {
		for (std::vector<TH1*>::iterator it=clearHistos.begin(); it!=clearHistos.end(); ++it)
			delete *it;
	}

	void Process(const MassBin& massBin, const bool twoMc) {
		TCanvas* c = GetCanvas();
		c->Clear();
		c->Divide(2, 2);

		size_t pos = templateName.find("Mc");
		const char* replacement[4] = { "Mc", "Data", "Diff", "RelDiff" };
		TH1* hist[4];

		for (unsigned int i=0; i<2; i++) {
			std::string histName(templateName);
			size_t len = (i==1&&twoMc) ? 5 : 2;
			histName.replace(pos, len, replacement[i]);

			massBin.GetDir()->GetObject(histName.c_str(), hist[i]);

			if (hist[i] == NULL) {
				std::cerr << "Could not find histogram '" << histName << "' in '" << massBin.GetDir()->GetName() << "'. The bookies created might be wrong." << std::endl;
			}
		}

		if (hist[0] == NULL || hist[1] == NULL) {
			return;
		}

		hist[0]->SetMaximum();
		hist[0]->SetMinimum();
		hist[1]->SetMaximum();
		hist[1]->SetMinimum();
		const double max = std::max(hist[0]->GetMaximum(), hist[1]->GetMaximum());
		hist[0]->SetMaximum(max);
		hist[1]->SetMaximum(max);

		// skip mass bins without any entries
		if (max == 0.) {
			return;
		}

		for (unsigned int i=2; i<4; i++) {
			std::string histName(templateName);
			size_t len = (i==1&&twoMc) ? 5 : 2;
			histName.replace(pos, len, replacement[i]);

			hist[i] = dynamic_cast<TH1*>(hist[0]->Clone(histName.c_str()));
			clearHistos.push_back(hist[i]);
		}

		hist[2]->Add(hist[1], -1.);

		hist[3]->Add(hist[1], -1.);
		hist[3]->Divide(hist[1]);

		const double max2 = std::max(std::abs(hist[2]->GetMaximum()), std::abs(hist[2]->GetMinimum()));
		hist[2]->SetMaximum( 1.1 * max2);
		hist[2]->SetMinimum(-1.1 * max2);

		// force histogram with relative differences to a range
		hist[3]->SetMaximum( 1.);
		hist[3]->SetMinimum(-1.);

		for (unsigned int i=0; i<4; i++) {
			c->cd(i+1);

			hist[i]->SetStats(false);

			hist[i]->SetContour(NUMBER_CONTOURS);
			hist[i]->Draw("AXIS");

			TExec* ex(NULL);
			if (i < 2) {
				ex = new TExec("ex", "gStyle->SetPalette(1);");
			} else {
				std::ostringstream os;
				os << "setDiffColorStyle(" << NUMBER_CONTOURS << ");";
				ex = new TExec("ex", os.str().c_str());
			}
			ex->Draw();

			hist[i]->Draw("COLZ SAME");

			c->cd(i+1)->RedrawAxis();
		}

		std::ostringstream os;
		os << "#scale[0.5]{Mass Bin [" << massBin.GetMassBinLower() << " - " << massBin.GetMassBinUpper() << "] MeV/c^{2}}";
		std::string labelText = os.str();

		c->cd();
		TLatex label(0.5, 1.0, labelText.c_str());
		label.SetNDC(true);
		label.SetTextAlign(23);
		label.Draw();

		c->Print((GetName() + ".ps").c_str());
	}

	std::vector<TCanvas*> Finalize() {
		std::vector<TCanvas*> ret;
		return ret;
	}

private:
	std::string templateName;

	std::vector<TH1*> clearHistos;
};

class AnglesComparison : public BookyDefinition {
public:
	AnglesComparison(const std::string& b, const std::string& s) : BookyDefinition(b, 480, 640), suffix(s) {
	}

	~AnglesComparison() {
		for (std::vector<TH1*>::iterator it=clearHistos.begin(); it!=clearHistos.end(); ++it)
			delete *it;
	}

	void Process(const MassBin& massBin, const bool twoMc) {
		TCanvas* c = GetCanvas();
		c->Clear();
		c->Divide(2, 3);

		const char* prefix[6] = { "hMIsobar", NULL, "hGJ", "hTY", "hHTheta", "hHPhi" };

		for (unsigned int i=0; i<6; i++) {
			c->cd(i+1);

			if (prefix[i] == NULL) {
				continue;
			}

			std::string histMcName;
			if (twoMc) {
				histMcName = GetHistogramName(prefix[i], "McAcc", suffix);
			} else {
				histMcName = GetHistogramName(prefix[i], "Mc", suffix);
			}

			TH1* histMc;
			massBin.GetDir()->GetObject(histMcName.c_str(), histMc);

			if (histMc == NULL) {
				std::cerr << "Could not find histogram '" << histMcName << "' in '" << massBin.GetDir()->GetName() << "'. The bookies created might be wrong." << std::endl;
				continue;
			}

			const std::string histDataName = GetHistogramName(prefix[i], "Data", suffix);

			TH1* histData;
			massBin.GetDir()->GetObject(histDataName.c_str(), histData);

			if (histData == NULL) {
				std::cerr << "Could not find histogram '" << histDataName << "' in '" << massBin.GetDir()->GetName() << "'. The bookies created might be wrong." << std::endl;
				continue;
			}

			double max = std::max(histData->GetMaximum(), histMc->GetMaximum());

			// skip mass bins without any entries
			if (max == 0.) {
				return;
			}

			double spacePlot = 0.7;
			double spaceDiff = 0.3;
			double spaceAcc = 0.;
			if (twoMc) {
				spacePlot = 0.6;
				spaceDiff = 0.2;
				spaceAcc = 0.2;
			}

			TH1* histFirst(NULL);
			TH1* histMcPsp(NULL);
			if (twoMc) {
				const std::string histMcPspName = GetHistogramName(prefix[i], "McPsp", suffix);
				massBin.GetDir()->GetObject(histMcPspName.c_str(), histMcPsp);

				if (histMcPsp == NULL) {
					std::cerr << "Could not find histogram '" << histMcPspName << "' in '" << massBin.GetDir()->GetName() << "'. The bookies created might be wrong." << std::endl;
					continue;
				}

				TH1* histMcPspScaled = dynamic_cast<TH1*>(histMcPsp->Clone((histMcPspName + "Scaled").c_str()));
				clearHistos.push_back(histMcPspScaled);

				double scalePsp = 1.;
				if (histMcPsp->Integral() != 0.)
					scalePsp = histData->Integral() / histMcPsp->Integral();
				histMcPspScaled->Scale(scalePsp);

				max = std::max(max, histMcPspScaled->GetMaximum());

				histMcPspScaled->SetFillColor(kBlue);
				histMcPspScaled->SetLineColor(kBlue);
				histMcPspScaled->SetMarkerColor(kBlue);
				histMcPspScaled->SetStats(false);
				histMcPspScaled->Draw("A E2");

				histFirst = histMcPspScaled;
			} else {
				histFirst = histMc;
			}

			histMc->SetFillColor(kRed);
			histMc->SetLineColor(kRed);
			histMc->SetMarkerColor(kRed);
			histMc->SetStats(false);
			if (twoMc) {
				histMc->Draw("E2 SAME");
			} else {
				histMc->Draw("A E2");
			}

			histData->SetStats(false);
			histData->Draw("SAME");

			const double boundUpper = 1.1*max;
			const double boundLower = 1.1*max * (spacePlot - 1.)/spacePlot;
			const double boundWidth = 1.1*max / spacePlot;

			histFirst->SetMaximum(boundUpper);
			histFirst->SetMinimum(boundLower);

			const std::string histDiffName = GetHistogramName(prefix[i], "Diff", suffix);
			TH1* histDiff = dynamic_cast<TH1*>(histMc->Clone(histDiffName.c_str()));
			histDiff->Add(histData, -1.);
			clearHistos.push_back(histDiff);

			histDiff->SetMaximum();
			histDiff->SetMinimum();
			const double maxDiff = std::max(std::abs(histDiff->GetMaximum()), std::abs(histDiff->GetMinimum()));

			const std::string histRelDiffName = GetHistogramName(prefix[i], "RelDiff", suffix);
			TH1* histRelDiff = dynamic_cast<TH1*>(histDiff->Clone(histRelDiffName.c_str()));
			histRelDiff->Divide(histData);
			clearHistos.push_back(histRelDiff);

			histRelDiff->SetMaximum();
			histRelDiff->SetMinimum();
			const double maxRelDiff = std::max(std::abs(histRelDiff->GetMaximum()), std::abs(histRelDiff->GetMinimum()));

			TransformHistogram(histDiff, -1.1*maxDiff, 1.1*maxDiff, boundWidth*spaceAcc + boundLower, boundWidth*(spaceAcc+spaceDiff) + boundLower);

			histDiff->SetFillColor(kWhite);
			histDiff->SetLineColor(kBlack);
			histDiff->SetStats(false);
			histDiff->Draw("SAME");

			TransformHistogram(histRelDiff, -1.1*maxRelDiff, 1.1*maxRelDiff, boundWidth*spaceAcc + boundLower, boundWidth*(spaceAcc+spaceDiff) + boundLower);

			histRelDiff->SetFillColor(kWhite);
			histRelDiff->SetLineColor(kBlue);
			histRelDiff->SetStats(false);
			histRelDiff->Draw("SAME");

			double maxAcceptance(0.);
			if (twoMc) {
				const std::string histAcceptanceName = GetHistogramName(prefix[i], "Acceptance", suffix);
				TH1* histAcceptance = dynamic_cast<TH1*>(histMc->Clone(histAcceptanceName.c_str()));
				histAcceptance->Divide(histMcPsp);
				clearHistos.push_back(histAcceptance);

				histAcceptance->SetMaximum();
				maxAcceptance = histAcceptance->GetMaximum();

				TransformHistogram(histAcceptance, 0., 1.1*maxAcceptance, boundLower, boundWidth*spaceAcc + boundLower);

				histAcceptance->SetFillColor(kWhite);
				histAcceptance->SetLineColor(kBlack);
				histAcceptance->SetStats(false);
				histAcceptance->Draw("SAME");
			}

			// draw axis
			const double xAxisMax = histFirst->GetXaxis()->GetXmax();
			const double xAxisMin = histFirst->GetXaxis()->GetXmin();
			const char* xAxisTitle = histFirst->GetXaxis()->GetTitle();

			TGaxis* xAxisBottom = new TGaxis(xAxisMin, boundLower, xAxisMax, boundLower, xAxisMin, xAxisMax, 510, "+");
			xAxisBottom->SetTitle(xAxisTitle);
			xAxisBottom->Draw();
			TGaxis* xAxisPlot = new TGaxis(xAxisMin, 0., xAxisMax, 0., xAxisMin, xAxisMax, 10, "-+SU");
			xAxisPlot->SetTickSize(0.01);
			xAxisPlot->Draw();
			if (twoMc) {
				TGaxis* xAxisDiff = new TGaxis(xAxisMin, boundLower + spaceAcc*boundWidth, xAxisMax, boundLower + spaceAcc*boundWidth, xAxisMin, xAxisMax, 10, "-+SU");
				xAxisDiff->SetTickSize(0.01);
				xAxisDiff->Draw();
			}

			const double yAxisPlotMax = 1.1*max;
			const double yAxisPlotMin = 0.;
			const char* yAxisPlotTitle = histFirst->GetYaxis()->GetTitle();
			const double yAxisDiffMax =  1.1*maxDiff;
			const double yAxisDiffMin = -1.1*maxDiff;
			const char* yAxisDiffTitle = "(M-D)";
			const double yAxisRelDiffMax =  1.1*maxRelDiff;
			const double yAxisRelDiffMin = -1.1*maxRelDiff;
			const char* yAxisRelDiffTitle = "(M-D)/D";

			TGaxis* yAxisPlot = new TGaxis(xAxisMin, 0., xAxisMin, boundUpper, yAxisPlotMin, yAxisPlotMax, 508, "-");
			yAxisPlot->SetTitle(yAxisPlotTitle);
			yAxisPlot->Draw();
			TGaxis* yAxisDiff = new TGaxis(xAxisMin, boundLower + spaceAcc*boundWidth, xAxisMin, 0., yAxisDiffMin, yAxisDiffMax, 3, "-");
			yAxisDiff->SetTitle(yAxisDiffTitle);
			yAxisDiff->Draw();
			TGaxis* yAxisRelDiff = new TGaxis(xAxisMax, boundLower + spaceAcc*boundWidth, xAxisMax, 0., yAxisRelDiffMin, yAxisRelDiffMax, 3, "L+");
			yAxisRelDiff->SetLabelColor(kBlue);
			yAxisRelDiff->SetLineColor(kBlue);
			yAxisRelDiff->SetTitle(yAxisRelDiffTitle);
			yAxisRelDiff->SetTitleColor(kBlue);
			yAxisRelDiff->Draw();
			if (twoMc) {
				const double yAxisAcceptanceMax = 1.1*maxAcceptance;
				const double yAxisAcceptanceMin = 0.;
				const char* yAxisAcceptanceTitle = "acc";

				TGaxis* yAxisAcceptance = new TGaxis(xAxisMin, boundLower, xAxisMin, boundLower + spaceAcc*boundWidth, yAxisAcceptanceMin, yAxisAcceptanceMax, 3, "-");
				yAxisAcceptance->SetTitle(yAxisAcceptanceTitle);
				yAxisAcceptance->Draw();
			}
		}

		std::ostringstream os;
		os << "#scale[0.5]{Mass Bin [" << massBin.GetMassBinLower() << " - " << massBin.GetMassBinUpper() << "] MeV/c^{2}}";
		std::string labelText = os.str();

		c->cd();
		TLatex label(0.5, 1.0, labelText.c_str());
		label.SetNDC(true);
		label.SetTextAlign(23);
		label.Draw();

		c->Print((GetName() + ".ps").c_str());
	}

	std::vector<TCanvas*> Finalize() {
		std::vector<TCanvas*> ret;
		return ret;
	}

private:
	std::string GetHistogramName(const std::string& p, const std::string& i, const std::string& s) const {
		return p + i + "_" + s;
	}

	void TransformHistogram(TH1* hist, const double oldMin, const double oldMax, const double newMin, const double newMax) const {
		for (Int_t i=1; i<=hist->GetNbinsX(); i++) {
			const double oldValue = hist->GetBinContent(i);
			const double newValue = (oldValue - oldMin) / (oldMax - oldMin) * (newMax - newMin) + newMin;
			hist->SetBinContent(i, newValue);
			//			std::cout << oldValue << " " << newValue << " " << oldMin << " " << oldMax << " " << newMin << " " << newMax << std::endl;
		}
	}

	std::string suffix;

	std::vector<TH1*> clearHistos;
};

class DiffVsMass : public BookyDefinition {
public:
	DiffVsMass(const std::string& b, TDirectory* o, const unsigned int mb, const double ml, const double mu) : BookyDefinition(b), out(o), massBins(mb), massBinsLower(ml), massBinsUpper(mu) {};

	void Process(const MassBin& massBin, const bool twoMc) {
		TIter histiter(massBin.GetDir()->GetListOfKeys());
		TKey* key;
		while ((key = dynamic_cast<TKey*>(histiter()))) {
			if (!TClass::GetClass(key->GetClassName())->InheritsFrom("TH1D"))
				continue;

			const std::string histNameMc(key->GetName());
			std::string histNameData(histNameMc);
			std::string histNameDiff(histNameMc);
			if (twoMc) {
				if ((histNameMc.length() >= 5 && histNameMc.substr(histNameMc.length()-5, 5) == "McAcc") || histNameMc.find("McAcc_") != std::string::npos) {
					size_t pos = histNameMc.find("McAcc");

					histNameData.erase(pos, 5);
					histNameData.insert(pos, "Data");

					histNameDiff.erase(pos, 5);
					histNameDiff.insert(pos, "Diff");
				} else {
					continue;
				}
			} else {
				if ((histNameMc.length() >= 2 && histNameMc.substr(histNameMc.length()-2, 2) == "Mc") || histNameMc.find("Mc_") != std::string::npos) {
					size_t pos = histNameMc.find("Mc");

					histNameData.erase(pos, 2);
					histNameData.insert(pos, "Data");

					histNameDiff.erase(pos, 2);
					histNameDiff.insert(pos, "Diff");
				} else {
					continue;
				}
			}

			TH1D* histBinMc;
			massBin.GetDir()->GetObject(histNameMc.c_str(), histBinMc);

			if (histBinMc == NULL) {
				continue;
			}

			TH1D* histBinData;
			massBin.GetDir()->GetObject(histNameData.c_str(), histBinData);

			if (histBinData == NULL) {
				continue;
			}

			const std::string n(histNameDiff + "VsMass");

			TH2D* histVsMass;
			out->GetObject(n.c_str(), histVsMass);
			if (histVsMass == NULL) {
				out->cd();
				histVsMass = new TH2D(n.c_str(), (std::string("difference of ") + histBinMc->GetTitle() + " vs. mass").c_str(),
				                      massBins, massBinsLower, massBinsUpper,
				                      histBinMc->GetNbinsX(), histBinMc->GetXaxis()->GetXmin(), histBinMc->GetXaxis()->GetXmax());
				histVsMass->SetOption("COLZ");
				histVsMass->SetXTitle("Resonance Mass [GeV/c^{2}]");
				histVsMass->SetYTitle(histBinMc->GetXaxis()->GetTitle());
			}

			for (Int_t i=1; i<=histBinMc->GetNbinsX(); i++) {
				const double x = massBin.GetMassBinCenter();
				const double y = histBinMc->GetBinCenter(i);
				const double z = histBinMc->GetBinContent(i) - histBinData->GetBinContent(i);

				histVsMass->Fill(x, y, z);
			}
		}
	}

	std::vector<TCanvas*> Finalize() {
		std::vector<TCanvas*> ret;

		TIter histiter(out->GetList());
		TObject *obj;
		while ((obj = dynamic_cast<TObject*>(histiter()))) {
			TH2D* histVsMass = dynamic_cast<TH2D*>(obj);
			if (histVsMass == NULL)
				continue;

			const std::string s(obj->GetName());
			if (s.find("VsMass") == std::string::npos) continue;

			const double max = std::max(std::abs(histVsMass->GetMaximum()), std::abs(histVsMass->GetMinimum()));
			histVsMass->SetMaximum( 1.1 * max);
			histVsMass->SetMinimum(-1.1 * max);

			histVsMass->SetStats(false);

			const std::string n("c" + s);

			TCanvas* c= new TCanvas(n.c_str(), n.c_str());
			c->cd();

			histVsMass->SetContour(NUMBER_CONTOURS);
			histVsMass->Draw("AXIS");

			std::ostringstream os;
			os << "setDiffColorStyle(" << NUMBER_CONTOURS << ");";
			TExec* ex = new TExec("ex", os.str().c_str());
			ex->Draw();

			histVsMass->Draw("COLZ SAME");

			c->cd()->RedrawAxis();

			c->Print((GetName() + ".ps").c_str());

			ret.push_back(c);
		}

		return ret;
	}

private:
	TDirectory* out;
	unsigned int massBins;
	double massBinsLower;
	double massBinsUpper;
};

std::map<std::string, BookyDefinition*>
setupBookies(const std::vector<MassBin>& massBins,
             TDirectory* outGlobal,
             const bool twoMc) {
	std::map<std::string, BookyDefinition*> definition;

	// add definitions for comparison of GJ and helicity angles for the
	// neutral and charged isobars
	definition.insert(std::pair<std::string, BookyDefinition*>("BookyNeutralIsobar", new AnglesComparison("BookyNeutralIsobar", "Neutral")));
	definition.insert(std::pair<std::string, BookyDefinition*>("BookyChargedIsobar", new AnglesComparison("BookyChargedIsobar", "Charged")));
	definition.insert(std::pair<std::string, BookyDefinition*>("BookyRhoIsobar", new AnglesComparison("BookyRhoIsobar", "ChargedRho")));

	// difference vs mass
	definition.insert(std::pair<std::string, BookyDefinition*>("BookyDiffVsMass", new DiffVsMass("BookyDiffVsMass", outGlobal, massBins.size(), massBins.front().GetMassBinLower(), massBins.back().GetMassBinUpper())));

	// loop over all mass bins
	for (std::vector<MassBin>::const_iterator massBin=massBins.begin(); massBin!=massBins.end(); massBin++) {
		TIter histiter(massBin->GetDir()->GetListOfKeys());
		TKey* key;
		while ((key = dynamic_cast<TKey*>(histiter()))) {
			if (!TClass::GetClass(key->GetClassName())->InheritsFrom("TH2"))
				continue;

			std::string s(key->GetName());
			if (twoMc) {
				if ((s.length() >= 5 && s.substr(s.length()-5, 5) == "McAcc") || s.find("McAcc_") != std::string::npos) {
					s.erase(s.find("McAcc"), 5);
				} else {
					continue;
				}
			} else {
				if ((s.length() >= 2 && s.substr(s.length()-2, 2) == "Mc") || s.find("Mc_") != std::string::npos) {
					s.erase(s.find("Mc"), 2);
				} else {
					continue;
				}
			}

			const std::string bookyName = "Booky" + s;
			if (definition.count(bookyName) == 0) {
				definition.insert(std::pair<std::string, BookyDefinition*>(bookyName, new Default2DComparison(bookyName, key->GetName())));
			}
		}
	}

	return definition;
}

void
setDiffColorStyle(const unsigned int NCont) {
	static Int_t *colors =0;
	static Bool_t initialized = kFALSE;

	Double_t Red[3]    = { 0.0, 1.0, 1.0 };
	Double_t Green[3]  = { 0.0, 1.0, 0.0 };
	Double_t Blue[3]   = { 1.0, 1.0, 0.0 };
	Double_t Length[3] = { 0.0, 0.50, 1.0 };

	if(!initialized){
		colors = new Int_t[NCont];
		Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,NCont);
		for (unsigned int i=0; i<NCont; i++) colors[i] = FI+i;
		initialized = kTRUE;
	}
	gStyle->SetPalette(NCont,colors);
}

bool
fileContainsTwoMc(TFile* in) {
	// guess if we are running on data that has been generated from one or
	// two sets of phasespace events (phasespace or accepted phasespace
	// only, or the combination of those two)
	TIter diriter(in->GetListOfKeys());
	TKey* keyO;

	bool twoMc(false);
	while ((keyO = dynamic_cast<TKey*>(diriter()))) {
		// check if keyO points to a TDirectory
		if (!TClass::GetClass(keyO->GetClassName())->InheritsFrom("TDirectory"))
			continue;

		// check if this directory is a mass bin dir
		const std::string nameDir = keyO->GetName();
		const size_t pointpos = nameDir.find(".");
		if (pointpos == 0 || pointpos == nameDir.size())
			continue;

		TDirectory* dir;
		in->GetObject(keyO->GetName(), dir);
		assert(dir != NULL);

		TIter histiter(dir->GetListOfKeys());
		TKey* keyI;
		while ((keyI = dynamic_cast<TKey*>(histiter()))) {
			const std::string s(keyI->GetName());
			if ((s.length() >= 5 && s.substr(s.length()-5, 5) == "McPsp") ||
			    (s.length() >= 5 && s.substr(s.length()-5, 5) == "McAcc") ||
			    s.find("McPsp_") != std::string::npos ||
			    s.find("McAcc_") != std::string::npos) {
				twoMc = true;
				break;
			}
		}

		if (twoMc) {
			break;
		}
	}

	return twoMc;
}

std::vector<MassBin>
getMassBins(TFile* in) {
	std::vector<MassBin> massBins;

	TIter diriter(in->GetListOfKeys());
	TKey* key;

	while ((key = dynamic_cast<TKey*>(diriter()))) {
		// check if keyO points to a TDirectory
		if (!TClass::GetClass(key->GetClassName())->InheritsFrom("TDirectory"))
			continue;

		// check if this directory is a mass bin dir
		const std::string nameDir = key->GetName();
		const size_t pointpos = nameDir.find(".");
		if (pointpos == 0 || pointpos == nameDir.size())
			continue;

		// extract upper and lower limits of massbin
		double massBinLower; std::istringstream(nameDir.substr(0, pointpos)) >> massBinLower;
		double massBinUpper; std::istringstream(nameDir.substr(pointpos+1)) >> massBinUpper;
		double massBinCenter = (massBinLower + massBinUpper) / 2.;

		TDirectory* inDir;
		in->GetObject(key->GetName(), inDir);
		assert(inDir != NULL);

		massBins.push_back(MassBin(inDir, massBinCenter, massBinLower, massBinUpper));
	}

	std::sort(massBins.begin(), massBins.end());
	return massBins;
}

void
addHistograms(TDirectory* out, const MassBin& massBin) {
	TIter histiter(massBin.GetDir()->GetListOfKeys());
	TKey* key;
	while ((key = dynamic_cast<TKey*>(histiter()))) {
		const std::string s(key->GetName());

		if (s.find("_Stack_") != std::string::npos) continue;
		if (s.find("Acceptance") != std::string::npos) continue;
		if (s.find("RelDiff") != std::string::npos) continue;

		TH1* histBin;
		massBin.GetDir()->GetObject(key->GetName(), histBin);

		TH1* histSum;
		out->GetObject(key->GetName(), histSum);
		if (histSum == NULL) {
			out->cd();
			histSum = dynamic_cast<TH1*>(histBin->Clone());
			histSum->SetMaximum();
			histSum->SetMinimum();
		} else {
			histSum->Add(histBin);
		}
	}
}

void
freeMemory(const MassBin& massBin) {
	TIter histiter(massBin.GetDir()->GetListOfKeys());
	TKey* key;
	while ((key = dynamic_cast<TKey*>(histiter()))) {
		const std::string s(key->GetName());

		TH1* histBin;
		massBin.GetDir()->GetObject(key->GetName(), histBin);
		delete histBin;
	}
}

void
plotGlobalWeightedEvts_3pin(const std::string& inFileName,
                            const std::string& outFileName) {
	TFile* inFile = TFile::Open(inFileName.c_str(), "READ");
	if (inFile == NULL) {
		std::cerr << "Input file '" << inFileName << "' could not be opened." << std::endl;
		return;
	}
	TFile* outFile = TFile::Open(outFileName.c_str(), "RECREATE");
	if (outFile == NULL) {
		std::cerr << "Output file '" << outFileName << "' could not be opened." << std::endl;
		return;
	}

	TDirectory* outGlobal = outFile->mkdir("global");
	TDirectory* outSum = outFile->mkdir("sum");

	// guess if we are running on data that has been generated from one or
	// two sets of phasespace events
	const bool twoMc = fileContainsTwoMc(inFile);

	// get list of mass bins from the file
	std::vector<MassBin> massBins = getMassBins(inFile);

	// get list of bookies to create
	std::map<std::string, BookyDefinition*> bookies = setupBookies(massBins, outGlobal, twoMc);

	for (std::vector<MassBin>::const_iterator it=massBins.begin(); it!=massBins.end(); it++) {
		MassBin massBin = *it;

		TDirectory* outDir = outFile->GetDirectory(it->GetDir()->GetName());
		if (outDir == NULL) {
			outDir = outFile->mkdir(it->GetDir()->GetName());
		}
		outDir->cd();

		for (std::map<std::string, BookyDefinition*>::const_iterator itBookies=bookies.begin(); itBookies!=bookies.end(); itBookies++) {
			const std::string bookyName = itBookies->first;
			BookyDefinition* bookyPlots = itBookies->second;

			bookyPlots->Process(massBin, twoMc);
			outDir->cd();
			bookyPlots->GetCanvas()->Write(NULL, TObject::kOverwrite);
		}

		addHistograms(outSum, massBin);

		freeMemory(massBin);
	}

	for (std::map<std::string, BookyDefinition*>::const_iterator itBookies=bookies.begin(); itBookies!=bookies.end(); itBookies++) {
		BookyDefinition* bookyPlots = itBookies->second;

		bookyPlots->Process(MassBin(outSum, -1., -1., -1.), twoMc);
		outGlobal->cd();
		bookyPlots->GetCanvas()->Write(NULL, TObject::kOverwrite);

		std::vector<TCanvas*> cs = bookyPlots->Finalize();
		for (std::vector<TCanvas*>::iterator c=cs.begin(); c!=cs.end(); c++) {
			outGlobal->cd();
			(*c)->Write(NULL, TObject::kOverwrite);
			delete *c;
		}
		delete itBookies->second;
	}

	std::cout<< "saving to disk..." <<std::endl;
	outFile->Write();
	outFile->Close();
	delete outFile;

	inFile->Close();
	delete inFile;
	std::cout<< "done!" <<std::endl;
}
