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

#include "TFile.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"
#include "TString.h"
#include "TObject.h"
#include "TCollection.h"
#include "TCanvas.h"
#include "TExec.h"
#include <iostream>
#include <string>
#include <map>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <algorithm>


int massbinwidth = 0;
const unsigned int numbercolors = 31;

struct booky_page {
    TCanvas *page_canvas; // canvas that will be print on that page
    double mass; // mass value with can be used for sorting purposes
    TString page_title; // title displayed on top of page
};

booky_page createBookyPage(TCanvas *c, double mass = 0.0, TString title = "") {
  booky_page p;
  p.page_canvas = c;
  p.mass = mass;
  p.page_title = title;
  return p;
}

// the make booky method will run through this map and create 1 booky for each map entry
// all canvases in the vector will be appended to this booky ps file
// key is the filename of the booky (without .ps), vector<TCanvas*> is the canvases that will be attached to this booky
// we save pairs in the vector instead, to be able to sort the canvases with ascending mass order (the double is the mass)
std::map<TString, std::vector<booky_page> > booky_map;

std::map<TString, std::vector<TString> > booky_setup_map;
// here just fill a vector with Histogram names and add this to the booky_setup_map with the name of the booky
void setupBookies() {
  // lets create a booky for the the 5 1D projections (3 FS particle case)
  std::vector<TString> histnames_neutral;
  std::vector<TString> histnames_charged;
  // neutral isobar decay
  histnames_neutral.push_back(TString("hMIsobarMC_Neutral")); // Canvas spot 1
  histnames_neutral.push_back(TString("spacer")); // Canvas spot 2
  histnames_neutral.push_back(TString("hGJMC_Neutral")); // Canvas spot 3
  histnames_neutral.push_back(TString("hTYMC_Neutral")); // Canvas spot 4
  histnames_neutral.push_back(TString("hHThetaMC_Neutral")); // Canvas spot 5
  histnames_neutral.push_back(TString("hHPhiMC_Neutral")); // Canvas spot 6
  // charged isobar decay
  histnames_charged.push_back(TString("hMIsobarMC_Charged")); // Canvas spot 1
  histnames_charged.push_back(TString("spacer")); // Canvas spot 2
  histnames_charged.push_back(TString("hGJMC_Charged")); // Canvas spot 3
  histnames_charged.push_back(TString("hTYMC_Charged")); // Canvas spot 4
  histnames_charged.push_back(TString("hHThetaMC_Charged")); // Canvas spot 5
  histnames_charged.push_back(TString("hHPhiMC_Charged")); // Canvas spot 6

  booky_setup_map.insert(std::pair<TString, std::vector<TString> >("Booky_neutral_isobar", histnames_neutral));
  booky_setup_map.insert(std::pair<TString, std::vector<TString> >("Booky_charged_isobar", histnames_charged));
}

bool comparePairs (const booky_page &i, const booky_page &j) {
  return (i.mass < j.mass);
}

void makeBookies() {
  // loop through map
  std::map<TString, std::vector<booky_page> >::iterator it;
  for(it = booky_map.begin(); it != booky_map.end(); it++) {
    TString psFileName(it->first);
    psFileName.Append(".ps");
    TCanvas dummyCanv("dummy", "dummy");
    dummyCanv.Print((psFileName + "["));
    // sort vector
    std::vector<booky_page> vec = it->second;
    std::sort(vec.begin(), vec.end(), comparePairs);
    // loop through vector and append canvases
    for (unsigned int i = 0; i < vec.size(); i++) {
      TString textlabel;
      if (vec[i].page_title.Length() == 0) {
        textlabel = "#scale[0.5]{Mass Bin ";
        textlabel += i + 1;
        textlabel += " [";
        textlabel += (int) (vec[i].mass * 1000 - massbinwidth / 2);
        textlabel += " - ";
        textlabel += ((int) (vec[i].mass * 1000) + massbinwidth / 2);
        textlabel += "] MeV/c^{2}}";
      }
      else {
        textlabel = vec[i].page_title;
      }
      TLatex label(0.37, 1.0, textlabel);
      label.SetNDC(true);
      label.SetTextAlign(23);
      vec[i].page_canvas->cd();
      label.Draw();

      vec[i].page_canvas->Print(psFileName);
    }
    dummyCanv.Print((psFileName + "]"));
  }
}

void fillDiffvsMassPlot(const TH1D* hist, std::string dirname, int massbins, double mass_start, double mass_end,
    std::map<std::string, std::pair<double, std::pair<double, double> > > diffbounds, TFile* ofile) {
  // first change into global directory of outfile
  ofile->cd("global");
  // now parse mass from dirname
  int pos = dirname.find(".");
  std::string masslow = dirname.substr(0, pos+1);
  std::string masshigh = dirname.substr(pos+1);
  double mass = atof(masslow.c_str()) + atof(masshigh.c_str());
  mass = mass/2/1000;
  // then get or create 2d hist with number of bins in x equal to number of massbins
  std::string s(hist->GetName());
  std::map<std::string, std::pair<double, std::pair<double, double> > >::iterator iter = diffbounds.find(s);
  s.insert(s.size(), "vsMass");
  TH2D* hdiffvsmass = (TH2D*)ofile->Get((std::string("global/")+s).c_str());
  if(!hdiffvsmass) {
    TCanvas *c = new TCanvas(s.c_str(), s.c_str());
    hdiffvsmass = new TH2D(s.c_str(), hist->GetTitle(), massbins, mass_start, mass_end,
        hist->GetNbinsX(), iter->second.second.first, iter->second.second.second);
    hdiffvsmass->SetContour(numbercolors);
    hdiffvsmass->SetXTitle("Resonance Mass [GeV/c^{2}]");
    hdiffvsmass->SetYTitle(hist->GetXaxis()->GetTitle());
    hdiffvsmass->SetMaximum(iter->second.first);
    hdiffvsmass->SetMinimum(-iter->second.first);
    hdiffvsmass->SetStats(0);
    hdiffvsmass->Draw("colz");

    // now lets add this canvas to its corresponding booky
    // first check if our booky_map already has an open booky for this type and get the vector
    std::vector<booky_page>& tempvec = booky_map["Booky_mass_overview"];
    // create new entry for this vector
    tempvec.push_back(createBookyPage(c, mass, TString(hist->GetName())));
  }
  // then loop over 1d hist bins -> corresponding value + get content and fill 2d hist
  for(int i = 0; i < hist->GetNbinsX(); i++) {
    double diffpos = hist->GetBinCenter(i);
    double diffval = hist->GetBinContent(i);
    if(diffval != 0.0)
      hdiffvsmass->Fill(mass, diffpos, diffval);
  }
}


void make2DOverviewCanvas(TH2D *mchist, TH2D *datahist, TH2D *diffhist, TH2D *reldiffhist, double mass) {
  if (mchist && datahist && diffhist && reldiffhist) {
    double scale = datahist->Integral();
    scale = scale / (mchist->Integral());
    mchist->Scale(scale);

    double max = mchist->GetMaximum();
    if(datahist->GetMaximum() > max)
      mchist->SetMaximum(datahist->GetMaximum());
    else
      datahist->SetMaximum(max);


    mchist->GetXaxis()->SetRangeUser(0.0, pow(mass,2));
    mchist->GetYaxis()->SetRangeUser(0.0, pow(mass,2));
    datahist->GetXaxis()->SetRangeUser(0.0, pow(mass,2));
    datahist->GetYaxis()->SetRangeUser(0.0, pow(mass,2));
    diffhist->GetXaxis()->SetRangeUser(0.0, pow(mass,2));
    diffhist->GetYaxis()->SetRangeUser(0.0, pow(mass,2));
    reldiffhist->GetXaxis()->SetRangeUser(0.0, pow(mass,2));
    reldiffhist->GetYaxis()->SetRangeUser(0.0, pow(mass,2));
    reldiffhist->SetMaximum(1.0);
    reldiffhist->SetMinimum(-1.0);


    TString s("setDiffColorStyle(");
    s += numbercolors;
    s += ");";

    char name[200];
    sprintf(name, "%s_%f", mchist->GetName(), mass);
    TCanvas *c = new TCanvas(name, mchist->GetTitle());
    c->Divide(2,2);
    c->cd(1);
    mchist->Draw("colz");
    mchist->SetContour(numbercolors);
    //TExec *ex1 = new TExec("ex1", "setNormalColorStyle();");
    TExec *ex1 = new TExec("ex1", "gStyle->SetPalette(1);");
    ex1->Draw();
    mchist->Draw("colz same");
    c->cd(2);
    datahist->Draw("colz");
    datahist->SetContour(numbercolors);
    //TExec *ex2 = new TExec("ex2", "setNormalColorStyle();");
    TExec *ex2 = new TExec("ex2", "gStyle->SetPalette(1);");
    ex2->Draw();
    datahist->Draw("colz same");
    c->cd(3);
    diffhist->Draw("colz");
    diffhist->SetContour(numbercolors);
    TExec *ex3 = new TExec("ex3", s);
    ex3->Draw();
    diffhist->Draw("colz same");
    c->cd(4);
    reldiffhist->Draw("colz");
    reldiffhist->SetContour(numbercolors);
    TExec *ex4 = new TExec("ex4", s);
    ex4->Draw();
    reldiffhist->Draw("colz same");
    c->Update();
    c->Write();

    // now lets add this canvas to its corresponding booky
    // first check if our booky_map already has an open booky for this type and get the vector
    std::vector<booky_page>& tempvec = booky_map[mchist->GetName()];
    // create new entry for this vector
    tempvec.push_back(createBookyPage(c, mass));
  }
}


void make1DOverviewCanvas(TFile *infile, TFile *outfile, TList *mclist, std::string dirname) {
  double mass = 0.0;
  int massstart = 0;
  int massend = 0;
  // check if directory is mass bin dir
  unsigned int pointpos = dirname.find(".");
  if (!(pointpos == 0 || pointpos == dirname.size())) {
    std::string masslow = dirname.substr(0, pointpos + 1);
    std::string masshigh = dirname.substr(pointpos + 1);
    massstart = atoi(masslow.c_str());
    massend = atoi(masshigh.c_str());
    mass = 1.0 * (massstart + massend) / 2 / 1000;
  }

  std::map<TString, std::vector<TString> >::iterator it;
  for (it = booky_setup_map.begin(); it != booky_setup_map.end(); it++) {
    TString name(it->first.Data());
    name += dirname.c_str();

    TCanvas *c = new TCanvas(name, "");
    c->Divide(2,3);
    std::vector<TString> histlist = it->second;

    for (unsigned int i = 0; i < histlist.size(); i++) {
      // CompareTo returns 0 if its a match....
      if (histlist[i].CompareTo("spacer")) {
        TIter histiter = TIter(mclist);
        TH1D *reldiffhist, *diffhist, *mchist, *datahist;
        // generate difference histograms
        std::string hnamemc(histlist[i].Data());
        // create new string with MC exchanged for Diff
        std::string hnamediff(hnamemc);
        int pos = hnamemc.find("MC");
        hnamediff.erase(pos, 2);
        hnamediff.insert(pos, "Diff");
        // create new string with MC exchanged for RelDiff
        std::string hnamereldiff(hnamemc);
        hnamereldiff.erase(pos, 2);
        hnamereldiff.insert(pos, "RelDiff");
        // create new string with MC exchanged for Data
        std::string hnamedata(hnamemc);
        hnamedata.erase(pos, 2);
        hnamedata.insert(pos, "Data");

        infile->GetObject((dirname + "/" + hnamereldiff).c_str(), reldiffhist);
        infile->GetObject((dirname + "/" + hnamediff).c_str(), diffhist);
        infile->GetObject((dirname + "/" + hnamedata).c_str(), datahist);
        infile->GetObject((dirname + "/" + hnamemc).c_str(), mchist);

        outfile->cd(dirname.c_str());
        if (mchist && datahist) {
          c->cd(i + 1);

          // std::cout<<i<<std::endl;
          double scale = datahist->Integral();
          scale = scale / (mchist->Integral());
          mchist->Scale(scale);

          mchist->SetLineColor(kRed);
          mchist->SetFillColor(kRed);
          mchist->Draw("E4");
          datahist->Draw("same");
          if (diffhist) {
            TLine* line = new TLine(mchist->GetXaxis()->GetXmin(), 0, mchist->GetXaxis()->GetXmax(), 0);
            line->SetLineStyle(3);
            line->Draw();
            diffhist->SetLineColor(kOrange - 3);
            diffhist->Draw("same");
          }
          double max = mchist->GetMaximum();
          double min = mchist->GetMinimum();
          if (max < datahist->GetMaximum())
            max = datahist->GetMaximum();
          if (diffhist)
            min = diffhist->GetMinimum();
          mchist->GetYaxis()->SetRangeUser(diffhist->GetMinimum() * 1.5, max * 1.2);
          c->Update();
        }
      }
    }
    c->Write();
    // now lets add this canvas to its corresponding booky
    // first check if our booky_map already has an open booky for this type and get the vector
    std::vector<booky_page>& tempvec = booky_map[it->first];
    // create new entry for this vector
    tempvec.push_back(createBookyPage(c, mass));
  }
}

/*
 * this script takes 2 TStrings as root filenames as a parameters
 * basic functionality:
 * loop through all directories (the mass bins) in the root file
 * -> create difference plots
 * -> create global plots
 * -> create 2D diff vs mass plots
 * -> etc...
 */
void plotGlobalWeightedEvts_3pin(TString input_filename, TString output_filename) {
  setupBookies();

  int massbins =0;
  double mass= 0.0, massstart =1000.0, massend=0.0;
  std::map<std::string, std::pair<double, std::pair<double, double> > > diffbounds;

  TFile* infile = TFile::Open(input_filename, "READ");
  TFile* outfile = new TFile(output_filename, "RECREATE");
  outfile->mkdir("global");

  TList *dirlist = infile->GetListOfKeys();
  massbins = dirlist->GetSize();
  infile->cd();
  TIter diriter(dirlist);
  TDirectory *dir;

  std::cout<< "scanning directories and creating overview canvases..." <<std::endl;
  while ((dir = (TDirectory *)diriter())) {
    std::string dirname = dir->GetName();
    // check if directory is mass bin dir
    unsigned int pointpos = dirname.find(".");
    if(pointpos == 0 || pointpos == dirname.size()) continue;
    std::string masslow = dirname.substr(0, pointpos+1);
    std::string masshigh = dirname.substr(pointpos+1);
    double massstarttemp = atof(masslow.c_str())/1000;
    double massendtemp = atof(masshigh.c_str())/1000;
    if((int)(massendtemp - massstarttemp) != massbinwidth)
      massbinwidth = (int)(massendtemp - massstarttemp);
    mass = (massstarttemp + massendtemp)/2;
    if(massstart > massstarttemp) massstart = massstarttemp;
    if(massend < massendtemp) massend = massendtemp;

    outfile->cd();
    outfile->mkdir(dir->GetName());
    infile->cd(dir->GetName());

    // make list of MC Histograms
    TList mclist;
    TList *histlist = gDirectory->GetListOfKeys();
    TIter histiter(histlist);
    TObject *obj;
    while ((obj = histiter())) {
      TString s(obj->GetName());
      if(s.EndsWith("MC"))
        mclist.Add(obj);
      else if(s.Contains("MC_"))
        mclist.Add(obj);
    }
    make1DOverviewCanvas(infile, outfile, &mclist, dirname);
    histiter = TIter(&mclist);
    TH1D *diffhist, *mchist;
    while ((mchist = (TH1D*)histiter())) {
      // generate difference histograms
      std::string hnamemc(mchist->GetName());
      // create new string with MC exchanged for Diff
      std::string hnamediff(hnamemc);
      int pos = hnamemc.find("MC");
      hnamediff.erase(pos, 2);
      hnamediff.insert(pos, "Diff");

      infile->GetObject((std::string(dir->GetName())+"/"+hnamediff).c_str(), diffhist);

      if (diffhist) {
        // get diff min max values
        std::pair<double, std::pair<double, double> > p;

        bool change =false;
        double maxdiff = diffhist->GetMaximum();
        double maxdifftemp = diffhist->GetMinimum();
        if(abs(maxdifftemp) > maxdiff) maxdiff = maxdifftemp;

        double diffmintemp = diffhist->GetXaxis()->GetXmin();
        double diffmaxtemp = diffhist->GetXaxis()->GetXmax();
        std::map<std::string, std::pair<double, std::pair<double, double> > >::iterator iter = diffbounds.find(
            diffhist->GetName());
        if (iter != diffbounds.end()) {
          p.first = iter->second.first;
          p.second.first = iter->second.second.first;
          p.second.second = iter->second.second.second;

          if (iter->second.first < maxdiff) {
            change = true;
            p.first = maxdiff;
          }
          if (iter->second.second.first > diffmintemp) {
            change = true;
            p.second.first = diffmintemp;
          }
          if (iter->second.second.second < diffmaxtemp) {
            change = true;
            p.second.second = diffmaxtemp;
          }

          if (change) {
            diffbounds[diffhist->GetName()] = p;
          }
        }
        else {
          p.first = maxdiff;
          p.second.first = diffmintemp;
          p.second.second = diffmaxtemp;
          diffbounds.insert(std::pair<std::string, std::pair<double, std::pair<double, double> > >(diffhist->GetName(),
              p));
        }
      }
    }
    histiter = TIter(&mclist);
    TH2D *reldiffhist2d, *diffhist2d, *mchist2d, *datahist2d;
    while ((mchist2d = (TH2D*) histiter())) {
      // generate difference histograms
      std::string hnamemc(mchist2d->GetName());
      // create new string with MC exchanged for Diff
      std::string hnamediff(hnamemc);
      int pos = hnamemc.find("MC");
      hnamediff.erase(pos, 2);
      hnamediff.insert(pos, "Diff");
      // create new string with MC exchanged for RelDiff
      std::string hnamereldiff(hnamemc);
      hnamereldiff.erase(pos, 2);
      hnamereldiff.insert(pos, "RelDiff");
      // create new string with MC exchanged for Data
      std::string hnamedata(hnamemc);
      hnamedata.erase(pos, 2);
      hnamedata.insert(pos, "Data");

      infile->GetObject((std::string(dir->GetName()) + "/" + hnamereldiff).c_str(), reldiffhist2d);
      infile->GetObject((std::string(dir->GetName()) + "/" + hnamediff).c_str(), diffhist2d);
      infile->GetObject((std::string(dir->GetName()) + "/" + hnamedata).c_str(), datahist2d);
      infile->GetObject((std::string(dir->GetName()) + "/" + hnamemc).c_str(), mchist2d);

      outfile->cd(dir->GetName());
      make2DOverviewCanvas(mchist2d, datahist2d, diffhist2d, reldiffhist2d, mass);
    }
  }

  dirlist = infile->GetListOfKeys();
  infile->cd();
  diriter = TIter(dirlist);

  std::cout << "creating global histograms and 2D diff vs mass plots..." << std::endl;
  while ((dir = (TDirectory *) diriter())) {
    std::string dirname = dir->GetName();
    // check if directory is mass bin dir
    unsigned int pointpos = dirname.find(".");
    if (pointpos == 0 || pointpos == dirname.size())
      continue;

    infile->cd(dir->GetName());

    // make list of MC Histograms
    TList mclist;
    TList *histlist = gDirectory->GetListOfKeys();
    TIter histiter(histlist);
    TObject *obj;
    while ((obj = histiter())) {
      TString s(obj->GetName());
      if (s.EndsWith("MC"))
        mclist.Add(obj);
      else if (s.Contains("MC_"))
        mclist.Add(obj);
    }
    histiter = TIter(&mclist);
    TH1D *hist;
    TH1D *diffhist, *mchist, *datahist;
    while ((hist = (TH1D*) histiter())) {
      // generate difference histograms
      std::string hnamemc(hist->GetName());
      // create new string with MC exchanged for Diff
      std::string hname(hnamemc);
      int pos = hnamemc.find("MC");
      hname.erase(pos, 2);
      hname.insert(pos, "Diff");
      // create new string with MC exchanged for Data
      std::string hnamedata(hnamemc);
      hnamedata.erase(pos, 2);
      hnamedata.insert(pos, "Data");
      // create new string for MC Global Histogram
      std::string hnamemcglob(hnamemc);
      hnamemcglob.insert(pos + 2, "Global");
      // create new string for MC Global Histogram
      std::string hnamedataglob(hnamemc);
      hnamedataglob.erase(pos, 2);
      hnamedataglob.insert(pos, "DataGlobal");

      infile->GetObject((std::string(dir->GetName()) + "/" + hname + ";1").c_str(), diffhist);
      infile->GetObject((std::string(dir->GetName()) + "/" + hnamedata + ";1").c_str(), datahist);
      infile->GetObject((std::string(dir->GetName()) + "/" + hnamemc + ";1").c_str(), mchist);
      if (datahist) {
        // make global histograms in global folder
        outfile->cd("global");
        TH1D* hmcglob = (TH1D*) outfile->Get(std::string("global/"+hnamemcglob).c_str());
        if (hmcglob == NULL)
          hmcglob = new TH1D(hnamemcglob.c_str(), mchist->GetTitle(), mchist->GetNbinsX(),
              mchist->GetXaxis()->GetXmin(), mchist->GetXaxis()->GetXmax());
        hmcglob->Add(mchist);
        TH1D* hdataglob = (TH1D*) outfile->Get(std::string("global/"+hnamedataglob).c_str());
        if (hdataglob == NULL)
          hdataglob = new TH1D(hnamedataglob.c_str(), datahist->GetTitle(), datahist->GetNbinsX(),
              datahist->GetXaxis()->GetXmin(), datahist->GetXaxis()->GetXmax());
        hdataglob->Add(datahist);

        // make diff vs. mass plots
        fillDiffvsMassPlot(diffhist, dir->GetName(), massbins, massstart, massend, diffbounds, outfile);
      }
    }
    histiter = TIter(&mclist);
    TH2D *mchist2d, *datahist2d;
    while ((mchist2d = (TH2D*) histiter())) {
      // generate difference histograms
      std::string hnamemc(mchist2d->GetName());
      // create new string with MC exchanged for Diff
      std::string hnamediff(hnamemc);
      int pos = hnamemc.find("MC");
      hnamediff.erase(pos, 2);
      hnamediff.insert(pos, "Diff");
      // create new string with MC exchanged for Data
      std::string hnamedata(hnamemc);
      hnamedata.erase(pos, 2);
      hnamedata.insert(pos, "Data");
      // create new string for MC Global Histogram
      std::string hnamemcglob(hnamemc);
      hnamemcglob.insert(pos + 2, "Global");
      // create new string for MC Global Histogram
      std::string hnamedataglob(hnamemc);
      hnamedataglob.erase(pos, 2);
      hnamedataglob.insert(pos, "DataGlobal");

      infile->GetObject((std::string(dir->GetName()) + "/" + hnamedata + ";1").c_str(), datahist2d);
      infile->GetObject((std::string(dir->GetName()) + "/" + hnamemc + ";1").c_str(), mchist2d);
      if (datahist2d) {
        // make global histograms in global folder
        outfile->cd("global");
        TH2D* hmcglob = (TH2D*) outfile->Get(std::string("global/" + hnamemcglob).c_str());
        if (hmcglob == NULL) {
          hmcglob = new TH2D(hnamemcglob.c_str(), mchist2d->GetTitle(), mchist->GetNbinsX(),
              mchist2d->GetXaxis()->GetXmin(), mchist2d->GetXaxis()->GetXmax(), mchist2d->GetNbinsY(),
              mchist2d->GetYaxis()->GetXmin(), mchist2d->GetYaxis()->GetXmax());
          hmcglob->SetXTitle(mchist2d->GetXaxis()->GetTitle());
          hmcglob->SetYTitle(mchist2d->GetYaxis()->GetTitle());
        }
        hmcglob->Add(mchist2d);
        TH2D* hdataglob = (TH2D*) outfile->Get(std::string("global/" + hnamedataglob).c_str());
        if (hdataglob == NULL) {
          hdataglob = new TH2D(hnamedataglob.c_str(), datahist2d->GetTitle(), datahist2d->GetNbinsX(),
              datahist2d->GetXaxis()->GetXmin(), datahist2d->GetXaxis()->GetXmax(), datahist2d->GetNbinsY(),
              datahist2d->GetYaxis()->GetXmin(), datahist2d->GetYaxis()->GetXmax());
          hdataglob->SetXTitle(datahist2d->GetXaxis()->GetTitle());
          hdataglob->SetYTitle(datahist2d->GetYaxis()->GetTitle());
        }
        hdataglob->Add(datahist2d);
      }
    }
  }

  makeBookies();

  std::cout<< "saving to disk..." <<std::endl;
  outfile->Write();
  std::cout<< "done!" <<std::endl;

  /*// ok lets make a canvas and plug some plots into it -> add to booky
  TCanvas* c = new TCanvas("KineValidate" + massbin, "Weighted Events", 10, 10, 600, 800);
  c->Divide(4, 4);

  // first column contains neutral isobar histograms
  c->cd(1);

  double totMC = GJHB_neutral_isobar.isobar_mass[0]->Integral();
  double totDATA = GJHB_neutral_isobar.isobar_mass[1]->Integral();

  if (totMC != 0)
    GJHB_neutral_isobar.isobar_mass[0]->Scale(totDATA / totMC);
  GJHB_neutral_isobar.isobar_mass[0]->SetLineColor(kRed);
  GJHB_neutral_isobar.isobar_mass[0]->SetFillColor(kRed);
  GJHB_neutral_isobar.isobar_mass[0]->Draw("E4");

  TH1D* hDiffMIsobar = new TH1D(*GJHB_neutral_isobar.isobar_mass[0]);
  hDiffMIsobar->Add(GJHB_neutral_isobar.isobar_mass[1], -1.);
  GJHB_neutral_isobar.isobar_mass[1]->Draw("same");
  hDiffMIsobar->SetLineColor(kOrange - 3);
  hDiffMIsobar->Draw("same");

  GJHB_neutral_isobar.isobar_mass[0]->GetYaxis()->SetRangeUser(hDiffMIsobar->GetMinimum() * 1.5,
      GJHB_neutral_isobar.isobar_mass[0]->GetMaximum() * 1.2);
  gPad->Update();

  c->cd(5);
  totMC = GJHB_neutral_isobar.costheta_GJF[0]->Integral();
  totDATA = GJHB_neutral_isobar.costheta_GJF[1]->Integral();
  if (totMC != 0)
    GJHB_neutral_isobar.costheta_GJF[0]->Scale(totDATA / totMC);
  GJHB_neutral_isobar.costheta_GJF[0]->SetLineColor(kRed);
  GJHB_neutral_isobar.costheta_GJF[0]->SetFillColor(kRed);
  GJHB_neutral_isobar.costheta_GJF[0]->Draw("E4");

  GJHB_neutral_isobar.costheta_GJF[1]->Draw("same E");

  totMC = GJHB_neutral_isobar.costheta_GJF_MC_raw->Integral();
  GJHB_neutral_isobar.costheta_GJF_MC_raw->Scale(totDATA / totMC);
  GJHB_neutral_isobar.costheta_GJF_MC_raw->Draw("same");

  GJHB_neutral_isobar.costheta_GJF[0]->GetYaxis()->SetRangeUser(0, GJHB_neutral_isobar.costheta_GJF[0]->GetMaximum() * 1.5);
  gPad->Update();

  c->cd(9);
  totMC = GJHB_neutral_isobar.phi_GJF[0]->Integral();
  totDATA = GJHB_neutral_isobar.phi_GJF[1]->Integral();
  GJHB_neutral_isobar.phi_GJF[0]->Sumw2();
  if (totMC != 0)
    GJHB_neutral_isobar.phi_GJF[0]->Scale(totDATA / totMC);
  GJHB_neutral_isobar.phi_GJF[0]->SetLineColor(kRed);
  GJHB_neutral_isobar.phi_GJF[0]->SetFillColor(kRed);
  GJHB_neutral_isobar.phi_GJF[0]->Draw("E4");

  GJHB_neutral_isobar.phi_GJF[1]->Draw("same E");
  GJHB_neutral_isobar.phi_GJF[0]->GetYaxis()->SetRangeUser(0, GJHB_neutral_isobar.phi_GJF[0]->GetMaximum() * 1.5);
  gPad->Update();

  c->cd(13);
  totMC = HHB_neutral_isobar.costheta_HF[0]->Integral();
  totDATA = HHB_neutral_isobar.costheta_HF[1]->Integral();
  HHB_neutral_isobar.costheta_HF[0]->Sumw2();
  if (totMC != 0)
    HHB_neutral_isobar.costheta_HF[0]->Scale(totDATA / totMC);
  HHB_neutral_isobar.costheta_HF[0]->SetLineColor(kRed);
  HHB_neutral_isobar.costheta_HF[0]->SetFillColor(kRed);
  HHB_neutral_isobar.costheta_HF[0]->Draw("E4");

  HHB_neutral_isobar.costheta_HF[1]->Draw("same E");
  HHB_neutral_isobar.costheta_HF[0]->GetYaxis()->SetRangeUser(0,
      HHB_neutral_isobar.costheta_HF[0]->GetMaximum() * 1.5);
  gPad->Update();

  c->cd(15);
  totMC = HHB_neutral_isobar.phi_HF[0]->Integral();
  totDATA = HHB_neutral_isobar.phi_HF[1]->Integral();
  HHB_neutral_isobar.phi_HF[0]->Sumw2();
  if (totMC != 0)
    HHB_neutral_isobar.phi_HF[0]->Scale(totDATA / totMC);
  HHB_neutral_isobar.phi_HF[0]->SetLineColor(kRed);
  HHB_neutral_isobar.phi_HF[0]->SetFillColor(kRed);
  HHB_neutral_isobar.phi_HF[0]->Draw("E4");

  HHB_neutral_isobar.phi_HF[1]->Draw("same E");
  HHB_neutral_isobar.phi_HF[0]->GetYaxis()->SetRangeUser(0,
      HHB_neutral_isobar.phi_HF[0]->GetMaximum() * 1.5);
  gPad->Update();

  c->cd(2);

  totMC = GJHB_charged_isobar.isobar_mass[0]->Integral();
  totDATA = GJHB_charged_isobar.isobar_mass[1]->Integral();

  if (totMC != 0)
    GJHB_charged_isobar.isobar_mass[0]->Scale(totDATA / totMC);
  GJHB_charged_isobar.isobar_mass[0]->SetLineColor(kRed);
  GJHB_charged_isobar.isobar_mass[0]->SetFillColor(kRed);
  GJHB_charged_isobar.isobar_mass[0]->Draw("E4");

  TH1D* hDiffMIsobar2 = new TH1D(*GJHB_charged_isobar.isobar_mass[0]);
  hDiffMIsobar2->Add(GJHB_charged_isobar.isobar_mass[1], -1.);
  GJHB_charged_isobar.isobar_mass[1]->Draw("same");
  hDiffMIsobar2->SetLineColor(kOrange - 3);
  hDiffMIsobar2->Draw("same");

  GJHB_charged_isobar.isobar_mass[0]->GetYaxis()->SetRangeUser(hDiffMIsobar->GetMinimum() * 1.5,
      GJHB_charged_isobar.isobar_mass[0]->GetMaximum() * 1.2);
  gPad->Update();

  c->cd(6);
  totMC = GJHB_charged_isobar.costheta_GJF[0]->Integral();
  totDATA = GJHB_charged_isobar.costheta_GJF[1]->Integral();
  //hGJ[0]->Sumw2();
  if (totMC != 0)
    GJHB_charged_isobar.costheta_GJF[0]->Scale(totDATA / totMC);
  GJHB_charged_isobar.costheta_GJF[0]->SetLineColor(kRed);
  GJHB_charged_isobar.costheta_GJF[0]->SetFillColor(kRed);
  GJHB_charged_isobar.costheta_GJF[0]->Draw("E4");

  GJHB_charged_isobar.costheta_GJF[1]->Draw("same E");

  totMC = GJHB_charged_isobar.costheta_GJF_MC_raw->Integral();
  GJHB_charged_isobar.costheta_GJF_MC_raw->Scale(totDATA / totMC);
  GJHB_charged_isobar.costheta_GJF_MC_raw->Draw("same");

  GJHB_charged_isobar.costheta_GJF[0]->GetYaxis()->SetRangeUser(0, GJHB_charged_isobar.costheta_GJF[0]->GetMaximum() * 1.5);
  gPad->Update();

  c->cd(10);
  totMC = GJHB_charged_isobar.phi_GJF[0]->Integral();
  totDATA = GJHB_charged_isobar.phi_GJF[1]->Integral();
  GJHB_charged_isobar.phi_GJF[0]->Sumw2();
  if (totMC != 0)
    GJHB_charged_isobar.phi_GJF[0]->Scale(totDATA / totMC);
  GJHB_charged_isobar.phi_GJF[0]->SetLineColor(kRed);
  GJHB_charged_isobar.phi_GJF[0]->SetFillColor(kRed);
  GJHB_charged_isobar.phi_GJF[0]->Draw("E4");

  GJHB_charged_isobar.phi_GJF[1]->Draw("same E");
  GJHB_charged_isobar.phi_GJF[0]->GetYaxis()->SetRangeUser(0, GJHB_charged_isobar.phi_GJF[0]->GetMaximum() * 1.5);
  gPad->Update();

  c->cd(14);
  totMC = HHB_charged_isobar.costheta_HF[0]->Integral();
  totDATA = HHB_charged_isobar.costheta_HF[1]->Integral();
  HHB_charged_isobar.costheta_HF[0]->Sumw2();
  if (totMC != 0)
    HHB_charged_isobar.costheta_HF[0]->Scale(totDATA / totMC);
  HHB_charged_isobar.costheta_HF[0]->SetLineColor(kRed);
  HHB_charged_isobar.costheta_HF[0]->SetFillColor(kRed);
  HHB_charged_isobar.costheta_HF[0]->Draw("E4");

  HHB_charged_isobar.costheta_HF[1]->Draw("same E");
  HHB_charged_isobar.costheta_HF[0]->GetYaxis()->SetRangeUser(0,
      HHB_charged_isobar.costheta_HF[0]->GetMaximum() * 1.5);
  gPad->Update();

  c->cd(16);
  totMC = HHB_charged_isobar.phi_HF[0]->Integral();
  totDATA = HHB_charged_isobar.phi_HF[1]->Integral();
  HHB_charged_isobar.phi_HF[0]->Sumw2();
  if (totMC != 0)
    HHB_charged_isobar.phi_HF[0]->Scale(totDATA / totMC);
  HHB_charged_isobar.phi_HF[0]->SetLineColor(kRed);
  HHB_charged_isobar.phi_HF[0]->SetFillColor(kRed);
  HHB_charged_isobar.phi_HF[0]->Draw("E4");

  HHB_charged_isobar.phi_HF[1]->Draw("same E");
  HHB_charged_isobar.phi_HF[0]->GetYaxis()->SetRangeUser(0,
      HHB_charged_isobar.phi_HF[0]->GetMaximum() * 1.5);
  gPad->Update();

  // add global diagrams

  c->cd(1);
  TLatex* Label = new TLatex();
  Label->PaintLatex(2, GJHB_neutral_isobar.isobar_mass[0]->GetMaximum() * 0.8, 0, 0.1, massbin.Data());

  c->Update();

  TList* Hlist = gDirectory->GetList();
  Hlist->Remove(mctr);
  Hlist->Remove(datatr);
  //Hlist->Remove("hWeights");

  TFile* outfile = TFile::Open(outfilename, "UPDATE");
  TString psFileName = outfilename;
  psFileName.ReplaceAll(".root", massbin);
  psFileName.Append(".ps");

  if (totMC != 0) {
    // get mass-integrated plots (or create them if not there)

    // neutral isobar
    TH1D* hMIsobarMCGlob = (TH1D*) outfile->Get("hMIsobarMCGlob");
    if (hMIsobarMCGlob == NULL)
      hMIsobarMCGlob = new TH1D("hMIsobarMCGlob", "2#pi neutral Isobar Mass (MC)", nbninsm, 0.2, 2.5);
    hMIsobarMCGlob->Add(GJHB_neutral_isobar.isobar_mass[0]);
    hMIsobarMCGlob->Write();
    TH1D* hMIsobarDataGlob = (TH1D*) outfile->Get("hMIsobarDataGlob");
    if (hMIsobarDataGlob == NULL)
      hMIsobarDataGlob = new TH1D("hMIsobarDataGlob", "2#pi neutral Isobar Mass (DATA)", nbninsm, 0.2, 2.5);
    hMIsobarDataGlob->Add(GJHB_neutral_isobar.isobar_mass[1]);
    hMIsobarDataGlob->Write();

    TH1D* hGJMCGlob = (TH1D*) outfile->Get("hGJMCGlob");
    if (hGJMCGlob == NULL)
      hGJMCGlob = new TH1D("hGJMCGlob", "Gottfried-Jackson Theta (MC)", nbinsang, -1, 1);
    hGJMCGlob->Add(GJHB_neutral_isobar.costheta_GJF[0]);
    hGJMCGlob->Write();
    TH1D* hGJDataGlob = (TH1D*) outfile->Get("hGJDataGlob");
    if (hGJDataGlob == NULL)
      hGJDataGlob = new TH1D("hGJDataGlob", "Gottfried-Jackson Theta (Data)", nbinsang, -1, 1);
    hGJDataGlob->Add(GJHB_neutral_isobar.costheta_GJF[1]);
    hGJDataGlob->Write();

    TH1D* hTYMCGlob = (TH1D*) outfile->Get("hTYMCGlob");
    if (hTYMCGlob == NULL)
      hTYMCGlob = new TH1D("hTYMCGlob", "Treiman-Yang Phi (MC)", nbinsang, -TMath::Pi(),TMath::Pi());
    hTYMCGlob->Add(GJHB_neutral_isobar.phi_GJF[0]);
    hTYMCGlob->Write();
    TH1D* hTYDataGlob = (TH1D*) outfile->Get("hTYDataGlob");
    if (hTYDataGlob == NULL)
      hTYDataGlob = new TH1D("hTYDataGlob", "Treiman-Yang Phi (Data)", nbinsang, -TMath::Pi(),TMath::Pi());
    hTYDataGlob->Add(GJHB_neutral_isobar.phi_GJF[1]);
    hTYDataGlob->Write();

    c->cd(3);
    hMIsobarMCGlob->SetLineColor(kRed);
    hMIsobarMCGlob->SetFillColor(kRed);
    hMIsobarMCGlob->Draw("E4");
    hMIsobarDataGlob->Draw("E SAME");

    c->cd(7);
    hGJMCGlob->SetLineColor(kRed);
    hGJMCGlob->SetFillColor(kRed);
    hGJMCGlob->Draw("E4");
    hGJDataGlob->Draw("E SAME");

    c->cd(11);
    hTYMCGlob->SetLineColor(kRed);
    hTYMCGlob->SetFillColor(kRed);
    hTYMCGlob->Draw("E4");
    hTYDataGlob->Draw("E SAME");


    // charged isobar
    TH1D* hMIsobarMCGlob2 = (TH1D*) outfile->Get("hMIsobarMCGlob2");
    if (hMIsobarMCGlob2 == NULL)
      hMIsobarMCGlob2 = new TH1D("hMIsobarMCGlob2", "2#pi charged Isobar Mass (MC)", nbninsm, 0.2,
          2.5);
    hMIsobarMCGlob2->Add(GJHB_charged_isobar.isobar_mass[0]);
    hMIsobarMCGlob2->Write();
    TH1D* hMIsobarDataGlob2 = (TH1D*) outfile->Get("hMIsobarDataGlob2");
    if (hMIsobarDataGlob2 == NULL)
      hMIsobarDataGlob2 = new TH1D("hMIsobarDataGlob2", "2#pi charged Isobar mass (DATA)", nbninsm,
          0.2, 2.5);
    hMIsobarDataGlob2->Add(GJHB_charged_isobar.isobar_mass[1]);
    hMIsobarDataGlob2->Write();

    TH1D* hGJMCGlob2 = (TH1D*) outfile->Get("hGJMCGlob2");
    if (hGJMCGlob2 == NULL)
      hGJMCGlob2 = new TH1D("hGJMCGlob2", "Gottfried-Jackson Theta (MC)", nbinsang, -1, 1);
    hGJMCGlob2->Add(GJHB_charged_isobar.costheta_GJF[0]);
    hGJMCGlob2->Write();
    TH1D* hGJDataGlob2 = (TH1D*) outfile->Get("hGJDataGlob2");
    if (hGJDataGlob2 == NULL)
      hGJDataGlob2 = new TH1D("hGJDataGlob2", "Gottfried-Jackson Theta (Data)", nbinsang, -1, 1);
    hGJDataGlob2->Add(GJHB_charged_isobar.costheta_GJF[1]);
    hGJDataGlob2->Write();

    TH1D* hTYMCGlob2 = (TH1D*) outfile->Get("hTYMCGlob2");
    if (hTYMCGlob2 == NULL)
      hTYMCGlob2 = new TH1D("hTYMCGlob2", "Treiman-Yang Phi (MC)", nbinsang, -TMath::Pi(),TMath::Pi());
    hTYMCGlob2->Add(GJHB_charged_isobar.phi_GJF[0]);
    hTYMCGlob2->Write();
    TH1D* hTYDataGlob2 = (TH1D*) outfile->Get("hTYDataGlob2");
    if (hTYDataGlob2 == NULL)
      hTYDataGlob2 = new TH1D("hTYDataGlob2", "Treiman-Yang Phi (Data)", nbinsang, -TMath::Pi(),TMath::Pi());
    hTYDataGlob2->Add(GJHB_charged_isobar.phi_GJF[1]);
    hTYDataGlob2->Write();

    c->cd(4);
    hMIsobarMCGlob2->SetLineColor(kRed);
    hMIsobarMCGlob2->SetFillColor(kRed);
    hMIsobarMCGlob2->Draw("E4");
    hMIsobarDataGlob2->Draw("E SAME");

    c->cd(8);
    hGJMCGlob2->SetLineColor(kRed);
    hGJMCGlob2->SetFillColor(kRed);
    hGJMCGlob2->Draw("E4");
    hGJDataGlob2->Draw("E SAME");

    c->cd(12);
    hTYMCGlob2->SetLineColor(kRed);
    hTYMCGlob2->SetFillColor(kRed);
    hTYMCGlob2->Draw("E4");
    hTYDataGlob2->Draw("E SAME");
  }

  c->Write();
  TCanvas dummyCanv("dummy", "dummy");
  c->Print((psFileName));*/
}
