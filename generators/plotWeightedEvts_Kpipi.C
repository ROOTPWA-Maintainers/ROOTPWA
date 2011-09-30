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

#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TLatex.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TLorentzRotation.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TClonesArray.h"
#include "TROOT.h"
#include "TMath.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TString.h"
#include <iostream>
#include <vector>
#include "NParticleEvent.h"
#include "TPRegexp.h"

using namespace std;

TString massbin;
TString mass;

int nbninsm = 150;
int nbinsang = 80;



struct GJHistBunch {
    // base histograms
    std::vector<TH1D*> isobar_mass;
    std::vector<TH1D*> costheta_GJF;
    std::vector<TH1D*> phi_GJF;
    TH1D* costheta_GJF_MC_raw;
    std::vector<TH2D*> costheta_GJF_tprime;

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

TString stripUnits(TString s) {
  //std::cout<<"stripunits-in: "<<s<<std::endl;
  TPRegexp re("^(.*)\\s*(\\[.*])\\s*$");
  re.Substitute(s,"$1");
  //std::cout<<"stripunits-out: "<<s<<std::endl;
  return s;
}

TString getUnits(TString s) {
  //std::cout << "getunits-in: " << s << std::endl;
  TPRegexp re("^.*(\\[.*])\\s*$");
  if(re.MatchB(s) > 0)
    re.Substitute(s,"$1");
  else
    s = "";
  //std::cout << "getunits-out: " << s << std::endl;
  return s;
}

GJHistBunch GJHistBunchFactory(TString name_prefix) {
  GJHistBunch temp;
  TH1D* hMIsobarMC = new TH1D("hMIsobarMC_" + name_prefix, name_prefix + " Isobar Mass (MC)", nbninsm, 0.0, 3.0);
  hMIsobarMC->SetXTitle("isobar mass [GeV]");
  hMIsobarMC->SetYTitle("# of events");
  temp.isobar_mass.push_back(hMIsobarMC);
  TH1D* hMIsobarData = new TH1D("hMIsobarData_" + name_prefix, name_prefix + " Isobar Mass (DATA)", nbninsm,
      0.0, 3.0);
  hMIsobarData->SetXTitle("isobar mass [GeV]");
  hMIsobarData->SetYTitle("# of events");
  temp.isobar_mass.push_back(hMIsobarData);
  temp.isobar_mass[0]->Sumw2();

  TH1D* hGJMC = new TH1D("hGJMC_" + name_prefix, name_prefix + " Isobar Cos Gottfried-Jackson Theta (MC)", nbinsang, -1, 1);
  hGJMC->SetXTitle("isobar cos(#theta_{GJ})");
  hGJMC->SetYTitle("# of events");
  temp.costheta_GJF.push_back(hGJMC);
  TH1D* hGJData = new TH1D("hGJData_" + name_prefix, name_prefix + " Isobar Cos Gottfried-Jackson Theta (DATA)", nbinsang, -1,
      1);
  hGJData->SetXTitle("isobar cos(#theta_{GJ})");
  hGJData->SetYTitle("# of events");
  temp.costheta_GJF.push_back(hGJData);
  temp.costheta_GJF[0]->Sumw2();

  temp.costheta_GJF_MC_raw = new TH1D("hGJMC_raw" + name_prefix,
      "Cos Gottfried-Jackson Theta (unweighted MC)", nbinsang, -1, 1);
  temp.costheta_GJF_MC_raw->SetXTitle("isobar cos(#theta_{GJ})");
  temp.costheta_GJF_MC_raw->SetYTitle("# of events");

  TH2D* hGJtMC = new TH2D("hGJtMC_" + name_prefix, name_prefix + " Isobar Cos GJ Theta vs t' (MC)", nbinsang, -1, 1, 20, 0,
      0.01);
  hGJtMC->SetXTitle("isobar cos(#theta_{GJ})");
  hGJtMC->SetYTitle("t' [GeV]");
  temp.costheta_GJF_tprime.push_back(hGJtMC);
  TH2D* hGJtData = new TH2D("hGJtData_" + name_prefix, name_prefix + " Isobar Cos GJ Theta vs t' (DATA)", nbinsang, -1, 1, 20,
      0, 0.01);
  hGJData->SetXTitle("isobar cos(#theta_{GJ})");
  hGJData->SetYTitle("t' [GeV]");
  temp.costheta_GJF_tprime.push_back(hGJtData);

  TH1D* hTYMC = new TH1D("hTYMC_" + name_prefix, name_prefix + " Isobar Treiman-Yang Phi (MC)", nbinsang, -TMath::Pi(),
      TMath::Pi());
  hTYMC->SetXTitle("isobar #phi_{TY} [rad]");
  hTYMC->SetYTitle("# of events");
  TH1D* hTYData = new TH1D("hTYData_" + name_prefix, name_prefix + " Isobar Treiman-Yang Phi (DATA)", nbinsang, -TMath::Pi(),
      TMath::Pi());
  hTYData->SetXTitle("isobar #phi_{TY} [rad]");
  hTYData->SetYTitle("# of events");
  temp.phi_GJF.push_back(hTYMC);
  temp.phi_GJF.push_back(hTYData);
  temp.phi_GJF[0]->Sumw2();

  return temp;
}

HelicityHistBunch HelicityHistBunchFactory(TString name_prefix) {
  HelicityHistBunch temp;
  TH1D* hHThetaMC = new TH1D("hHThetaMC_" + name_prefix, name_prefix + " Isobar Helicity Cos Theta (MC)", nbinsang, -1, 1);
  hHThetaMC->SetXTitle("cos(#theta_{hel} of negative particle from isobar)");
  hHThetaMC->SetYTitle("# of events");
  temp.costheta_HF.push_back(hHThetaMC);
  TH1D* hHThetaData = new TH1D("hHThetaData_" + name_prefix, name_prefix + " Isobar Helicity Cos Theta (DATA)", nbinsang, -1,
      1);
  hHThetaData->SetXTitle("cos(#theta_{hel}) of negative particle from isobar");
  hHThetaData->SetYTitle("# of events");
  temp.costheta_HF.push_back(hHThetaData);
  temp.costheta_HF[0]->Sumw2();

  TH1D* hHPhiMC = new TH1D("hHPhiMC_" + name_prefix, name_prefix + " Isobar Helicity Phi (MC)", nbinsang, -TMath::Pi(),
      TMath::Pi());
  hHPhiMC->SetXTitle("#phi_{hel} [rad] of negative particle from isobar");
  hHPhiMC->SetYTitle("# of events");
  TH1D* hHPhiData = new TH1D("hHPhiData_" + name_prefix, name_prefix + " Isobar Helicity Phi (DATA)", nbinsang, -TMath::Pi(),
      TMath::Pi());
  hHPhiData->SetXTitle("#phi_{hel} [rad] of negative particle from isobar");
  hHPhiData->SetYTitle("# of events");
  temp.phi_HF.push_back(hHPhiMC);
  temp.phi_HF.push_back(hHPhiData);
  temp.phi_HF[0]->Sumw2();

  return temp;
}

void fillWeightedHelicityAnglePlots(const HelicityAngles &ha, double weight, unsigned int tree_index, HelicityHistBunch &hhb) {
  hhb.costheta_HF[tree_index]->Fill(ha.cosTheta, weight);
  hhb.phi_HF[tree_index]->Fill(ha.phi, weight);
}

void fillWeightedGJAnglePlots(const TLorentzVector &isobar, double weight, double tprime, unsigned int tree_index, GJHistBunch &hBunch) {
  hBunch.costheta_GJF[tree_index]->Fill(isobar.CosTheta(), weight);
  if (tree_index == 0)
    hBunch.costheta_GJF_MC_raw->Fill(isobar.CosTheta());
  hBunch.costheta_GJF_tprime[tree_index]->Fill(isobar.CosTheta(), tprime, weight);
  hBunch.phi_GJF[tree_index]->Fill(isobar.Phi(), weight);
  hBunch.isobar_mass[tree_index]->Fill(isobar.M(), weight);
}

HelicityAngles calculateHelicityAngles(const NParticleState &isobar, TLorentzVector *beam = NULL) {
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
  if(isobar.getParticle(0)->q() == -1)
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


void makeDifferencePlots(TFile *outfile) {
  outfile->cd();
  TList *dirlist = outfile->GetListOfKeys();
  TIter diriter(dirlist);
  TDirectory *dir;

  std::cout << "scanning directories and creating diff histograms..." << std::endl;
  while ((dir = (TDirectory *) diriter())) {
    std::string dirname = dir->GetName();
    // check if directory is mass bin dir
    unsigned int pointpos = dirname.find(".");
    if (pointpos == 0 || pointpos == dirname.size())
      continue;

    outfile->cd(dir->GetName());

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
    TH1D *diffhist, *reldiffhist, *mchist, *datahist;
    double scale;
    while ((mchist = (TH1D*) histiter())) {
      // generate difference histograms
      std::string hnamemc(mchist->GetName());
      int pos = hnamemc.find("MC");
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
      hnamedata.erase(pos, 2);
      hnamedata.insert(pos, "Data");

      outfile->GetObject((std::string(dir->GetName()) + "/" + hnamediff).c_str(), diffhist);
      outfile->GetObject((std::string(dir->GetName()) + "/" + hnamereldiff).c_str(), reldiffhist);
      outfile->GetObject((std::string(dir->GetName()) + "/" + hnamedata).c_str(), datahist);
      outfile->GetObject((std::string(dir->GetName()) + "/" + hnamemc).c_str(), mchist);
      if (!diffhist) {
        if (datahist) {
          outfile->cd(dir->GetName());
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
      int pos = hnamemc.find("MC");
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
      hnamedata.erase(pos, 2);
      hnamedata.insert(pos, "Data");

      outfile->GetObject((std::string(dir->GetName()) + "/" + hnamediff).c_str(), diffhist2d);
      outfile->GetObject((std::string(dir->GetName()) + "/" + hnamereldiff).c_str(), reldiffhist2d);
      outfile->GetObject((std::string(dir->GetName()) + "/" + hnamedata).c_str(), datahist2d);
      outfile->GetObject((std::string(dir->GetName()) + "/" + hnamemc).c_str(), mchist2d);
      if (!diffhist2d) {
        if (datahist2d) {
          outfile->cd(dir->GetName());
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
}

TH2D* createDalitzHistogram(TString name, TString title, double mass, unsigned int treeentries) {
  double rangelow = 0.0;
  double rangehigh = pow(mass-0.1, 2);
  unsigned int nbins = 0;
  nbins = 0.6*TMath::Sqrt(treeentries)*(rangehigh-rangelow);
  std::cout<<"nbins for mass "<<mass<<": "<<nbins<<std::endl;
  // lower bound is 20 bins
  if(nbins < 20) nbins = 20;
  // upper bound is 130 bins
  if(nbins > 130) nbins = 130;
  TH2D* phist = new TH2D(name, title, nbins, rangelow, rangehigh, nbins, rangelow, rangehigh);
  return phist;
}


void plotWeightedEvts_Kpipi(TTree* mctr, TTree* datatr, TString outfilename = "kineplots.root",
    TString mass_ = "000") {
  cout << " running histo creator " << endl;
  double massval = 0.0;
  unsigned int datatreeentries = 0;
  mass = mass_;
  if (mass != "000") {
    massbin = ("_m") + mass;
    massbin.ReplaceAll(" ", "");
  }

  std::string binname(mass.Data());
  unsigned int pointpos = binname.find(".");
  if(pointpos == 0 || pointpos == binname.size())
    std::cout<<"Warning: Bad massbin name!"<<std::endl;
  std::string masshigh = binname.substr(pointpos+1);
  massval = atof(masshigh.c_str());
  massval /=1000;
  datatreeentries = datatr->GetEntries();


  gROOT->SetStyle("Plain");
  TFile* outfile = TFile::Open(outfilename, "UPDATE");
  outfile->cd();
  gDirectory->cd();
  if (!gDirectory->GetDirectory(mass)) {
    gDirectory->mkdir(mass);
  }
  gDirectory->cd(mass);

  // --------------- global diagrams
  std::vector<TH1D*> hM;
  TH1D* hMMC = new TH1D("hResMassMC", "Mass (MC)", 60, 0.5, 2.5);
  hM.push_back(hMMC);
  TH1D* hMData = new TH1D("hResMassData", "Mass (DATA)", 60, 0.5, 2.5);
  hM.push_back(hMData);

  vector<TH1D*> hThetaLab;
  TH1D* hThetaLabMC = new TH1D("hThetaLabMC", "Cos Theta Lab (MC)", nbninsm, 0.997, 1);
  hThetaLab.push_back(hThetaLabMC);
  TH1D* hThetaLabData = new TH1D("hThetaLabData", "Cos Theta Lab (Data)", nbninsm, 0.997, 1);
  hThetaLab.push_back(hThetaLabData);
  hThetaLab[0]->Sumw2();

  // Dalitz plots
  std::vector<TH2D*> dalitz_neutral;
  /*  TH2D* dalitz = createDalitzHistogram("hDalitzMC",
      "Dalitz Plot #pi^{0}#pi^{0} vs. #pi^{-}#pi^{0} (MC)", massval, datatreeentries);
  dalitz->SetXTitle("mass^{2}(#pi^{0}#pi^{0}) [GeV^{2}/c^{4}]");
  dalitz->SetYTitle("mass^{2}(#pi^{-}#pi^{0}) [GeV^{2}/c^{4}]");
  dalitz->SetStats(0);
  dalitz_neutral.push_back(dalitz);
  dalitz = createDalitzHistogram("hDalitzData",
      "Dalitz Plot #pi^{0}#pi^{0} vs. #pi^{-}#pi^{0} (Data)", massval, datatreeentries);
  dalitz->SetXTitle("mass^{2}(#pi^{0}#pi^{0}) [GeV^{2}/c^{4}]");
  dalitz->SetYTitle("mass^{2}(#pi^{-}#pi^{0}) [GeV^{2}/c^{4}]");
  dalitz->SetStats(0);
  dalitz_neutral.push_back(dalitz);*/
TH2D* dalitz = createDalitzHistogram("hDalitzMC",
      "Dalitz Plot K^{-}#pi^{+} vs. #pi^{-}#pi^{+} (MC)", massval, datatreeentries);
  dalitz->SetXTitle("mass^{2}(#pi^{-}#pi^{+}) [GeV^{2}/c^{4}]");
  dalitz->SetYTitle("mass^{2}(K^{-}#pi^{+}) [GeV^{2}/c^{4}]");
  dalitz->SetStats(0);
  dalitz_neutral.push_back(dalitz);
  dalitz = createDalitzHistogram("hDalitzData",
      "Dalitz Plot K^{-}#pi^{+} vs. #pi^{-}#pi^{+} (MC)", massval, datatreeentries);
  dalitz->SetXTitle("mass^{2}(#pi^{-}#pi^{+}) [GeV^{2}/c^{4}]");
  dalitz->SetYTitle("mass^{2}(K^{-}#pi^{+}) [GeV^{2}/c^{4}]");
  dalitz->SetStats(0);
  dalitz_neutral.push_back(dalitz);
  /*std::vector<TH2D*> dalitz_charged;
  dalitz = new TH2D("hDalitzChargedMC", "Dalitz Plot #pi^{-}#pi^{0} vs. #pi^{-}#pi^{0} (MC)", nbninsm, 0.0,
      2.2, nbninsm, 0.0, 2.2);
  dalitz_charged.push_back(dalitz);
  dalitz = new TH2D("hDalitzChargedData", "Dalitz Plot #pi^{-}#pi^{0} vs. #pi^{-}#pi^{0} (Data)", nbninsm,
      0.0, 2.2, nbninsm, 0.0, 2.2);
  dalitz_charged.push_back(dalitz);*/

  // --------------- generate histogram bunches
  // K-pi+ isobar
  // GJ Histogram Bunch
  GJHistBunch GJHB_Kpi_isobar = GJHistBunchFactory("Kpi");
  // Helicity Histogram Bunch
  HelicityHistBunch HHB_Kpi_isobar = HelicityHistBunchFactory("Kpi");

  // pi-pi+ isobar
  // GJ Histogram Bunch MC
  GJHistBunch GJHB_pipi_isobar = GJHistBunchFactory("pipi");
  // Helicity Histogram Bunch
  HelicityHistBunch HHB_pipi_isobar = HelicityHistBunchFactory("pipi");

  double avweight = 1;

  /*TH1D* hCosTheta[2] = {new TH1D("hCosThetaMC", "Cos Theta (Data)", nbninsm, -1, 1), new TH1D("hCosThetaData", "Cos Theta (MC)", nbninsm, -1, 1) };
  TH1D* foohTY[2] = {new TH1D("foohTYMC_", " Isobar Treiman-Yang Phi (MC)", nbinsang, -TMath::Pi(),
      TMath::Pi()), new TH1D("foohTYData_", " Isobar Treiman-Yang Phi (DATA)", nbinsang, -TMath::Pi(),
      TMath::Pi()) };
  hCosTheta[0]->Sumw2();
  foohTY[0]->Sumw2();

  TH1D* hhCosTheta[2] = {new TH1D("hhCosThetaMC", "Cos Theta (Data)", nbninsm, -1, 1), new TH1D("hhCosThetaData", "Cos Theta (MC)", nbninsm, -1, 1) };
    TH1D* hfoohTY[2] = {new TH1D("hfoohTYMC_", " Isobar Treiman-Yang Phi (MC)", nbinsang, -TMath::Pi(),
        TMath::Pi()), new TH1D("hfoohTYData_", " Isobar Treiman-Yang Phi (DATA)", nbinsang, -TMath::Pi(),
        TMath::Pi()) };
    hhCosTheta[0]->Sumw2();
    hfoohTY[0]->Sumw2();*/

  //Loop both over data and mc tree
  // itree = 0: mc tree
  // itree = 1: data tree
  for (unsigned int itree = 0; itree < 2; ++itree) {
    TTree* tr = (itree == 0) ? mctr : datatr;
    if (tr == NULL)
      continue;

    double weight = 1;
    double impweight = 1;
    double maxweight = 0;
    TClonesArray* p = new TClonesArray("TLorentzVector");
    TLorentzVector* beam = NULL;
    UInt_t qbeam = -1;
    std::vector<int>* q = NULL;
    if (itree == 0)
      tr->SetBranchAddress("weight", &weight);
    if (itree == 0)
      tr->SetBranchAddress("impweight", &impweight);
    tr->SetBranchAddress("p", &p);
    tr->SetBranchAddress("beam", &beam);
    tr->SetBranchAddress("qbeam", &qbeam);
    tr->SetBranchAddress("q", &q);

    TVector3 vertex;

    int _qbeam = -1;
    NParticleEvent event(p, q, beam, &_qbeam, &vertex);

    // loop over tree entries
    unsigned int nevt = tr->GetEntries();
    for (unsigned int i = 0; i < nevt-1; ++i) {
      //int percent = i*100/nevt;
      //	 if (percent % 10 == 0)
      //	 	std::cout << i << "% " << std::endl; 
      tr->GetEntry(i);
      // in case its data tree (itree=1) put weights to 1
      if (itree == 1) {
        weight = 1;
        impweight = 1;
      }
      if (weight == 0)
        continue;

      _qbeam = -1;// (int) qbeam;
      // this builds all subsystems
      event.refresh();
      //cout << i << " of " << nevt << endl;
      if (impweight != 0)
        weight /= impweight;

      if (weight > maxweight)
        maxweight = weight;
      if (itree == 0)
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
      /*
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
      */
// now fill corresponding dalitz plots
      // this part is pi-pi0pi0 specific (old)
	  // K-pi+pi- (new)

      // put K- pi+ and then pi- pi+ combination to front of vector to make combination picks easier later on
      for(unsigned int i = 0; i < comb.size(); i++) {
		// first pass put the charged combination to the end
        if(comb[i].first.first < 0 && comb[i].first.second < 0) {
		  if (i==2) break;
		  else {
			std::pair<std::pair<float, float>, double> temp = comb[2];
            comb[2] = comb[i];
            comb[i] = temp;
		  }
		}
      }
	  // second pass swap the havier combination to the front	
      if (fabs(comb[1].first.first)+fabs(comb[1].first.second) > fabs(comb[0].first.first)+fabs(comb[0].first.second)){
	    std::pair<std::pair<float, float>, double> temp = comb[0];
	    comb[0] = comb[1];
	    comb[1] = temp;
	  }
      // now actually fill the histograms
      dalitz_neutral[itree]->Fill(comb[1].second, comb[0].second, weight);

     /* unsigned int index = 0;
      TLorentzVector fsParticles[3];
      TLorentzVector beam = *event.beam();
      TLorentzVector X(0, 0, 0, 0);
      for (unsigned int is = 0; is < 3; ++is) {
        if (event.getParticle(is).q() == -1)
          index = is;
        fsParticles[is] = event.getParticle(is).p();
        X += fsParticles[is];
      }
      TVector3 xboost = X.BoostVector();
      for (unsigned int is = 0; is < 3; ++is) {
        fsParticles[is].Boost(-xboost);
      }
      beam.Boost(-xboost);
      TVector3 zaxis = beam.Vect().Unit();
      TVector3 yaxis = beam.Vect().Cross(X.Vect()).Unit();
      TVector3 xaxis = yaxis.Cross(zaxis);

      double phi = 0.0;
      double theta = zaxis.Angle(fsParticles[index].Vect());
      double cosTheta = TMath::Cos(theta);
      double fX = xaxis.Dot(fsParticles[index].Vect());
      double fY = yaxis.Dot(fsParticles[index].Vect());
      if (!(fX == 0.0 && fY == 0.0))
        phi = TMath::ATan2(fY, fX);

      hCosTheta[itree]->Fill(cosTheta, weight);
      foohTY[itree]->Fill(phi, weight);

      TLorentzVector isobar(0, 0, 0, 0);
      unsigned int hindex = 0;
      for (unsigned int is = 0; is < 3; ++is) {
        if (is == index)
          continue;
        isobar += fsParticles[is];
        hindex = is;
      }

      // create helicity frame coordinate system
      TVector3 hzaxis = isobar.Vect().Unit();
      TVector3 hyaxis = zaxis.Cross(hzaxis).Unit();
      TVector3 hxaxis = hyaxis.Cross(hzaxis);

      // boost NParticleState into isobar rest frame
      const TVector3 hboost = isobar.BoostVector();
      fsParticles[hindex].Boost(-hboost);
      // theta
      double htheta = hzaxis.Angle(fsParticles[hindex].Vect());
      double hcosTheta = TMath::Cos(htheta);
      double hfX = hxaxis.Dot(fsParticles[hindex].Vect());
      double hfY = hyaxis.Dot(fsParticles[hindex].Vect());
      double hphi = 0.0;
      if (!(hfX == 0.0 && hfY == 0.0))
        hphi = TMath::ATan2(hfY, hfX);

      hhCosTheta[itree]->Fill(hcosTheta, weight);
      hfoohTY[itree]->Fill(hphi, weight);*/

      // transform into GJ
      event.toGJ();
      // build again to get particles in GJF
      event.build();

      /*TLorentzVector temp(0,0,0,0);
      //cout<<"number of particles in event: "<<event.nParticles()<<endl;
      for (unsigned int is = 0; is < 3; ++is) {
        temp += event.getParticle(is).p();
      }
      temp.Dump();*/

      // loop over all states that contain n-1 final state particles
      // and plot angles
      unsigned int nstates = event.nStates();
      //cout<<"number of substates: "<<event.nStates()<<endl;
      for (unsigned int is = 0; is < nstates; ++is) {

        const NParticleState& state = event.getState(is);
        if (state.n() == npart) {
          hM[itree]->Fill(state.p().M());
        
	}
	if (state.n() == 2 && state.q() == 0){
	  float sum_mass = state.getParticle(0)->p().M()+state.getParticle(1)->p().M();
	  if (fabs(sum_mass-0.633) < 0.010) {
	    // this is an isobar decaying into K-pi+ final state particles
	    fillWeightedGJAnglePlots(state.p(), weight, tprime, itree, GJHB_Kpi_isobar);
	    fillWeightedHelicityAnglePlots(calculateHelicityAngles(state), weight, itree, HHB_Kpi_isobar);
	  }
	  else { 
	    if (fabs(sum_mass-0.279) < 0.010) {
	      // this is an isobar decaying into pi-pi+ final state particles
	      fillWeightedGJAnglePlots(state.p(), weight, tprime, itree, GJHB_pipi_isobar);
	      fillWeightedHelicityAnglePlots(calculateHelicityAngles(state), weight, itree, HHB_pipi_isobar);
	    }
	  }
	}
      }

    }// end loop over events
    if (itree == 0)
      avweight /= (double) nevt;
    cout << "Maxweight=" << maxweight << endl;
    cout << "Average weight=" << avweight << endl;
  }// end loop over trees

  outfile->Write();
  makeDifferencePlots(outfile);

  TList* Hlist = gDirectory->GetList();
  Hlist->Remove(mctr);
  Hlist->Remove(datatr);
  //Hlist->Remove("hWeights");
  int nobj = Hlist->GetEntries();
  std::cout << "Found " << nobj << " Objects in HList" << std::endl;

  outfile->Close();

  gROOT->cd();

}
