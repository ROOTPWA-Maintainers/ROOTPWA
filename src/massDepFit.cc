///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
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
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      fitting program for massdependent fit rootpwa
//      minimizes massDepFitLikeli function
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <complex>
#include <cassert>
#include <time.h>

#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "TComplex.h"
#include "TRandom3.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

#include "reportingUtils.hpp"
#include "fitResult.h"
#include "pwacomponent.h"
#include "massDepFitLikeli.h"


using namespace std;
using namespace ROOT::Math;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{
  cerr << "usage:" << endl
       << progName
       << " -l # -u # -i inputfile [-o outfile -n normfile"
       << "  -M minimizer [-m algorithm] -t # -q -h]" << endl
       << "    where:" << endl
       << "        -l #       lower edge of mass bin [MeV/c^2]" << endl
       << "        -u #       upper edge of mass bin [MeV/c^2]" << endl
       << "        -i file    path to input file" << endl
       << "        -d dir     path to directory with decay amplitude files (default: '.')" << endl
       << "        -o file    path to output file (default: 'fitresult.root')" << endl
       << "        -n file    path to normalization integral file (default: 'norm.int')" << endl
    //       << "        -r #       rank of spin density matrix (default: 1)" << endl
       << "        -M name    minimizer (default: Minuit2)" << endl
       << "        -m name    minimization algorithm (optional, default: Migrad)" << endl
       << "                   available minimizers: Minuit:      Migrad, Simplex, Minimize, Migrad_imp" << endl
       << "                                         Minuit2:     Migrad, Simplex, Combined, Scan, Fumili" << endl
       << "                                         GSLMultiMin: ConjugateFR, ConjugatePR, BFGS, BFGS2, SteepestDescent" << endl
       << "                                         GSLMultiFit: -" << endl
       << "                                         GSLSimAn:    -" << endl
       << "                                         Linear:      Robust" << endl
       << "                                         Fumili:      -" << endl
       << "        -t #       minimizer tolerance (default: 1e-10)" << endl
       << "        -q         run quietly (default: false)" << endl
       << "        -h         print help" << endl
       << endl;
  exit(errCode);
}

//  function loops through fitResults and puts phasespace values into a graph for interpolation
TGraph*
getPhaseSpace(TTree* tree, const std::string& wave){
  unsigned int n=tree->GetEntries();
  TGraph* graph=new TGraph(n);
  fitResult* res=0;
  tree->SetBranchAddress("fitResult_v2",&res);
  for(unsigned int i=0; i<n; ++i){
    tree->GetEntry(i);
    double m=res->massBinCenter();
    double ps=res->phaseSpaceIntegral(wave);
    graph->SetPoint(i,m,ps);
  }
  return graph;
}


int
main(int    argc,
     char** argv)
{
  // --------------------------------------------------------------------------
   // internal parameters
  const string       valTreeName         = "pwa";
  const string       valBranchName       = "fitResult_v2";
  // double             defaultStartValue   = 0.01;
//   bool               useFixedStartValues = false;
//   double             startValStep        = 0.0005;
  const unsigned int maxNmbOfIterations  = 20000;
  const bool         runHesse            = true;
  const bool         runMinos            = false;
  //unsigned int maxParNameLength = 20;       // maximum length of parameter names
//   int                startValSeed        = 1234567;
 // parse command line options
  const string progName           = argv[0];
  double       massBinMin         = 0;                      // [MeV/c^2]
  double       massBinMax         = 5000;                      // [MeV/c^2]
  
  string       inFileName         = "fitresult.root";       // input filename
  string       outFileName        = "mDep.result.root";       // output filename
  string       normIntFileName    = "";                     // file with normalization integrals
  string       minimizerType[2]   = {"Minuit2", "Migrad"};  // minimizer, minimization algorithm
  double       minimizerTolerance = 1e-10;                  // minimizer tolerance
  bool         quiet              = false;
  
extern char* optarg;
  // extern int optind;
  int ca;
  while ((ca = getopt(argc, argv, "l:u:i:o:n:r:M:m:t:qh")) != -1)
    switch (ca) {
    case 'l':
      massBinMin = atof(optarg);
      break;
    case 'u':
      massBinMax = atof(optarg);
      break;
     case 'o':
      outFileName = optarg;
      break;
     case 'i':
      inFileName = optarg;
      break;
    case 'n':
      normIntFileName = optarg;
      break;
    case 'M':
      minimizerType[0] = optarg;
      break;
    case 'm':
      minimizerType[1] = optarg;
      break;
    case 't':
      minimizerTolerance = atof(optarg);
      break;
    case 'q':
      quiet = true;
      break;
    case 'h':
      usage(progName);
      break;
    }


  printInfo << "creating and setting up likelihood function" << endl;
  

  // open input file and get results tree
  TFile* infile=TFile::Open(inFileName.c_str());
  if(infile==NULL){
    cerr << "Input file " << inFileName <<" not found."<< endl;
    return 1;
  }
  TTree* tree=(TTree*)infile->Get(valTreeName.c_str());
  if(tree==NULL){
    cerr << "Input tree " << valTreeName <<" not found."<< endl;
    return 1;
  }



  std::map<std::string,pwachannel > channels;
  std::string ch1="1-2-+0+rho770_02_a21320=pi-_2_rho770.amp";
  std::string ch2="1-2-+0+pi-_02_f21270=pi-+_1_a11269=pi+-_0_rho770.amp";
  std::string ch3="1-2-+0+rho31690=rho770_03_f21270_13_pi-.amp";
  std::string ch4="1-2-+0+rho770_02_a11269=pi-_0_rho770.amp";
  channels[ch1]=pwachannel(complex<double>(50,0),getPhaseSpace(tree,ch1));
  channels[ch2]=pwachannel(complex<double>(50,0),getPhaseSpace(tree,ch2));
  //channels[ch3]=pwachannel(complex<double>(5,0),getPhaseSpace(tree,ch3));
  channels[ch4]=pwachannel(complex<double>(50,0),getPhaseSpace(tree,ch4));
  pwacomponent comp1("pi2(1880)",1880,150,channels);
  comp1.setLimits(1870,1900,50,300);
  comp1.setFixed(0,0);

  channels["1-2-+0+rho770_02_a21320=pi-_2_rho770.amp"].setCoupling(complex<double>(1,0));
  channels["1-2-+0+pi-_02_f21270=pi-+_1_a11269=pi+-_0_rho770.amp"].setCoupling(complex<double>(5,0));
  pwacomponent comp2("pi2(2300)",2300,500,channels);
  comp2.setLimits(2250,2500,100,800);
  comp2.setFixed(1,0);

  channels["1-2-+0+rho770_02_a21320=pi-_2_rho770.amp"].setCoupling(complex<double>(50,0));
  channels["1-2-+0+pi-_02_f21270=pi-+_1_a11269=pi+-_0_rho770.amp"].setCoupling(complex<double>(50,0));
  pwacomponent comp3("pi2(1670)",1672,260,channels);
  comp3.setLimits(1600,1700,200,400);
  comp3.setFixed(1,1);

  channels["1-2-+0+rho770_02_a21320=pi-_2_rho770.amp"].setCoupling(complex<double>(20,0));
  channels["1-2-+0+pi-_02_f21270=pi-+_1_a11269=pi+-_0_rho770.amp"].setCoupling(complex<double>(5,0));
  pwacomponent comp4("pi2(2100)",2090,200,channels);
  comp4.setLimits(2000,2200,100,800);
  comp4.setFixed(1,0);


  std::map<std::string,pwachannel > channels0mp;
  std::string ch5="1-0-+0+pi-_00_f01500=rho770_00_rho770.amp";
  std::string ch6="1-0-+0+rho770_00_a11269=pi-_0_rho770.amp";
  
  
  channels0mp[ch5]=pwachannel(complex<double>(50,0),getPhaseSpace(tree,ch5));
  channels0mp[ch6]=pwachannel(complex<double>(50,0),getPhaseSpace(tree,ch5));
  channels0mp[ch6].setCoupling(std::complex<double>(0.1,0));
  pwacomponent comp5("pi(1800)",1840,208,channels0mp);
  comp5.setLimits(1750,1850,100,300);
  comp5.setFixed(0,0);
  //channels0mp.clear();
  channels0mp[ch6]=pwachannel(complex<double>(50,0),getPhaseSpace(tree,ch5));
  channels0mp[ch5].setCoupling(std::complex<double>(0.5,0));
  channels0mp[ch6].setCoupling(std::complex<double>(50,0));
  pwacomponent comp6("pi(1700)",1700,200,channels0mp);
  comp6.setLimits(1600,1750,100,500);
  comp6.setFixed(0,0);
  pwacomponent comp13("pi(2100)",2200,200,channels0mp);
  comp13.setLimits(2000,2300,100,500);
  comp13.setFixed(1,0);

  channels0mp.clear();
  channels0mp[ch6]=pwachannel(complex<double>(50,0),getPhaseSpace(tree,ch6));
  pwabkg bkg1=pwabkg("0-+Bkg",0,0.0001,channels0mp);
  bkg1.setIsobars(770,1269);
  bkg1.setLimits(0,0,0,1);
  bkg1.setFixed(1,0);

  std::map<std::string,pwachannel> channels2mp;
  channels2mp[ch4]=pwachannel(complex<double>(50,0),getPhaseSpace(tree,ch4));
  pwabkg bkg2=pwabkg("rhoa1bkg",0,0.0001,channels2mp);
  bkg2.setIsobars(770,1269);
  bkg2.setLimits(0,0,0,1);
  bkg2.setFixed(1,0);

  channels2mp.clear();
  channels2mp[ch3]=pwachannel(complex<double>(50,0),getPhaseSpace(tree,ch3));
  pwabkg bkg3=pwabkg("pirho3",0,0.0001,channels2mp);
  bkg3.setIsobars(139,1690);
  bkg3.setLimits(0,0,0,1);
  bkg3.setFixed(1,0);

  std::map<std::string,pwachannel> channels1pp;
  std::string anchorwave("1-1++0+sigma_01_a11269=pi-_0_rho770.amp");
  string ch1ppsa1("1-1++0+sigma_01_a11269=pi-_0_rho770.amp");
  string ch1ppra11("1-1++0+rho770_11_a11269=pi-_0_rho770.amp");
  string ch1ppra12("1-1++0+rho770_12_a11269=pi-_0_rho770.amp");
  string ch1ppf1pi("1-1++0+pi-_11_f11285=pi-+_11_a11269=pi+-_0_rho770.amp");
		  
  //  channels1pp[ch1ppsa1]=pwachannel(complex<double>(50,0),
  //				   getPhaseSpace(tree,ch1ppsa1));
//channels1pp[ch1ppra11]=pwachannel(complex<double>(50,0),
  //			   getPhaseSpace(tree,ch1ppra11));
// channels1pp[ch1ppra12]=pwachannel(complex<double>(50,0),
  //			   getPhaseSpace(tree,ch1ppra12));
  channels1pp[ch1ppf1pi]=pwachannel(complex<double>(50,0),
				   getPhaseSpace(tree,ch1ppf1pi));

  pwacomponent comp7("a1(1640)",1647,254,channels1pp);
  comp7.setLimits(1600,1700,100,500);
  comp7.setFixed(0,0);
  pwacomponent comp8("a1(1930)",1930,155,channels1pp);
  comp8.setLimits(1900,1980,100,300);
  comp8.setFixed(1,1);
  pwacomponent comp9("a1(2095)",2069,451,channels1pp);
  comp9.setLimits(1950,2150,100,300);
  comp9.setFixed(0,0);
  
  channels1pp.clear();
  channels1pp[ch1ppsa1]=pwachannel(complex<double>(50,0),
				   getPhaseSpace(tree,ch1ppsa1));
  pwabkg bkg4=pwabkg("sigmaa1",0,0.0001,channels1pp);
  bkg4.setIsobars(600,1269);
  bkg4.setLimits(0,0,0,1);
  bkg4.setFixed(1,0);

  channels1pp.clear();
  channels1pp[ch1ppra11]=pwachannel(complex<double>(50,0),
				   getPhaseSpace(tree,ch1ppra11));
  pwabkg bkg5=pwabkg("1pprhoa1bkg",0,0.0001,channels1pp);
  bkg5.setIsobars(770,1269);
  bkg5.setLimits(0,0,0,1);
  bkg5.setFixed(1,0);

  pwacompset compset;
  compset.add(&comp1); // pi2(1880)
  //compset.add(&bkg2);  
  //compset.add(&comp2); // pi2(2300)
  compset.add(&comp3); // pi2(1670)
  //compset.add(&bkg3);
  compset.add(&comp4); // pi2(2100)
  compset.add(&comp5); // pi(1800)
  compset.add(&comp6); // pi(1700)
 compset.add(&comp13); // pi(2100)
  compset.add(&bkg1);
  //compset.add(&comp7);
  //compset.add(&comp8);
  //compset.add(&comp9);
  //compset.add(&bkg4);
  


  massDepFitLikeli L;
  L.init(tree,&compset,massBinMin,massBinMax);

  const unsigned int nmbPar  = L.NDim();
  
  printInfo << nmbPar << " Parameters in fit" << endl;

  // ---------------------------------------------------------------------------
  // setup minimizer
  printInfo << "creating and setting up minimizer " << minimizerType[0] << " using algorithm " << minimizerType[1] << endl;
  Minimizer* minimizer = Factory::CreateMinimizer(minimizerType[0], minimizerType[1]);
  if (!minimizer) { 
    printErr << "could not create minimizer! exiting!" << endl;
    throw;
  }
  minimizer->SetFunction(L);
  minimizer->SetPrintLevel((quiet) ? 0 : 3);

  // ---------------------------------------------------------------------------

  // Set startvalues
  unsigned int parcount=0;
  for(unsigned int ic=0;ic<compset.n();++ic){
    const pwacomponent& comp=*compset[ic];
    TString name(comp.name());
    double mmin,mmax,gmin,gmax;
    comp.getLimits(mmin,mmax,gmin,gmax);
    if(comp.fixM())minimizer->SetFixedVariable(parcount++,
					       (name+"_M").Data() ,
					       comp.m0());
    else minimizer->SetLimitedVariable(parcount++, 
				       (name+"_M").Data(), 
				       comp.m0(), 
				       0.10,
				       mmin,mmax);
    if(comp.fixGamma())minimizer->SetFixedVariable(parcount++,
						   (name+"_Gamma").Data() ,
						   comp.gamma());
    else minimizer->SetLimitedVariable(parcount++, 
				       (name+"_Gamma").Data(), 
				       comp.gamma(), 
				       0.01,
				       gmin,gmax);
    std::map<std::string,pwachannel >::const_iterator it=comp.channels().begin();
    while(it!=comp.channels().end()){
      minimizer->SetVariable(parcount++,(name + "_ReC" + it->first).Data() , it->second.C().real(), 0.10);
      
      // fix one phase
      if(it->first==anchorwave)minimizer->SetFixedVariable(parcount++,(name + "_ImC" + it->first).Data() , 0.0);
	
	else {minimizer->SetVariable(parcount++,(name + "_ImC" + it->first).Data() , it->second.C().imag(), 0.10);}
      
      ++it;
    } // end loop over channels

  }
 



  // find minimum of likelihood function
  printInfo << "performing minimization." << endl;
  {
    minimizer->SetMaxIterations(maxNmbOfIterations);
    minimizer->SetMaxFunctionCalls(maxNmbOfIterations*5);
    minimizer->SetTolerance    (minimizerTolerance);
    bool success = minimizer->Minimize();
    if (success)
      printInfo << "minimization finished successfully." << endl;
    else
      printWarn << "minimization failed." << endl;
    if (runHesse) {
      printInfo << "calculating Hessian matrix." << endl;
      success = minimizer->Hesse();  // comes only with ROOT 5.24+
      if (!success)
	printWarn << "calculation of Hessian matrix failed." << endl;
    }
    printInfo << "minimization stopped after " << minimizer->NCalls() << " function calls. minimizer status summary:" << endl
	      << "    total number of parameters .......................... " << minimizer->NDim()             << endl
	      << "    number of free parameters ........................... " << minimizer->NFree()            << endl
	      << "    maximum allowed number of iterations ................ " << minimizer->MaxIterations()    << endl
	      << "    maximum allowed number of function calls ............ " << minimizer->MaxFunctionCalls() << endl
	      << "    minimizer status .................................... " << minimizer->Status()           << endl
	      << "    minimizer provides error and error matrix ........... " << minimizer->ProvidesError()    << endl
	      << "    minimizer has performed detailed error validation ... " << minimizer->IsValidError()     << endl
	      << "    estimated distance to minimum ....................... " << minimizer->Edm()              << endl
	      << "    statistical scale used for error calculation ........ " << minimizer->ErrorDef()         << endl
	      << "    minimizer strategy .................................. " << minimizer->Strategy()         << endl
	      << "    absolute tolerance .................................. " << minimizer->Tolerance()        << endl;
  }

  // ---------------------------------------------------------------------------
  // print results
  printInfo << "minimization result:" << endl;
  for (unsigned int i = 0; i< nmbPar; ++i) {
    cout << "    parameter [" << setw(3) << i << "] ";
      //	 << setw(maxParNameLength); //<< L.parName(i) << " = ";
    //if (parIsFixed[i])
    //  cout << minimizer->X()[i] << " (fixed)" << endl;
    //else {
      cout << setw(12) << maxPrecisionAlign(minimizer->X()[i]) << " +- "
	   << setw(12) << maxPrecisionAlign(minimizer->Errors()[i]);
      if (runMinos && (i == 156)) {  // does not work for all parameters
	double minosErrLow = 0;
	double minosErrUp  = 0;
	const bool success = minimizer->GetMinosError(i, minosErrLow, minosErrUp);
	if (success)
	  cout << "    Minos: " << "[" << minosErrLow << ", +" << minosErrUp << "]" << endl;
      } else
	cout << endl;
  }

  const double* par=minimizer->X();
  // matrix of couplings

  compset.setPar(par);

  cerr << compset << endl;

  cerr << "Fitting finished... Start building graphs ... " << endl;

   std::vector<std::string> wl=compset.wavelist();
   std::map<std::string, unsigned int> wmap;
   unsigned int ndatabins=tree->GetEntries();

   std::vector<TGraphErrors*> datagraphs;
   std::vector<TMultiGraph*> graphs;
   for(unsigned int iw=0; iw<wl.size();++iw){
     wmap[wl[iw]]=iw;
     graphs.push_back(new TMultiGraph);
     graphs[iw]->SetName(wl[iw].c_str());
     graphs[iw]->SetTitle(wl[iw].c_str());
     graphs[iw]->SetDrawOption("AP");
     datagraphs.push_back(new TGraphErrors(ndatabins));
     string name("data_");name.append(wl[iw]);
     datagraphs[iw]->SetName(name.c_str());
     datagraphs[iw]->SetTitle(name.c_str());
     datagraphs[iw]->SetDrawOption("AP");
     //datagraphs[iw]->SetLineColor(kRed);
     //datagraphs[iw]->SetMarkerColor(kRed);

     graphs[iw]->Add(datagraphs[iw],"P");
   }

 


     // build fitgraphs
   unsigned int nbins=ndatabins;//200;
   //double mmin=1200.;
   //double md=10.;
   std::vector<TGraph*> fitgraphs;
  

   for(unsigned int iw=0; iw<wl.size();++iw){
     fitgraphs.push_back(new TGraph(nbins));
     string name("fit_");name.append(wl[iw]);
     fitgraphs[iw]->SetName(name.c_str());
     fitgraphs[iw]->SetTitle(name.c_str());
     fitgraphs[iw]->SetLineColor(kRed);
     fitgraphs[iw]->SetLineWidth(2);
     fitgraphs[iw]->SetMarkerColor(kRed);
     fitgraphs[iw]->SetDrawOption("AP");
     //fitgraphs[iw]->SetMarkerStyle(22);
     graphs[iw]->Add(fitgraphs[iw],"cp");
     graphs[iw]->Add(getPhaseSpace(tree,wl[iw]));
   }

   std::vector<TGraph*> compgraphs; // individual components
   // loop over components and build graphs
     for(unsigned int ic=0;ic<compset.n();++ic){
       const pwacomponent* c=compset[ic];
       std::map<std::string,pwachannel >::const_iterator it=c->channels().begin();
       while(it!=c->channels().end()){
	 string name=c->name();name.append("__");
	 name.append(it->first);
	 TGraph* gcomp=new TGraph(nbins);
	 gcomp->SetName(name.c_str());
	 gcomp->SetTitle(name.c_str());
	 unsigned int color=kBlue;
	 if(dynamic_cast<const pwabkg*>(c)!=NULL)color=kMagenta;
	 gcomp->SetLineColor(color);
	 gcomp->SetMarkerColor(color);
	 
	 compgraphs.push_back(gcomp);
	 graphs[wmap[it->first]]->Add(gcomp,"cp");
	 ++it;
       }// end loop over channels
      
     }// end loop over components

   std::vector<TGraphErrors*> phasedatagraphs;
   std::vector<TMultiGraph*> phasegraphs;
   std::vector<TGraph*> phasefitgraphs;
   unsigned int c=0;

  for(unsigned int iw=0; iw<wl.size();++iw){
     for(unsigned int iw2=iw+1; iw2<wl.size();++iw2){
       phasegraphs.push_back(new TMultiGraph);
       string name("dPhi_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       phasegraphs[c]->SetName(name.c_str());
       phasegraphs[c]->SetTitle(name.c_str());
       phasedatagraphs.push_back(new TGraphErrors(3*nbins));
       name=("dPhi_data_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       phasedatagraphs[c]->SetName(name.c_str());
       phasedatagraphs[c]->SetTitle(name.c_str());
       phasegraphs[c]->Add(phasedatagraphs[c],"p");
       ++c;
     }
   }

   c=0;
   for(unsigned int iw=0; iw<wl.size();++iw){
     for(unsigned int iw2=iw+1; iw2<wl.size();++iw2){
       phasefitgraphs.push_back(new TGraph(nbins));
       string name("dPhi_fit_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       phasefitgraphs[c]->SetName(name.c_str());
       phasefitgraphs[c]->SetTitle(name.c_str());
       phasefitgraphs[c]->SetLineColor(kRed);
       phasefitgraphs[c]->SetMarkerColor(kRed);
       phasefitgraphs[c]->SetDrawOption("AP");
       phasefitgraphs[c]->SetMarkerStyle(22);
       phasegraphs[c]->Add(phasefitgraphs[c],"cp");
       ++c;
     }
   }
   
   std::vector<TGraph2D*> phase2d;
   c=0;
   for(unsigned int iw=0; iw<wl.size();++iw){
     for(unsigned int iw2=iw+1; iw2<wl.size();++iw2){
       phase2d.push_back(new TGraph2D(nbins));
       string name("dPhi_2d_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       phase2d[c]->SetName(name.c_str());
       phase2d[c]->SetTitle(name.c_str());
       phase2d[c]->SetLineColor(kRed);
       phase2d[c]->SetMarkerColor(kRed);
       //phasegraphs[c]->Add(phasefitgraphs[c],"cp");
       ++c;
     }
   }




  // get data
   fitResult* rho=0;
   tree->SetBranchAddress(valBranchName.c_str(),&rho);
   vector<double> prevps(wl.size());
   double mprev=0;
   double binwidth=30; // half binwidth
   //double w=2*30/10;
   for(unsigned int i=0;i<ndatabins;++i){
     tree->GetEntry(i);
     double m=rho->massBinCenter();
     unsigned int c=0;
     for(unsigned int iw=0; iw<wl.size();++iw){
       double ps=rho->phaseSpaceIntegral(wl[iw].c_str());
       datagraphs[iw]->SetPoint(i,m,rho->intensity(wl[iw].c_str()));
       datagraphs[iw]->SetPointError(i,binwidth,rho->intensityErr(wl[iw].c_str()));
       fitgraphs[iw]->SetPoint(i,m,compset.intensity(wl[iw],m)*ps*ps);           
       // second loop to get phase differences
       for(unsigned int iw2=iw+1; iw2<wl.size();++iw2){
	 double ps2=rho->phaseSpaceIntegral(wl[iw2].c_str());
	 
	 phasedatagraphs[c]->SetPoint(i,m,rho->phase(wl[iw].c_str(),
	 					     wl[iw2].c_str()));
	 phasedatagraphs[c]->SetPoint(i+ndatabins,m,rho->phase(wl[iw].c_str(),
	 				     wl[iw2].c_str())-360);
	 phasedatagraphs[c]->SetPoint(i+2*ndatabins,m,rho->phase(wl[iw].c_str(),
	 					     wl[iw2].c_str())+360);
	 TVector2 v;v.SetMagPhi(1,rho->phase(wl[iw].c_str(),
					     wl[iw2].c_str())/TMath::RadToDeg());

	 //phasedatagraphs[c]->SetPoint(i,m*v.X(),m*v.Y());

	 phase2d[c]->SetPoint(i,v.X(),v.Y(),m);

	 
	 phasedatagraphs[c]->SetPointError(i,binwidth,
					   rho->phaseErr(wl[iw].c_str(),
							 wl[iw2].c_str()));
	 phasefitgraphs[c]->SetPoint(i,m,compset.phase(wl[iw],ps,
						       wl[iw2],ps2,m)*TMath::RadToDeg());
	 c++;
       }// end inner loop over waves

       prevps[iw]=ps;
       
        
     } // end loop over waves

     // loop over components to fill individual graphs
     unsigned int compcount=0;
       for(unsigned int ic=0;ic<compset.n();++ic){
	 const pwacomponent* c=compset[ic];
	 std::map<std::string,pwachannel >::const_iterator it=c->channels().begin();
	 while(it!=c->channels().end()){
	   double I=norm(c->val(m)*it->second.C())*it->second.ps(m)*it->second.ps(m);
	   compgraphs[compcount]->SetPoint(i,m,I);
	   ++compcount;
	   ++it;
	 } // end loop over channels
       }// end loop over components

     mprev=m;
   }




//    for(unsigned int im=0;im<nbins;++im){ // fine loop in masses -> fits
     
//      double m=mmin+im*md;
//      for(unsigned int iw=0; iw<wl.size();++iw){
//        fitgraphs[iw]->SetPoint(im,m,compset.intensity(wl[iw],m));
//      } // end loop over waves
     
     
//    }// end loop over mass bins
   
   



   TFile* outfile=TFile::Open(outFileName.c_str(),"RECREATE");
   for(unsigned int iw=0; iw<wl.size();++iw){
     TGraph* g=(TGraph*)graphs[iw]->GetListOfGraphs()->At(2);
     for(unsigned int ib=0;ib<nbins;++ib){
       double m,ps;g->GetPoint(ib,m,ps);
       g->SetPoint(ib,m,ps*2000);
     }
     
     graphs[iw]->Write();
   }
   for(unsigned int iw=0; iw<phasegraphs.size();++iw){

/// rectivfy phase graphs
   
   // unsigned int refbin=4;
//    double m;
//    double predph;
//    phasedatagraphs[iw]->GetPoint(refbin,m,predph);
//    double prefph;
//    phasefitgraphs[iw]->GetPoint(refbin,m,prefph);
//    for(unsigned int ib=refbin+1;ib<nbins;++ib){
//      double dph; phasedatagraphs[iw]->GetPoint(ib,m,dph);
//      double fph; phasefitgraphs[iw]->GetPoint(ib,m,fph);
//      double dp,dm;dp=dph+360;dm=dph-360;
//      double fp,fm;fp=fph+360;fm=fph-360;
//      if(0){
//        if(fabs(dp-predph)<fabs(dph-predph) && fabs(dp-predph)<fabs(dm-predph))
// 	 phasedatagraphs[iw]->SetPoint(ib,m,dp);
//        else if(fabs(dm-predph)<fabs(dph-predph) && fabs(dm-predph)<fabs(dp-predph))
// 	 phasedatagraphs[iw]->SetPoint(ib,m,dm);
       
//        if(fabs(fp-prefph)<fabs(fph-prefph) && fabs(fp-prefph)<fabs(fm-prefph))
// 	 phasefitgraphs[iw]->SetPoint(ib,m,fp);
//        else if(fabs(fm-prefph)<fabs(fph-prefph) && fabs(fm-prefph)<fabs(fp-prefph))
// 	 phasefitgraphs[iw]->SetPoint(ib,m,fm);
       
//        phasedatagraphs[iw]->GetPoint(ib,m,predph);
//        phasefitgraphs[iw]->GetPoint(ib,m,prefph);
//      }
//    }
   
     phasegraphs[iw]->Write();
     //phase2d[iw]->Write();
   }
   outfile->Close();
   
   return 0;
   
}
