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
//      BE CAREFULL: Apart from the dynamic width the sqrts of the phase space factors are needed
//                   So at most places in this code ps acctually is the sqrt!!! 
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>
#include <complex>
#include <cassert>
#include <time.h>

#include "TTree.h"
#include "TF1.h"
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

#include "libconfig.h++"

using namespace std;
using namespace libconfig;
using namespace ROOT::Math;
using namespace rpwa;


void
usage(const string& progName,
      const int     errCode = 0)
{
  cerr << "usage:" << endl
       << progName
       << " -c configfile -i inputfile [-o outfile -l # -u #"
       << "  -M minimizer [-m algorithm] -t # -q -h] [-S fitResultFiles]" << endl
       << "    where:" << endl
       << "        -c file    path to config File" << endl
       << "        -i file    path to input file" << endl
       << "        -o file    path to output file (default: 'mDep.result.root')" << endl
    //       << "        -r #       rank of spin density matrix (default: 1)" << endl
       << "        -l # -u #  lower and upper mass range used for fit" << endl
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
       << "        -P         plotting only - no fit" << endl
       << "        -S files   Systematic error plotting. give list of files" << endl

       << endl;
  exit(errCode);
}

//  function loops through fitResults and puts phasespace values into a graph for interpolation
//  THIS CONTAINS NOW THE RIGHT VALUE NOT!!! THE SQRT!!!
TGraph*
getPhaseSpace(TTree* tree, TF1* fsps,const std::string& wave){
  unsigned int n=tree->GetEntries();
  TGraph* graph=new TGraph(n);
  fitResult* res=0;
  tree->SetBranchAddress("fitResult_v2",&res);
  for(unsigned int i=0; i<n; ++i){
    tree->GetEntry(i);
    double m=res->massBinCenter();
    double ps=res->phaseSpaceIntegral(wave);
    ps*=ps; // remember that phaseSpaceIntegral returns sqrt of integral!!! 
    //ps*=fsps->Eval(m);
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
  bool               onlyPlotting        = false;
  bool               sysPlotting         = false;

  //unsigned int maxParNameLength = 20;       // maximum length of parameter names
//   int                startValSeed        = 1234567;
 // parse command line options
  const string progName           = argv[0];
  double       massBinMin         = 0;                      // [MeV/c^2]
  double       massBinMax         = 5000;                      // [MeV/c^2]
  
  string       inFileName         = "fitresult.root";       // input filename
  string       outFileName        = "mDep.result.root";       // output filename
  //string       normIntFileName    = "";                     // file with normalization integrals
  string       minimizerType[2]   = {"Minuit2", "Migrad"};  // minimizer, minimization algorithm
  double       minimizerTolerance = 1e-10;                  // minimizer tolerance
  bool         quiet              = false;
  
  string       configFile;        // configuration file

extern char* optarg;
 extern int optind;
  // extern int optind;
  int ca;
  while ((ca = getopt(argc, argv, "c:i:o:u:l:M:m:t:qhPS")) != -1)
    switch (ca) {
    case 'c':
      configFile = optarg;
      break;
    case 'o':
      outFileName = optarg;
      break;
     case 'i':
      inFileName = optarg;
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
   case 'l':
      massBinMin = atof(optarg);
      break;
   case 'u':
      massBinMax = atof(optarg);
      break;
    case 'q':
      quiet = true;
      break;
    case 'h':
      usage(progName);
      break;
    case 'P':
      onlyPlotting=true;
      break;
    case 'S':
      sysPlotting=true;
      break;
    }


  TF1* fPS=new TF1("fps","[3] * (x-[0])*exp( (x-[0])*([1]+(x-[0])*[2]) )",900,3000);
  fPS->SetParameter(0,698); //5pi threshold
  fPS->SetParLimits(0,698,698);
  fPS->SetParameter(1,6.27068E-3); // slope
  fPS->SetParLimits(1,6.27068E-3,6.27068E-3); // slope
  //fPS->SetParameter(2,-0.939708E-6); // correction
  fPS->SetParameter(2,-2.5E-6); // correction
  fPS->SetParLimits(2,-3E-6,0);
  fPS->SetParameter(3,1E-6); // Normalization
  fPS->SetParLimits(3,1E-6,1E-6); // Normalization
  //fPS->SetParameter(3,10); // Normalization
  //TF1* fPS=new TF1("fps","1.",900,4000);


  vector<TTree*> sysTrees;
  if(sysPlotting){
    // open files with fits  
    for(int i=optind;i<argc;++i){
      // open input file and get results tree
      TFile* infile=TFile::Open(argv[i]);
      if(infile==NULL){
	cerr << "Systematics Input file " << inFileName <<" not found."<< endl;
	return 1;
      }
      TTree* tree=(TTree*)infile->Get(valTreeName.c_str());
      if(tree==NULL){
	cerr << "Input tree " << valTreeName <<" not found."<< endl;
	return 1;
      }
      sysTrees.push_back(tree);
    }
    printInfo << sysTrees.size() << " files for systematics found " << endl;
  }// end if sysPlotting

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

  // Setup Component Set (Resonances + Background)
  pwacompset compset;
  Config Conf;
  Conf.readFile(configFile.c_str());
  const Setting& root = Conf.getRoot();
  bool check=true;
  // Resonances
  if(Conf.exists("components.resonances")){
    const Setting &bws = root["components"]["resonances"];
    // loop through breitwigners
    int nbw=bws.getLength();
    printInfo << "found " << nbw << " Resonances in config" << endl;
    for(int ibw = 0; ibw < nbw; ++ibw) {
      const Setting &bw = bws[ibw];
      string jpc;
      string name;
      double mass=-1;double ml,mu;int mfix; 
      double width=-1;double wl,wu;int wfix;int wdyn;
     
      check&=bw.lookupValue("name",     name);
      check&=bw.lookupValue("jpc",       jpc);
      const Setting &massSet = bw["mass"];
      check&=massSet.lookupValue("val",        mass);
      check&=massSet.lookupValue("lower",        ml);
      check&=massSet.lookupValue("upper",        mu);
      check&=massSet.lookupValue("fix",        mfix);
      const Setting &widthSet = bw["width"];
      check&=widthSet.lookupValue("val",       width);
      check&=widthSet.lookupValue("lower",       wl);
      check&=widthSet.lookupValue("upper",       wu);
      check&=widthSet.lookupValue("fix",       wfix);
      bool checkdyn=widthSet.lookupValue("dyn",       wdyn);
      if(!checkdyn)wdyn=0;
      cout << "---------------------------------------------------------------------" << endl;
      cout << name << "    JPC = " << jpc << endl;
      cout << "mass(limits)  = " << mass <<" ("<<ml<<","<<mu<<") MeV/c^2";
      if(mfix==1)cout<<"  -- FIXED";
      cout<< endl;
      cout << "width(limits) = " << width <<" ("<<wl<<","<<wu<<") MeV/c^2";
      if(wfix==1)cout<<"  -- FIXED";
      if(wdyn!=0)cout<<"  -- DYNAMIC WIDTH";
      else cout<<"  -- CONST WIDTH";
      cout<< endl;
      const Setting &channelSet = bw["decaychannels"];
      unsigned int nCh=channelSet.getLength();
      cout << "Decaychannels (coupling):" << endl;
      std::map<std::string,pwachannel > channels;
      for(unsigned int iCh=0;iCh<nCh;++iCh){
	const Setting &ch = channelSet[iCh];
	string amp;
	double cRe=0;
	double cIm=0;
	check&=ch.lookupValue("amp",amp);
	check&=ch.lookupValue("coupling_Re",cRe);
	check&=ch.lookupValue("coupling_Im",cIm);
	complex<double> C(cRe,cIm);
	cout << "   " << amp << "  " << C << endl;
	channels[amp]=pwachannel(C,getPhaseSpace(tree,fPS,amp));
      }// end loop over channels
      if(!check){
	printErr << "Bad config value lookup! Check your config file!" << endl;
	return 1;
      }
      pwacomponent* comp1=new pwacomponent(name,mass,width,channels);
      comp1->setLimits(ml,mu,wl,wu);
      comp1->setFixed(mfix,wfix);
      if(wdyn==0)comp1->setConstWidth();
      compset.add(comp1);
      cout << "CHECK val(m0)="<< comp1->val(mass) << endl;
    }// end loop over resonances
  }
  cout << endl;
  // Background components
  if(Conf.exists("components.background")){
    const Setting &bws = root["components"]["background"];
    // loop through breitwigners
    int nbw=bws.getLength();
    printInfo << "found " << nbw << " Background components in config" << endl;
    for(int ibw = 0; ibw < nbw; ++ibw) {
      const Setting &bw = bws[ibw];
      string name;
      double mass=-1;double ml,mu;int mfix; 
      double width=-1;double wl,wu;int wfix;
     
      check&=bw.lookupValue("name",     name);
      const Setting &massSet = bw["m0"];
      check&=massSet.lookupValue("val",        mass);
      check&=massSet.lookupValue("lower",        ml);
      check&=massSet.lookupValue("upper",        mu);
      check&=massSet.lookupValue("fix",        mfix);
      const Setting &widthSet = bw["g"];
      check&=widthSet.lookupValue("val",       width);
      check&=widthSet.lookupValue("lower",       wl);
      check&=widthSet.lookupValue("upper",       wu);
      check&=widthSet.lookupValue("fix",       wfix);
      cout << "---------------------------------------------------------------------" << endl;
      cout << name << endl;
      cout << "mass-offset(limits)  = " << mass <<" ("<<ml<<","<<mu<<") MeV/c^2";
      if(mfix==1)cout<<"  -- FIXED";
      cout<< endl;
      cout << "g(limits)            = " << width <<" ("<<wl<<","<<wu<<") MeV/c^2";
      if(wfix==1)cout<<"  -- FIXED";
      cout<< endl;
      std::map<std::string,pwachannel > channels;
      string amp;
      double cRe=0;
      double cIm=0;
      double mIso1=0;
      double mIso2=0;
      check&=bw.lookupValue("amp",amp);
      check&=bw.lookupValue("coupling_Re",cRe);
      check&=bw.lookupValue("coupling_Im",cIm);
      check&=bw.lookupValue("mIsobar1",mIso1);
      check&=bw.lookupValue("mIsobar2",mIso2);
      complex<double> C(cRe,cIm);
      cout << "Decaychannel (coupling):" << endl;
      cout << "   " << amp << "  " << C << endl;
      cout << "   Isobar masses: " << mIso1<<"  "<< mIso2<< endl;
      channels[amp]=pwachannel(C,getPhaseSpace(tree,fPS,amp));

      if(!check){
	printErr << "Bad config value lookup! Check your config file!" << endl;
	return 1;
      }
      pwabkg* bkg=new pwabkg(name,mass,width,channels);
      bkg->setIsobars(mIso1,mIso2);
      bkg->setLimits(ml,mu,wl,wu);
      bkg->setFixed(mfix,wfix);
      compset.add(bkg);
    }// end loop over background
  }// endif

  // add phase space
  compset.setPS(fPS);


 cout << "---------------------------------------------------------------------" << endl << endl;

  // set anchorwave
  vector<string> anchorwave_channel;
  vector<string> anchorwave_reso;
  if(Conf.exists("components.anchorwave")){
    const Setting &anc = root["components"]["anchorwave"];
    // loop through breitwigners
    unsigned int nanc=anc.getLength();
    for(unsigned int ianc=0;ianc<nanc;++ianc){
      string ch,re;
      const Setting &anco = anc[ianc];
      anco.lookupValue("channel",ch);
      anco.lookupValue("resonance",re);
      cout << "Ancorwave: "<< endl;
      cout << "    " << re << endl;
      cout << "    " << ch << endl;
      anchorwave_channel.push_back(ch);
      anchorwave_reso.push_back(re);
    }
  }


    cout << "---------------------------------------------------------------------" << endl << endl;
 
 

  massDepFitLikeli L;
  L.init(tree,fPS,&compset,massBinMin,massBinMax);

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
      if(find(anchorwave_reso.begin(),anchorwave_reso.end(),name)!=anchorwave_reso.end() && find(anchorwave_channel.begin(),anchorwave_channel.end(),it->first)!=anchorwave_channel.end()){
	minimizer->SetFixedVariable(parcount++,(name + "_ImC" + it->first).Data() , 0.0);
      }
      else {minimizer->SetVariable(parcount++,(name + "_ImC" + it->first).Data() , it->second.C().imag(), 0.10);}
      
      ++it;
    } // end loop over channels

  }
 
  const unsigned int nfree=minimizer->NFree();
  printInfo <<  nfree  << " Free Parameters in fit" << endl;


  // find minimum of likelihood function
  double chi2=0;
  if(onlyPlotting) printInfo << "Plotting mode, skipping minimzation!" << endl;
  else {
    printInfo << "performing minimization." << endl;
    
    minimizer->SetMaxIterations(maxNmbOfIterations);
    minimizer->SetMaxFunctionCalls(maxNmbOfIterations*5);
    minimizer->SetTolerance    (minimizerTolerance);
    bool success = minimizer->Minimize();
    const double* par=minimizer->X();
    compset.setPar(par);
    cerr << compset << endl;
    if (success){
      printInfo << "minimization finished successfully." << endl;
      chi2=minimizer->MinValue();
    }
    else
      printWarn << "minimization failed." << endl;
    if (runHesse) {
      printInfo << "calculating Hessian matrix." << endl;
      success = minimizer->Hesse();  // comes only with ROOT 5.24+
      if (!success)
	printWarn << "calculation of Hessian matrix failed." << endl;
    }
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
  

  // ---------------------------------------------------------------------------
  // print results
  //map<TString, double> errormap;
  printInfo << "minimization result:" << endl;
  for (unsigned int i = 0; i< nmbPar; ++i) {
    cout << "    parameter [" << setw(3) << i << "] ";
    cout << minimizer->VariableName(i) << " " ;
      //	 << setw(maxParNameLength); //<< L.parName(i) << " = ";
    //if (parIsFixed[i])
    //  cout << minimizer->X()[i] << " (fixed)" << endl;
    //else {
      cout << setw(12) << maxPrecisionAlign(minimizer->X()[i]) << " +- "
	   << setw(12) << maxPrecisionAlign(minimizer->Errors()[i]);
      //errormap[minimizer]=minimizer->Errors()[i];
      

      if (runMinos && (i == 156)) {  // does not work for all parameters
	double minosErrLow = 0;
	double minosErrUp  = 0;
	const bool success = minimizer->GetMinosError(i, minosErrLow, minosErrUp);
	if (success)
	  cout << "    Minos: " << "[" << minosErrLow << ", +" << minosErrUp << "]" << endl;
      } else
	cout << endl;
  }

 cout << "---------------------------------------------------------------------" << endl;
 // Reduced chi2
 
 printInfo << chi2 << " chi2" << endl;
 unsigned int numdata=L.NDataPoints();
 // numDOF
 unsigned int numDOF=numdata-nfree;
 printInfo << numDOF << " degrees of freedom" << endl;
 double redChi2 = chi2/(double)numDOF;
 printInfo << redChi2 << " chi2/nDF" << endl;
 cout << "---------------------------------------------------------------------" << endl;


  // write out results
  // Likelihood and such
 const Setting& fitqualS= root["fitquality"];
 Setting& chi2S=fitqualS["chi2"];
 chi2S=chi2;
 Setting& ndfS=fitqualS["ndf"];
 ndfS=(int)numDOF;
 Setting& redchi2S=fitqualS["redchi2"];
 redchi2S=redChi2;

  // Setup Component Set (Resonances + Background)
  const Setting& bws= root["components"]["resonances"];
  const Setting& bkgs= root["components"]["background"];
  unsigned int nbws=bws.getLength();
  unsigned int nbkgs=bkgs.getLength();
  // loop over components
  unsigned int nc=compset.n();
  for(unsigned int ic=0;ic<nc;++ic){
    const pwacomponent* comp=compset[ic];
    string name=comp->name();
    // search corresponding setting
   
    string sname;
    bool found=false;
    for(unsigned int is=0;is<nbws;++is){
      const Setting& bw = bws[is];
      bw.lookupValue("name",     sname);
      if(sname==name){
	found=true;
	// set values to this setting
	Setting& sm = bw["mass"];
	Setting& smval = sm["val"];
	smval = comp->m0();
	Setting& smerr = sm["error"];
	TString merrname=name+"_M";
	smerr=minimizer->Errors()[minimizer->VariableIndex(merrname.Data())];
	
	Setting& sw = bw["width"];
	Setting& swval = sw["val"];
	swval = comp->gamma();

	Setting& swerr = sw["error"];
	TString werrname=name+"_Gamma";
	swerr=minimizer->Errors()[minimizer->VariableIndex(werrname.Data())];

	cout << name 
	     << "   mass="<<double(smval)<<" +- "<<double(smerr)
	     << "   width="<<double(swval)<<" +- "<<double(swerr)<< endl;

	// loop through channel and fix couplings
	const Setting& sChn=bw["decaychannels"];
	unsigned int nCh=sChn.getLength();
	const std::map<std::string,pwachannel >& ch=comp->channels();
	std::map<std::string,pwachannel>::const_iterator it=ch.begin();
	for(;it!=ch.end();++it){
	  std::complex<double> c= it->second.C();
	  string ampname=it->first;
	  // loop through channels in setting
	  for(unsigned int isc=0;isc<nCh;++isc){
	    Setting& sCh=sChn[isc];
	    string amp; sCh.lookupValue("amp",amp);
	    if(amp==ampname){
	      Setting& sRe=sCh["coupling_Re"];
	      sRe=c.real();
	      Setting& sIm=sCh["coupling_Im"];
	       sIm=c.imag();
	      break;
	    } // endif
	  } // end loop through cannels in setting
	  
	} // end loop through channels of component

	break;
      }
    }

    // loop over background settings
    if(!found){
      for(unsigned int is=0;is<nbkgs;++is){
	const Setting& bw = bkgs[is];
	bw.lookupValue("name",     sname);
	if(sname==name){
	  Setting& sm = bw["m0"];
	  Setting& smval = sm["val"];
	  smval = comp->m0();
	  Setting& sw = bw["g"];
	  Setting& swval = sw["val"];
	  swval = comp->gamma();

	  const pwachannel& ch=comp->channels().begin()->second;
	  std::complex<double> c=ch.C();
	  Setting& sRe=bw["coupling_Re"];
	  sRe=c.real();
	  Setting& sIm=bw["coupling_Im"];
	  sIm=c.imag();
	  break;
	}
      }
    }
    

    

  }




 //  bool check=true;
//   // Resonances
//   if(Conf.exists("components.resonances")){
//     const Setting &bws = root["components"]["resonances"];
//     // loop through breitwigners
//     int nbw=bws.getLength();
//     printInfo << "found " << nbw << " Resonances in config" << endl;
//     for(int ibw = 0; ibw < nbw; ++ibw) {
//       const Setting &bw = bws[ibw];
//       string jpc;
//       string name;
//       double mass=-1;double ml,mu;int mfix; 
//       double width=-1;double wl,wu;int wfix;
     
//       check&=bw.lookupValue("name",     name);
//       check&=bw.lookupValue("jpc",       jpc);
//       const Setting &massSet = bw["mass"];
//       check&=massSet.lookupValue("val",        mass);
//       check&=massSet.lookupValue("lower",        ml);
//       check&=massSet.lookupValue("upper",        mu);
//       check&=massSet.lookupValue("fix",        mfix);
//       const Setting &widthSet = bw["width"];
//       check&=widthSet.lookupValue("val",       width);
//       check&=widthSet.lookupValue("lower",       wl);
//       check&=widthSet.lookupValue("upper",       wu);
//       check&=widthSet.lookupValue("fix",       wfix);
//       cout << "---------------------------------------------------------------------" << endl;
//       cout << name << "    JPC = " << jpc << endl;
//       cout << "mass(limits)  = " << mass <<" ("<<ml<<","<<mu<<") MeV/c^2";
//       if(mfix==1)cout<<"  -- FIXED";
//       cout<< endl;
//       cout << "width(limits) = " << width <<" ("<<wl<<","<<wu<<") MeV/c^2";
//       if(wfix==1)cout<<"  -- FIXED";
//       cout<< endl;
//       const Setting &channelSet = bw["decaychannels"];
//       unsigned int nCh=channelSet.getLength();
//       cout << "Decaychannels (coupling):" << endl;
//       std::map<std::string,pwachannel > channels;
//       for(unsigned int iCh=0;iCh<nCh;++iCh){
// 	const Setting &ch = channelSet[iCh];
// 	string amp;
// 	double cRe=0;
// 	double cIm=0;
// 	check&=ch.lookupValue("amp",amp);
// 	check&=ch.lookupValue("coupling_Re",cRe);
// 	check&=ch.lookupValue("coupling_Im",cIm);
// 	complex<double> C(cRe,cIm);
// 	cout << "   " << amp << "  " << C << endl;
// 	channels[amp]=pwachannel(C,getPhaseSpace(tree,amp));
//       }// end loop over channels
//       if(!check){
// 	printErr << "Bad config value lookup! Check your config file!" << endl;
// 	return 1;
//       }
//       pwacomponent* comp1=new pwacomponent(name,mass,width,channels);
//       comp1->setLimits(ml,mu,wl,wu);
//       comp1->setFixed(mfix,wfix);
//       compset.add(comp1);
//     }// end loop over resonances
//   }
//   cout << endl;
//   // Background components
//   if(Conf.exists("components.background")){
//     const Setting &bws = root["components"]["background"];
//     // loop through breitwigners
//     int nbw=bws.getLength();
//     printInfo << "found " << nbw << " Background components in config" << endl;
//     for(int ibw = 0; ibw < nbw; ++ibw) {
//       const Setting &bw = bws[ibw];
//       string name;
//       double mass=-1;double ml,mu;int mfix; 
//       double width=-1;double wl,wu;int wfix;
     
//       check&=bw.lookupValue("name",     name);
//       const Setting &massSet = bw["m0"];
//       check&=massSet.lookupValue("val",        mass);
//       check&=massSet.lookupValue("lower",        ml);
//       check&=massSet.lookupValue("upper",        mu);
//       check&=massSet.lookupValue("fix",        mfix);
//       const Setting &widthSet = bw["g"];
//       check&=widthSet.lookupValue("val",       width);
//       check&=widthSet.lookupValue("lower",       wl);
//       check&=widthSet.lookupValue("upper",       wu);
//       check&=widthSet.lookupValue("fix",       wfix);
//       cout << "---------------------------------------------------------------------" << endl;
//       cout << name << endl;
//       cout << "mass-offset(limits)  = " << mass <<" ("<<ml<<","<<mu<<") MeV/c^2";
//       if(mfix==1)cout<<"  -- FIXED";
//       cout<< endl;
//       cout << "g(limits)            = " << width <<" ("<<wl<<","<<wu<<") MeV/c^2";
//       if(wfix==1)cout<<"  -- FIXED";
//       cout<< endl;
//       std::map<std::string,pwachannel > channels;
//       string amp;
//       double cRe=0;
//       double cIm=0;
//       double mIso1=0;
//       double mIso2=0;
//       check&=bw.lookupValue("amp",amp);
//       check&=bw.lookupValue("coupling_Re",cRe);
//       check&=bw.lookupValue("coupling_Im",cIm);
//       check&=bw.lookupValue("mIsobar1",mIso1);
//       check&=bw.lookupValue("mIsobar2",mIso2);
//       complex<double> C(cRe,cIm);
//       cout << "Decaychannel (coupling):" << endl;
//       cout << "   " << amp << "  " << C << endl;
//       cout << "   Isobar masses: " << mIso1<<"  "<< mIso2<< endl;
//       channels[amp]=pwachannel(C,getPhaseSpace(tree,amp));

//       if(!check){
// 	printErr << "Bad config value lookup! Check your config file!" << endl;
// 	return 1;
//       }
//       pwabkg* bkg=new pwabkg(name,mass,width,channels);
//       bkg->setIsobars(mIso1,mIso2);
//       bkg->setLimits(ml,mu,wl,wu);
//       bkg->setFixed(mfix,wfix);
//       compset.add(bkg);
//     }// end loop over background
//   }// endif
  
  string outconfig(outFileName);
  outconfig.append(".conf");
  Conf.writeFile(outconfig.c_str());

  cerr << "Fitting finished... Start building graphs ... " << endl;

  int syscolor=kAzure-9;

   std::vector<std::string> wl=compset.wavelist();
   std::map<std::string, unsigned int> wmap;
   unsigned int ndatabins=tree->GetEntries();

   std::vector<TGraphErrors*> datagraphs;
   std::vector<TGraphErrors*> intenssysgraphs;

   std::vector<TMultiGraph*> graphs;

   for(unsigned int iw=0; iw<wl.size();++iw){
     wmap[wl[iw]]=iw;
     graphs.push_back(new TMultiGraph);

     intenssysgraphs.push_back(new TGraphErrors(ndatabins));
     string name=("sys_");name.append(wl[iw]);
     intenssysgraphs[iw]->SetName(name.c_str());
     intenssysgraphs[iw]->SetTitle(name.c_str());
     intenssysgraphs[iw]->SetLineColor(syscolor);
     intenssysgraphs[iw]->SetFillColor(syscolor);
     intenssysgraphs[iw]->SetDrawOption("2");
     graphs[iw]->Add(intenssysgraphs[iw],"2");



     graphs[iw]->SetName(wl[iw].c_str());
     graphs[iw]->SetTitle(wl[iw].c_str());
     graphs[iw]->SetDrawOption("AP");
     datagraphs.push_back(new TGraphErrors(ndatabins));
     name="data_";name.append(wl[iw]);
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
   std::vector<TGraph*> absphasegraphs;
  

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
     graphs[iw]->Add(getPhaseSpace(tree,fPS,wl[iw]));

     absphasegraphs.push_back(new TGraph(nbins));
     name="absphase_";name.append(wl[iw]);
     absphasegraphs[iw]->SetName(name.c_str());
     absphasegraphs[iw]->SetTitle(name.c_str());
     absphasegraphs[iw]->SetLineColor(kRed);
     absphasegraphs[iw]->SetLineWidth(2);
     absphasegraphs[iw]->SetMarkerColor(kRed);
     absphasegraphs[iw]->SetDrawOption("AP");



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
   std::vector<TGraphErrors*> phasesysgraphs;
   

   std::vector<TGraphErrors*> realdatagraphs;
   std::vector<TGraphErrors*> realsysgraphs;
   std::vector<TGraphErrors*> imagdatagraphs;
   std::vector<TGraphErrors*> imagsysgraphs;

   std::vector<TGraph*> realfitgraphs;
   std::vector<TGraph*> imagfitgraphs;
   
    std::vector<TMultiGraph*> phasegraphs;
   std::vector<TMultiGraph*> overlapRegraphs;
   std::vector<TMultiGraph*> overlapImgraphs;

   std::vector<TGraph*> phasefitgraphs;
   unsigned int c=0;
  

  for(unsigned int iw=0; iw<wl.size();++iw){
     for(unsigned int iw2=iw+1; iw2<wl.size();++iw2){

  

       phasegraphs.push_back(new TMultiGraph);
       overlapImgraphs.push_back(new TMultiGraph);
       overlapRegraphs.push_back(new TMultiGraph);
       string name("dPhi_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       phasegraphs[c]->SetName(name.c_str());
       phasegraphs[c]->SetTitle(name.c_str());
       name="Re_";name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       overlapRegraphs[c]->SetName(name.c_str());
       overlapRegraphs[c]->SetTitle(name.c_str());
       name="Im_";name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       overlapImgraphs[c]->SetName(name.c_str());
       overlapImgraphs[c]->SetTitle(name.c_str());

       phasesysgraphs.push_back(new TGraphErrors(nbins));
       name=("dPhi_sys_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       phasesysgraphs[c]->SetName(name.c_str());
       phasesysgraphs[c]->SetTitle(name.c_str());
       phasesysgraphs[c]->SetLineColor(syscolor);
       phasesysgraphs[c]->SetFillColor(syscolor);
       phasesysgraphs[c]->SetDrawOption("2");
       //phasesysgraphs[c]->SetLineWidth(1);
       //phasesysgraphs[c]->SetFillStyle(1001);
       
       phasegraphs[c]->Add(phasesysgraphs[c],"2");


       phasedatagraphs.push_back(new TGraphErrors(nbins));
       name=("dPhi_data_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       phasedatagraphs[c]->SetName(name.c_str());
       phasedatagraphs[c]->SetTitle(name.c_str());
       phasegraphs[c]->Add(phasedatagraphs[c],"p");

    
       realsysgraphs.push_back(new TGraphErrors(nbins));
       name=("RE_sys_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       realsysgraphs[c]->SetName(name.c_str());
       realsysgraphs[c]->SetTitle(name.c_str());
       realsysgraphs[c]->SetLineColor(syscolor);
       realsysgraphs[c]->SetFillColor(syscolor);
       realsysgraphs[c]->SetDrawOption("2");
       overlapRegraphs[c]->Add(realsysgraphs[c],"2");

       imagsysgraphs.push_back(new TGraphErrors(nbins));
       name=("IM_sys_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       imagsysgraphs[c]->SetName(name.c_str());
       imagsysgraphs[c]->SetTitle(name.c_str());
       imagsysgraphs[c]->SetLineColor(syscolor);
       imagsysgraphs[c]->SetFillColor(syscolor);
       imagsysgraphs[c]->SetDrawOption("2");
       overlapImgraphs[c]->Add(imagsysgraphs[c],"2");


       realdatagraphs.push_back(new TGraphErrors(nbins));
       name=("RE_data_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       realdatagraphs[c]->SetName(name.c_str());
       realdatagraphs[c]->SetTitle(name.c_str());
       overlapRegraphs[c]->Add(realdatagraphs[c],"p");

       imagdatagraphs.push_back(new TGraphErrors(nbins));
       name=("IM_data_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       imagdatagraphs[c]->SetName(name.c_str());
       imagdatagraphs[c]->SetTitle(name.c_str());
       //imagdatagraphs[c]->SetLineStyle(2);
       overlapImgraphs[c]->Add(imagdatagraphs[c],"p");

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
       phasefitgraphs[c]->SetDrawOption("P");
       phasefitgraphs[c]->SetMarkerStyle(24);
       phasefitgraphs[c]->SetMarkerSize(0.2);
       phasegraphs[c]->Add(phasefitgraphs[c],"cp");

       realfitgraphs.push_back(new TGraph(nbins));
       name=("Re_fit_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       realfitgraphs[c]->SetName(name.c_str());
       realfitgraphs[c]->SetTitle(name.c_str());
       realfitgraphs[c]->SetLineColor(kRed);
       realfitgraphs[c]->SetMarkerColor(kRed);
       realfitgraphs[c]->SetDrawOption("AP");
       realfitgraphs[c]->SetMarkerStyle(24);
       realfitgraphs[c]->SetMarkerSize(0.2);
       overlapRegraphs[c]->Add(realfitgraphs[c],"cp");

       imagfitgraphs.push_back(new TGraph(nbins));
       name=("Im_fit_");name.append(wl[iw]);
       name.append("---");name.append(wl[iw2]);
       imagfitgraphs[c]->SetName(name.c_str());
       imagfitgraphs[c]->SetTitle(name.c_str());
       imagfitgraphs[c]->SetLineColor(kRed);
       //imagfitgraphs[c]->SetLineStyle(2);
       imagfitgraphs[c]->SetMarkerColor(kRed);
       imagfitgraphs[c]->SetDrawOption("AP");
       imagfitgraphs[c]->SetMarkerStyle(24);
       imagfitgraphs[c]->SetMarkerSize(0.2);
       overlapImgraphs[c]->Add(imagfitgraphs[c],"cp");

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
   vector<double> prevphase(wl.size());
   double binwidth=30; // half binwidth
   //double w=2*30/10;
   for(unsigned int i=0;i<ndatabins;++i){
     tree->GetEntry(i);
     double m=rho->massBinCenter();
     //cout << "MASS: "<<m << endl;
     double fsps=sqrt(fPS->Eval(m));
     unsigned int c=0;
     for(unsigned int iw=0; iw<wl.size();++iw){

       //cout << wl[iw] << endl;

       double ps=rho->phaseSpaceIntegral(wl[iw].c_str())*fsps;
    
       datagraphs[iw]->SetPoint(i,m,rho->intensity(wl[iw].c_str()));
       datagraphs[iw]->SetPointError(i,binwidth,rho->intensityErr(wl[iw].c_str()));
      fitgraphs[iw]->SetPoint(i,m,compset.intensity(wl[iw],m));      
      double absphase=compset.phase(wl[iw],m)*TMath::RadToDeg();
      if(i>0){
	double absp=absphase+360;
	double absm=absphase-360;
	if(fabs(absphase-prevphase[iw])>fabs(absp-prevphase[iw])){
	  absphase=absp;
	}
	else if(fabs(absphase-prevphase[iw])>fabs(absm-prevphase[iw])){
	  absphase=absm;
	}
      }
      prevphase[iw]=absphase;
      absphasegraphs[iw]->SetPoint(i,m,absphase);      
      if(sysPlotting){
	double maxIntens=-10000000;
	double minIntens=10000000;
	for(unsigned int iSys=0;iSys<sysTrees.size();++iSys){
	  // get data
	  fitResult* rhoSys=0;
	  sysTrees[iSys]->SetBranchAddress(valBranchName.c_str(),&rhoSys);
	  sysTrees[iSys]->GetEntry(i);
	  // check if waves are in fit
	  if(rhoSys->waveIndex(wl[iw])==-1)continue;
	  double myI=rhoSys->intensity(wl[iw].c_str());
	  if(maxIntens<myI)maxIntens=myI;
	  if(minIntens>myI)minIntens=myI;
	  delete rhoSys;
	} // end loop over systematic trees
	
	intenssysgraphs[iw]->SetPoint(i,m,(maxIntens+minIntens)*0.5);
	intenssysgraphs[iw]->SetPointError(i,binwidth,(maxIntens-minIntens)*0.5);
      }
      
      //cout << "getting phases" << endl;

       // second loop to get phase differences
       unsigned int wi1=rho->waveIndex(wl[iw].c_str());
       
       for(unsigned int iw2=iw+1; iw2<wl.size();++iw2){
	 //double ps2=rho->phaseSpaceIntegral(wl[iw2].c_str())*fsps;
	  
	 unsigned int wi2=rho->waveIndex(wl[iw2].c_str());
	 complex<double> r=rho->spinDensityMatrixElem(wi1,wi2);
	 TMatrixT<double> rCov=rho->spinDensityMatrixElemCov(wi1,wi2);

	 realdatagraphs[c]->SetPoint(i,
				     rho->massBinCenter(),
				     r.real());
	 realdatagraphs[c]->SetPointError(i,
			    binwidth,
			    sqrt(rCov[0][0]));
	 imagdatagraphs[c]->SetPoint(i,
		       rho->massBinCenter(),
		       r.imag());
     
	 imagdatagraphs[c]->SetPointError(i,
			    binwidth,
			    sqrt(rCov[1][1]));


	 double dataphi=rho->phase(wl[iw].c_str(),wl[iw2].c_str());

	 phasedatagraphs[c]->SetPoint(i,m,dataphi);
	   
	 TVector2 v;v.SetMagPhi(1,rho->phase(wl[iw].c_str(),
					     wl[iw2].c_str())/TMath::RadToDeg());
	 phase2d[c]->SetPoint(i,v.X(),v.Y(),m);

	 phasedatagraphs[c]->SetPointError(i,binwidth,
					   rho->phaseErr(wl[iw].c_str(),
							 wl[iw2].c_str()));
	 double fitphase=compset.phase(wl[iw],wl[iw2],m)*TMath::RadToDeg();

	 if(sysPlotting){
	   //cout << "start sysplotting" << endl;
	   // loop over systematics files
	   double maxPhase=-10000;
	   double minPhase=10000;
	   double maxRe=-10000000;
	   double minRe=10000000;
	   double maxIm=-10000000;
	   double minIm=10000000;
	  
	   for(unsigned int iSys=0;iSys<sysTrees.size();++iSys){
	     // cerr << iSys;
	   // get data
	     fitResult* rhoSys=0;
	     sysTrees[iSys]->SetBranchAddress(valBranchName.c_str(),&rhoSys);
	     sysTrees[iSys]->GetEntry(i);
	        

	     // check if waves are in fit
	     if(rhoSys==NULL || rhoSys->waveIndex(wl[iw])==-1 || rhoSys->waveIndex(wl[iw2]) ==-1)continue;
	     // get correct wave indices!!!
	     unsigned int wi1Sys=rhoSys->waveIndex(wl[iw].c_str());
	     unsigned int wi2Sys=rhoSys->waveIndex(wl[iw2].c_str());

	     double myphi=rhoSys->phase(wl[iw].c_str(),wl[iw2].c_str());
	     double myphiplus=myphi+360;
	     double myphiminus=myphi-360;
	     // translate by 2pi to get closest solution to fit
	     if(fabs(myphiplus-dataphi)<fabs(myphi-dataphi)){
		 myphi=myphiplus;
		 //cout << "myphiminus" << endl;
	     }
	     if(fabs(myphiminus-dataphi)<fabs(myphi-dataphi)){
	       myphi=myphiminus;
	       //cout << "myphiplus" << endl;
	     }

	     if(myphi>maxPhase)maxPhase=myphi;
	     if(myphi<minPhase)minPhase=myphi;

	     // real and imag part:
	     complex<double> r=rhoSys->spinDensityMatrixElem(wi1Sys,wi2Sys);
	     if(maxRe<r.real())maxRe=r.real();
	     if(minRe>r.real())minRe=r.real();
	     if(maxIm<r.imag())maxIm=r.imag();
	     if(minIm>r.imag())minIm=r.imag();
	     delete rhoSys; rhoSys=0;
	   }// end loop over sys trees
	   //cerr << "loop over systrees finished" << endl;

	   phasesysgraphs[c]->SetPoint(i,m,(maxPhase+minPhase)*0.5);
	   phasesysgraphs[c]->SetPointError(i,binwidth,(maxPhase-minPhase)*0.5);
	   
	   realsysgraphs[c]->SetPoint(i,m,(maxRe+minRe)*0.5);
	   realsysgraphs[c]->SetPointError(i,binwidth,(maxRe-minRe)*0.5);
	   imagsysgraphs[c]->SetPoint(i,m,(maxIm+minIm)*0.5);
	   imagsysgraphs[c]->SetPointError(i,binwidth,(maxIm-minIm)*0.5);
	

	 }// end if sysplotting
	 //cout << "sysplotting finished" << endl;

	 phasefitgraphs[c]->SetPoint(i,m,fitphase);

	 complex<double> fitval=compset.overlap(wl[iw],wl[iw2],m);

	 realfitgraphs[c]->SetPoint(i,m,fitval.real());
	 imagfitgraphs[c]->SetPoint(i,m,fitval.imag());


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
	   double I=norm(c->val(m)*it->second.C())*it->second.ps(m);
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
     TGraph* g=(TGraph*)graphs[iw]->GetListOfGraphs()->At(3);
      for(unsigned int ib=0;ib<nbins;++ib){
        double m,ps;g->GetPoint(ib,m,ps);
        g->SetPoint(ib,m,ps*500.);
       }
     
     graphs[iw]->Write();
     //absphasegraphs[iw]->Write();
   }


 for(unsigned int iw=0; iw<phasegraphs.size();++iw){

/// rectivfy phase graphs
   
   unsigned int refbin=6;
   double m;
   double predph;
   phasedatagraphs[iw]->GetPoint(refbin,m,predph);
   double prefph;
   phasefitgraphs[iw]->GetPoint(refbin,m,prefph);
   for(unsigned int ib=refbin+1;ib<nbins;++ib){
     double dph; phasedatagraphs[iw]->GetPoint(ib,m,dph);
     double fph; phasefitgraphs[iw]->GetPoint(ib,m,fph);
     double dp,dm;dp=dph+360;dm=dph-360;
     double fp,fm;fp=fph+360;fm=fph-360;
     if(1){
       if(fabs(fp-prefph)<fabs(fph-prefph) && fabs(fp-prefph)<fabs(fm-prefph))
   	 phasefitgraphs[iw]->SetPoint(ib,m,fp);
       else if(fabs(fm-prefph)<fabs(fph-prefph) && fabs(fm-prefph)<fabs(fp-prefph))
   	 phasefitgraphs[iw]->SetPoint(ib,m,fm);
       phasefitgraphs[iw]->GetPoint(ib,m,prefph);

       if(fabs(dp-prefph)<fabs(dph-prefph) && fabs(dp-prefph)<fabs(dm-prefph))
   	 phasedatagraphs[iw]->SetPoint(ib,m,dp);
       else if(fabs(dm-prefph)<fabs(dph-prefph) && fabs(dm-prefph)<fabs(dp-prefph))
   	 phasedatagraphs[iw]->SetPoint(ib,m,dm);
       
     
       
       phasedatagraphs[iw]->GetPoint(ib,m,predph);
     


   // put systematic error closest to fit/data
       double sph; phasesysgraphs[iw]->GetPoint(ib,m,sph);
       double sp,sm;sp=sph+360;sm=sph-360;
       if(fabs(sp-prefph)<fabs(sph-prefph) && fabs(sp-prefph)<fabs(sm-prefph))
	 phasesysgraphs[iw]->SetPoint(ib,m,sp);
       else if(fabs(sm-prefph)<fabs(sph-prefph) && fabs(sm-prefph)<fabs(sp-prefph))
	 phasesysgraphs[iw]->SetPoint(ib,m,sm);


     }
   }
   // backward:
   phasedatagraphs[iw]->GetPoint(refbin,m,predph);
   phasefitgraphs[iw]->GetPoint(refbin,m,prefph);
   for(unsigned int i=0;i<refbin;++i){
       unsigned int ib=refbin-i-1;
       double dph; phasedatagraphs[iw]->GetPoint(ib,m,dph);
     double fph; phasefitgraphs[iw]->GetPoint(ib,m,fph);
     double dp,dm;dp=dph+360;dm=dph-360;
     double fp,fm;fp=fph+360;fm=fph-360;
     if(1){

       if(fabs(fp-prefph)<fabs(fph-prefph) && fabs(fp-prefph)<fabs(fm-prefph)){
   	 phasefitgraphs[iw]->SetPoint(ib,m,fp);
	 
       }
       else if(fabs(fm-prefph)<fabs(fph-prefph) && fabs(fm-prefph)<fabs(fp-prefph)){
   	 phasefitgraphs[iw]->SetPoint(ib,m,fm);
	 
       }
       
       phasefitgraphs[iw]->GetPoint(ib,m,prefph);
       
       
       
       if(fabs(dp-prefph)<fabs(dph-prefph) && fabs(dp-prefph)<fabs(dm-prefph)){
   	 phasedatagraphs[iw]->SetPoint(ib,m,dp);
       }
       
       else if(fabs(dm-prefph)<fabs(dph-prefph) && fabs(dm-prefph)<fabs(dp-prefph)){
   	 phasedatagraphs[iw]->SetPoint(ib,m,dm);
       }
       
       phasedatagraphs[iw]->GetPoint(ib,m,predph);
       
       
       // put systematic error closest to fit
       double sph; phasesysgraphs[iw]->GetPoint(ib,m,sph);
       double sp,sm;sp=sph+360;sm=sph-360;
       if(fabs(sp-predph)<fabs(sph-predph) && fabs(sp-predph)<fabs(sm-predph))
	 phasesysgraphs[iw]->SetPoint(ib,m,sp);
       else if(fabs(sm-predph)<fabs(sph-predph) && fabs(sm-predph)<fabs(sp-predph))
	 phasesysgraphs[iw]->SetPoint(ib,m,sm);

      
     }
   }

     
     phasegraphs[iw]->Write();
     phasesysgraphs[iw]->Write();
     overlapRegraphs[iw]->Write();
     overlapImgraphs[iw]->Write();
     //phase2d[iw]->Write();
 } // end loop over waves

outfile->Close();
   
   return 0;
   
}
