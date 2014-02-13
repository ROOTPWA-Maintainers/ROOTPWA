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

#include <boost/assign/std/vector.hpp>

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
#include "TStopwatch.h"
#include "libconfig.h++"
#include <libConfigUtils.hpp>

#define MASSSCALE 0.001


using namespace std;
using namespace libconfig;
using namespace ROOT::Math;
using namespace rpwa;
using namespace boost;
using namespace boost::assign;


void
usage(const string& progName,
      const int     errCode = 0)
{
  cerr << "usage:" << endl
       << progName
       << " -c configfile -i inputfile [-o outfile"
       << "  -M minimizer [-m algorithm] -t # -q -h -P -R -C] [-S fitResultFiles]" << endl
       << "    where:" << endl
       << "        -c file    path to config File" << endl
       << "        -i file    path to input file" << endl
       << "        -o file    path to output file (default: 'mDep.result.root')" << endl
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
       << "        -P         plotting only - no fit" << endl
       << "        -R         plot in fit range only" << endl
       << "        -C         switch OFF covariances between real and imag part" << endl
       << "        -S files   Systematic error plotting. give list of files" << endl

       << endl;
  exit(errCode);
}

//  function loops through fitResults and puts phasespace values into a graph for interpolation
//  THIS CONTAINS NOW THE RIGHT VALUE NOT!!! THE SQRT!!!
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
    ps*=ps; // remember that phaseSpaceIntegral returns sqrt of integral!!! 
    graph->SetPoint(i,m,ps);
  }
  return graph;
}

// changes status of variables (fixed/released)
// fixed values from config remain fixed
// parameters are taken from current status of fitter
// level
// 0 = release only couplings
// 1 = release couplings and masses
// 2 = release couplings, masses and widths
void releasePars(Minimizer* minimizer, const pwacompset& compset,
		 const vector<string>& anchorwave_reso,
		 const vector<string>& anchorwave_channel,
		 int level){
  // copy state
  unsigned int npar=minimizer->NDim();
  double par[npar];
  for(unsigned int i=0;i<npar;++i)par[i]=minimizer->X()[i];
  minimizer->Clear();

  unsigned int parcount=0;
  for(unsigned int ic=0;ic<compset.n();++ic){
    const pwacomponent& comp=*compset[ic];
    TString name(comp.name());
    double mmin,mmax,gmin,gmax;
    comp.getLimits(mmin,mmax,gmin,gmax);
    if(comp.fixM() || level==0)minimizer->SetFixedVariable(parcount,
					       (name+"_M").Data() ,
					       par[parcount]);
    else minimizer->SetLimitedVariable(parcount,
				       (name+"_M").Data(),
				       par[parcount],
				       5.0,
				       mmin,mmax);
    if(level==0 && !comp.fixM()) printInfo << minimizer->VariableName(parcount)
			   << " fixed to " << par[parcount] << endl;
    ++parcount;
    if(comp.fixGamma() || level < 2)minimizer->SetFixedVariable(parcount,
						   (name+"_Gamma").Data() ,
						    par[parcount]);
    else minimizer->SetLimitedVariable(parcount,
				       (name+"_Gamma").Data(),
				       par[parcount],
				       5.0,
				       gmin,gmax);
    if(level<2 && !comp.fixGamma()) printInfo << minimizer->VariableName(parcount)
			  << " fixed to " << par[parcount] << endl;
    ++parcount;

    std::map<std::string,pwachannel >::const_iterator it=comp.channels().begin();
    while(it!=comp.channels().end()){
      minimizer->SetVariable(parcount,(name + "_ReC" + it->first).Data() , par[parcount], 10.0);
      ++parcount;
      // fix one phase
      if(find(anchorwave_reso.begin(),anchorwave_reso.end(),name)!=anchorwave_reso.end() && find(anchorwave_channel.begin(),anchorwave_channel.end(),it->first)!=anchorwave_channel.end()){
	minimizer->SetFixedVariable(parcount,(name + "_ImC" + it->first).Data() , 0.0);
      }
      else {minimizer->SetVariable(parcount,(name + "_ImC" + it->first).Data() , par[parcount], 0.10);}
      ++parcount;
      ++it;
    } // end loop over channels
  }// end loop over components
  // set phase space
  unsigned int nfreeFsmd=compset.nFreeFsmdPar();
  for(unsigned int ifreeFsmd=0;ifreeFsmd<nfreeFsmd;++ifreeFsmd){
    double val,lower,upper;
    val=par[parcount];
    compset.getFreeFsmdLimits(ifreeFsmd,lower,upper);
    TString name("PSP_"); name+=+ifreeFsmd;
    minimizer->SetLimitedVariable(parcount, 
				  name.Data(), 
				  val, 0.0001 ,lower,upper);
  }



  const unsigned int nfree=minimizer->NFree();
  printInfo <<  nfree  << " Free Parameters in fit" << endl;


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
  bool               rangePlotting        = false;
  bool               sysPlotting         = false;
  bool               doCov               = true;

  //unsigned int maxParNameLength = 20;       // maximum length of parameter names
//   int                startValSeed        = 1234567;
 // parse command line options
  const string progName           = argv[0];
  
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
  while ((ca = getopt(argc, argv, "c:i:o:M:m:t:qhPRCS")) != -1)
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
    case 'q':
      quiet = true;
      break;
    case 'h':
      usage(progName);
      break;
    case 'P':
      onlyPlotting=true;
      break;
    case 'R':
      rangePlotting=true;
      break;
    case 'C':
      doCov=false;
      break;
    case 'S':
      sysPlotting=true;
      break;
    }


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


  vector<TTree*> sysTrees;
  if(sysPlotting){
    // add this fit
    sysTrees.push_back(tree);
    // open files with fits
    for(int i=optind;i<argc;++i){
      // open input file and get results tree
      TFile* infile=TFile::Open(argv[i]);
      if(infile==NULL){
	cerr << "Systematics Input file " << inFileName <<" not found."<< endl;
	return 1;
      }
      TTree* systree=(TTree*)infile->Get(valTreeName.c_str());
      if(systree==NULL){
	cerr << "Input tree " << valTreeName <<" not found."<< endl;
	return 1;
      }
      sysTrees.push_back(systree);
    }
    printInfo << sysTrees.size() << " files for systematics found " << endl;
  }// end if sysPlotting

  printInfo << "creating and setting up likelihood function" << endl;
  printInfo << "doCovariances = " << doCov << endl;



  // Setup Component Set (Resonances + Background)
  pwacompset compset;
  Config Conf;
  Conf.readFile(configFile.c_str());
  const Setting& root = Conf.getRoot();
  bool check=true;

	// input section
	const Setting* configInput = findLibConfigGroup(root, "input");
	if(not configInput) {
		printErr << "'input' section in configuration file does not exist." << endl;
		exit(1);
	}

	// get information of waves to be used in the fit
	const Setting* configInputWaves = findLibConfigList(*configInput, "waves");
	if(not configInputWaves) {
		printErr << "'waves' list does not exist in section 'input' in configuration file." << endl;
		exit(1);
	}

	vector<string> waveNames;
	vector<pair<double, double> > waveMassLimits;
	{
		const int nrWaves = configInputWaves->getLength();
		printInfo << "going to read information of " << nrWaves << " waves to be used in the fit." << endl;

		for(int idxWave=0; idxWave<nrWaves; ++idxWave) {
			const Setting* configInputWave = &((*configInputWaves)[idxWave]);

			map<string, Setting::Type> mandatoryArguments;
			insert(mandatoryArguments)
			    ("name", Setting::TypeString);
			if(not checkIfAllVariablesAreThere(configInputWave, mandatoryArguments)) {
				printErr << "'waves' list in 'input' section in configuration file contains errors." << endl;
				exit(1);
			}

			string name;
			configInputWave->lookupValue("name", name);

			double massLower;
			if(not configInputWave->lookupValue("massLower", massLower)) {
				massLower = -1.;
			}
			double massUpper;
			if(not configInputWave->lookupValue("massUpper", massUpper)) {
				massUpper = -1.;
			}

			// check that wave does not yet exist
			if(find(waveNames.begin(), waveNames.end(), name) != waveNames.end()) {
				printErr << "wave '" << name << "' defined twice." << endl;
				exit(1);
			}


			waveNames.push_back(name);
			waveMassLimits.push_back(make_pair(massLower, massUpper));

			printInfo << idxWave << ": " << name << " (mass range: " << massLower << "-" << massUpper << " MeV/c^2)" << endl;
		}

		if(waveNames.size()!= nrWaves || waveMassLimits.size()!= nrWaves) {
			printErr << "detected a messed up internal array." << endl;
			exit(1);
		}

		compset.setWaveList(waveNames);
	}


  // overall final-state mass dependence
  TF1* fPS = NULL;
  const Setting* configFsmd = findLibConfigGroup(root, "finalStateMassDependence", false);
  if(configFsmd){
    map<string, Setting::Type> mandatoryArguments;
    insert(mandatoryArguments)
        ("formula", Setting::TypeString)
        ("val", Setting::TypeArray)
        ("lower", Setting::TypeArray)
        ("upper", Setting::TypeArray)
        ("error", Setting::TypeArray)
        ("fix", Setting::TypeArray);
    if(not checkIfAllVariablesAreThere(configFsmd, mandatoryArguments)) {
      printErr << "'finalStateMassDependence' section in configuration file contains errors." << endl;
      exit(1);
    }
    
    string formula;
    configFsmd->lookupValue("formula", formula);

    fPS = new TF1("fps", formula.c_str(), 900, 3000);
    const unsigned int nrPar = fPS->GetNpar();

    const Setting& configFsmdValue = (*configFsmd)["val"];
    if(configFsmdValue.getLength() != nrPar) {
      printErr << "'val' in 'finalStateMassDependence' has to have a length of " << nrPar << "." << endl;
      return false;
    }
    if(configFsmdValue.getLength() > 0 and not configFsmdValue[0].isNumber()) {
      printErr << "'val' in 'finalStateMassDependence' has to be made up of numbers." << endl;
      return false;
    }

    const Setting& configFsmdLower = (*configFsmd)["lower"];
    if(configFsmdLower.getLength() != nrPar) {
      printErr << "'lower' in 'finalStateMassDependence' has to have a length of " << nrPar << "." << endl;
      return false;
    }
    if(configFsmdLower.getLength() > 0 and not configFsmdLower[0].isNumber()) {
      printErr << "'lower' in 'finalStateMassDependence' has to be made up of numbers." << endl;
      return false;
    }

    const Setting& configFsmdUpper = (*configFsmd)["upper"];
    if(configFsmdUpper.getLength() != nrPar) {
      printErr << "'upper' in 'finalStateMassDependence' has to have a length of " << nrPar << "." << endl;
      return false;
    }
    if(configFsmdUpper.getLength() > 0 and not configFsmdUpper[0].isNumber()) {
      printErr << "'upper' in 'finalStateMassDependence' has to be made up of numbers." << endl;
      return false;
    }

    const Setting& configFsmdError = (*configFsmd)["error"];
    if(configFsmdError.getLength() != nrPar) {
      printErr << "'error' in 'finalStateMassDependence' has to have a length of " << nrPar << "." << endl;
      return false;
    }
    if(configFsmdError.getLength() > 0 and not configFsmdError[0].isNumber()) {
      printErr << "'error' in 'finalStateMassDependence' has to be made up of numbers." << endl;
      return false;
    }

    const Setting& configFsmdFix = (*configFsmd)["fix"];
    if(configFsmdFix.getLength() != nrPar) {
      printErr << "'fix' in 'finalStateMassDependence' has to have a length of " << nrPar << "." << endl;
      return false;
    }
    if(configFsmdFix.getLength() > 0 and not configFsmdFix[0].isNumber()) {
      printErr << "'fix' in 'finalStateMassDependence' has to be made up of numbers." << endl;
      return false;
    }

    for (unsigned int i=0; i<nrPar; ++i) {
      fPS->SetParameter(i, configFsmdValue[i]);
      fPS->SetParError(i, configFsmdError[i]);
      fPS->SetParLimits(i, configFsmdLower[i], configFsmdUpper[i]);

      if (((int)configFsmdFix[i]) != 0) {
        fPS->FixParameter(i, configFsmdValue[i]);
      }
    }
    printInfo << "using final-state mass dependence as defined in the configuration file." << endl;
  } else {
    printInfo << "not using final-state mass dependence." << endl;
  }
  compset.setFuncFsmd(fPS);


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
	channels[amp]=pwachannel(C,getPhaseSpace(tree,amp));
      }// end loop over channels
      if(!check){
	printErr << "Bad config value lookup! Check your config file!" << endl;
	return 1;
      }
      pwacomponent* comp1=new pwacomponent(name,mass,width,channels);
      cerr << "created component" << endl;
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
      channels[amp]=pwachannel(C,getPhaseSpace(tree,amp));

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


 cout << "---------------------------------------------------------------------" << endl << endl;

	if(not compset.doMapping()) {
		printErr << "error while mapping the waves to the decay channels and components." << endl;
		exit(1);
	}

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
	if(not L.init(tree,&compset,waveNames,waveMassLimits,doCov)) {
		printErr << "error while initializing likelihood function." << endl;
		exit(1);
	}

   const unsigned int nmbPar  = L.NDim();
  // double par[nmbPar];
  // for(unsigned int ip=0;ip<nmbPar;++ip)par[ip]=1.4;


  // TStopwatch watch;
  // L.DoEval(par);
  // watch.Stop();


  //printInfo << "TESTCALL TO LIKELIHOOD takes " <<  maxPrecisionAlign(watch.CpuTime()) << " s" << endl;

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
  }// end loop over components
  // set phase space
  unsigned int nfreeFsmd=compset.nFreeFsmdPar();
  for(unsigned int ifreeFsmd=0;ifreeFsmd<nfreeFsmd;++ifreeFsmd){
    double val,lower,upper;
    val=compset.getFreeFsmdPar(ifreeFsmd);
    compset.getFreeFsmdLimits(ifreeFsmd,lower,upper);
    TString name("PSP_"); name+=+ifreeFsmd;
    minimizer->SetLimitedVariable(parcount++, 
				  name.Data(), 
				  val, 0.0001 ,lower,upper);
  }



  const unsigned int nfree=minimizer->NFree();
  printInfo <<  nfree  << " Free Parameters in fit" << endl;


  // find minimum of likelihood function
  double chi2=0;
  if(onlyPlotting) printInfo << "Plotting mode, skipping minimzation!" << endl;
  else {
    printInfo << "performing minimization. MASSES AND WIDTHS FIXED" << endl;

    minimizer->SetMaxIterations(maxNmbOfIterations);
    minimizer->SetMaxFunctionCalls(maxNmbOfIterations*5);
    minimizer->SetTolerance    (minimizerTolerance);
    // only do couplings
    TStopwatch fitW;
    // releasePars(minimizer,compset,anchorwave_reso,anchorwave_channel,0);
    bool success = minimizer->Minimize();
    if(!success)printWarn << "minimization failed." << endl;
    else printInfo << "minimization successful." << endl;
    printInfo << "Minimization took " <<  maxPrecisionAlign(fitW.CpuTime()) << " s" << endl;
    //release masses
    releasePars(minimizer,compset,anchorwave_reso,anchorwave_channel,1);
    printInfo << "performing minimization. MASSES RELEASED" << endl;
    fitW.Start();
    success &= minimizer->Minimize();
    if(!success)printWarn << "minimization failed." << endl;
    else printInfo << "minimization successful." << endl;
    printInfo << "Minimization took " <<  maxPrecisionAlign(fitW.CpuTime()) << " s" << endl;
    //release widths
    releasePars(minimizer,compset,anchorwave_reso,anchorwave_channel,2);
    printInfo << "performing minimization. ALL RELEASED" << endl;
    fitW.Start();
    success &= minimizer->Minimize();
    printInfo << "Minimization took " <<  maxPrecisionAlign(fitW.CpuTime()) << " s" << endl;

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

  if(configFsmd){
    const Setting& configFsmdValue = (*configFsmd)["val"];
    const Setting& configFsmdError = (*configFsmd)["error"];
    const Setting& configFsmdFix = (*configFsmd)["fix"];

    const unsigned int nrPar = fPS->GetNpar();
    unsigned int iPar = 0;
    for(unsigned int i=0; i<nrPar; ++i) {
      if(((int)configFsmdFix[i]) == 0) {
        ostringstream sName;
        sName << "PSP_" << iPar;
        configFsmdValue[i] = compset.getFreeFsmdPar(iPar);
        configFsmdError[i] = minimizer->Errors()[minimizer->VariableIndex(sName.str())];
        ++iPar;
      }
    }
    assert(compset.nFreeFsmdPar() == iPar);
  }
 
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


  
  string outconfig(outFileName);
  outconfig.append(".conf");
  Conf.writeFile(outconfig.c_str());

  cerr << "Fitting finished... Start building graphs ... " << endl;

  int syscolor=kAzure-9;

   std::vector<std::string> wl=compset.getWaveList();
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
     graphs[iw]->Add(getPhaseSpace(tree,wl[iw]));

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
   vector<double> prevphase(wl.size());
   double binwidth=MASSSCALE*30; // half binwidth
   //double w=2*30/10;
   for(unsigned int i=0;i<ndatabins;++i){
     tree->GetEntry(i);
     double m=rho->massBinCenter();
     double mScaled=m*MASSSCALE;
     //cout << "MASS: "<<m << endl;
     unsigned int c=0;
     for(unsigned int iw=0; iw<wl.size();++iw){
		// check that this mass bin should be taken into account for this
		// combination of waves
		if(rangePlotting && (i < L.getMassBinLimits()[iw][iw].first || i > L.getMassBinLimits()[iw][iw].second)) {
			continue;
		}

       datagraphs[iw]->SetPoint(i,mScaled,rho->intensity(wl[iw].c_str()));
       datagraphs[iw]->SetPointError(i,binwidth,rho->intensityErr(wl[iw].c_str()));
      fitgraphs[iw]->SetPoint(i,mScaled,compset.intensity(wl[iw],m));
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
      absphasegraphs[iw]->SetPoint(i,mScaled,absphase);
      if(sysPlotting){
	double maxIntens=-10000000;
	double minIntens=10000000;
	for(unsigned int iSys=0;iSys<sysTrees.size();++iSys){
	  // get data
	  fitResult* rhoSys=0;
	  if(iSys==0)rhoSys=rho;
	  else {
	    sysTrees[iSys]->SetBranchAddress(valBranchName.c_str(),&rhoSys);
	    sysTrees[iSys]->GetEntry(i);
	  }
// rename one wave
	  string mywave1=wl[iw];


	  if(mywave1=="1-1++0+pi-_01_rho1700=pi-+_10_pi1300=pi+-_00_sigma.amp" && iSys>0)mywave1="1-1++0+pi-_01_eta11600=pi-+_10_pi1300=pi+-_00_sigma.amp";

	  // check if waves are in fit
	  if(rhoSys->waveIndex(mywave1)==-1){
	    delete rhoSys;
	    continue;
	  }
	  double myI=rhoSys->intensity(mywave1.c_str());
	  if(maxIntens<myI)maxIntens=myI;
	  if(minIntens>myI)minIntens=myI;
	  if(iSys>0)delete rhoSys;
	} // end loop over systematic trees

	intenssysgraphs[iw]->SetPoint(i,mScaled,(maxIntens+minIntens)*0.5);
	intenssysgraphs[iw]->SetPointError(i,binwidth,(maxIntens-minIntens)*0.5);
      }

      //cout << "getting phases" << endl;

       // second loop to get phase differences
       unsigned int wi1=rho->waveIndex(wl[iw].c_str());

       for(unsigned int iw2=iw+1; iw2<wl.size();++iw2){
			// check that this mass bin should be taken into account for this
			// combination of waves
			if(rangePlotting && (i < L.getMassBinLimits()[iw][iw2].first || i > L.getMassBinLimits()[iw][iw2].second)) {
				continue;
			}

	 //double ps2=rho->phaseSpaceIntegral(wl[iw2].c_str())*fsps;

	 unsigned int wi2=rho->waveIndex(wl[iw2].c_str());
	 complex<double> r=rho->spinDensityMatrixElem(wi1,wi2);
	 TMatrixT<double> rCov=rho->spinDensityMatrixElemCov(wi1,wi2);

	 realdatagraphs[c]->SetPoint(i,
				     mScaled,
				     r.real());
	 realdatagraphs[c]->SetPointError(i,
			    binwidth,
			    sqrt(rCov[0][0]));
	 imagdatagraphs[c]->SetPoint(i,
		       mScaled,
		       r.imag());

	 imagdatagraphs[c]->SetPointError(i,
			    binwidth,
			    sqrt(rCov[1][1]));


	 double dataphi=rho->phase(wl[iw].c_str(),wl[iw2].c_str());

	 phasedatagraphs[c]->SetPoint(i,mScaled,dataphi);

	 TVector2 v;v.SetMagPhi(1,rho->phase(wl[iw].c_str(),
					     wl[iw2].c_str())/TMath::RadToDeg());
	 phase2d[c]->SetPoint(i,v.X(),v.Y(),mScaled);

	 phasedatagraphs[c]->SetPointError(i,binwidth,
					   rho->phaseErr(wl[iw].c_str(),
							 wl[iw2].c_str()));
	 double fitphase=compset.phase(wl[iw],wl[iw2],m)*TMath::RadToDeg();

	 if(sysPlotting){
	   //cerr << "start sysplotting" << endl;
	   // loop over systematics files
	   double maxPhase=-10000;
	   double minPhase=10000;
	   double maxRe=-10000000;
	   double minRe=10000000;
	   double maxIm=-10000000;
	   double minIm=10000000;

	   for(unsigned int iSys=0;iSys<sysTrees.size();++iSys){
	     //cerr << iSys;
	   // get data
	     fitResult* rhoSys=0;
	     if(iSys==0)rhoSys=rho;
	     else {
	       sysTrees[iSys]->SetBranchAddress(valBranchName.c_str(),&rhoSys);
	       if(i<sysTrees[iSys]->GetEntries()) sysTrees[iSys]->GetEntry(i);
	       else{
		 delete rhoSys;
		 continue;
	       }
	     }
	     // rename one wave
	     string mywave1=wl[iw];
	     string mywave2=wl[iw2];

if(mywave1=="1-1++0+pi-_01_rho1700=pi-+_10_pi1300=pi+-_00_sigma.amp" && iSys>0)mywave1="1-1++0+pi-_01_eta11600=pi-+_10_pi1300=pi+-_00_sigma.amp";
if(mywave2=="1-1++0+pi-_01_rho1700=pi-+_10_pi1300=pi+-_00_sigma.amp" && iSys>0)mywave2="1-1++0+pi-_01_eta11600=pi-+_10_pi1300=pi+-_00_sigma.amp";

	     // check if waves are in fit
 if(rhoSys==NULL || rhoSys->waveIndex(mywave1.c_str())==-1 || rhoSys->waveIndex(mywave2.c_str()) ==-1){
	       delete rhoSys;
	       continue;
	     }
	     // get correct wave indices!!!
	     unsigned int wi1Sys=rhoSys->waveIndex(mywave1.c_str());
	     unsigned int wi2Sys=rhoSys->waveIndex(mywave2.c_str());

	     double myphi=rhoSys->phase(mywave1.c_str(),mywave2.c_str());
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
	     if(iSys>0)delete rhoSys;
	     rhoSys=0;
	   }// end loop over sys trees
	   // cerr << "loop over systrees finished" << endl;

	   phasesysgraphs[c]->SetPoint(i,mScaled,(maxPhase+minPhase)*0.5);
	   phasesysgraphs[c]->SetPointError(i,binwidth,(maxPhase-minPhase)*0.5);

	   realsysgraphs[c]->SetPoint(i,mScaled,(maxRe+minRe)*0.5);
	   realsysgraphs[c]->SetPointError(i,binwidth,(maxRe-minRe)*0.5);
	   imagsysgraphs[c]->SetPoint(i,mScaled,(maxIm+minIm)*0.5);
	   imagsysgraphs[c]->SetPointError(i,binwidth,(maxIm-minIm)*0.5);


	 }// end if sysplotting
	 //cerr << "sysplotting finished" << endl;

	 phasefitgraphs[c]->SetPoint(i,mScaled,fitphase);

	 complex<double> fitval=compset.overlap(wl[iw],wl[iw2],m);

	 realfitgraphs[c]->SetPoint(i,mScaled,fitval.real());
	 imagfitgraphs[c]->SetPoint(i,mScaled,fitval.imag());


	 c++;
       }// end inner loop over waves

     } // end loop over waves
     //cerr << "outer loop over waves finished" << endl;
     // loop over components to fill individual graphs
     unsigned int compcount=0;
       for(unsigned int ic=0;ic<compset.n();++ic){
	 const pwacomponent* c=compset[ic];
	 std::map<std::string,pwachannel >::const_iterator it=c->channels().begin();
	 while(it!=c->channels().end()){
	   double I=norm(c->val(m)*it->second.C())*it->second.ps(m)*compset.calcFsmd(m);
	   compgraphs[compcount]->SetPoint(i,mScaled,I);
	   ++compcount;
	   ++it;
	 } // end loop over channels
       }// end loop over components

       cerr << "Finished plotting mass-bin " << m << endl;
       //mprev=m;
   }
   cerr << "Finished Loop Over DataBins" << endl;



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
        g->SetPoint(ib,m*MASSSCALE,ps*500.);
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

 if (fPS != NULL) {
   fPS->Write("funcFinalStateMassDependence");
 }

outfile->Close();

   return 0;

}
