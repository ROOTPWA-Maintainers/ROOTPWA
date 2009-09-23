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
// Hamiltonian Markov Chain Monte Carlo Simlation of PWALikelihood
//
//    Copyright Sebastian Neubert, 2009
//
//    This file is part of rootpwa
//
//    Foobar is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Foobar is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
//

#include <iostream>
#include <list>
#include <vector>
#include <string>
#include <map>
#include <complex>
#include <assert.h>
#include <time.h>
#include "TBranch.h"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TComplex.h"
#include "TRandom3.h"
#include "TFitBin.h"
#include "TStopwatch.h"
#include "TPWALikelihood.h"
#include "TLogMultiGaus.h"
#include "TMCMCMeta.h"

//#include "TMinuitMinimizer.h"
#include "Minuit2/Minuit2Minimizer.h"
using namespace std;
using namespace ROOT::Minuit2;

int lineno = 1; // global variables needed for lex (not understood)
char *progname;


int main(int argc, char** argv){

  gRandom->SetSeed(time(NULL));

  //if(argc<2)return 1;
  char* progname=argv[0];

  // set options
  unsigned int Nit=5000; // number of iterations
  double low_mass=0;
  double high_mass=0;
  string wavefile_name; // wavelist filename
  TString ofile_name="MCMCresult.root"; // output filename
  TString start_file_name; // file with starting values
  double step=0.1;
  string norm_file_name; // file with normalization integrals
  bool quiet=false;
  unsigned int rank=1; // rank
  double Ttotal=14000; // max runtime <4 hours


  extern char *optarg;
  // extern int optind;
  int c;
  while ((c = getopt(argc, argv, "l:u:i:w:o:s:S:n:r:qh")) != -1)
    switch (c) {
    case 'i':
      Nit = atoi(optarg);
      break;
    case 'l':
      low_mass = atof(optarg);
      break;
    case 'u':
      high_mass = atof(optarg);
      break;
    case 'o':
      ofile_name = optarg;
      break;
    case 's':
      step = atof(optarg);
      break;
    case 'S':
       start_file_name= optarg;
      break;
    case 'n':
      norm_file_name = optarg;
      break;
    case 'w':
      wavefile_name = optarg;
      break;
    case 'r':
      rank = atoi(optarg);
      break;

    case 'q':
      quiet = true;
      break;
    case 'h':
      cerr << "usage:" << endl;
      cerr << progname << " -l # -u # [-i -h -q] [-n normfile -s start values -r rank] -w wavelist -o outfile " << endl;
      cerr << "    where:" << endl;
      cerr << "        -i #:     iterations" << endl;
      cerr << "        -t #:     maximal runtime in sec (default: 14000)" << endl;
      cerr << "        -l #:     lower edge of bin" << endl;
      cerr << "        -u #:     upper edge of bin" << endl;
      cerr << "        -w file:  set wavelist file" << endl;
      cerr << "        -o file:  set output file" << endl;
      cerr << "        -n file:  set normalization file" << endl;
      cerr << "        -s #:     set stepsize" <<endl;
      cerr << "        -r #:     set rank" <<endl;
      cerr << "        -q:       run quietly" << endl;
      cerr << "        -h:       print help" << endl;
      exit(0);
      break;
    }
  

  if(norm_file_name.length()<=1)norm_file_name="norm.int";
  if(wavefile_name.length()<=1){
    cerr << "No wavelist specified! Aborting!" << endl;
    return 1;
  }
  
  
  TPWALikelihood<double> L;
  if(quiet)L.SetQuiet();
  L.SetWavelist(wavefile_name);
  L.SetRank(rank);
  cout<<L.NDim()<<endl;
  L.LoadAmplitudes();
  L.LoadIntegrals(norm_file_name,norm_file_name);
  
  
  

  // TLogMultiGaus L;
  //vector<double> prec(2);
  //prec[1]=0.5;
  //prec[0]=1;
  //L.Set(prec);

  cout << "Installing parameters:" << endl;
  unsigned int ndim=L.NDim();
  double dL[ndim];   // gradient of loglikelihood
  double dLlast[ndim]; 
  double E;        // likelihood-value
  double H;        // Hamiltonian
  double X[ndim];  // position in phasespace
  double Xlast[ndim]; // position in PS for last MCMC step
  double R[ndim];  // vector of velocities in PS
  double Rlast[ndim];

  bool hasstart=false;
  if(start_file_name.Length()>2){
    TFile* file=TFile::Open(start_file_name,"READ");
    TFitBin* start=new TFitBin();
    // check if tree exists
    TTree* tree=(TTree*)file->Get("pwa");
    if(tree==NULL){
      cout << "pwa Tree not found for starting values." << endl;
    }
    else {
      tree->SetBranchAddress("fitbin",&start);
      // search for entry with same mass bin
      for(int i=0;i<tree->GetEntriesFast();++i){      
	tree->GetEntry(i);
	if(start->mass()<=high_mass && start->mass()>=low_mass){
	  cout << "Start values: found suiting mass bin: " << start->mass() <<endl;
	  break;
	}
      }
      // if no suiting mass bin is found use last one
      start->getParameters(X);
      hasstart=true;
    }
    file->Close();
  } // endif (start_file_name.Length()>2)
  
  for(unsigned int i=0;i<L.NDim();++i){
    if(!hasstart)X[i]=0.1; // set default values
    cout << L.parname(i) << " " << X[i] << endl;
    R[i]=gRandom->Gaus(0,1); // starting values of velocities
  }
  


  // Output tree with Markov Chain samples
  TFile* file=TFile::Open(ofile_name,"RECREATE");

  

  TTree* tree=new TTree("MarkovChain","MarkovChain");
  tree->Branch("N",&ndim,"N/I");
  tree->Branch("dL",&dL,"dL[N]/D");
  tree->Branch("X",&X,"X[N]/D");
  tree->Branch("L",&E,"L/D");
  tree->Branch("H",&H,"H/D");
  tree->Branch("R",&R,"R[N]/D");
  tree->Branch("Rlast",&Rlast,"Rlast[N]/D");

   

  unsigned int Nacc=0;
  unsigned int Nleap=5;  // max number of leap frog steps per iteration
  double eps=step; // time step
  double eps2=0.5*eps;
  double m=0.1; // virtual mass

  // write meta info ***********************************
  TMCMCMeta* meta=new TMCMCMeta();
  for(unsigned int i=0;i<L.NDim();++i){
    meta->parnames.push_back(L.parname(i));    
  }
  meta->stepsize=eps;
  meta->virtmass=m;
  meta->Nleap=Nleap;
  meta->low_mass=low_mass;
  meta->high_mass=high_mass;
  meta->NEvents=L.nevents();
  L.getIntCMatrix(meta->Norm,meta->Norm); // BE CAREFULL HERE!!!
  meta->Write("MetaInfo");


  TStopwatch watch;
  watch.Start();
  double Tmean=0; // mean time for one simulation step;
  double Tacc=0;   // accumulated time
  
  // Start itertions
  for(unsigned it=0;it<Nit;++it){ // iteration loop
    double tstart=watch.CpuTime();watch.Continue();
    cout << "Iteration" << it << endl;
    
    // decide which way to go:
    //double sign= gRandom->Uniform()>0.5 ? +1 : -1;
    //eps*=sign;
    //eps2*=sign;

    // Calculate gradient at starting point
    L.FdF(X,E,dL);
    // Save Starting values;
    if(it==0)tree->Fill();

    double Hlast=0;  // Hamilton Function at last step
    // calculate Hamiltonian:
    // add kinetical term
    for(unsigned int ip=0; ip<ndim;++ip){ // loop over velocities
      Hlast+=R[ip]*R[ip];
    } // end loop over velocities
    Hlast*=0.5/m;
    Hlast+=E;


    // Shift velocity update by eps/2:
    for(unsigned int ip=0; ip<ndim;++ip){ // loop over parameters
      // store previous 
      Xlast[ip]=X[ip];
      dLlast[ip]=dL[ip];
      Rlast[ip]=R[ip];
      R[ip]=R[ip]-eps2*dL[ip]; // update velocities at middle of step
    }// end loop over parameters

    unsigned int nl=gRandom->Integer(Nleap);
    for(unsigned il=0; il<nl;++il){ // leap frog loop
      for(unsigned int ip=0; ip<ndim;++ip){ // loop over postions
	X[ip]=X[ip]+eps/m*R[ip]; // update position
      }// end loop over positions
      // Calculate Likelihood update with gradient
      
      L.FdF(X,E,dL);
      // check for last step which has to be a half step again!
      if(il<nl-1){
	for(unsigned int ip=0; ip<ndim;++ip){ // loop over velocities
	  R[ip]=R[ip]-eps*dL[ip];
	} // end loop over velocities
      }
      else {
	for(unsigned int ip=0; ip<ndim;++ip){ // loop over velocities
	  R[ip]=R[ip]-eps2*dL[ip];
	} // end loop over velocities
      }
    } // end Hamilton integration through Nleap steps
    // Do Metropolis Hastings Step
    // calculate Hamiltonian:
    H=0;
    // add kinetical term (and perturb velocities):
    double P=0;
    for(unsigned int ip=0; ip<ndim;++ip){ // loop over velocities
	H+=R[ip]*R[ip];
	R[ip]=gRandom->Gaus(0,1); // Gibbs step for velocity
	P+=R[ip]*R[ip];
    } // end loop over velocities
    // Normalize Updated Velocity vector
    double plength=gRandom->Gaus(0,1);
    plength/=sqrt(P);
    for(unsigned int ip=0; ip<ndim;++ip)R[ip]*=plength;


    H*=0.5/m;
    H+=E;
    // in principle H should be constant! Appart from numerical errors
    cout << ">>>> (Hlast-H)/H = " << (Hlast-H)/H << endl;
    if(H>Hlast){
      tree->Fill();  
      ++Nacc;
    }
    else{ 
      if(gRandom->Uniform()<TMath::Exp(H-Hlast)){
	tree->Fill();  
	++Nacc;
      }
      else {
	// Restore old values
	//for(unsigned int ip=0; ip<ndim;++ip){ // loop over postions
	//  // store previous 
	//  X[ip]=Xlast[ip];
	//  dL[ip]=dLlast[ip];
	//}// end loop over positions
	//tree->Fill();
      }
    }
      
    cout << "Acceptance Rate: " << (double)Nacc/(double)(it+1) << endl;
    double tend=watch.CpuTime();watch.Continue();
    
    // recalculate mean tiem for step;
    Tacc+=tend-tstart;
    Tmean=Tacc/(double)(it+1);
    
    if(it%10==0){
      cout << "Mean time per iteration: "<< Tmean << "sec"<<endl;
      cout << "Max time left for this job: " 
	   << Ttotal-Tacc << "sec ("<< Tacc/Ttotal*100 << "%)" << endl;
    }

    // Stop if we are close to time running out
    if(Tacc+5*Tmean>Ttotal)break;
  } // end iteration loop

  
  tree->Write();
  file->Close();
  cout << "Acceptance Rate: " << (double)Nacc/(double)Nit << endl;


  return 0;
/*
  // write out result ------------------------------------------------------
  // ----------------------------------------------------------
  TFile* outfile=new TFile(ofile_name,"RECREATE");
  TFitBin* result=new TFitBin();

  
  double mass=0.5*(low_mass+high_mass);
  vector<TString> wavenames;
  const vector<string>& wnames=L.wavenames();
  for(int i=0;i<wnames.size();++i){
    wavenames.push_back(TString(wnames[i].c_str()));
  }

  double logli=mini.MinValue();
  vector<TComplex> amps;
  vector<complex<double> > V;
  L.buildCAmps(mini.X(),V,true);
  // convert to TComplex;
  for(int i=0;i<V.size();++i)amps.push_back(TComplex(V[i].real(),V[i].imag()));

  
 // error matrix;
  int npar=L.NDim();
  TMatrixD errm(npar,npar);
  for(int i=0;i<npar;++i){
    for(int j=0;j<npar;++j){
      errm[i][j]=mini.CovMatrix(i,j);
    }
  }

  int n=wavenames.size();
  // indices of real and imag part of parameter in error matrix;
  int parcount=0;
  vector<pair<int,int> > indices(npar); 
  for(int ir=0;ir<rank;++ir){
    for(int ia=ir;ia<2*n;++ia){
      if(ia==ir){
	indices[parcount]=pair<int,int>(parcount,-1);
	++parcount;
      }
      else {
	indices[parcount]=pair<int,int>(parcount,parcount+1);
	parcount+=2;
      }
    }
  }
  assert(parcount=n+1);

  // get normalization integral
  integral norm=L.normInt();
  TCMatrix integr;
  integr.ResizeTo(n,n);
   for(int i=0;i<n;++i){// begin outer loop
     TString w1=wavenames[i];
     for(int j=0;j<n;++j){
       TString w2=wavenames[j];
       if(w1=="flat" || w2=="flat"){
	 if(w1==w2)integr.set(i,j,complex<double>(1,0));
	 else integr.set(i,j,complex<double>(0,0));
       }
       else integr.set(i,j,norm.val(w1.Data(),w2.Data()));
     }// end inner loop
   } // end outere loop over parameters;

   cout << "filling TFitBin" << endl;

  result->fill(amps,
	       indices,
	       wavenames,
	       L.nevents(),
	       mass,
	       integr,
	       errm,
	       logli,
	       rank);
    
  TString binname="fitbin";binname+=low_mass;binname+="_";binname+=high_mass;
  binname.ReplaceAll(" ","");
  result->Write(binname);
  */
  return 0;
}

// dummy function needed since we link to but do not use minuit.o
int mnparm(int, string, double, double, double, double) {
    cerr << "this is impossible" << endl;
    throw "aFit";
    return 0;
}
