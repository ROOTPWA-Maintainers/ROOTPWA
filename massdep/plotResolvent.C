#include <vector>
#include <iostream>
#include <fstream>
#include <complex>

#include "TMath.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TString.h"
#include "TAxis.h"

using namespace std;

double const mpi=0.13957018;
double const mpi2=mpi*mpi;
double const pi=TMath::Pi();

typedef complex<double> cnum;

// analytic function J eqn 39
cnum J(cnum z){
  if(z.real()==0)return cnum(1,0);
  cnum zm=sqrt(z-4*mpi2);
  cnum term1=0.5*zm/sqrt(z);
  cnum num=sqrt(z)+zm;
  cnum denom=sqrt(z)-zm;
  cnum term2=log(num/denom) - cnum(0,pi);
  return term1*term2;
}

cnum u(unsigned int a, cnum t){
  if(a==0)return J(t);
  cnum tt=t-4*mpi2;
  if(a==1) return sqrt(tt);
  else if(a==2) return tt;
  else if(a==3) return tt*tt;
  return 0;
}

int plotResolvent(char* name){

  unsigned int nIso=2;
  unsigned int nBase=4;

  ifstream infile(name);
  complex<double> s;
  
  vector<complex<double> > IA(nBase);
  vector<TString> IsoNames(nIso*nIso);

  vector<TMultiGraph*> graphsRe(nIso*nIso);
  for(unsigned int i=0;i<nIso*nIso;++i)graphsRe[i]=new TMultiGraph();
  vector<TMultiGraph*> graphsIm(nIso*nIso);
  for(unsigned int i=0;i<nIso*nIso;++i)graphsIm[i]=new TMultiGraph();

 
  double dt=0.05;
  double t0=0.001;
  unsigned int count=0;
  // read values
  while(infile.good()){
    infile >> s;
    // calculate t-range
    double tHat=sqrt(s.real()-mpi)*sqrt(s.real()-mpi);
    unsigned int nt=(unsigned int)floor(tHat/(double)dt);
    string name;
    for(unsigned int i=0;i<nIso;++i){
      for(unsigned int j=0;j<nIso;++j){
	// create graphs;
	TGraph* gRe=new TGraph(nt);
	TGraph* gIm=new TGraph(nt);
	graphsRe[i*nIso+j]->Add(gRe);
	graphsIm[i*nIso+j]->Add(gIm);
	infile >> name;
	if(count==0)IsoNames[i*nIso+j]=name.c_str();
	TString label("RE-");label+=name.c_str();
	gRe->SetTitle(label);
	label.ReplaceAll("RE","IM");
	gIm->SetTitle(label);
	

	for(unsigned int a=0;a<nBase;++a){
	  infile >> IA[a];
	}// end reading loop over basis functions 
	
	// loop over t and draw graphs
	for(unsigned it=0;it<nt;++it){
	  double t=t0+it*dt;
	  complex<double> I=0;
	  for(unsigned int a=0;a<nBase;++a){
	    I+=IA[a]*u(a,t);
	  } // end loop over basis 
	  // fill graph
	  gRe->SetPoint(it,t,I.real());
	  gIm->SetPoint(it,t,I.imag());
	} // end loop over t

      } // end loop over isobars
    } // end loop over isobars
    ++count;
  } // end loop over s
  
  // Draw results

  TCanvas* c=new TCanvas("c","c",10,10,1000,1000);
  c->Divide(nIso*2,nIso);
  // loop over Isobars
   for(unsigned int i=0;i<nIso;++i){
      for(unsigned int j=0;j<nIso;++j){
	c->cd(1+i*2*nIso+j*2);
	
	graphsRe[i*nIso+j]->Draw("AC");
	graphsRe[i*nIso+j]->GetXaxis()->SetTitle("(m_{2#pi})^{2}");
	graphsRe[i*nIso+j]->GetYaxis()->SetTitle(IsoNames[i*nIso+j]+"-RE");
	graphsRe[i*nIso+j]->GetYaxis()->SetTitleOffset(1.2);
	c->cd(1+i*2*nIso+j*2+1);
	graphsIm[i*nIso+j]->Draw("AC");
	graphsIm[i*nIso+j]->GetXaxis()->SetTitle("(m_{2#pi})^{2}");
	graphsIm[i*nIso+j]->GetYaxis()->SetTitle(IsoNames[i*nIso+j]+"-IM");
	graphsIm[i*nIso+j]->GetYaxis()->SetTitleOffset(1.4);
      }
   } // end loop over isobars


  return 0;
  
}
