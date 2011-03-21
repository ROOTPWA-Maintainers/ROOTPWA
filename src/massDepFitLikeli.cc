#include "massDepFitLikeli.h"
#include "TTree.h"
#include "TMatrixD.h"
#include "fitResult.h"

#include <vector>
#include <complex>
#include <iostream>

using namespace std;

unsigned int
rpwa::massDepFitLikeli::NDim() const {return _compset->numPar();}

unsigned int
rpwa::massDepFitLikeli::NDataPoints() const {
  // loop through input data an check how many bins are in massrange
  unsigned int nbins=_tree->GetEntries();
  cout << nbins << " mass bins in input data." << endl;
  unsigned int nbinsInRange=0;
  // loop over mass-bins
  for(unsigned im=0;im<nbins;++im){
    _tree->GetEntry(im);
    double mass=_rhom->massBinCenter();
    if(mass>=_mmin && mass<=_mmax)++nbinsInRange;
  }
  cout << nbinsInRange << " mass bins in allowed mass range {"
       << _mmin<<","<<_mmax<<"}"<< endl;
  
  // calculate data points, remember (Re,Im)=> factor 2:
  // a complex, symmetric matrix with real diagonal has n^2 elements:
  return nbinsInRange*_wlist.size()*_wlist.size();
}

void
rpwa::massDepFitLikeli::init(TTree* sp, pwacompset* compset,
			     double mmin, double mmax){
  _mmin=mmin;
  _mmax=mmax;
  _compset=compset;
  _tree=sp;
  _rhom=NULL;//new fitResult();
  if(_tree->SetBranchAddress("fitResult_v2",&_rhom))cerr<<"Branch not found!"<<endl;

  _tree->GetEntry(0);
  
  // build waveindex map;
  _wlist.clear();
  _wlist=_compset->wavelist();
  _index.clear();

  cerr << "Number of components: "<<_compset->n()<< endl;
  cerr << (*_compset) << endl;

  cerr << "Number of waves in fit: "<< _wlist.size() << endl;
 for(unsigned int i=0;i<_wlist.size();++i){
   _index.push_back(_rhom->waveIndex(_wlist[i]));
   cerr <<  _wlist[i] << endl;
 }

 

}


double
rpwa::massDepFitLikeli::DoEval(const double* par) const {
  // rank==1 mass dependent fit
  
  // input values: m
  unsigned int nwaves=_index.size();
     
  _compset->setPar(par);

  // Breitwigner parameters
  //masses;
  //widths;
  
  double chi2=0;
  unsigned int nbins=_tree->GetEntries();
  
  

  // loop over mass-bins
  for(unsigned im=0;im<nbins;++im){
    _tree->GetEntry(im);
    double mass=_rhom->massBinCenter();
    if(mass<_mmin)continue;
    if(mass>_mmax)continue;
    //if(mass>2000)continue;
    //cout << "Mass=" << mass << endl;
    // inpu values: measured spin density matrix
    
       
    // sum over the contributions to chi2 -> rho_ij
      for(unsigned int i=0; i<nwaves; ++i){
	for(unsigned int j=i; j<nwaves; ++j){
	  // calculate target spin density matrix element
	  complex<double> rho;
	  // loop over all waves
	  complex<double> f1;
	  complex<double> f2;
	  const string w1=_wlist[i];
	  const string w2=_wlist[j];
	  // loop over components and look which contains these waves
	  
	  for(unsigned int k=0;k<_compset->n();++k){
	    if((*_compset)[k]->channels().count(w1)>0)
	      f1+=(*_compset)[k]->val(mass)* ((*_compset)[k]->channels()).find(w1)->second.C();
	    if((*_compset)[k]->channels().count(w2)>0)
	      f2+=(*_compset)[k]->val(mass)* ((*_compset)[k]->channels()).find(w2)->second.C(); 
	  }
	  f1*=_rhom->phaseSpaceIntegral(_index[i]);
	  f2*=_rhom->phaseSpaceIntegral(_index[j]);
	  rho=f1*conj(f2);
	  // compare to measured spin density matrix element
	  complex<double> rhom=rho-_rhom->spinDensityMatrixElem(_index[i],_index[j]);
	    
	  double dchi=norm(rhom);
	  TMatrixT<double> cov= _rhom->spinDensityMatrixElemCov(_index[i],
								_index[j]);
	  //cov.Print();
	  if(i==j)dchi/=cov[0][0];
	  else {
	    
	    dchi=rhom.real()*rhom.real()*cov[1][1]-(2*cov[0][1])*rhom.real()*rhom.imag()+cov[0][0]*rhom.imag()*rhom.imag();
	    dchi/=cov[0][0]*cov[1][1]-cov[0][1]*cov[0][1];

	    //dchi=rhom.real()*rhom.real()/cov[0][0]+rhom.imag()*rhom.imag()/cov[1][1];

	  }
	  
	  //cerr << "d-rho("<<i<<","<<j<<")=" << dchi<<endl;;
	  //cerr << "sigma ="<< sigma << endl;
	  //if(i==j)
	  chi2+=dchi;
	}
      }// end loop over i
  } // end loop over mass-bins
  
  return chi2;
} // end DoEval
