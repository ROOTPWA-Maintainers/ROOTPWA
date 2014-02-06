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
//      Data storage class for PWA fit result of one mass bin
//
// Environment:
//      Software developed for the COMPASS experiment at CERN
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#include <algorithm>

#include "reportingUtils.hpp"
#include "TFitBin.h"


using namespace std;
using namespace rpwa;


ClassImp(TFitBin);


TFitBin::TFitBin()
  : _hasErrors(true)
{ }


TFitBin::~TFitBin()
{
  Reset();
}


void
TFitBin::Reset()
{
  _amps.clear();
  _wavenames.clear();
}


void 
TFitBin::fill(const std::vector<TComplex>&             prodAmps,
	      const std::vector<std::pair<int, int> >& indices,
	      const std::vector<TString>&              prodAmpNames,
	      const int                                nevt,
	      const unsigned int                       nmbEvents,
	      const double                             massBinCenter,
	      const complexMatrix&                     normIntegral,
	      const TMatrixD&                          fitParCovMatrix,
	      const double                             logLikelihood,
	      const int                                rank)
{
// !!! add some consistency checks
  _int.ResizeTo (normIntegral.nRows(),       normIntegral.nCols());
  _errm.ResizeTo(fitParCovMatrix.GetNrows(), fitParCovMatrix.GetNcols());
  _amps      = prodAmps;
  _indices   = indices;
  _wavenames = prodAmpNames;
  _nevt      = nevt;
  _rawevents = nmbEvents;
  _mass      = massBinCenter;
  _int       = normIntegral;
  if (_errm.GetNrows() == 0)  // check whether there really is an error matrix
    _hasErrors = false;
  _errm  = fitParCovMatrix;
  _logli = logLikelihood;
  _rank  = rank;
  buildWaveMap();
}


void
TFitBin::buildWaveMap() {
  int n=_wavenames.size();
  for(int i=0;i<n;++i){
    // strip rank
    TString title=wavetitle(i);
    //cout << "title=="<<title<<endl;
    if(find(_wavetitles.begin(),_wavetitles.end(),title)==_wavetitles.end()){
      //cout<<"Wave: "<<title<<endl;
      _wavetitles.push_back(title);
    }
    
    // look for index of first occurence
    int j;
    for(j=0;j<n;++j){
      if(_wavenames[j].Contains(title))break;
    }
    _wavemap[i]=j;
  }
}

// take care of flat wave which has no integral!
TComplex
TFitBin::getInt(int i, int j) const {
  // TODO: idea - prestore index of flat wave!
  bool flat1=_wavenames[i].Contains("flat");
  bool flat2=_wavenames[j].Contains("flat");
  if(flat1 && flat2) return 1;
  else if(flat1 || flat2) return 0;
  else {
    map<int, int>::const_iterator indexA = _wavemap.find(i);
    map<int, int>::const_iterator indexB = _wavemap.find(j);
    if ((indexA == _wavemap.end()) || (indexB == _wavemap.end())) {
      printWarn << "Amplitude index " << i << " or " << j << " is out of bound." << endl;
      return 0;
    }
    return _int(indexA->second, indexB->second);
  }
}

Int_t
TFitBin::rankofwave(int i) const {
  if(!_wavenames[i].Contains("flat")){
    return TString(_wavenames[i](1,1)).Atoi();
  }
  else return -1;
}

// OBSERVABLES CALCULATIONS:

double 
TFitBin::intens() const {  // total intensity
  return intens("");
}

double
TFitBin::intens(const char* tag) const {
  return norm(tag)*_nevt;
}

double 
TFitBin::norm(const char* tag) const // added intensity of waves containing tag
{
  TComplex c(0,0);
  int n=_wavenames.size();
  // loop over rank:
  for(int r=-1; r<_rank; ++r){
    // double loop
    for(int i=0;i<n;++i){// outer loop
      if(!(_wavenames[i].Contains(tag) && rankofwave(i)==r))continue;
      for(int j=0;j<n;++j){
	if(!(_wavenames[j].Contains(tag) && rankofwave(j)==r))continue;
	c+=_amps[i]*TComplex::Conjugate(_amps[j])*getInt(i,j);
	//cout<< _wavenames[i] << " | " << _wavenames[j] << endl;
      }
    }
  }
  return c.Re();
}

double
TFitBin::getParameter(const char* name) const {
  // check if real or imag
  TString na(name);
  bool re=false;
  if(na.Contains("flat")){
    na="flat";
    re=true;
  }
  else re=na.Contains("RE");
  // find wave -> remove _IM/_RE
  if(re)na.ReplaceAll("_RE","");
  else na.ReplaceAll("_IM","");
  unsigned int n=_wavenames.size();
  for(unsigned int i=0;i<n;++i){// outer loop
    if(na==_wavenames[i]){
      if(re)return _amps[i].Re();
      else return _amps[i].Im();
    }
  }
  // not found -> return standard
  return 0;
}


double
TFitBin::intens(int i) const {
  //cout<<wavename(i)<<"    r="<<rankofwave(i)<<endl;
  //cout<<wavetitle(i)<<endl;
  return norm(wavetitle(i).Data())*_nevt;
}

TComplex
TFitBin::spinDens(int i, int j) const {
  // account for find all waves with same title
  TString w1=wavetitle(i);
  TString w2=wavetitle(j);
  int n=_wavenames.size();
  TComplex result(0,0); // list for ranks
  TComplex a1(0,0);
  TComplex a2(0,0);
  for(int r=0; r<_rank; ++r){
    for(int k=0;k<n;++k){
      if(rankofwave(k)!=r)continue;
      if(_wavenames[k].Contains(w1))a1=_amps[k];
      if(_wavenames[k].Contains(w2))a2=_amps[k];
    }
    result+=a1*TComplex::Conjugate(a2);
  }
  
  return result;
}


double 
TFitBin::err(const char* tag) const {
  if(!_hasErrors)return 0;

  // double loop
  TMatrixD C;
  vector<int> cpar;
  getCov(tag,C,cpar);
  if(cpar.size()==0)return 0;
  //C.Print();

  // Now build jacobian of intensity
  TMatrixD J(2,cpar.size()*2);
  // number of sub-jacobians
  int nsub=cpar.size();
  for(int i=0;i<nsub;++i){//loop over subjacobians
    TMatrixD Ji(2,2);Ji[0][0]=0;Ji[0][1]=0;Ji[1][0]=0;Ji[1][1]=0;
    // which rank am I in?
    int currentrank=rankofwave(cpar[i]);
    //cout << "current rank ="<<currentrank<<endl;
    for(int j=0;j<nsub;++j){ // loop over j
      // only consider same rank contributions
      if(rankofwave(cpar[j])!=currentrank){
	//cout << " exclude this rank from Ji" << endl;
	continue;
      }
      //cout << "vpsi"<<cpar[j]<<" , "<<cpar[i]<<endl;
      TComplex vpsi=_amps[cpar[j]]*getInt(j,i);
      Ji[0][0]+=vpsi.Re();
      Ji[0][1]+=vpsi.Im();
    }
    // set submatrix
    J.SetSub(0,2*i,Ji);
  
  }// end loop over subjacobians
  //J.Print();
  J*=2;
  TMatrixD JT(TMatrixD::kTransposed,J);
  TMatrixD CJT=C*JT;
  TMatrixD x=J*CJT;
  return _nevt*TMath::Sqrt(x[0][0]);
}

double 
TFitBin::err(int i) const {
  //cout<<_wavenames[i]<<endl;
  return err(wavetitle(i).Data());
}


double
TFitBin::phase(int i, int j) const // phase difference between wave i and j
{ 
  return spinDens(i,j).Theta()/TMath::Pi()*180.;
}

double
TFitBin::phaseErr(int i, int j) const 
{
  if(!_hasErrors)return 0;
  TMatrixD C=spinDensErr(i,j);
  TComplex rho=spinDens(i,j);
  TMatrixD J(1,2);
  J[0][0]=1./(rho.Re()+rho.Im()*rho.Im()/rho.Re());
  J[0][1]=-rho.Im()/(rho.Re()*rho.Re()+rho.Im()*rho.Im());
  TMatrixD JT(TMatrixD::kTransposed,J);
  TMatrixD CJT=C*JT;
  TMatrixD x=J*CJT;
  return sqrt(x[0][0])/TMath::Pi()*180.;
}

TMatrixD
TFitBin::spinDensErr(int i, int j) const {
  TMatrixD C;
  getCov(i,j,C);
  TMatrixD J(2,4);
  TMatrixD Ji=M(TComplex::Conjugate(_amps[j]));
  J.SetSub(0,0,Ji);
  Ji=M(_amps[i]);
  TMatrixD Jc(2,2);Jc[0][0]=1;Jc[0][1]=0;Jc[1][0]=0;Jc[1][1]=-1;
  Ji*=Jc;
  J.SetSub(0,2,Ji);
  TMatrixD JT(TMatrixD::kTransposed,J);
  TMatrixD CJT=C*JT;
  TMatrixD x=J*CJT;
  return x;
}

TMatrixD
TFitBin::M(const TComplex& c) const {
  TMatrixD m(2,2);
  m[0][0]=c.Re();
  m[0][1]=-c.Im();
  m[1][0]=c.Im();
  m[1][1]=c.Re();
  return m;
}


double
TFitBin::coh(int i, int j) const // coherence between wave i and j
{ return 0;}


void
TFitBin::listwaves() const {
  for(unsigned int i=0; i<_wavenames.size();++i)cout<<i<<"  "<<_wavenames[i]<<endl;
}


TMatrixD 
TFitBin::getErr(pair<int,int> k) const {
  TMatrixD m(2,2);
  m[0][0]=_errm[k.first][k.first];
  m[0][1]=_errm[k.first][k.second];
  m[1][0]=_errm[k.second][k.first];
  m[1][1]=_errm[k.second][k.second];
  return m;
}


void
TFitBin::getCov(const char* tag, TMatrixD& C, vector<int>& cpar) const {
  // loop over wavenames to get list of paramters
  cpar.clear();
  int n=_wavenames.size();
  vector<int> par; // vector of (real parameters)
  for(int i=0;i<n;++i){// outer loop
    if(_wavenames[i].Contains(tag)){
      par.push_back(_indices[i].first);
      par.push_back(_indices[i].second);
      cpar.push_back(i); // vector of complex parameter indices
      //cout << _wavenames[i] << endl;
    }
  }
  // now build up matrix:
  int n2=par.size();
  C.ResizeTo(n2,n2);
  for(int i=0;i<n2;++i){
    for(int j=0;j<n2;++j){
      //cout<<"par="<<par[i]<<"   "<<par[j]<<endl;
      if(par[i]<0 || par[j]<0)C[i][j]=0;
      else C[i][j]=_errm[par[i]][par[j]];
    }
  }
}


void
TFitBin::getCov(int i, int j, TMatrixD& C) const {
  // now build up matrix:
  int n2=4;
  C.ResizeTo(n2,n2);
  vector<int> par; // vector of (real parameters)
  par.push_back(_indices[i].first);
  par.push_back(_indices[i].second);
  par.push_back(_indices[j].first);
  par.push_back(_indices[j].second);
  for(int i=0;i<n2;++i){
    for(int j=0;j<n2;++j){
      //cout<<"par="<<par[i]<<"   "<<par[j]<<endl;
      if(par[i]<0 || par[j]<0)C[i][j]=0;
      else C[i][j]=_errm[par[i]][par[j]];
    }
  }
}


void
TFitBin::PrintParameters() const {
  for(unsigned int i=0;i<_amps.size();++i){
    cout << _wavenames[i] << "  "  << _amps[i] << " +- ("
	 << sqrt( _errm[_indices[i].first][_indices[i].first]  )
	 << "," ;
    if(_indices[i].second>=0){
      cout << sqrt( _errm[_indices[i].second][_indices[i].second]  )
	   << ")" << endl;
    }
    else {
      cout << "0)" << endl;
    }
  } 

}

void 
TFitBin::printAmpsGenPW(ostream& s) const {
 for(unsigned int i=0;i<_amps.size();++i){
    s << _wavenames[i] << " "  
      << _amps[i].Re() << " "
      << _amps[i].Im() << endl;
 }
}




void 
TFitBin::getParameters(double* par) const {
  unsigned int c=0;
  for(unsigned int i=0;i<_amps.size();++i){
    par[c++]=_amps[i].Re();
    if(_indices[i].second>=0){ // check if this one exists
      par[c++]=_amps[i].Im();
    }
  }
}
