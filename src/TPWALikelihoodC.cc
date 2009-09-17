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
//      Implementation of class TPWALikelihoodC
//      see TPWALikelihoodC.hh for details
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

// This Class' Header ------------------
#include "TPWALikelihoodC.h"

// C/C++ Headers ----------------------
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <list>
#include <algorithm>
using namespace std;
using std::cout;
using std::endl;
using std::ifstream;
using std::pair;
using std::stringstream;
#include <assert.h>

// Collaborating Class Headers --------
#include "TString.h"
#include "TMath.h"
#include "TSystem.h"
#include "TCMatrix.h"
#include "TStopwatch.h"
#include "TPWARealConstraint.h"
#include "TPWANullConstraint.h"
#include "TPWAPhaseConstraint.h"

// Class Member definitions ----------


TPWALikelihoodC::TPWALikelihoodC() : _rank(1),_dim(0),_nposrefl(0),_nnegrefl(0), _nevents(0),_nwaves(0),_ncalls(0), _Ltime(0), _Ntime(0), _quiet(false),_useNorm(true),_Vflat("flat__+",0),_NumNullConstraints(0)
{
  _Vflat.setConstraint(new TPWARealConstraint());
}

TPWALikelihoodC::~TPWALikelihoodC()
{
  ClearCache();
}

void
TPWALikelihoodC::ClearCache(){
  for(unsigned int i=0;i<_data.size();++i){
    if(_data[i]!=NULL){
      _data[i]->clear();
      delete (_data[i]);
      _data[i]=0;
    }
  }
  _data.clear();
}


void 
TPWALikelihoodC::Init(const TString& wavelist, unsigned int rank,
		      const TString& norm,
		      const TString& acceptance,
		      int scaleevents,
		      const TString constraints){
  cout<<"Reading Integrals "<<endl;
  ifstream file(norm.Data());
  _I.scan(file);
  //if(!_quiet)_I.print(cout);
  file.close();
  _mat=_I.mat();
  list<string> normwaves=_I.files();
  // load acceptance integral
  file.open(acceptance.Data());
  _Acc.scan(file);
  //if(!_quiet)_Acc.print(cout);
  file.close();
  _Acc.events(scaleevents);
  _accmat=_Acc.mat();
  list<string> accwaves=_Acc.files();

  cout<<"Reading amplitude names from wavelist "<<wavelist<<endl;
  // clear production amplitude vector
  _V.clear();
  _V.reserve(_rank);
  _parnames.clear();
  _rank=rank;
  // setup vectors
  for(unsigned int ir=0;ir<rank;++ir){
    _V.push_back(new vector<TPWAAmp>);
  }
  assert(_V.size()==_rank);
  // read wavelist and fill wavelist vector
  unsigned int npos=0; // number of positive reflectivity waves
  unsigned int nneg=0; // number of negative reflectivity waves

  file.open(wavelist.Data());
  string line;
  while(file.good()){
    getline(file,line);
    unsigned int pos=line.find(" ");
    string name=line.substr(0,pos);
    double thres;
    if(pos<line.length()){
      thres=atof(line.substr(pos,line.length()).c_str());
    } else thres=0;

    if(!_quiet){
      cout << name << " thresh=" << thres << endl;
    }
    if(line.length()>1){
      // check normalization and acceptance integral index
      if(find(normwaves.begin(),normwaves.end(),name)==normwaves.end()){
	cout << "Wave " << name 
	     << " not in normalization integrals! Aborting!" << endl;
	throw;
      }
      unsigned int normindex=_I.index(name);
      if(find(accwaves.begin(),accwaves.end(),name)==accwaves.end()){
	cout << "Wave " << name 
	     << " not in acceptance integrals! Aborting!" << endl;
	throw;
      }
      unsigned int accindex=_Acc.index(name);

      _wavenames.push_back(name);

      // assign to rank
      for(unsigned int ir=0;ir<_rank;++ir){
	// positivity constraints:
	TPWAAmp amp(name,ir,thres,normindex,accindex);
	unsigned int nrefl;
	if(amp.reflectivity()==1){
	  nrefl=npos;
	  if(ir==_rank-1)++npos;
	}
	else{
	  nrefl=nneg;
	  if(ir==_rank-1)++nneg;
	}
	if(nrefl<ir){
	  amp.setConstraint(new TPWANullConstraint());
	  ++_NumNullConstraints; // number of null-constraint waves
	}
	if(nrefl==ir)amp.setConstraint(new TPWARealConstraint());
	_V[ir]->push_back(amp);

	if(!_quiet){
	  cout << amp.name() << "  rank" << ir << "   " << amp.type() << endl;
    }
      } // end loop over rank
    } // end checking line length
  } // end reading waveset file

  // IMPORTANT TO PUT THIS HERE!!!
  _nwaves=_V[0]->size();

  // apply constraints
  if(constraints.Length()>1){
    cout << "Loading Constraint Config: "<< constraints << endl;
    addConstraints(constraints);
  }
  // calculate dimension of parameter vector
  _dim=0;
  for(unsigned int iamp=0;iamp<_nwaves; ++iamp){
    for(unsigned int ir=0;ir<_rank;++ir){
      _dim+=(*_V[ir])[iamp].npar();
      for(int ipar=0;ipar<(*_V[ir])[iamp].npar();++ipar){
	_parnames.push_back((*_V[ir])[iamp].parname(ipar));
	_parthresholds.push_back((*_V[ir])[iamp].threshold());
      }
    }
  }

 
  
  // add flat
   _dim+=_Vflat.npar();
   _parnames.push_back("V_flat_RE");

   assert(_dim==_parnames.size());
  _parcache.resize(_dim);
  _dLcache.resize(_dim);
  // Number of waves from first rank
  

  // prepare some memory:
  for(unsigned int ir=0;ir<_rank;++ir){
    _dL.push_back(new vector<complex<double> >(_nwaves,complex<double>(0,0)));
    _dl.push_back(new vector<complex<double> >(_nwaves,complex<double>(0,0)));
    _dLN.push_back(new vector<complex<double> >(_nwaves,complex<double>(0,0)));
  }
}

// reads in all amplitude files and stores values in memory
void 
TPWALikelihoodC::LoadAmplitudes(){
  // for normalization the integrals need to be loaded first!
  if( _mat.nrows()==0 ) {
    cout << "Normalization integrals have to be loaded before amplitudes! "
	 << " ABORTING!" << endl;
    throw;
  }
  
  ProcInfo_t info1, info2;
  gSystem->GetProcInfo(&info1);
  if(!_quiet){
    cout << "Memory before: " << info1.fMemResident << endl;
    cout<<"Loading amplitudes"<<endl;
  }
  
  // first clear cache
  ClearCache();

  // loop over amplitudes and read in data:
  
  for(unsigned int i=0;i<_nwaves;++i){
    ifstream ampfile(_V[0]->at(i).name().c_str());
    if(!ampfile.good()){
      cout<<" File not found! Aborting!"<<endl;
      throw;
    }
    // get integral
    unsigned int index=_V[0]->at(i).normindex();
    complex<double> I=_mat.el(index,index);
    // create cache:
    vector<complex<double> >* amps=new vector<complex<double> >();
    amps->reserve(_nevents);                
    _data.push_back(amps);
    // read file
    complex<double> a1;
    while (ampfile.read((char*) &a1,sizeof(complex<double>))){
      //cout << a1 << endl;
      // if option is switched on -> renormalize:
      if(_useNorm){
	// rescale decay amplitude
	a1/=sqrt(I.real());
      }
      amps->push_back(a1);
    }
    if(!_quiet)cout<<amps->size()<<" events"<<endl;
    _nevents=amps->size(); 
    // TODO: cross check that all files have the same number of events!
  }
  
  if(!_quiet){
  gSystem->GetProcInfo(&info2);
  cout << "Memory after: " << info2.fMemResident << endl;
  cout << "Memory used for amplitudes: " 
       << info2.fMemResident-info1.fMemResident << endl;
  }// endif quiet

  // ***********************************************************
  // rescale integrals if necessary
  // remember: we have two copies of the matrix! -> change both!
  // so far only one is changed!!!!!!!!!!!!!!!
  if(_useNorm){
    if(!_quiet)cout << "Renormalizing Integrals ... ";
    // Normalization:
    for(unsigned int i=0;i<_nwaves;++i){// outer loop over waves
      unsigned int index1=_V[0]->at(i).normindex();
      double Ni=sqrt(_mat.el(index1,index1).real());
      for(unsigned int j=0;j<_nwaves;++j){
	if(i==j)continue; // do diagonal terms later!
	unsigned int index2=_V[0]->at(j).normindex();
	double Nj=sqrt(_mat.el(index2,index2).real());
	_mat.el(index1,index2)/=(Ni*Nj);
      } // end inner loop
    }// end out loop
    if(!_quiet)cout << "Renormalizing Acceptance Integrals ... " << endl;
    if(!_quiet)cout << "Note: Diagonal Terms of Renorm.Integrals are not = 1 here!";
    // Acceptance:
    // Remember that the matrices _mat and _accmat are already normalized
    // to the number of montecarlo events!
    for(unsigned int i=0;i<_nwaves;++i){// outer loop over waves
      unsigned int index1=_V[0]->at(i).normindex();
      unsigned int accindex1=_V[0]->at(i).accindex();
      double Ni=sqrt(_mat.el(index1,index1).real());
      for(unsigned int j=0;j<_nwaves;++j){
	unsigned int index2=_V[0]->at(j).normindex();
	unsigned int accindex2=_V[0]->at(j).accindex();
	double Nj=sqrt(_mat.el(index2,index2).real());
	_accmat.el(accindex1,accindex2)/=(Ni*Nj);
	if(!_quiet) cout <<" Norm("<<i<<j<<"): "
	     <<_mat.el(index1,index2) 
	     << "   Acc("<<i<<j<<"): "
	     <<_accmat.el(accindex1,accindex2)<< endl;
      } // end inner loop
    }// end out loop
    
    // do diagonal terms of normalization only after all other things have 
    // been done
    if(!_quiet)cout << "Renormalizing Diagonal terms ... " << endl;
    for(unsigned int i=0;i<_nwaves;++i){
      unsigned int index1=_V[0]->at(i).normindex();
      _mat.el(index1,index1)=1;
      if(!_quiet) cout <<" Norm("<<i<<i<<"): "
		       << _mat.el(index1,index1) << endl;
    }
    if(!_quiet)cout << "...finished." << endl;
  }// end if useNorm
} // end load amplitudes



void
TPWALikelihoodC::partoamp(const double* x) const {
  double par[2];
  unsigned int k=0;
  
  for(unsigned int r=0;r<_rank;++r){
    vector<TPWAAmp>* V=_V[r];
    unsigned int nwaves=V->size();
    for(unsigned int i=0;i<nwaves;++i){
      unsigned int npar=(*V)[i].npar();
      for(unsigned int ipar=0;ipar<npar;++ipar){
	par[ipar]=x[k++];
      }
      (*V)[i].setPar(par);
    } // loop over waves
  } // end loop over rank
  // loop again to update contraints
  for(unsigned int r=0;r<_rank;++r){
    for(unsigned int i=0;i<_nwaves;++i){
      (*_V[r])[i].updateAmp();
    } // loop over waves
  } // end loop over rank
  // do background:
  par[0]=x[k];par[1]=0;
  assert(k==_dim-1);
  _Vflat.setPar(par);
}

void
TPWALikelihoodC::amptopar(double* x) const {
  unsigned int k=0;
  for(unsigned int r=0;r<_rank;++r){
    unsigned int nwaves=_V[r]->size();
    for(unsigned int i=0;i<nwaves;++i){
      unsigned int npar=(*_V[r])[i].npar();
      for(unsigned int ipar=0;ipar<npar;++ipar){
	x[k++]=(*_V[r])[i].par(ipar);
      }
    } // loop over waves
  } // end loop over rank
  x[k]=_Vflat.par(0);
}

void
TPWALikelihoodC::Print() const {
  for(unsigned int r=0;r<_rank;++r){
    unsigned int nwaves=_V[r]->size();
    for(unsigned int i=0;i<nwaves;++i){
      cout << (*_V[r])[i] << endl;
    }
  }
  cout << _Vflat << endl;
}


// Likelyhood and first derivatives in one go: *****************************

void 
TPWALikelihoodC::FdF(const double* x, double& f, double* df) const {
  if(!_quiet)cout << " Calling FdF !!!!!!!!!!! " << endl;
  
  ++_ncalls;

  // check timing
  TStopwatch timer;
  timer.Start(true);

  // parameter cache:
  for(unsigned int ix=0;ix<_dim;++ix)(*const_cast<vector<double>*>(&_parcache))[ix]=x[ix];

  // build complex numbers from parameters
  // remember rank restrictions!
  // vector<complex<double> > V(_rank*_nwaves);

  // load parameters into amplitudes
  partoamp(x);

   // loglikelyhood:
  double L=0;
  // reset derivatives
  for(unsigned int ir=0;ir<_rank;++ir){
    for(unsigned int ia=0;ia<_nwaves;++ia){
      (*_dL[ir])[ia]=complex<double>(0,0);
      (*_dl[ir])[ia]=complex<double>(0,0);
      (*_dLN[ir])[ia]=complex<double>(0,0);
    }
  }

  double dLdflat=0;
  // amplitude:
  // two contributions for two reflectivities
  complex<double> Aplus(0,0);
  complex<double> Aminus(0,0); 
  // loop over events and calculate log likelihood
  // and calculate derivatives with respect to parameters
  for(unsigned int ievt=0;ievt<_nevents;++ievt){
    //if(ievt<3 || _nevents-ievt<4)cout << "event"<< ievt << "_____________________" << endl;
    double l=0; // likelihood contribution of this event
    // incoherent sum over ranks: loop
    for(unsigned int r=0;r<_rank;++r){
      for(unsigned int iamp=0;iamp<_nwaves;++iamp){ // loop over waves
	// contribution to likelihood:
	complex<double> a=(*_V[r])[iamp] * ((*_data[iamp])[ievt]);
	int thisrefl=(*_V[r])[iamp].reflectivity();
	if(thisrefl==-1)Aminus+=a;
	else Aplus+=a;
	//if(ievt<3 || _nevents-ievt<4)cout << "refl="<<thisrefl<< a << endl;
	// contributions to derivatives:
	for(unsigned int k=0;k<_nwaves; ++k){ // loop over derivatives 
	  // loop (only inside current rank) and only for waves with same refl
	  if(thisrefl==(*_V[r])[k].reflectivity())(*_dl[r])[k]+=a;
	} // end loop over derivatives
	
      } // end loop over waves
      l+=std::norm(Aplus);
      l+=std::norm(Aminus);
      Aplus.real()=0;Aplus.imag()=0;Aminus.real()=0;Aminus.imag()=0;
      assert(l>=0);
      for(unsigned int k=0;k<_nwaves; ++k){ // again loop over derivatives
	// loop (only inside current rank)
	(*_dl[r])[k]*=std::conj((*_data[k])[ievt]);
      } // end loop over derivatives
    
    }// end loop over rank
    l+=std::norm(_Vflat.amp());
    L-=TMath::Log(l);
    //loop over derivatives to incorporate factor 2/sigma
    double g=2./l;
    for(unsigned int ir=0;ir<_rank;++ir){
      for(unsigned int id=0;id<_nwaves; ++id){
	(*_dL[ir])[id]-=(*_dl[ir])[id]*g;
	(*_dl[ir])[id].real()=0;
	(*_dl[ir])[id].imag()=0;
      }
    }
    dLdflat-=_Vflat.amp().real()*g;
  }// end loop over events for calculation of intensity term
 
   
  double t1=timer.RealTime();
  timer.Start(true);

  // add normalization integrals ?? TODO: can we exploit symmetry here?
  complex<double> N(0,0);
  unsigned int outcount=0; // output array counter
  
  // event number normalization
  double nevt;
  double n2;
  if(!_useNorm){
    nevt=(double)_nevents;
    n2=2*nevt;
  }
  else {
    nevt=1;
    n2=2;
  }

  for(unsigned int r=0; r<_rank;++r){ // loop over rank **********************
    //cout <<"/n Rank"<<r<<endl;
    for(unsigned int iamp=0;iamp<_nwaves;++iamp){// outer loop *************
      TPWAAmp& ampi=(*_V[r])[iamp];
      for(unsigned int jamp=0;jamp<_nwaves;++jamp){ // inner loop ********
	TPWAAmp& ampj=(*_V[r])[jamp];
	// check if i and j are of the same reflectivity
	if(ampi.reflectivity()!=ampj.reflectivity())continue;
	complex<double> p=ampi.amp()*std::conj(ampj.amp());
	//cout<<_wavenames[iamp]<<"   "<<_wavenames[jamp]<<endl;
	complex<double> I=_accmat.el(ampi.accindex(),ampj.accindex());
	//cout<<"I="<<I<<endl;
	//cout<<"p="<<p<<endl;
	p*=I;
	N+=p;
	// calculate contribution to derivative
	(*_dLN[r])[iamp]+=ampj.amp()*std::conj(I);
      }// end inner loop over waves **************************************
      // account for 2nevents:
      (*_dL[r])[iamp]+=(*_dLN[r])[iamp]*n2;
      // ---------------------------------------
      // sort results for derivatives into output array and cache:
      unsigned int npar=ampi.npar();
      for(unsigned int ipar=0;ipar<npar;++ipar){
	// dL/dpar = dL/dRE*dR/dpar + dL/dIM*dIM/dpar
	std::complex<double> deriv=ampi.dampdpar(ipar);
	_dLcache[outcount]=(*_dL[r])[iamp].real()*deriv.real()+(*_dL[r])[iamp].imag()*deriv.imag();
	df[outcount]=_dLcache[outcount];
	++outcount;
      } // end loop over parameters (given by possibly contraint amp)

      //cout << "df(r="<<r<<",i="<<iamp<<")="<<dL[_offset[r]+iamp]<<endl;
    }// end outer loop over waves ******************************************
  } // end loop over rank ****************************************************
  // take care of derivative for flat:
  _dLcache[outcount]=dLdflat+n2*_Vflat.amp().real();
  df[outcount]=dLdflat+n2*_Vflat.amp().real();

  
  N.real()+=std::norm(_Vflat.amp()); 

  double t2=timer.RealTime();
  timer.Stop();

  _Ltime+=t1;
  _Ntime+=t2;
  
  if(!_quiet){
    cout << "LikelihoodC: "<<L<<endl;
    cout << "Normalization: "<< N.real() << endl;
    cout << "Normalized: "<<L + nevt*N.real() << endl;
    cout << "Time for LikelihoodC: "<< t1 <<endl;
    cout << "Time for Normalization: "<< t2 <<endl;
  }
  f=L + nevt*N.real(); // return value of likelihood


  


  return;
}


//************************************************************************

double 
TPWALikelihoodC::DoEval(const double* x) const
{

  // call FdF
  double df[_dim];
  double L;
  FdF(x, L,df);
  return L;
  /*
 
  // build complex numbers from parameters
  // remember rank restrictions!
  ++(*const_cast<unsigned int*>(&_ncalls));

  // check timing
  TStopwatch timer;
  timer.Start(true);


  vector<complex<double> > V(_rank*_nwaves);
  double re,im;
  unsigned int k=0;
  for(unsigned int r=0;r<_rank;++r){
    //cout<<"Rank "<<r<<"    offset="<<_offset[r]<<endl;
    for(unsigned int i=0;i<_nwaves;++i){
      if(i<r){re=0;im=0;}
      else if(i==r){re=x[k++];im=0;} // real parameter
      else {
	re=x[k++];
	im=x[k++];
      }
      V[_offset[r]+i]=complex<double>(re,im);
      //cout<<"Wave"<<i<<"="<<V[_offset[r]+i]<<endl;
    }
  } // end loop over rank
  double flat=x[k];
  double flat2=flat*flat;

  double L=0;
  // loop over events and calculate log likelyhood
  for(unsigned int ievt=0;ievt<_nevents;++ievt){
    double l=0;
    // incoherent sum over ranks: loop
    for(unsigned int r=0;r<_rank;++r){
      complex<double> A(0,0);
      for(unsigned int iamp=0;iamp<_nwaves;++iamp){
	A+=(V[_offset[r]+iamp]*((*_data[iamp])[ievt]));
      } // end loop over waves
      //cout << "A2=" << std::norm(A) << endl;
      l+=std::norm(A);
      assert(l>=0);
      //cout << "l="<< l << endl;
    }// end loop over rank
    l+=flat2;
    L-=TMath::Log(l);
    //cout << "L="<< L << endl;
  }// end loop over events
  
   double t1=timer.RealTime();
  timer.Start(true);
  // add normalization integrals ?? TODO: can we exploit symmetry here?
  complex<double> N(0,0);
  for(unsigned int r=0; r<_rank;++r){
    //cout <<"/n Rank"<<r<<endl;
    for(unsigned int iamp=r;iamp<_nwaves;++iamp){
      for(unsigned int jamp=r;jamp<_nwaves;++jamp){
	complex<double> p=V[_offset[r]+iamp]*std::conj(V[_offset[r]+jamp]);
	//cout<<_wavenames[iamp]<<"   "<<_wavenames[jamp]<<endl;
	complex<double> I=(const_cast<intmat*>(&_mat)->el(_wmap[iamp],_wmap[jamp]));
	//cout<<"I="<<I<<endl;
	//cout<<"p="<<p<<endl;
	p*=I;
	N+=p;
      }// end inner loop over waves
    }// end outer loop over waves
  } // end loop over rank
  N.real()+=flat2; 

  double t2=timer.RealTime();
  timer.Stop();

  (*const_cast<double*>(&_Ltime))+=t1;
  (*const_cast<double*>(&_Ntime))+=t2;

  cout << "LikelihoodC: "<<L<<endl;
  cout << "Normalization: "<< N.real() << endl;
  cout << "Normalized: "<<L + (double)_nevents*N.real() << endl;
  cout << "Time for LikelihoodC: "<< t1 <<endl;
  cout << "Time for Normalization: "<< t2 <<endl;

  return L + (double)_nevents*N.real();
  */
}

double 
TPWALikelihoodC::DoDerivative(const double* x, unsigned int idf) const
{
  //cout << " Calling Derivative !!!!!!!!!!! " << endl;

  // check if we have cached this
  bool same_x=true;
  for(unsigned int i=0; i<_dim; ++i){
    same_x &= (_parcache[i]==x[i]);
  }
  if(same_x){
    //cout << "using cached derivative! " << endl;
    return _dLcache[idf];
  }

  // call FdF
  double L;
  double df[_dim];
  FdF(x, L,df);
 
  return df[idf];
}

unsigned int 
TPWALikelihoodC::NDim() const {
  return _dim;
}

// Complex valued Amplitudes and
// mapping of real and imaginary part of amplitudes
// in error matrix (needed by TFitBin)
// takes into account real-valued parameters
void 
TPWALikelihoodC::buildCAmps(const double* x, vector<complex<double> >& V, 
			    vector<pair<int,int> >& indices, 
			    vector<string>& names, 
			    const TMatrixD& errpar,
			    TMatrixD& erramp,
			    bool withFlat){
  // build complex numbers from parameters
  // remember rank restrictions!
  V.clear();
  indices.clear();
  names.clear();
  erramp.Clear();
  //unsigned int namp=_rank*_nwaves;
  //if(withFlat)namp+=1;


  partoamp(x);
  unsigned int m=2*(_nwaves*_rank);
  if(withFlat)m+=2;
  erramp.ResizeTo(m,m);
  // transformation matrix from parameters to amps
  // (linearized version of partoamp)
  TMatrixD T(m,NDim());


  unsigned int parcount=0;
  for(unsigned int r=0;r<_rank;++r){
    //cout<<"Rank "<<r<<"    offset="<<_offset[r]<<endl;
    for(unsigned int i=0;i<_nwaves;++i){
      string wavename=(*_V[r])[i].name();
      V.push_back((*_V[r])[i].amp());
      stringstream name;
      name << "V" << r << "_" << wavename;
      names.push_back(name.str());
      // now produce error matrix for amplitudes
      // loop over amplitudes to get all the entries for this row
      unsigned int npar=(*_V[r])[i].npar();
      for(unsigned int ipar=0;ipar<npar;++ipar){
	complex<double> dAdp=(*_V[r])[i].dampdpar(ipar);
	T[r*2*_nwaves+2*i][parcount]=dAdp.real(); // d(real)/dpar
	T[r*2*_nwaves+2*i+1][parcount]=dAdp.imag(); // d(im)/dpar
	++parcount;
      }
      indices.push_back(pair<int,int>(r*2*_nwaves+2*i,r*2*_nwaves+2*i+1));
    } // end loop over waves
  } // end loop over rank

  if(withFlat){
    V.push_back(complex<double>(x[parcount],0));
    T[_rank*2*_nwaves][parcount]=1;
    T[_rank*2*_nwaves+1][parcount]=0;
    indices.push_back(pair<int,int>(_rank*2*_nwaves,_rank*2*_nwaves+1));
    names.push_back("V_flat");
  }

  TMatrixD TT(TMatrixD::kTransposed,T);
  TMatrixD CTT=errpar*TT;
  erramp=T*CTT;

}



// depends on naming convention for waves!!!
// VR_IGJPCMEIso....
int
TPWALikelihoodC::geteps(TString wavename){
  int eps=0;
  unsigned int i=6;
  if(wavename[0] == 'V') {
    i=9; // check if we have a parameter name
    //cout << "geteps from parameter name" << endl;
  }
  if(wavename[i] == '-') eps= -1;
  else if(wavename[i] == '+') eps= 1;
  else {
    cout<<"Wrong wavename format. "
	<<"Could not determine reflectivity! Aborting!"<<endl;
    throw;
  }
  //cout << "eps=" << eps << endl;
  return eps;
}

void
TPWALikelihoodC::getIntCMatrix(TCMatrix& integr, TCMatrix& acceptance){
  //integr.ResizeTo(_nwaves,_nwaves
   for(unsigned int i=0;i<_nwaves;++i){// begin outer loop
     for(unsigned int j=0;j<_nwaves;++j){
        integr.set(i,j,_mat.el((*_V[0])[i].normindex(),(*_V[0])[j].normindex()));
	acceptance.set(i,j,_mat.el((*_V[0])[i].accindex(),(*_V[0])[j].accindex()));
     }// end inner loop
   } // end outere loop over parameters;
   // add flat
   integr.set(_nwaves,_nwaves,1);
   acceptance.set(_nwaves,_nwaves,1);
}

// A simple factory for creating contraint objects
// Supports: TPWAPhaseContsraint
// config file format:
// 
// <ConstraintType> {
// ConstraintWave
// Parameter2
// ...
// ParameterN
// }
//
// For example:
// Phase {
// ConstraintWave
// Masterwave
// RelativePhase
// }

unsigned int
TPWALikelihoodC::addConstraints(TString ConstraintsFile){
  ifstream cfile(ConstraintsFile.Data());
  unsigned int ccounter=0;
  char line[256];
  if(!cfile.good()){
    cout << "Invalid constraint file " << ConstraintsFile << endl;
    return 0;
  }
  
  while(cfile.good()){
    cfile.getline(line,256);
    // check if we have a valid token
    if(strchr(line,'{')!=NULL){
      // process constraint:
      // which contraint?
      if(strstr(line,"Phase")){
	// get 3 paramters:
	cfile.getline(line,256);
	string slavewave(line);
	cfile.getline(line,256);
	string masterwave(line);
	cfile.getline(line,256);
	double phase=atof(line);
	// check if master and slave are in our waveset
	unsigned int imaster=_nwaves+1;
	unsigned int islave=_nwaves+1;
	for(unsigned int iamp=0;iamp<_nwaves; ++iamp){
	  if((*_V[0])[iamp].name()==masterwave)imaster=iamp;
	  else if((*_V[0])[iamp].name()==slavewave)islave=iamp;
	} // end loop over waves
	if(imaster>=_nwaves){
	  cout << "Master wave " << masterwave << " not in waveset. Skipping phase constraint of " << slavewave << endl;
	}
	if(islave>=_nwaves){
	  cout << "Constraint wave " << slavewave << " not in waveset. Skipping phase constraint to " << masterwave << endl;
	}
	else { // checks succeed
	  if(!_quiet) cout << "Constraining " << slavewave 
			   << " to " << masterwave 
			   << " with relative phase dPhi=" << phase << endl;
	  // INSTALL CONSTRAINT IN ALL RENKS
	  for(unsigned int ir=0; ir<_rank; ++ir){
	    (*_V[ir])[islave].setConstraint(new TPWAPhaseConstraint(phase,&(*_V[0])[imaster]));
	    ++ccounter;
	  } // end loop over rank
	} // end checks if waves are found soucceed 
      } // end PHASE CONSTRAINT --------------------------------------
      else {
	cout << "Unknown constraint" << line << " Aborting." << endl;
	throw;
      } // end Unknown constraint -------------------------------------
      cfile.getline(line,256); // get trailing line.
    } // end constraint block
  } // end loop through constraint config file
  cfile.close();
  if(!_quiet)cout << "Installed "<<ccounter<< " contraints."<<endl;
  return ccounter;
}
