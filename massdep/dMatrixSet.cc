//-----------------------------------------------------------
//
// Description:
//      Implementation of class pwacomponent
//      see dMatrixSet.h for details
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

// This Class' Header ------------------
#include "dMatrixSet.h"


// C/C++ Headers ----------------------
#include <algorithm>
#include <iostream>

// Collaborating Class Headers --------
#include "TF1.h"

// Class Member definitions -----------

using namespace std;



vector<string> 
dMatrixSet::wavelist() const {
   vector<string> wl;
//   for(unsigned int i=0;i<n();++i){
//     const map<string,pwachannel >& channellist=_comp[i]->channels();
//     map<string,pwachannel >::const_iterator it=channellist.begin();
//     while(it!=channellist.end()){
//       if(find(wl.begin(),wl.end(),it->first)==wl.end())
// 	 wl.push_back(it->first);
//       ++it;
//     }
//   }
  return wl;
}


void
dMatrixSet::setPS(TF1* fPS){
  _phasespace=fPS;
  // check if there are free parameters in the phase space that should be fitted
  unsigned int nparPS=_phasespace->GetNpar();
  // loop over parameters and check limits
  // remember which parameters to let float
  _freePSpar.clear();
  for(unsigned int i=0;i<nparPS;++i){
    double min,max;
    _phasespace->GetParLimits(i,min,max);
    if(min!=max){
      _freePSpar.push_back(i);
      cout << "PS parameter "<< i << " floating in ["
	   << min  << "," << max << "]" << endl;
    }
  }// end loop over parameters
  _numpar+=_freePSpar.size();
}

double 
dMatrixSet::getFreePSPar(unsigned int i){
  if(i<_freePSpar.size())
    return _phasespace->GetParameter(_freePSpar[i]);
  else return 0;
  return 0;
}


void 
dMatrixSet::getFreePSLimits(unsigned int i, double& lower, double& upper){
  if(i<_freePSpar.size()){
    _phasespace->GetParLimits(_freePSpar[i],lower,upper);
  }
}

void
dMatrixSet::setPar(const double* par){ // set parameters
  // unsigned int parcount=0;
//   // components
//   for(unsigned int i=0;i<n();++i){
//     _comp[i]->setPar(par[parcount],par[parcount+1]);
//     parcount+=2;
//     _comp[i]->setCouplings(&par[parcount]);
//     parcount+=_comp[i]->numChannels()*2; // RE and Im for each channel
//   } // end loop over components
//   // phase space
//   unsigned int nfreepar=_freePSpar.size();
//   for(unsigned int ipar=0;ipar<nfreepar;++ipar){
//     _phasespace->SetParameter(_freePSpar[ipar],par[parcount]);
//     ++parcount;
//   }
}


void 
dMatrixSet::getPar(double* par){       // return parameters 
  // unsigned int parcount=0;
//   // components
//   for(unsigned int i=0;i<n();++i){
//     par[parcount]=_comp[i]->m0();
//     par[parcount+1]=_comp[i]->gamma();
//     parcount+=2;
//     _comp[i]->getCouplings(&par[parcount]);
//     parcount+=_comp[i]->numChannels()*2; // RE and Im for each channel
//   }
//  // phase space
//   unsigned int nfreepar=_freePSpar.size();
//   for(unsigned int ipar=0;ipar<nfreepar;++ipar){
//     par[parcount]=_phasespace->GetParameter(_freePSpar[ipar]);
//     ++parcount;
//   }
}

unsigned int 
dMatrixSet::numPar() const {
  return 0;
}

double 
dMatrixSet::ps(double m){return 0;}//_phasespace->Eval(m);}

double 
dMatrixSet::intensity(const std::string& wave, double m){
  // loop over all components and pick up those that contribute to this channel
  // complex<double> rho(0,0);
//   for(unsigned int ic=0;ic<n();++ic){
//     if(_comp[ic]->channels().count(wave)==0)continue;
//     else {
//       rho+=_comp[ic]->val(m)*_comp[ic]->channels().find(wave)->second.C()*sqrt(_comp[ic]->channels().find(wave)->second.ps(m));
//     }

//   }
//   return norm(rho)*_phasespace->Eval(m);
  return 0;
}

double 
dMatrixSet::phase(const std::string& wave, double m){
  // loop over all components and pick up those that contribute to this channel
 //  complex<double> rho(0,0);
//   for(unsigned int ic=0;ic<n();++ic){
//     if(_comp[ic]->channels().count(wave)==0)continue;
//     else {
//       rho+=_comp[ic]->val(m)*_comp[ic]->channels().find(wave)->second.C();
//     }

//   }
//   return arg(rho);
  return 0;
}


double 
dMatrixSet::phase(const std::string& wave1,
			const std::string& wave2,
			double m){
  return arg(overlap(wave1,wave2,m));
}

std::complex<double>
dMatrixSet::overlap(const std::string& wave1,
			  const std::string& wave2,
			  double m){
//     // loop over all components and pick up those that contribute to this channel
   complex<double> rho1(0,0);
//   complex<double> rho2(0,0);

//   for(unsigned int ic=0;ic<n();++ic){
//     if(_comp[ic]->channels().count(wave1)!=0){
//       rho1+=_comp[ic]->val(m)*_comp[ic]->channels().find(wave1)->second.C()*sqrt(_comp[ic]->channels().find(wave1)->second.ps(m));
//     }
//     if(_comp[ic]->channels().count(wave2)!=0){
//       rho2+=_comp[ic]->val(m)*_comp[ic]->channels().find(wave2)->second.C()*sqrt(_comp[ic]->channels().find(wave2)->second.ps(m));
//     }
//   }
//   return rho1*conj(rho2)*_phasespace->Eval(m);
   return rho1;
}


std::ostream& operator<< (std::ostream& out,const dMatrixSet& cs){
 
  return out;
}

///////////////////////////////////////////////////////////////////////////

qntag::qntag(const TString& tag){
  _I=atoi(tag(0,1).Data());
  _G= tag[1]=='+' ? +1 : -1;
  _J=atoi(tag(2,1).Data());
  _P= tag[2]=='+' ? +1 : -1;
  _C= tag[3]=='+' ? +1 : -1;
}

bool operator== (const qntag& lhs,const qntag& rhs){
  return (lhs._I==rhs._I) && (lhs._G==rhs._G) &&(lhs._J==rhs._J) &&(lhs._P==rhs._P) &&(lhs._C==rhs._C);
}

bool operator< (const qntag& lhs,const qntag& rhs){
  if(lhs._I==rhs._I){
    if(lhs._G==rhs._G){
      if(lhs._J==rhs._J){
	if(lhs._P==rhs._P){
	  if(lhs._C==rhs._C)return false;
	  else return (lhs._C<rhs._C);}
	else return (lhs._P<rhs._P);}
      else return (lhs._J<rhs._J);}
    else return (lhs._G<rhs._G);}
  else return (lhs._I<rhs._I);
}

std::ostream& operator<< (std::ostream& o,const qntag& tag){
  o << tag._I << tag._G << tag._J << tag._P << tag._C;
  return o;
}
