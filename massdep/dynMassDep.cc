//-----------------------------------------------------------
//
// Description:
//      Implementation of class dynMassDep
//      see dynMassDep.h for details
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------



// This Class' Header ------------------
#include "dynMassDep.h"

// C/C++ Headers ----------------------
#include <vector>
#include "TMath.h"

// Collaborating Class Headers --------
#include "cauchyIntegral.h"

// Class Member definitions -----------

using namespace std;


rpwa::dynMassDep::dynMassDep(double M, double width,
			     unsigned int nparticles, double* masses) 
  : mS(M*M), mM(M), mWidth(width){
  ps=new mcPhaseSpace(nparticles,masses,0,40,400,50000);
}



void
rpwa::dynMassDep::addDecayChannel(absDecayChannel* ch){
  _channels.push_back(ch);
  ps->addDecayChannel(ch);
}


rpwa::cd 
rpwa::dynMassDep::val_bnl(double m, unsigned int i){

  const double m0     = mM;
  const double Gamma0 = mWidth;
  
  // test rho rho
  const double s1=0.77549*.77549;;
  const double s=m*m;
  const double s0=m0*m0;
  
  const double lam=s*s+2*s1*s1 - 2*(s*s1+s1*s1+s1*s);
  const double lam0=s0*s0+2*s1*s1 - 2*(s0*s1+s1*s1+s1*s0);
  

  const double q      = sqrt(fabs(lam)/4/s);
  const double q0     = sqrt(fabs(lam0)/4/s0);
 
  cerr << "m ="<<m 
       << "   q ="<< q 
       << "   q0="<< q0 << endl; 

  //const int    l      = 0;

  const double    GammaV = Gamma0 * (m0 / m) * (q / q0);
  complex<double> ret    = (m0 * Gamma0) / (m0 * m0 - m * m - complex<double>(0, 1) * m0 * GammaV);
  return ret;


}


rpwa::cd
rpwa::dynMassDep::val(double m,unsigned int i){
  if(m<ps->thres()){
     cout << "Requested m=" << m << " which is below threshold "<< ps->thres() << endl;
    //throw;
    return 0;
  }
  double br=1;
  if(_channels.size()>0)br=_channels[i]->branching();
  double lambda=br*ps->eval(m,i)/get_rho0(i);
  double s=m*m;
  double ms=0;
  double r=0;
  for(unsigned int ich=0;ich<ps->nChannels();++ich){
    ms+=get_ms(s,ich);
    double br=1;
    if(_channels.size()>0)br=_channels[i]->branching();
    r+=br*ps->eval(m,ich)/get_rho0(ich);
  }
  cd N(mWidth*mM*lambda,0);
  cd D(mS-s-ms,-mM*mWidth*r);
  return N/D;
}

rpwa::cd
rpwa::dynMassDep::val_static(double m,unsigned int i){
 if(m<ps->thres()){
     cout << "Requested m=" << m << " which is below threshold "<< ps->thres() << endl;
    //throw;
    return 0;
  }
  double s=m*m;
  //double ms=get_ms(s);
  cd N(mWidth*mM,0);
  cd D(mS-s,-mM*mWidth);
  return N/D;
}


rpwa::cd
rpwa::dynMassDep::val_nodisperse(double m, unsigned int i){
 if(m<ps->thres()){
     cout << "Requested m=" << m << " which is below threshold "<< ps->thres() << endl;
    //throw;
    return 0;
  }
double br=1;
 if(_channels.size()>0)br=_channels[i]->branching();
 double lambda=br*ps->eval(m,i)/get_rho0(i);
 double s=m*m;
 double r=0;
 for(unsigned int ich=0;ich<ps->nChannels();++ich){
   double br=1;
   if(_channels.size()>0)br=_channels[i]->branching();
    r+=br*ps->eval(m,ich)/get_rho0(ich);
  }
 cd N(mWidth*mM*lambda,0);
 cd D(mS-s,-mM*mWidth*r);
  return N/D;
}


double
rpwa::dynMassDep::get_rho0(unsigned int i)const{
  if(!ps->hasCached())throw;
  return ps->eval(mM,i);
}

double
rpwa::dynMassDep::disperse(double* x, double* par) const{
  // x[0] is s' (which will be integrated over in the cauchyIntegral)
  // par[0] is s (which the dispersion integral depends on)
  // par[1] is the channel;
  unsigned int i=(unsigned int)par[1];
  //cout << "Calling disperse" << endl;
  double nom=mM*mWidth*get_rho(sqrt(x[0]),i)/get_rho0(i);
  double denom=(x[0]-par[0])*(x[0]-mS);
  return nom/denom;

}

double
rpwa::dynMassDep::calc_ms(double s, unsigned int i) const {
  if(s==mS)return 0;
  double cutoffs=1600;
  TF1* f=new TF1("ms",this,&dynMassDep::disperse,ps->thres()*ps->thres(),cutoffs,2,"dynMassDep","disperse");
  f->SetParameter(0,s);
  f->SetParameter(1,i);
  vector<realPole> p;
  double residual_s=mM*mWidth*get_rho(sqrt(s),i)/get_rho0(i)/(s-mS);
  double residual_mS=mM*mWidth*get_rho(sqrt(mS),i)/get_rho0(i)/(mS-s);
  p.push_back(realPole(s,residual_s));
  p.push_back(realPole(mS,residual_mS));
  double low=ps->thres();low*=low;
  cauchyIntegral cauchy(f,p,low,cutoffs);
  double I=cauchy.eval_Hunter(4);
  //cout << "I= " << setprecision(12);
  //cout << I << endl; 
  delete f;
  double br=1;
  if(_channels.size()>0)br=_channels[i]->branching();
  return br*(s-mS)/TMath::Pi()*I;
}


void
rpwa::dynMassDep::store_ms(double maxM, unsigned int steps){
  // loop over decay channels
  unsigned int nch=ps->nChannels();
  for(unsigned int ich=0;ich<nch;++ich){
    _graphms.push_back(new TGraph(steps));
    double m0=ps->thres();
    double step=(maxM-m0)/(double)steps;
    for(unsigned int i=0;i<steps;++i){
      double x=m0+((double)i+0.5)*step;
      double s=x*x;
      double valms=calc_ms(s,ich);
      if(valms!=valms){
	cerr<<"Nan in ms!"<<endl;
	cerr<<"m="<<x<<endl;
	throw;
      }
      _graphms[ich]->SetPoint(i,s,valms);
    }
  }// end loopmover decay channels
}
  
double
rpwa::dynMassDep::get_ms(double s, unsigned int i) const{
  if(_graphms.size()==0)throw;
  else return _graphms[i]->Eval(s);
}
