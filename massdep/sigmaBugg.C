// David Buggs sigma parameterization
// From: J.Phys.G34 (2007) 151

#include <complex>
#include <vector>

#include <iostream>
#include <iomanip>
#include "TMath.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "cauchyIntegral.h"
#include "Math/IFunction.h"
#include "Math/GSLIntegrator.h"

//#include "SpecFuncMathMore.h"

//using namespace ROOT::Math;
using namespace std;

double const mpi=0.13957018;
double const mpi2=mpi*mpi;
double const mK=0.493677;
double const mK2=mK*mK;
double const mEta=0.547853;
double const mEta2=mEta*mEta;

double const pi=TMath::Pi();

// sigma mass
double const M=0.900;//1.3150;
double const M2=M*M;

double const b1=3.728;//1.302;
double const b2=0.092;//0.340;
double const A=2.882;//2.426;
double const sa=0.41*mpi2;
double const gKK=0.6;
double const gEtaEta=0.2;
double const g4=0.012;
double const Lambda=3.39;
double const s0=3.238;


typedef complex<double> comp;

comp const ic(0,1);

// phase space factor analytically continued below threshold
comp rho(double s, double m2){
  if(s>4*m2)return sqrt(1-4*m2/s);
  else return ic*sqrt(4*m2/s-1);
}

comp dampedRho(double s, double m2){
  double alpha=5.2;
  if(s-4*m2<0)alpha=8.4;
  return exp(-alpha*fabs(s-4*m2));
}

// 4pi rho parametr
comp rho4pi(double s){
  if(s<16*mpi2)return 0;
  //return sqrt(1-16*mpi2/s)/(1+exp(Lambda*(s0-s)));//
  return sqrt(1-16*mpi2/s)/(1+exp(Lambda*(s0-s)))*exp(-2*(s-1.45*1.45));
}


// empirical factor from BES analysis
comp g1(double s){
  return M*(b1+b2*s)*exp(-(s-M2)/A);
}

comp adler(double s){
  return g1(s)*(s-sa)/(M2-sa);
}

// pipi width
comp gamma1(double s){
  return adler(s)*rho(s,mpi2);
}



// KK
comp gamma2(double s){
  return gKK*g1(s)*dampedRho(s,mK2)*rho(s,mK2);
}

// etaeta
comp gamma3(double s){
  return gEtaEta*g1(s)*dampedRho(s,mEta2)*rho(s,mEta2);
}

// 4pi
comp gamma4(double s){
  if(s<16*mpi2)return 0;
  else return M*g4*rho4pi(s)/rho4pi(M2);
}

comp gammaTot(double s){
  return gamma1(s)+gamma2(s)+gamma3(s)+gamma4(s);
}


// dispersion term (only needed for pipi -> r=real
comp myj1(double s){
  double r=rho(s,mpi2).real();
  return 1./pi*(2 + r*TMath::Log((1-r)/(1+r)));
}

comp z1(double s){
  return myj1(s)-myj1(M2);
}


comp msParam(double s){
  comp result;
  double a1=17.051;
  double s1=3.0533;double ds1=(s-s1)*(s-s1);
  double w1=1.0448;double w12=w1*w1;
  double dM1=(M2-s1)*(M2-s1);
  double a2=-0.0536;
  double s2=-0.0975;double ds2=(s-s2)*(s-s2);
  double w2=0.2801;double w22=w2*w2;
  double dM2=(M2-s2)*(M2-s2);

 //  double a1=3.5320;
//   double s1=2.9876;double ds1=(s-s1)*(s-s1);
//   double w1=0.8804;double w12=w1*w1;
//   double dM1=(M2-s1)*(M2-s1);
//   double a2=-0.0427;
//   double s2=-0.4619;double ds2=(s-s2)*(s-s2);
//   double w2=-0.0036;double w22=w2*w2;
//   double dM2=(M2-s2)*(M2-s2);
  result=a1/(ds1+w12)-a1/(dM1+w12)+a2/(ds2+w22)-a2/(dM2+w22);


  return result;

}

// dispersion function for pipi
double disperse(double* x, double* par){
  // x[0] is s' (which will be integrated over in the cauchyIntegral)
  // par[0] is s (which the dispersion integral depends on)
  // par[1] is M2 (pole mass of sigma squared)
  //cout << "Calling disperse" << endl;
  double nom=gamma4(x[0]).real();
  double denom=(x[0]-par[0])*(x[0]-par[1]);
  return nom/denom;
}



// the function may contain only one pole, which has to ly outside of integration region. The other pole will be added by the principal value integrator. For normal integration both poles can be added
class dispInt: public ROOT::Math::IBaseFunctionOneDim {
public:
  dispInt(double pole):_pole(pole),_pole2(0){};
  double DoEval(double x) const;
  ROOT::Math::IBaseFunctionOneDim* Clone() const {return new dispInt(*this);}
  
  void setPole(double s){_pole=s;_pole2=0;}
  void setPoles(double s, double s2){_pole=s;_pole2=s2;}
  

private:
  double _pole;
  double _pole2;

};

double 
dispInt::DoEval(double x) const{
   double nom=gamma4(x).real();
  double denom=(x-_pole);
  if(_pole2!=0){denom*=(x-_pole2);

    //cout << "2poles" << endl;
  }
  return nom/denom;
}

// even component of dispersion integrand
// Used for Longman method of integration
class dispIntEven: public ROOT::Math::IBaseFunctionOneDim {
public:
  dispIntEven(double pole1, double pole2, double center):_pole(pole1),_pole2(pole2),_center(center){};
  double DoEval(double x) const;
  ROOT::Math::IBaseFunctionOneDim* Clone() const {return new dispIntEven(*this);}
  
  // center should coincide with one of the poles
  void setPoles(double s, double s2, double center){_pole=s;_pole2=s2;_center=center;}
  
private:
  double _pole;
  double _pole2;
  double _center;

};

double 
dispIntEven::DoEval(double x) const{
  double xplus=_center+x;
  double xminus=_center-x;
  double nom1=gamma4(xplus).real();
  double nom2=gamma4(xminus).real();  
 
  double denom1=(xplus-_pole)*(xplus-_pole2);
  double denom2=(xminus-_pole)*(xminus-_pole2);

  return 2*(nom1/denom1+nom2/denom2);
}




comp ms(double s, double cutoff=10){
  if(s==M2)return 0;
  double W=1.315; // mass
  double W2=W*W;
  double cutoffs=s*2;//cutoff;
  if(s<W2)cutoffs=W2*2;
  TF1* f=new TF1("ms",&disperse,16*mpi2,cutoff,2,"disperse");
  f->SetParameter(0,s);
  f->SetParameter(1,W2);
  vector<realPole> p;
  double residual_s=gamma4(s).real()/(s-W2);
  double residual_mS=gamma4(W2).real()/(W2-s);
  double low=16*mpi2;
  if(s>low)p.push_back(realPole(s,residual_s)); // this pole only exists above threshold!
  p.push_back(realPole(W2,residual_mS));

  cauchyIntegral cauchy(f,p,low,cutoffs);
  double I=cauchy.eval_Hunter(4);
  // calculate the rest of the integral cutoffs-> inifinty

  double I1=f->Integral(cutoffs,cutoff);
  //cout << cutoff << "   " << I1 << endl;


  delete f;
  
  
  return (s-W2)/TMath::Pi()*(I+I1);
}


comp msGSL(double s){
  ROOT::Math::GSLIntegrator Integrator;
 dispInt d(0);
 
 double mym2=0.9*0.9;
  
 // split integral in three parts containting at most one pole each
 double mid=0.5*(s+mym2);
 double top= s > mym2 ? s : mym2;
 double bottom= s > mym2 ? mym2 : s;
 
 d.setPole(top);
 Integrator.SetFunction(d);
 double val1=Integrator.IntegralCauchy(16*mpi2,mid,bottom);
 d.setPole(bottom);
 Integrator.SetFunction(d);
 double val2=Integrator.IntegralCauchy(d,mid,top+1,top);
 d.setPoles(bottom,top);
 Integrator.SetFunction(d);
 double val3=Integrator.IntegralUp(d,top+1);
 double msgsl=(s-mym2)/TMath::Pi()*(val1+val2+val3);
 return msgsl;


}

// using Longmans decomposition method
comp msLong(double s){
  ROOT::Math::GSLIntegrator Integrator(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);

 double mym2=0.9*0.9;
  
 // split integral into 5 parts containting at most one pole each
 // such that the pole lies in the center of the interval

 double mid=0.5*(s+mym2);
 double diff1=0.2*fabs(s-mym2);
 double top= s > mym2 ? s : mym2;
 double bottom= s > mym2 ? mym2 : s;
 double diff2=fabs(bottom-4*mpi2);
 double diff=diff1 > diff2 ? diff2 : diff1;
 // intervals:
 // [16mpi2, bottom-diff] no pole
 // [bottom-diff, bottom+diff] bottom pole
 // [bottom+diff, top-diff] no pole  
 // [top -diff, top+diff] top pole
 // [top+diff, infty] no pole
 dispInt d(0);
 dispIntEven dE(0,0,0);

 // [16mpi2, bottom-diff] no pole
 d.setPoles(bottom,top);
 Integrator.SetFunction(d);
 double val1=Integrator.Integral(16*mpi2,bottom-diff);
 // [bottom-diff, bottom+diff] bottom pole
 dE.setPoles(bottom,top,bottom);
 Integrator.SetFunction(dE);
 double val2=2*Integrator.Integral(0,+diff);
 // [bottom+diff, top-diff] no pole  
 Integrator.SetFunction(d);
 double val3=Integrator.Integral(bottom+diff,top-diff);
 // [top -diff, top+diff] top pole
 dE.setPoles(bottom,top,top);
 Integrator.SetFunction(dE);
 double val4=2*Integrator.Integral(0,+diff);
// [top+diff, infty] no pole
 Integrator.SetFunction(d);
 double val5=Integrator.IntegralUp(top+diff);

 double msgsl=(s-mym2)/TMath::Pi()*(val1+val2+val3+val4+val5);
 return msgsl;


}




comp D(double s){
  
  return M2 - s - adler(s)*z1(s)  - ic*gammaTot(s) - msGSL(s); 


}


comp tscat(double m){
  double s=m*m;
  return M*gamma1(s)/D(s); 

}

comp tprod(double m){
  double s=m*m;
  return 1./D(s); 

}


void sigmaBugg(){


 unsigned int npoints=240;
  TGraph* gRe=new TGraph(npoints);
  TGraph* gIm=new TGraph(npoints);
  TGraph* gIntense=new TGraph(npoints);
  TGraph* gPhase=new TGraph(npoints);
  TGraph* g4piPS=new TGraph(npoints);
  TGraph* gMS=new TGraph(npoints);
  TGraph* gMSGSL=new TGraph(npoints);
  TGraph* gMSLong=new TGraph(npoints);
 
  std::vector<TGraph*> vgMSCut;
  unsigned int ng=80;
  TMultiGraph* multi=new TMultiGraph();
  for(unsigned int ig=0;ig<ng;++ig){
    vgMSCut.push_back(new TGraph(npoints));
    multi->Add(vgMSCut[ig]);
  }
  TGraph* gMSParam=new TGraph(npoints);

  gMS->SetTitle("Dispersive Term m(s)");
  gRe->SetTitle("Scattering");
  gIm->SetTitle("Scattering Imag");
  gIntense->SetTitle("Scattering Intensity");
  gPhase->SetTitle("Scattering Phase");

  TGraph* gReProd=new TGraph(npoints);
  TGraph* gImProd=new TGraph(npoints);
  TGraph* gIntenseProd=new TGraph(npoints);
  TGraph* gPhaseProd=new TGraph(npoints);
  gReProd->SetTitle("Production");
  gImProd->SetTitle("Production Imag");
  gIntenseProd->SetTitle("Production Intensity");
  gPhaseProd->SetTitle("Production Phase");

  double step=0.02;
  double Mass=2*mpi+step;
  for(unsigned int i=0;i<npoints; ++i){
    complex<double> amp=tscat(Mass); 
    complex<double> pamp=tprod(Mass);
    double p=sqrt(Mass*Mass/4.-mpi2);
   
    g4piPS->SetPoint(i,Mass,rho4pi(Mass*Mass).real());
    //gMS->SetPoint(i,Mass,ms(Mass*Mass,100).real());
    gMSGSL->SetPoint(i,Mass,msGSL(Mass*Mass).real());
    gMSLong->SetPoint(i,Mass,msLong(Mass*Mass).real());
    gMSParam->SetPoint(i,Mass,msParam(Mass*Mass).real());

 //    for(unsigned int j=0;j<ng;++j){
//       double cut=8+(double)j*5;
//       vgMSCut[j]->SetPoint(i,Mass,ms(Mass*Mass,cut).real());
//     }


    gRe->SetPoint(i,Mass,amp.real()*Mass/(2*p));
    gIm->SetPoint(i,Mass,amp.imag()*Mass/(2*p));
    gIntense->SetPoint(i,Mass,norm(amp));
    gPhase->SetPoint(i,Mass,arg(amp));

    gReProd->SetPoint(i,Mass,pamp.real());
    gImProd->SetPoint(i,Mass,pamp.imag());
    gIntenseProd->SetPoint(i,Mass,norm(pamp));
    gPhaseProd->SetPoint(i,Mass,arg(pamp));

    Mass+=step;
  }
  cout << "finished" << endl;

  TCanvas* c=new TCanvas("c","c",10,10,1200,600);
  c->Divide(4,1);
  // c->cd(1);
  //gRe->Draw("APL");
  //gIm->Draw("PL same");
  //c->cd(2);
  //gIntense->Draw("APL");
  //c->cd(3);
  //gPhase->Draw("APL");
  
  c->cd(1);
  gReProd->Draw("APL");
  gImProd->Draw("PL same");
  c->cd(2);
  gIntenseProd->Draw("APL");
  //g4piPS->Draw("PL same");
   c->cd(3);
  gPhaseProd->Draw("APL");
  c->cd(4);
   gMSGSL->Draw("APL");
   gMSParam->SetLineColor(kRed);
   gMSParam->SetMarkerColor(kRed);
 gMSParam->Draw("same PL");
 g4piPS->SetLineColor(kBlue);
 g4piPS->SetMarkerColor(kBlue);
 g4piPS->Draw("same PL");
 gMSLong->SetLineColor(kMagenta);
 gMSLong->SetMarkerColor(kMagenta);
 gMSLong->Draw("same PL");




//   TF1* f=new TF1("f","0.5*exp(0.5*x)/cos(TMath::Pi()*0.5*x)",-2,2);
//  //TF1* f=new TF1("f","x*exp(-x)",0,20);
//   unsigned int n=400;
//  TGraph* gf=new TGraph(400);
//  //f->SetParameter(0,1);
//  vector<realPole> p;
//  p.push_back(realPole(-1,exp(-0.5)/TMath::Pi()));
//  for(unsigned int ip=1;ip<4;++ip){
//    cout << ip << "   ";
//    for(int sig=-1;sig<2;sig+=2){
   
//      double k=sig*(double)ip;
//      double res=pow(-1,k)*exp(k-0.5)/TMath::Pi();
//      if(fabs(k-0.5)<1){p.push_back(realPole(2*(k-0.5),res));
//        cout << sig << "   z=" <<k-0.5<<"   res="<< res << "   ";}
//    }
//    cout << endl;
//  }

 //p.push_back(realPole(4,36.63127778));
 
 //p.push_back(realPole(6,1));
 
 //  for(double k=0;k<=1;k+=1){
 //     double res=pow(-1.,k)*TMath::Exp(k-0.5)/TMath::Pi();
 //     p.push_back(realPole(k-0.5,res));
 
 //   }

 // cauchyIntegral cauchy(f,p,-2,2);
//  cout << "Calculating hunter: " << setprecision(12);
//  double I=cauchy.eval_Hunter(4); cout << I << endl;
 
 unsigned int n=200;
 double s=1.5*1.5;
 TGraph* gms=new TGraph(n);
 for(unsigned int i=0;i<n;++i){
   double m=(2+(double)i*0.01);
   double cut=m*m;
   gms->SetPoint(i,m,ms(s,cut).real());
 }

 TGraph* ggsl=new TGraph(n);
 ROOT::Math::GSLIntegrator Integrator;
 dispInt d(0);
 for(unsigned int i=0;i<n;++i){
   double m=(2*mpi+(double)i*0.01);
   double s=m*m;double mym2=1.35*1.35;
  
   double mid=0.5*(s+mym2);
   double top= s > mym2 ? s : mym2;
   double bottom= s > mym2 ? mym2 : s;
   
   cout << "s=" << s 
	<<"   bottom="<<bottom
	<<"   mid="<<mid
	<<"   top="<<top<<endl;

   d.setPole(top);
   Integrator.SetFunction(d);
   double val1=Integrator.IntegralCauchy(16*mpi2,mid,bottom);
   d.setPole(bottom);
   Integrator.SetFunction(d);
   double val2=Integrator.IntegralCauchy(d,mid,3*top,top);
   d.setPoles(bottom,top);
   Integrator.SetFunction(d);
   double val3=Integrator.IntegralUp(d,3*top);
   double msgsl=(s-M2)/TMath::Pi()*(val1+val2+val3);
   ggsl->SetPoint(i,m,msgsl);
 }



  TCanvas* c2=new TCanvas();
  //c2->Divide(1,2);
  //c2->cd(1);gms->Draw("APL");
  //c2->cd(2);gf->Draw("APL");
  //c2->cd(1);
  multi->Draw("APL");
  gMSParam->Draw("same");
  gMSGSL->SetLineColor(kBlue);
  gMSGSL->SetMarkerColor(kBlue);
  gMSGSL->Draw("same");
  
	       //c2->cd(2);ggsl->Draw("APL");

 //  double W=0.9128; // mass
//   double W2=W*W;
//   TF1* gf=new TF1("ms",&disperse,16*mpi2,100,2,"disperse");
//   double s1=2.5;
//   gf->SetParameter(0,s1);
//   gf->SetParameter(1,W2);
//   vector<realPole> p;
//   double residual_s=gamma4(s).real()/(s1-W2);
//   double residual_mS=gamma4(W2).real()/(W2-s1);
//   double low=16*mpi2;
//   if(s>low)p.push_back(realPole(s1,residual_s)); // this pole only exists above threshold!
//   p.push_back(realPole(W2,residual_mS));
  

   //c2->cd(2);gf->Draw("");
}









