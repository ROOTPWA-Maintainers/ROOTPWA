//-----------------------------------------------------------
// File and Version Information:
// $Id$
//

// This Class' Header ------------------
#include "mcPhaseSpace.h"

// C/C++ Headers ----------------------
#include <iostream>
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TGraph.h"

// Collaborating Class Headers --------
#include "nBodyPhaseSpaceGen.h"


using namespace std;
using namespace rpwa;

// Class Member definitions -----------



mcPhaseSpace::mcPhaseSpace(unsigned int n,
			   double* masses,
			   double minM,double maxM,
			   unsigned int nsteps,
			   unsigned int nsamples,
			   int seed)
  : _thres(minM),_nsamples(nsamples),_nparticles(n)
{
  // prepare generator

  _gen = new nBodyPhaseSpaceGen();
  _gen->setWeightType(nBodyPhaseSpaceGen::S_U_CHUNG);
  _gen->setKinematicsType(nBodyPhaseSpaceGen::BLOCK);
  _gen->setDecay(n,masses);
  _gen->setSeed(seed);
  _wmax=1.01 * _gen->estimateMaxWeight(maxM);
  cout << "Maximal weight=" << _wmax << endl;
  _gen->setMaxWeight(_wmax);

  gRandom->SetSeed(seed);
  

  
  // check minimum mass;
  double threshold=0;
  for(unsigned k=0;k<n;++k)threshold+=masses[k];
  if(_thres<threshold){
    cerr << "Setting minimum mass to threshold "<< threshold << endl;
    _thres=threshold;
  }
  

  double step=(maxM-_thres)/(double)nsteps;
  
  cout << "Prepcalculating PhaseSpace..." << endl;
  cout << "   nsteps   =  " << nsteps << endl;
  cout << "   minM     =  " << _thres << endl;
  cout << "   maxM     =  " << maxM << endl;
  cout << "   step     =  " << step << endl;
  cout << "   nSamples =  " << nsamples << endl;
  _graph=new TGraph(nsteps);
  for(unsigned int istep=0;istep<nsteps;++istep){
    double m=_thres+((double)istep+0.5)*step;
    _graph->SetPoint(istep,m,rho(m));
    
  }
  cout << " ... done " << endl;
  

}


double 
mcPhaseSpace::Evaluate(double *x, double *p)const {
  return _graph->Eval(x[0],0,"S");
}

double
mcPhaseSpace::eval(double m) const {
  return _graph->Eval(m,0,"S");
}


double
mcPhaseSpace::rho(double m)const {
  
  double ua=m;
  for(unsigned ip=0;ip<_nparticles;++ip){
    ua-=_gen->daughter(ip).M();
  }
  double f4=1./19.58525983E3;

  TLorentzVector mother(0,0,0,m);
  double I=0;
  //double wsum=0;
  for(unsigned int i=0;i<_nsamples;++i){
    double w=_gen->generateDecay(mother);
    //w+=1;
    // get event
    vector<TLorentzVector> p(_nparticles);
    for(unsigned ip=0;ip<_nparticles;++ip){
      p[ip]=_gen->daughter(ip);
    }
    // calculate contribution of this event
    // // for the moment say we combine particle 1&2 into an isobar
    TLorentzVector Iso1=p[0]+p[1];
    // go into isobar rest frame 
    TVector3 b=-Iso1.BoostVector();
    p[0].Boost(b);
    double k1=p[0].Vect().Mag();
    int l1=1;
    double rho1=pow(k1,2*l1+1);///Iso1.M();

    TLorentzVector Iso2=p[2]+p[3];
    // go into isobar rest frame 
    TVector3 b2=-Iso2.BoostVector();
    p[2].Boost(b2);
    p[3].Boost(b2);
    //cout << "p2+p3=" << (p[2].Vect()+p[3].Vect()).Mag() << endl;
    double k2=p[2].Vect().Mag();
    int l2=1;
    double rho2=pow(k2,2*l2+1);///Iso2.M();
    

    
    double gamma=0.1462;
    double mass=0.77549;
    //double gamma=0.6;
    //double mass=0.4;
    
    double BWwheight1=TMath::BreitWigner(Iso1.M(),mass,gamma);
    double BWwheight2=TMath::BreitWigner(Iso2.M(),mass,gamma);


  //   // 4pi=a1 pi=rhopi pi=pipipipi
//     TLorentzVector Iso1=p[0]+p[1];
//     TLorentzVector Iso2=p[0]+p[1]+p[2];
//     // go into isobar rest frame 
//     TVector3 b=-Iso1.BoostVector();
//     p[0].Boost(b);
//     double k1=p[0].Vect().Mag();
//     int l1=1;
//     double rho1=pow(k1,2*l1+1);///Iso1.M();
//     TVector3 b2=-Iso2.BoostVector();
//     // go into a1 restframe
//     p[2].Boost(b2);
//     double k2=p[2].Vect().Mag();
//     int l2=0;
//     double rho2=pow(k2,2*l2+1);///Iso1.M();
    
    
//     double BWwheight1=TMath::BreitWigner(Iso1.M(),0.77549,0.1462);
//     double BWwheight2=TMath::BreitWigner(Iso2.M(),1.230,0.280);



    // breakup momentum
    double k = Iso2.Vect().Mag();
    int l=0;
    double g= w* BWwheight1 * BWwheight2 *pow(k,2*l+1)/m *rho1 *rho2;//; * ua*ua/2;
    g+=w;
    I+= 1./w *  pow(k,2*l+1)/m*  BWwheight1 * BWwheight2 *rho1 *rho2 ;
    //I+= w* pow(k,2*l)/m *rho1 *rho2;
 
    //I+=w;
    
  }// end loop over samples
  I*=  TMath::Exp(-(m*m-2.1025))*ua*ua/2. * f4 / (double)_nsamples;


  return I;
}

