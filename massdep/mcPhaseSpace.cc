//-----------------------------------------------------------
// File and Version Information:
// $Id$
//

// This Class' Header ------------------
#include "mcPhaseSpace.h"

// C/C++ Headers ----------------------
#include <iostream>
#include <assert.h>
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TGraph.h"

// Collaborating Class Headers --------
#include "nBodyPhaseSpaceGen.h"
#include "absMassDep.h"

using namespace std;
using namespace rpwa;

// Class Member definitions -----------



rpwa::mcPhaseSpace::mcPhaseSpace(unsigned int n,
			   double* masses,
			   double minM,double maxM,
			   unsigned int nsteps,
			   unsigned int nsamples,
			   int seed)
  : _wmax(0),_thres(minM),_mMax(maxM),_nsteps(nsteps),_nsamples(nsamples),_nparticles(n),_opt(0),_graph(NULL)
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
  

}


void
rpwa::mcPhaseSpace::doCalc(int l) {
  double step=(_mMax-_thres)/(double)_nsteps;
  
  cout << "Prepcalculating PhaseSpace..." << endl;
  cout << "   nsteps   =  " << _nsteps << endl;
  cout << "   minM     =  " << _thres << endl;
  cout << "   maxM     =  " << _mMax << endl;
  cout << "   step     =  " << step << endl;
  cout << "   nSamples =  " << _nsamples << endl;
  _graph=new TGraph(_nsteps);
  int tenpercent=(int)(_nsteps*0.1);
 
  for(unsigned int istep=0;istep<_nsteps;++istep){
    // ************ Progress messages ************************
    if(istep%10==0){std::cout<<".";cout.flush();}
    if(istep%tenpercent==0){
      std::cout<<"["
               <<ceil((double)istep*100/(double)_nsteps)<<"%"
               <<"]";
      std::cout.flush();
    }
    // ******************************************************

    double m=_thres+((double)istep+0.5)*step;
    _graph->SetPoint(istep,m,rho(m,l));
    
  }
  cout << " ... done " << endl;
}



double 
rpwa::mcPhaseSpace::Evaluate(double *x, double *p) {
  if(!hasCached())throw;
  return _graph->Eval(x[0]);
}

double
rpwa::mcPhaseSpace::eval(double m) {
  if(!hasCached())throw;
  return _graph->Eval(m);
}


void 
rpwa::mcPhaseSpace::setSubSystems21(absMassDep* iso){ // 3-body decay (pp)p
  _opt=1;
  _isobars.push_back(iso);
}

void 
rpwa::mcPhaseSpace::setSubSystems22(absMassDep* iso1, 
		     absMassDep* iso2){  // 4-body decay (pp)(pp)
  _opt=2;
  _isobars.push_back(iso1);
  _isobars.push_back(iso2);
  
}

void 
rpwa::mcPhaseSpace::setSubSystems132(absMassDep* iso1, 
			absMassDep* iso2){ // 4-body decay  p((pp)p)
  _opt=3;
  _isobars.push_back(iso1);
  _isobars.push_back(iso2);
}



double
rpwa::mcPhaseSpace::rho(double m, int l)const {
  //cout << "rho::M=" << m << endl;
  double ua=m;
  for(unsigned ip=0;ip<_nparticles;++ip){
    ua-=_gen->daughter(ip).M();
  }
  double fn=1./(2*TMath::Power(2*TMath::Pi(),2*_nparticles-3));
  double u=1;
  if(_nparticles>2)u=TMath::Power(ua,_nparticles-2)/TMath::Factorial(_nparticles-2);

  

  TLorentzVector mother(0,0,0,m);
  double I=0;
  //double wsum=0;
  for(unsigned int i=0;i<_nsamples;++i){
    double w=_gen->generateDecay(mother);
    if(_nparticles==2)w=1;
    
    // get event
    vector<TLorentzVector> p(_nparticles);
    for(unsigned ip=0;ip<_nparticles;++ip){
      p[ip]=_gen->daughter(ip);
    }
    // calculate contribution of this event

    double decayprop=0;
    double breakup=0;

    switch(_opt){
      // 2body
    case 0:  
      {
	
	TLorentzVector Iso1=p[0]+p[1];
	// go into isobar rest frame 
	TVector3 b=-Iso1.BoostVector();
	p[0].Boost(b);
	breakup=p[0].Vect().Mag();
	decayprop=0.5;
	break;
      }
      // 3body -> (12)3

    case 1:
      {
	assert(_nparticles==3);
	TLorentzVector Iso1=p[0]+p[1];
	// go into isobar rest frame 
	TVector3 b=-Iso1.BoostVector();
	p[0].Boost(b);
	double k1=p[0].Vect().Mag();
	int l1=_isobars[0]->l();
	double rho1=pow(k1,2*l1+1);///Iso1.M();
	double BW=norm(_isobars[0]->val(Iso1.M()));
	decayprop=rho1*BW;
	breakup=p[2].Vect().Mag();
	break;
      }

	// 4body -> (12)(34)
    case 2: 
      {
	// double mass=0.77549;
        // double gamma=0.1462;

	assert(_nparticles==4);
	//cout << "(12)(34) phasespace" << endl;
	TLorentzVector Iso1=p[0]+p[1];
	// go into isobar rest frame 
	TVector3 b=-Iso1.BoostVector();
	p[0].Boost(b);
	double k1=p[0].Vect().Mag();
	int l1=_isobars[0]->l();
	double rho1=pow(k1,2*l1+1);///Iso1.M();
	double BW1=norm(_isobars[0]->val(Iso1.M()));
	
	TLorentzVector Iso2=p[2]+p[3];
	// go into isobar rest frame 
	TVector3 b2=-Iso2.BoostVector();
	p[2].Boost(b2);
	//cout << "p2+p3=" << (p[2].Vect()+p[3].Vect()).Mag() << endl;
	double k2=p[2].Vect().Mag();
	int l2=_isobars[1]->l();
	double rho2=pow(k2,2*l2+1);///Iso2.M();
	double BW2=norm(_isobars[1]->val(Iso2.M()));

	breakup=Iso1.Vect().Mag();
	decayprop=rho1*BW1*rho2*BW2;
	if(decayprop!=decayprop){
	  cerr << "mcPhaseSpace:: decayprop("<<m<<")==Nan"<<endl;
	  throw;
	}
	if(decayprop==0){
	  cerr << "mcPhaseSpace:: decayprop("<<m<<")==0"<<endl;
	  throw;
	}
	break;
      }
  
    case 3: 
      {

	//  (1(23))4
	// double mass=0.77549;
        // double gamma=0.1462;

	assert(_nparticles==4);
	//cout << "(12)3))4 phasespace" << endl;
	TLorentzVector Iso1=p[0]+p[1];
	// go into isobar rest frame 
	TVector3 b=-Iso1.BoostVector();
	p[0].Boost(b);
	double k1=p[0].Vect().Mag();
	int l1=_isobars[0]->l();
	double rho1=pow(k1,2*l1+1);///Iso1.M();
	double BW1=norm(_isobars[0]->val(Iso1.M()));
	
	TLorentzVector Iso2=Iso1+p[2];
	// go into isobar rest frame 
	TVector3 b2=-Iso2.BoostVector();
	p[2].Boost(b2);
	//cout << "p2+p3=" << (p[2].Vect()+p[3].Vect()).Mag() << endl;
	double k2=p[2].Vect().Mag();
	int l2=_isobars[1]->l();
	double rho2=pow(k2,2*l2+1);///Iso2.M();
	double BW2=norm(_isobars[1]->val(Iso2.M()));

	breakup=p[3].Vect().Mag();
	decayprop=rho1*BW1*rho2*BW2;
	if(decayprop!=decayprop){
	  cerr << "mcPhaseSpace:: decayprop("<<m<<")==Nan"<<endl;
	  throw;
	}
	if(decayprop==0){
	  cerr << "mcPhaseSpace:: decayprop("<<m<<")==0"<<endl;
	  throw;
	}
	break;
      }
    }

// double BWwheight1=TMath::BreitWigner(Iso1.M(),mass,gamma);
//     double BWwheight2=TMath::BreitWigner(Iso2.M(),mass,gamma);

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
    int l=0;

    //cout << "W=" <<w << endl;
    //cout << fn*u*pow(breakup,2*l+1)/m * decayprop << endl;

    I += 1./w * pow(breakup,2*l+1)/m * decayprop;
    
    //I+= w* pow(k,2*l)/m *rho1 *rho2;
 
    //I+=w;
    
  }// end loop over samples
  I*= fn * u / (double)_nsamples;
  if(I!=I){
    cerr << "mcPhaseSpace:: rho("<<m<<")==Nan"<<endl;
    throw;
  }
  //cerr << "mcPhaseSpace:: rho("<<m<<")="<<I<<endl;
  return I;
}

