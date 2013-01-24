//-----------------------------------------------------------
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
#include "absDecayChannel.h"


using namespace std;
using namespace rpwa;

// Class Member definitions -----------



rpwa::mcPhaseSpace::mcPhaseSpace(unsigned int n,
			   double* masses,
			   double minM,double maxM,
			   unsigned int nsteps,
			   unsigned int nsamples,
			   int seed)
  : _wmax(0),_thres(minM),_mMax(maxM),_nsteps(nsteps),_nsamples(nsamples),_nparticles(n)
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
rpwa::mcPhaseSpace::addDecayChannel(absDecayChannel* ch)
{
  _channels.push_back(ch);
}
 

unsigned int 
rpwa::mcPhaseSpace::nChannels()const{
  unsigned int nchan=_channels.size();
  if(nchan==0)nchan=1;
  return nchan;
}

void
rpwa::mcPhaseSpace::doCalc() {
  double step=(_mMax-_thres)/(double)_nsteps;
  
  cout << "Prepcalculating PhaseSpace..." << endl;
  cout << "   nsteps   =  " << _nsteps << endl;
  cout << "   minM     =  " << _thres << endl;
  cout << "   maxM     =  " << _mMax << endl;
  cout << "   step     =  " << step << endl;
  cout << "   nSamples =  " << _nsamples << endl;

  for(unsigned int ich=0;ich<nChannels();++ich){
    _graph.push_back(new TGraph(_nsteps));
  }
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
    // calculate phase space for all decay channels
    vector<double> rh(nChannels());
    rho(m,rh);
    for(unsigned int ich=0;ich<nChannels();++ich){
      _graph[ich]->SetPoint(istep,m,rh[ich]);
    }
    
  }
  cout << " ... done " << endl;
}



// double 
// rpwa::mcPhaseSpace::Evaluate(double *x, double *p) {
//   if(!hasCached())throw;
//   return _graph[0]->Eval(x[0]);
// }

double
rpwa::mcPhaseSpace::eval(double m,unsigned int i) {
  if(!hasCached())throw;
  return _graph[i]->Eval(m);
}


void
rpwa::mcPhaseSpace::rho(double m, std::vector<double>& results)const {
  //cout << "rho::M=" << m << endl;
  double ua=m;
  for(unsigned ip=0;ip<_nparticles;++ip){
    ua-=_gen->daughter(ip).M();
  }
  double fn=1./(2*TMath::Power(2*TMath::Pi(),2*_nparticles-3));
  double u=1;
  if(_nparticles>2)u=TMath::Power(ua,_nparticles-2)/TMath::Factorial(_nparticles-2);

  TLorentzVector mother(0,0,0,m);

  // calculate all channels in parallel
  results.clear();
  results.resize(nChannels());

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

    if(_channels.size()==0){
	TLorentzVector Iso1=p[0]+p[1];
	// go into isobar rest frame 
	TVector3 b=-Iso1.BoostVector();
	p[0].Boost(b);
	double breakup=p[0].Vect().Mag();
	double decayprop=0.5;
	results[0] += 1./w * breakup/m * decayprop;
    }
    else { // loop over decay channels
      unsigned int nch=_channels.size();
      for(unsigned int ich=0;ich<nch;++ich){
	double prob;double brkup;
	_channels[ich]->tau(p,prob,brkup);
	results[ich] += 1./w * pow(brkup,2*_channels[ich]->l()+1)/m * prob;
      }// end loop over decay channels
    }

     
  }// end loop over samples
  // loop over decay channels and renormalize
  for(unsigned int ich=0;ich<results.size();++ich){
    results[ich]*= fn * u / (double)_nsamples;
    if(results[ich]!=results[ich]){
      cerr << "mcPhaseSpace:: rho_"<<ich<<"("<<m<<")==Nan"<<endl;
      throw;
    }
  }
}

