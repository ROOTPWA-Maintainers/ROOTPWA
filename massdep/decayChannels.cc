//-----------------------------------------------------------
//
// Description:
//      Implementation of classes in decayChannels
//      see decayChannels.h for details
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

// This Class' Header ------------------
#include "decayChannels.h"

// C/C++ Headers ----------------------
#include<assert.h>
#include<vector>
#include "TVector3.h"
#include "TLorentzVector.h"


// Collaborating Class Headers --------
#include "absMassDep.h"
#include <iostream>

// Class Member definitions -----------

using namespace rpwa;
using namespace std;

void
rpwa::decay21::tau(std::vector<TLorentzVector>& p,
	     double& evtweight, // returns isobar part of phase space for event
	     double& breakup)    // returns breakup momentum of top vertex
{
  assert(p.size()==3);
  TLorentzVector Iso1=p[0]+p[1];
  // go into isobar rest frame 
  TVector3 b=-Iso1.BoostVector();
  p[0].Boost(b);
  double k1=p[0].Vect().Mag();
  int l1=isobar->l();
  double rho1=pow(k1,2*l1+1);///Iso1.M();
  double BW=norm(isobar->val(Iso1.M()));
  evtweight=rho1*BW;
  breakup=p[2].Vect().Mag();
  
}

void
rpwa::decay22::tau(std::vector<TLorentzVector>& p,
	     double& evtweight, // returns isobar part of phase space for event
	     double& breakup)    // returns breakup momentum of top vertex
{
	// 4body -> (12)(34)
  assert(p.size()==4);
  //cout << "(12)(34) phasespace" << endl;
  TLorentzVector Iso1=p[0]+p[1];
  // go into isobar rest frame 
  TVector3 b=-Iso1.BoostVector();
  p[0].Boost(b);
  double k1=p[0].Vect().Mag();
  int l1=isobar1->l();
  double rho1=pow(k1,2*l1+1);///Iso1.M();
  double BW1=norm(isobar1->val(Iso1.M()));
	
  TLorentzVector Iso2=p[2]+p[3];
  // go into isobar rest frame 
  TVector3 b2=-Iso2.BoostVector();
  p[2].Boost(b2);
  //cout << "p2+p3=" << (p[2].Vect()+p[3].Vect()).Mag() << endl;
  double k2=p[2].Vect().Mag();
  int l2=isobar2->l();
  double rho2=pow(k2,2*l2+1);///Iso2.M();
  double BW2=norm(isobar2->val(Iso2.M()));
  
  breakup=Iso1.Vect().Mag();
  evtweight=rho1*BW1*rho2*BW2;
  if(evtweight!=evtweight){
    cerr << "decay22:: decayprop==Nan"<<endl;
    throw;
  }
  if(evtweight==0){
    cerr << "decay22:: decayprop==0"<<endl;
    throw;
  }
}



void
rpwa::decay23::tau(std::vector<TLorentzVector>& p,
	     double& evtweight, // returns isobar part of phase space for event
	     double& breakup)    // returns breakup momentum of top vertex
{
  //  (1(23))4
  assert(p.size()==4);
  //cout << "(12)3))4 phasespace" << endl;
  TLorentzVector Iso1=p[0]+p[1];
  // go into isobar rest frame 
  TVector3 b=-Iso1.BoostVector();
  p[0].Boost(b);
  double k1=p[0].Vect().Mag();
  int l1=isobar2->l();
  double rho1=pow(k1,2*l1+1);///Iso1.M();
  double BW1=norm(isobar2->val(Iso1.M()));
  
  TLorentzVector Iso2=Iso1+p[2];
  // go into isobar rest frame 
  TVector3 b2=-Iso2.BoostVector();
  p[2].Boost(b2);
  //cout << "p2+p3=" << (p[2].Vect()+p[3].Vect()).Mag() << endl;
  double k2=p[2].Vect().Mag();
  int l2=isobar3->l();
  double rho2=pow(k2,2*l2+1);///Iso2.M();
  double BW2=norm(isobar3->val(Iso2.M()));
  
  breakup=p[3].Vect().Mag();
  evtweight=rho1*BW1*rho2*BW2;
  if(evtweight!=evtweight){
    cerr << "decay23:: decayprop==Nan"<<endl;
    throw;
  }
  if(evtweight==0){
    cerr << "decay23:: decayprop==0"<<endl;
    throw;
  }
}
