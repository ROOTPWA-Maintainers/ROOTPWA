//-----------------------------------------------------------
// Description:
//      Implementation of class NParticleEvent
//      see NParticleEvent.h for details
//
// Environment:
//      Software developed for the COMPASS experiment at CERN.
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

// This Class' Header ------------------
#include "NParticleEvent.h"

// C/C++ Headers ----------------------
#include "TVector3.h"
#include <iomanip>

// Collaborating Class Headers --------
#include "FSParticle.h"

// Class Member definitions -----------
NParticleEvent::NParticleEvent(TClonesArray* fs_momenta, 
			       std::vector<int>* fs_charges,
			       TLorentzVector* beam,
			       int* beam_charge,
			       TVector3* vertex)
  : _fsmomenta(fs_momenta), _fs_charges(fs_charges),
    _beam(beam),_qbeam(beam_charge),_vertex(vertex)
{}  



void
NParticleEvent::refresh(){
  // clear
  _NPStates.clear();
  _fsparticles.clear();

  // build FSParticles
  unsigned int nfs=_fsmomenta->GetEntries();
  for(unsigned int ifs=0;ifs<nfs;++ifs){
    _fsparticles.push_back(FSParticle(*(TLorentzVector*)_fsmomenta->At(ifs),
				      *_vertex,
				      _fs_charges->at(ifs)));

  }

  build();
}


TLorentzVector
NParticleEvent::p(){
  TLorentzVector result;
  unsigned int npart=_fsparticles.size();
  for(unsigned int i=0;i<npart;++i){
    result+=_fsparticles[i].p();
  }
  return result;
}


double
NParticleEvent::tprime(){
  TLorentzVector beam=*_beam;
  TLorentzVector p=this->p();
  // recalibrate beam -- assumes exclusivity!
  TVector3 dir=beam.Vect();
  double const mpi=beam.M();// 0.13957; not always a pion!
  double k=sqrt(p.E()*p.E()-mpi*mpi)/dir.Mag();
  dir*=k;
  beam.SetVectM(dir,mpi);
  return -(beam-p).M2();
}




void
NParticleEvent::toGJ(){
  TLorentzVector tempX=p();
  // rotate event into scattering plane
  // get normal vector
  TVector3 y(0,1,0);
  TVector3 N=_beam->Vect().Cross(tempX.Vect());
  TVector3 rot=N.Cross(y);
  TRotation t;
  double a=N.Angle(y);
  t.Rotate(a,rot);
  //t.SetXEulerAngles(N.Phi(),N.Theta()-TMath::Pi()*0.5,-TMath::Pi()*0.5);
  TLorentzRotation T(t);
  TLorentzRotation L1(T);
  tempX*=T;
  //tempX.Vect().Print();
  _beam->Transform(T);
  //_beamPi.p().Vect().Print();

  // boost to X rest frame
  TVector3 boost=-tempX.BoostVector();
  TLorentzRotation b;
  b.Boost(boost);
  tempX*=b;
  //tempX.Vect().Print();
  _beam->Transform(b);

  // put beam along z-axis
  TVector3 beamdir=_beam->Vect();
  //std::cout<<"beamDir before rotation:";beamdir.Print();
  a=beamdir.Angle(TVector3(0,0,1));
  //std::cout<<"angle="<<a<<std::endl;
  TRotation t2;
  t2.Rotate(a,TVector3(0,1,0));
  T=TLorentzRotation(t2);
  _beam->Transform(T);

  // transform pions
  unsigned int npart=_fsparticles.size();
  for(unsigned int i=0; i<npart; ++i){
    getParticle(i).Transform(L1);
    getParticle(i).Transform(b);
    getParticle(i).Transform(T);
  }
  
  // check basic things:
  //TLorentzVector p2=p();
  //std::cout<<"In GJ-Frame:  P=("<<p2.Px()
  //	   <<","<<p2.Py()<<","<<p2.Pz()<<","<<p2.E()<<std::endl;
  //std::cout<<"In GJ-Frame:  Ptemp=("<<tempX.Px()
  //	   <<","<<tempX.Py()<<","<<tempX.Pz()<<","<<tempX.E()<<std::endl;
  //std::cout<<"Beam direction:";
  //_beamPi.p().Vect().Print();

  // we need to rebuild if the event was build before...
 
}



unsigned int
NParticleEvent::build(){
  // clear previous builds
  _NPStates.clear();
  // build n-particle state
  unsigned int n=_fsparticles.size();
  //std::cout<<n<<" pions in event"<<std::endl;
  for(unsigned int np=1;np<=n;++np){// np = number of particles in (sub)state
    // select np particless
    int* permu=new int[np];
    NParticleEvent::permute(n,np,permu);
    delete[] permu;
  } // end loop npi
  return _NPStates.size();
}


// Recursive method to find all permutations of final state particles
void
NParticleEvent::permute(int n, int k, int* permu, int x, int i){
  // computes all (k out of n) combinations of pions
  // start values: x=-1 i=1;
  if(i>k)return;
  for(int y=x+1;y<n-k+i;++y){
    permu[i-1]=y;
    NParticleEvent::permute(n,k,permu,y,i+1);
    if(i==k){
      NParticleState s;
      s.setBeam(*_beam);
      // add particles to state
      bool ok=true;
      for(int j=0;j<k;++j){
	//std::cout<<"x["<<j<<"]="<<permu[j]<<"  ";
	FSParticle* part=&(_fsparticles.at(permu[j]));
	if(!s.addParticle(part)){// returns false if pion cannot 
	  //be added because doubel counting
	  ok=false;
	  break;
	}
      } // end loop over particles;
      if(ok)_NPStates.push_back(s);
    } // end if(i==k)
  } // end loop over remaining states
}

void
NParticleEvent::writeGAMP(ostream& out){
  // number of particles
  out<<nParticles()+1<<std::endl;
  int geantId=9;
  // beam:
  TLorentzVector p=*_beam;
  int q=*_qbeam>0 ? +1 : -1;
  out << geantId <<" "
      <<q<<" "<< std::setprecision(9)
      <<p.X()<<" " // px
      <<p.Y()<<" " // py
      <<p.Z()<<" " // pz
      <<p.E()<<std::endl; // E

  // pions:
  for(unsigned int i=0;i<nParticles();++i){
    FSParticle pi=_fsparticles[i];
    p=pi.p();
    q=pi.q()>0 ? +1 : -1;
    geantId= q>0 ? 8 : 9;
    out << geantId <<" "
        <<q<<" "
        <<p.X()<<" " // px
        <<p.Y()<<" " // py
        <<p.Z()<<" " // pz
        <<p.E()<<std::endl; // E
  }
}
