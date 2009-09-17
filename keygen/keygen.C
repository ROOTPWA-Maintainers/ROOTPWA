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

// This is a library macro

#include "TString.h"
#include "TMath.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

class partkey {
public:
  partkey(const TString& name, partkey* iso, partkey* batch, int l, int s=-1)
    : _name(name), _iso(iso), _batch(batch), _l(l),_s(s),_isFS(false)
  {
    _isofsp=_iso->fsp();
    _batchfsp=_batch->fsp();
  }

  partkey(const TString& name) 
    : _name(name), _iso(0), _batch(0), _l(0),_s(-1), _isFS(true), _id(0)
  {
    _isofsp=1;
    _batchfsp=0;
  }
    
  TString name()const {return _name;}
  partkey* iso(){return _iso;}
  partkey* batch(){return _batch;}
  int l() const {return _l;}
  int s() const {return _s;}

  void write(std::ostream& s, int offset=0);
  void setIds(std::vector<int>& ids);
  // returns number of final state particles in this branch:
  int fsp(){return _isofsp+_batchfsp;}
  void setS(int s){_s=s;}
  void setName(TString name){_name=name;}
  TString wavename();
  void getFSPattern(std::vector<int>& pattern); // fills pattern with + or -

private:
  TString _name;
  partkey* _iso;
  partkey* _batch;
  int _isofsp;   //number of final state particles in iso
  int _batchfsp; //number of final state particles in batch
  int _l; // angular momentum
  int _s; // spin (has to be setup separately)
  bool _isFS; // true if this particle is a final state particle
  int _id; // for final state particles this is the id
  
  void tab(std::ostream& s, int offset){for(int i=0;i<offset;++i)s<<" ";}
};


void
partkey::write(std::ostream& s, int o){
  tab(s,o);
  int len=_name.Length()+4;
  if(!_isFS){
    s << _name << "{" << std::endl;
    _iso->write(s,o+len);
    _batch->write(s,o+len);
    tab(s,o+len);
    s << "l="<<2*_l<<std::endl;
    if(_s!=-1){
      tab(s,o+len);
    s << "s="<<2*_s<<std::endl;
    }
    tab(s,o);
    s << "}" ;
    if(_name=="sigma") s<<" massdep=amp_ves";
    s<<std::endl;

  }
  else {
    s << _name << "["<<_id<<"]"<< std::endl;
  }
}


void 
partkey::setIds(std::vector<int>& ids){
  if(_isFS){
    if(ids.size()!=1)std::cout<<_name<<": inconsistent id distribution!"
			      <<std::endl;
    _id=ids[0];
  }
  else{
    // split the ids and pass them on to isobar and batchelor:
    std::vector<int> isoids;
    for(int i=0;i<_isofsp;++i)isoids.push_back(ids[i]);
    _iso->setIds(isoids);
    std::vector<int> batchids;
    for(int i=_isofsp;i<fsp();++i)batchids.push_back(ids[i]);
    _batch->setIds(batchids);
  }

}


TString 
partkey::wavename(){
  TString wn(_name);
  if(fsp()>2){
    wn+="=";
    wn+=_iso->wavename();
    wn+="_";wn+=_l;
    if(_s>=0)wn+=_s;
    wn+="_";
    wn+=_batch->wavename();
  }
  wn.ReplaceAll("(","");
  wn.ReplaceAll(")","");
  return wn;
}


// this is designed for pure charged pion final states
void 
partkey::getFSPattern(std::vector<int>& pattern){ // fills pattern with + or -
  if(_isFS){
    if(_name.Contains("pi-"))pattern.push_back(-1);
    else if(_name.Contains("pi+"))pattern.push_back(+1);
    else if(_name.Contains("pi0"))pattern.push_back(0);
  }
  else {
    _iso->getFSPattern(pattern);
    _batch->getFSPattern(pattern);
  }
}


// **************************************************************



class wavekey {
public:
  wavekey(int J, int P, int M, int EPS,  partkey* key);

  void buildName(int j, int p, int m);
  void writeSingleAmp(int m, std::ostream& s); // amp for a given m
  void write(std::ostream& s); // combine 2 amps to get into reflectivity base
  void write(); // creates file by himself
  TString wavename(bool file=false);
  void distributeIds(std::vector<int>& ids,
		     const std::vector<int>& pimid,  // piminus ids
		     const std::vector<int>& pipid,  // piplus ids
		     const std::vector<int>& pi0id,  // pizero ids
		     const std::vector<int>& pattern);
private:
  int _J;  // Total Angular momentum
  int _P;  // Parity
  int _M;  // Z-Component of J --> M>=0 in Reflectivity basis
  int _EPS; //Reflectivity
  TString _iso;
  int _l;
  int _s;
  TString _batch;
  partkey* _key;
};


wavekey::wavekey(int J, int P, int M, int EPS, partkey* key)
  : _J(J),_P(P),_M(M),_EPS(EPS), _key(key)
{
  if(_M<0){
    std::cout<<"*******************************************************"
	     <<std::endl;
    std::cout<<"m<0 not allowed for relfectivity basis! Setting m=-m>0!"
	     <<std::endl;
    _M*=-1;
  }

  _iso=key->iso()->name();
  _l=key->l();
  _s=key->s();
  _batch=key->batch()->name();
  buildName(_J,_P,_M);
  // build permutations
  // we use the convention, that in the decay chain the final state particles
  // are always given like that: pi+pi-pi+pi-pi-
  // the 12 permutations are then
  // (11223)(11232)(12213)(12231)(13212)(13221)
  // (21123)(21132)(22113)(22131)(23112)(23121)
  
}


void
wavekey::buildName(int j, int p, int m){ // build name of _key from quantum numbers
  TString name;
  name+="J = ";
  name+=j*2;
  name+=" P = ";
  name+=p;
  name+=" M = ";
  name+=m*2;
  name+=" ";
  _key->setName(name);
}

// to build an amplitude in the reflectivity base
// |ame>=[|am>-refl|a-m>]theta
// refl=eP(-)^(J-m)
void 
wavekey::writeSingleAmp(int refl, std::ostream& s){
  //if(sign_m<0)buildName(_J,_P,-_M);
  double theta= (_M==0 ? 0.5 : 0.7071067811865);
  
  // build permutations
  // check where the pi- and the pi plusses are:
  std::vector<int> pipattern;
  _key->getFSPattern(pipattern);

  // make this generall for N pi.
  // how many pi-minus?
  int npim;
  npim=std::count(pipattern.begin(),pipattern.end(),-1);
  int npip;
  npip=std::count(pipattern.begin(),pipattern.end(),+1);
  int npi0;
  npi0=std::count(pipattern.begin(),pipattern.end(),0);


  std::cout<<"Number of pi-:"<<npim<<std::endl;
  std::cout<<"Number of pi+:"<<npip<<std::endl;
  std::cout<<"Number of pi0:"<<npi0<<std::endl;

  // setup ids for the pions
  std::vector<int> pimid(npim);
  for(int ipi=0;ipi<npim;++ipi){
    pimid[ipi]=ipi+1;
  };
  std::vector<int> pipid(npip);
  for(int ipi=0;ipi<npip;++ipi){
    pipid[ipi]=ipi+1;
  };
  std::vector<int> pi0id(npi0);
  for(int ipi=0;ipi<npi0;++ipi){
    pi0id[ipi]=ipi+1;
  };

  
  // number of permutations
  double npermpip=TMath::Factorial(npip);
  double npermpim=TMath::Factorial(npim);
  double npermpi0=TMath::Factorial(npi0);


  double pipfactor=1./TMath::Sqrt(npermpip);
  double pimfactor=1./TMath::Sqrt(npermpim);
  double pi0factor=1./TMath::Sqrt(npermpi0);


  //::std::showpoint;
  s <<  std::showpoint << std::setprecision(12) << pi0factor << " * (" <<std::endl; // normalization for pi0 symmetrization
  
  int count0=0;
  do{
     int countp=0;
     s << pipfactor << " * ( "<< std::endl; // norm for pip symmetrization 
     do{// permute the pi+s
       int count=0;
       s << pimfactor << " * ( "<< std::endl; // norm for pi- symmetrization 
       do{// permute the pi-s
	 std::vector<int> ids(npip+npim+npi0);
	 distributeIds(ids,pimid,pipid,pi0id,pipattern);
	 _key->setIds(ids);
	 s << " "<<theta<<" * (" << std::endl;
	 buildName(_J,_P,_M);
	 _key->write(s,3);
	 s<< (refl>=0 ? " +" : " - ");
	 buildName(_J,_P,-_M);
	 _key->write(s,3);
	 s<< ")"<<std::endl;
	 if(count++<npermpim-1) s << " + ";
       }
       while(std::next_permutation(pimid.begin(),pimid.end()));
       s<< ")"<<std::endl;
       if(countp++<npermpip-1) s << " + ";
     }
     while(std::next_permutation(pipid.begin(),pipid.end()));
     if(count0++<npermpi0-1) s << ") + ";
  }
  while(std::next_permutation(pi0id.begin(),pi0id.end()));
  s << " )\n)";
     

  /*
  // permute pi+es
  pipid[1]=1;pipid[0]=2;
  count=0;
  do{
    std::vector<int> ids(5);
    distributeIds(ids,pimid,pipid,pipattern);
    _key->setIds(ids);
    s << " "<<theta<<" * (" << std::endl;
    buildName(_J,_P,_M);
    _key->write(s,3);
    s<< (refl>0 ? " + " : " - ");
    buildName(_J,_P,-_M);
    _key->write(s,3);
    s<< ")"<<std::endl;
    if(count++<5)s << "+";
  }
  while(std::next_permutation(pimid.begin(),pimid.end()));
  s << " )\n)";
  */


}


void
wavekey::write(std::ostream& s){
  s << "########################################################" << std::endl;
  s << "# 5 charged pion final state                           #" << std::endl;
  s << "# gamp key file                                        #" << std::endl;
  s << "########################################################" << std::endl;
  s << "# IG JPC Isobar1 [lsm] Isobar2                         #" << std::endl;
  s << "#                                                      #" << std::endl;
  s << "# "<<wavename()<< std::endl;
  s << "#                                                      #" << std::endl;
  s << "########################################################" << std::endl;
  s << "# from here on all ang.mom. in units of 1/2hbar !!!    #" << std::endl;
  s << "########################################################" << std::endl;
  s << std::endl;
  s << "debug = 0;" << std::endl;
  s << "channel = t;" << std::endl;
  s << "mode = binary;" << std::endl;
  s << std::endl;
  // build reflectivity amplitude
  double refl=-(double)_EPS*(double)_P*pow(-1.0,_J-_M);
  std::cerr << "reflectivity factor: "<< refl << std::endl;
  writeSingleAmp(refl,s);
  s<< ";" <<std::endl;
}


void 
wavekey::distributeIds(std::vector<int>& ids,
		       const std::vector<int>& pimid,  // piminus ids
		       const std::vector<int>& pipid,  // piplus ids
		       const std::vector<int>& pi0id,  // pizero ids
		       const std::vector<int>& pattern){
  assert(pattern.size()==pimid.size()+pipid.size()+pi0id.size());
  int pimcount=0;
  int pipcount=0;
  int pi0count=0;
  for(int i = 0; i<pattern.size();++i){
    if(pattern[i]==-1)ids[i]=pimid[pimcount++];
    else if(pattern[i]==+1)ids[i]=pipid[pipcount++];
    else ids[i]=pi0id[pi0count++];
  }
}


void
wavekey::write(){
  ofstream file(wavename(true).Data());
  std::cout<< "Creating keyfile "<<wavename(true)<<std::endl;
  write(file);
  file.close();

}

TString
wavekey::wavename(bool file){
  std::stringstream wv;
  TString res;
  wv<<"1-"<<_J<<(_P>0 ? "+" : "-")<<"+"
    <<_M<<(_EPS>=0 ? "+" : "-")
    <<_key->iso()->wavename()<<"_"<<_l<<_s<<"_"
    <<_key->batch()->wavename();
  if(file)wv << ".key";
  res=wv.str().c_str();
  return res;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

   
void keygen(){

  partkey p1("pi+");
  partkey p2("pi-");
  partkey p3("pi+");
  partkey p4("pi-");
  partkey p5("pi-");
  
  partkey rho1("rho(770)",&p1,&p2,1);

  partkey rho2("rho(770)",&p3,&p4,1);

  partkey a1("a1(1269)",&rho2,&p5,0);
  
  partkey X("X",&rho1,&a1,0,0);

  //            J,P,M,EPS  
  wavekey mykey(1,-1,1,1,&X);

  mykey.write(std::cout);
  
  std::vector<int> pat;
  X.getFSPattern(pat);
  for(int i=0;i<pat.size();++i)std::cout<<pat[i];
  std::cout<<std::endl;
  

}
