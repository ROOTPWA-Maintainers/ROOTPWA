///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
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
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      class that describes decay vertex of isobar into two particles
//      the isobar -> particle1 + particle 2 vertex has exactly one
//      incoming mother and two outgoing daughter particle
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "utilities.h"
#include "isobarDecayVertex.h"
#include "angMomCoupl.h"

	
using namespace std;
using namespace rpwa;


bool isobarDecayVertex::_debug = false;


isobarDecayVertex::isobarDecayVertex(particle&          mother,
				     particle&          daughter1,
				     particle&          daughter2,
				     const unsigned int L,
				     const unsigned int S)
  : interactionVertex(),
    _L               (L),
    _S               (S)
{
  if (_debug)
    printInfo << "contructing isobar decay vertex "
	      << mother.name()    << "  --->  "
	      << daughter1.name() << "  +  " << daughter2.name() << endl;
  interactionVertex::addInParticle (mother);
  interactionVertex::addOutParticle(daughter1);
  interactionVertex::addOutParticle(daughter2);
}


isobarDecayVertex::isobarDecayVertex(const isobarDecayVertex& vert)
{
  *this = vert;
}


isobarDecayVertex::~isobarDecayVertex()
{ }


isobarDecayVertex&
isobarDecayVertex::operator = (const isobarDecayVertex& vert)
{
  if (this != &vert) {
    interactionVertex::operator = (vert);
    _L = vert._L;
    _S = vert._S;
  }
  return *this;
}

    
bool
isobarDecayVertex::checkMultiplicativeQn(const int     mQn,
					 const int     d1Qn,
					 const int     d2Qn,
					 const string& qnName)
{
  if (mQn == d1Qn * d2Qn)
    return true;
  else {
    if (_debug)
      printWarn << ((qnName == "") ? "multiplicative quantum number": qnName) << " mismatch: "
		<< mother().name() << " = " << mQn << " != " << d1Qn * d2Qn  << " = "
		<< "(" << daughter1().name() << " = " << d1Qn << ") * "
		<< "(" << daughter2().name() << " = " << d2Qn << ")" << endl;
    return false;
  }
}


bool
isobarDecayVertex::checkAdditiveQn(const int     mQn,
				   const int     d1Qn,
				   const int     d2Qn,
				   const string& qnName)
{
  if (mQn == d1Qn + d2Qn)
    return true;
  else {
    if (_debug)
      printWarn << ((qnName == "") ? "additive quantum number": qnName) << " mismatch: "
		<< mother().name() << " = " << mQn << " != " << d1Qn + d2Qn  << " = "
		<< "(" << daughter1().name() << " = " << d1Qn << ") + "
		<< "(" << daughter2().name() << " = " << d2Qn << ")" << endl;
    return false;
  }
}


bool 
isobarDecayVertex::checkConsistency(){
  bool vertexConsistent = true;
  // check multiplicative quantum numbers
  // G-parity
  if (!checkMultiplicativeQn(mother().G(), daughter1().G(), daughter2().G(), "G-parity"))
    vertexConsistent = false;
  // C-parity
  if (!checkMultiplicativeQn(mother().C(), daughter1().C(), daughter2().C(), "C-parity"))
    vertexConsistent = false;
  // Parity
  const int angMomParity = (_L % 4 == 0) ? 1 : -1; // modulo 4 because L is in units of hbar/2
  if (mother().P() != daughter1().P() * daughter2().P() * angMomParity) {
    if (_debug)
      printWarn << "parity mismatch: "
		<< mother().name() << " = " << mother().P() << " != "
		<< daughter1().P() * daughter2().P() * angMomParity  << " = "
		<< "(" << daughter1().name() << " = " << daughter1().P() << ") * "
		<< "(" << daughter2().name() << " = " << daughter2().P() << ") * "
		<< "(ang. momentum = " << angMomParity << ")" << endl;
    vertexConsistent = false;
  }
  // check additive quantum numbers
  // charge
  if (!checkAdditiveQn(mother().charge(),      daughter1().charge(),      daughter2().charge(),      "charge"))
    vertexConsistent = false;
  // baryon number
  if (!checkAdditiveQn(mother().baryonNmb(),   daughter1().baryonNmb(),   daughter2().baryonNmb(),   "baryonNmb"))
    vertexConsistent = false;
  // strangeness
  if (!checkAdditiveQn(mother().strangeness(), daughter1().strangeness(), daughter2().strangeness(), "strangeness"))
    vertexConsistent = false;
  // charm
  if (!checkAdditiveQn(mother().charm(),       daughter1().charm(),       daughter2().charm(),       "charm"))
    vertexConsistent = false;
  // beautty
  if (!checkAdditiveQn(mother().beauty(),      daughter1().beauty(),      daughter2().beauty(),      "beauty"))
    vertexConsistent = false;
  // check angular momentum like quantum numbers
  // spin coupling: S in {|s1 - s2|, ..., s1 + s2}
  if (!angMomCoupl(daughter1().J(), daughter2().J()).inRange(_S)) {
    if(_debug)
      printWarn << "spins "
		<< "(" << daughter1().name() << " 2J = " << daughter1().J() << ") and "
		<< "(" << daughter2().name() << " 2J = " << daughter2().J() << ") "
		<< "cannot couple to total spin 2S = " << _S << endl;
    vertexConsistent = false;
  }
  // L-S coupling: J in {|L - S|, ..., L + S}
  if (!angMomCoupl(_L, _S).inRange(mother().J())) {
    if (_debug)
      printWarn << "orbital angular momentum 2L = " << _L << " and spin 2S = " << _S
		<< " cannot couple to angular momentum 2J = " << mother().J() << endl;
    vertexConsistent = false;
  }
  // isospin coupling: I in {|I_1 - I_2|, ..., I_1 + I_2}
  if (!angMomCoupl(daughter1().isospin(), daughter2().isospin()).inRange(mother().isospin())) {
    if (_debug)
      printWarn << "isospins "
		<< "(" << daughter1().name() << " 2I = " << daughter1().isospin() << ") and "
		<< "(" << daughter2().name() << " 2I = " << daughter2().isospin() << ") "
		<< "cannot couple to total isospin 2I = " << mother().isospin() << endl;
    vertexConsistent = false;
  }
  //!!! missing: spin projections
  if (!vertexConsistent && _debug)
    printWarn << "vertex data are inconsistent (see warings above):" << endl
	      << *this << flush;
  return vertexConsistent;
}


const TLorentzVector&
isobarDecayVertex::updateMotherLzVec()
{
  if (_debug)
    printInfo << "updating Lorentz-vector of particle " << mother().name()
	      << " p_before = " << mother().lzVec() << " GeV, " << flush;
  mother().setLzVec(daughter1().lzVec() + daughter2().lzVec());
  cout << "p_after = " << mother().lzVec() << " GeV" << endl;
  return mother().lzVec();
}


void
isobarDecayVertex::getListOfValidDecays(vector<isobarDecayVertex*>& d1list,
					vector<isobarDecayVertex*>& d2list,
					vector<isobarDecayVertex*>& outlist,
					int maxl,bool blockExotic){
  // save original daughters
  particle* d1orig=outParticles()[0];
  particle* d2orig=outParticles()[1];
 
  // here we exploit the fact that in the isobarmodel each vertex
  // has a unique association with its mother particle
  unsigned int n1=d1list.size();
  unsigned int n2=d2list.size();
  for(unsigned int i1=0;i1<n1;++i1){ // loop over daughters1
    if(d1list[i1]!=NULL)_outParticles[0]=&(d1list[i1]->mother());
    for(unsigned int i2=0;i2<n2;++i2){ // loop over daughters2
      if(d2list[i2]!=NULL)_outParticles[1]=&(d2list[i2]->mother());
      getListOfValidDecays(outlist,maxl,blockExotic);
    } // end loop over daughters2
  } // end loop over daughters1

  // restore original daughters
  _outParticles[0]=d1orig;
  _outParticles[1]=d2orig;

}



void
isobarDecayVertex::getListOfValidDecays(vector<isobarDecayVertex*>& outlist, 
					int maxl,
					bool blockExotic){
  
  int q=daughter1().charge()*daughter2().charge();
  int G=daughter1().G()*daughter2().G();
  int strangeness=daughter1().strangeness()+daughter2().strangeness();
  int charm=daughter1().charm()+daughter2().charm();
  int beauty=daughter1().beauty()+daughter2().beauty();
  
  // s-s coupling loop
  int s1=daughter1().J();
  int s2=daughter2().J();
  int I1=daughter1().isospin();
  int I2=daughter2().isospin();
  
  //vector<int> scoupl=angMomCoupl(s1,s2).getCombinations(); Is this worth it?
  for(int s=abs(s1-s2);s<=(s1+s2);s+=2){
    // l loop
    for(int l=0;l<=maxl;l+=2){
      int P=daughter1().P() * daughter2().P() * (l % 4 == 0 ? 1 : -1);
      // l-s coupling loop
      for(int J=abs(l-s);J<=l+s;J+=2){	
	for(int I=abs(I1-I2); I<= I1+I2 && I<=2; I+=2){

	  int C=G * (I % 4 == 0 ? 1 : -1);

	  cout <<"IG="<<I<<sign(G)<<"   Jpc="<<J<<sign(P)<<sign(C)
	       <<"   s="<<s<<"   l="<<l<<endl;
	  
	  if(blockExotic){
	    // Quark model boundary conditions:
	    // for P==C everything ok
	    // check P=(-1)^(J+1)
	    if(P!=C && (C != (J % 4 == 0 ? 1 : -1) || P != (J+2 % 4 == 0 ? 1 : -1)) ){
	      cout << "Blocking spin-exotic isobar!" << endl;
	      continue;
	    }
	  }
	

	  // build new mother & vertex
	  particle* newmo=new particle("Isobar",I,G,J,P,C,0);
	  newmo->setCharge(q);newmo->setStrangeness(strangeness);
	  newmo->setCharm(charm);newmo->setBeauty(beauty);
	  outlist.push_back(new isobarDecayVertex(*newmo,
						  daughter1(),
						  daughter2(),l,s));

	}// end isospin loop
	
      }// end l-s coupling loop
      
    }// end l loop
    
  }// end s-s coupling loop
  
  
  
}









ostream&
isobarDecayVertex::print(ostream& out) const
{
  out << "isobar decay vertex "
      << "(data are " << ((!dataAreValid()) ? "not " : "") << "valid):" << endl
      << "    mother "     << *(inParticles()[0])  << endl
      << "    daughter 1 " << *(outParticles()[0]) << endl
      << "    daughter 2 " << *(outParticles()[1]) << endl
      << "    2L = " << _L << ", 2S = " << _S << endl;
  return out;
}
