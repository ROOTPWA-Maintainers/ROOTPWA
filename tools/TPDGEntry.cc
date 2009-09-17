//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
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
// Description:
//      Implementation of class TPDGEntry
//      see TPDGEntry.h for details
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

// This Class' Header ------------------
#include "TPDGEntry.h"

// C/C++ Headers ----------------------
#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;

// Collaborating Class Headers --------


// Class Member definitions -----------

ClassImp(TPDGEntry)


bool operator== (const TPDGEntry& lhs, const TPDGEntry& rhs)
{
  return lhs.pdgID()==rhs.pdgID();
}




std::istream& operator>> (std::istream& s, TPDGEntry& me){
  char delim;
  char buffer[20];
  char** mark;
  s >> me._mass >>delim>>me. _mass_ep>> delim>> me._mass_en 
    >> delim>> me._width >> delim>> me._width_ep >>delim>>  me._width_en>>delim;
  // read I
  s.get(buffer,10,',');
  // analyse isospin 
  if(strpbrk(buffer,"</?")==NULL){ // we have a integer number
    me._I=atof(buffer);
  }
  else if(strpbrk(buffer,"<")!=NULL){ // photon
    me._I=atof(&(buffer[1]))-1;
  }
  else if(strpbrk(buffer,"/")!=NULL){ //
    double a=atof(buffer);
    double b=atof(&buffer[2]);
    if(b!=0)me._I=a/b;
    else me._I=-1;
  }
  else me._I=-1;
  s >> delim;

  // read g-Parity
  char buffer2[20];
  s.get(buffer2,10,',');
  // analyse isospin 
  if(strpbrk(buffer2,"+-")==NULL){ // we have a integer number
    me._G=0;
  }
  else if(strpbrk(buffer2,"+")!=NULL){ // photon
    me._G=1;
  }
  else me._G=-1;
  s >> delim;

  // read J
  char buffer3[20];
  s.get(buffer3,10,',');
  // analyse isospin 
  if(strpbrk(buffer3,"</?")==NULL){ // we have a integer number
    me._J=atof(buffer3);
  }
  else if(strpbrk(buffer3,"<")!=NULL){ // photon
    me._J=atof(&(buffer3[1]))-1;
  }
  else if(strpbrk(buffer3,"/")!=NULL){ //
    double a=atof(buffer3);
    double b=atof(&buffer3[2]);
    if(b!=0)me._J=a/b;
    else me._J=-1;
  }
  else me._J=-1;
  s >> delim;

  // read Parity
  char buffer4[20];
  s.get(buffer4,10,',');
  //cout<< "P:: "<<buffer<<"|"<<endl;
  // analyse isospin 
  if(strpbrk(buffer4,"+-")==NULL){ // we have a integer number
    me._P=0;
  }
  else if(strpbrk(buffer4,"+")!=NULL){ // photon
    me._P=1;
  }
  else me._P=-1;
  s >> delim;

// read C-Parity
  s.get(buffer,10,',');
  //cout<< "C:: "<<buffer<<"|"<<endl;
  // analyse isospin 
  if(strpbrk(buffer4,"+-")==NULL){ // we have a integer number
    me._C=0;
  }
  else if(strpbrk(buffer4,"+")!=NULL){ // photon
    me._C=1;
  }
  else me._C=-1;
  s >> delim;


  // read antiparticle flag
  char buffer5[20];
  s.get(buffer5,20,',');
  if(strpbrk(buffer5,"BF")==NULL){
    me._aflag = blank;
  }
  else if(strpbrk(buffer5,"B")!=NULL){
    me._aflag = B; 
  }
  else  me._aflag = F;
  s >> delim;
    
  // read pdg  
  char buffer6[20];
  s.get(buffer6,20,',');
  //cout<< "PDG:: "<<buffer<<"|"<<endl;
  me._pdgID=atoi(buffer6);
  s >> delim;
    
  
  // read Q
  char buffer7[20];
  s.get(buffer7,10,',');
  if(strpbrk(buffer7,"+-/")==NULL){ // we have a integer number
    me._q=atof(buffer7);
  }
  else if(strpbrk(buffer7,"/")!=NULL){ //
    double a=atof(buffer7);
    double b=atof(&buffer7[3]);
    me._q=a/b;
  }
  else if(strpbrk(buffer7,"+")!=NULL){ // photon
    me._q=1;
  }
  else me._q=-1;
  s >> delim ;

  char buffer10[20];
  s.get(buffer10,20,',');
  if(strpbrk(buffer10,"1234")==NULL){
    me._R = NoBaryon;
  }
  else if(strpbrk(buffer10,"4")!=NULL){
     me._R = Certain;
  }
  else if(strpbrk(buffer10,"3")!=NULL){
     me._R = Likely;
  }
  else if(strpbrk(buffer10,"2")!=NULL){
     me._R = Fair;
  }
  else  me._R = Poor;
  s >> delim;
  

  char buffer9[20];
  s.get(buffer9,20,',');
  if(strpbrk(buffer9,"RDSF")==NULL){
    me._status = NotDefined;
  }
  else if(strpbrk(buffer9,"R")!=NULL){
     me._status = Established;
  }
  else if(strpbrk(buffer9,"D")!=NULL){
     me._status = Omitted;
  }
  else if(strpbrk(buffer9,"S")!=NULL){
     me._status = NeedsConfirmation;
  }
  else  me._status = FurtherMeson;
  s >> delim;
 
 return s;
 
}

void
TPDGEntry::Print(){
  cout << "---------------------" << endl;
  cout << "PDG Entry: " << _name << "   ID: " << _pdgID << endl;
  cout << "Q=" << _q << endl;
  cout << "mass="<<_mass <<" +- ("<<_mass_ep<<","<< _mass_en << ")" << endl;
  cout << "width="<<_width <<" +- ("<<_width_ep<<","<< _width_en << ")" << endl;
  cout << "IGJPC="<<Istr()<<Gstr()<<Jstr()<<Pstr()<<Cstr()<<endl;
  cout << "Status="<<Statstr() << endl;
}


TString 
TPDGEntry::Istr(){
  if(_I<0)return "?";
  TString result;
  result+=_I;result.ReplaceAll(" ","");
  return result;
}

TString 
TPDGEntry::Gstr(){
  if(_G==0)return "?";
  else if(_G>0) return "+";
  else return "-";
}



TString 
TPDGEntry::Jstr(){
  if(_J<0)return "?";
  TString result;
  result+=_J;result.ReplaceAll(" ","");
  return result;
}

TString 
TPDGEntry::Pstr(){
  if(_P==0)return "?";
  else if(_P>0) return "+";
  else return "-";
}

TString 
TPDGEntry::Cstr(){
  if(_C==0)return "?";
  else if(_C>0) return "+";
  else return "-";
}

TString 
TPDGEntry::Statstr(){
  switch (_status){
  case Established : {return "Established";}
  case Omitted : {return "Omitted";}
  case NeedsConfirmation : {return "Needs Confirmation";}
  case FurtherMeson : {return "Further Meson";}
  default : {return "Undefined";}
  }

}
