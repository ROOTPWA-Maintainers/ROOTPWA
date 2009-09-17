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

// creates latex output from wavefile name
// reades from stdin

#include "TString.h"
#include "TObjString.h"
#include <iostream>
#include <map>
#include "TObjArray.h"

using namespace std;


int main(int argc, char** argv){
  
  // setup isobar dictionary key->tex
  map<TString,TString> isobars;
  isobars["pi-"]="\\pi^-";
  isobars["pi+"]="\\pi^+";
  isobars["pi+-"]="\\pi^\\pm";
  isobars["pi-+"]="\\pi^\\mp";
  isobars["sigma"]="\\sigma";
  isobars["rho770"]="\\rho(770)";
  isobars["a11269"]="a_1(1269)";
  isobars["a21320"]="a_2(1320)";
  isobars["rho1450"]="\\rho(1450)";
  isobars["rho1700"]="\\rho(1700)";
  isobars["pi1300"]="\\pi(1300)";
  isobars["pi1800"]="\\pi(1800)";
  isobars["pi21670"]="\\pi_2(1670)";
  isobars["f01370"]="f_0(1370)";
  isobars["f01500"]="f_0(1500)";
  isobars["f01700"]="f_0(1700)";
  isobars["f11285"]="f_1(1285)";
  isobars["f11420"]="f_1(1420)";
  isobars["b11235"]="b_1(1235)";
  isobars["b11800"]="b_1(1800)";
  isobars["b11500"]="b_1(1500)";
  isobars["f21270"]="f_2(1270)";
  isobars["f21950"]="f_2(1950)";
  isobars["f21565"]="f_2(1565)";
  isobars["f21270"]="f_2(1270)";
  isobars["f22010"]="f_2(2010)";
  isobars["f11420"]="f_1(1420)";
  isobars["eta1440"]="\\eta(1420)";
  isobars["eta21645"]="\\eta_2(1645)";
  isobars["rho31690"]="\\rho_3(1690)";
  isobars["a21320"]="a_2(1320)";


  // print latex header
  cout<<"\\documentclass[12pt,a4paper]{article}"<<endl;
  cout<<"\\usepackage{amsmath, amsthm, amssymb}"<<endl;
  cout<<"\\begin{document}"<<endl;
  cout << "\\begin{align*}" << endl;
  cout<<"\\begin{aligned}"<<endl;
  char line[300];
  int count=0;
  bool wasthr=false;
  while(!(cin>>line).eof()) { // begin event loop

    TString ls(line);
    // check if there is a threshold
    
    
    if(ls.IsAlnum()){
      cout << ls;
      continue;
    }
    

    if(count>0)cout << " \\\\" << endl;
    // make pagebreak every 15 waves
    if(count>0 && count%15==0){
       cout<<"\\end{aligned}"<<endl;
       cout<<"\\end{align*}"<<endl;
       cout<<"\\pagebreak"<<endl;
       cout << "\\begin{align*}" << endl;
       cout<<"\\begin{aligned}"<<endl;
    }

    // tokenize input
    ++count;
    // extract header
    TString head=ls(0,7);
    // tokenize header
    TString I=head(0,1); 
    TString g=head(1,1);
    TString J=head(2,1);
    TString p=head(3,1);
    TString c=head(4,1);
    TString m=head(5,1);
    TString eps=head(6,1);
    // remove file extension
    ls.Remove(ls.Length()-4);
    ls.Remove(0,7);
    //cout << I<<g<<J<<p<<c<<m<<eps <<" -- " << ls << endl;

    TObjArray* tokens=ls.Tokenize("_=");

    //cout << tokens->GetEntries() << " tokens found " << endl;

    int mode=0;

    cout << I<<"^"<<g<<J<<"^{"<<p<<c<<"}"<<m<<"^"<<eps<<"\\quad & ";

    for(int i=0;i<tokens->GetEntries();++i){
      TString tok=((TObjString*)tokens->At(i))->GetString();
      //cout << tok << endl;
      if(mode==0){ // isobar mode
	if(isobars[tok].Length()==0)cout << tok <<" ";
	else cout << isobars[tok]<<" " ;
	// check which mode to switch to 
	// depending whether we get _ or =
	ls.Remove(0,tok.Length());
	if(ls(0,1)=="_")mode=1;
	else mode=2;
      }
      else if(mode==1){ // ls mode
	if(tok.Length()==1){// only l
	  cout << "["<<tok<<"]"<<" ";
	}
	else {
	  cout << "\\left[\\begin{array}{c}"<<tok(0,1)<<"\\\\"
	       <<tok(1,1)<<"\\end{array}\\right]";
	}
	ls.Remove(0,tok.Length());
	mode=0;
      }
      else if(mode==2){
	if(isobars[tok].Length()==0)cout << "\\rightarrow " << tok <<" ";
	else cout << "\\rightarrow " << isobars[tok] <<" ";
	ls.Remove(0,tok.Length());
	if(ls(0,1)=="_")mode=1;
	else mode=2;
      }
      //cout << " mode = " << mode << endl;
      ls.Remove(0,1); // remove delimiter
    }
    cout << " & ";    

  } // input loop
  cout<<"\\end{aligned}"<<endl;
  cout<<"\\end{align*}"<<endl;
  cout<<"\\end{document}"<<endl;

};
