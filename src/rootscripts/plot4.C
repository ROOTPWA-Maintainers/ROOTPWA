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
void plot4(int i,int j,double mlow=1, double mup=4){
  gROOT->Reset();
  gROOT->ProcessLine(".L ./plotwave.C");
  gROOT->ProcessLine(".L ./plotphase.C");
  
  mlow*=1000;
  mup*=1000;

  TString sel=",\"_mass>=";sel+=mlow;sel+=" && ";sel+="_mass<=";sel+=mup;
  sel+="\")";

  TCanvas* c=new TCanvas("c","c",10,10,600,800);
  c->Divide(1,3);
 
  c->cd(1);
  TString com="plotwave(\"";com+=i;com+="\"";com+=sel;
  gROOT->ProcessLine(com);

  c->cd(2);
  com="plotwave(\"";com+=j;com+="\"";com+=sel;
  gROOT->ProcessLine(com);

  c->cd(3);
  com="plotphase(";com+=i;com+=",";com+=j;com+=sel;
  gROOT->ProcessLine(com);

  return;
}
 
