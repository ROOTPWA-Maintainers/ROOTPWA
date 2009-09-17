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
void plotmcmc(TString tag, const char* select=0, const char* title=0, const char* opt="COL", bool save=false){

  TString tit(tag);

  TString com="intens(";
  if(tag.IsDigit()){
    com+=tag;
    com+="):err(";
    com+=tag;
    com+="):_mass";
    if(title==NULL){
      TFitBin* bin=new TFitBin();
      pwa->SetBranchAddress("fitbin",&bin);
      pwa->GetEntry(0);
      tit=bin->wavename(atoi(tag.Data()));
    }

  }
  else{
    com+="\"";
    com+=tag;
    com+="\"):_mass";
  }


  TString sel;
  if(select!=0){
    sel=select;
  }

  cout << com << endl;
  
  pwa->Draw(com.Data(),sel.Data(),"COL");
  
  return;
}
