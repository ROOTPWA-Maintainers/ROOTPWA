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
void plotwave(TString tag, const char* select=0, const char* title=0, const char* opt="APZ", double norm=1, int color=kBlack, bool save=false){

  TString tit(tag);

  TString com="intens(";
  if(tag.IsDigit()){
    com+=tag;
    com+="):err(";
    com+=tag;
    com+="):_mass>>h";
    com+=tit;
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
    com+="\"):err(\"";
    com+=tag;
    com+="\"):_mass>>h";
    com+=tit;
  }


  TString sel;
  if(select!=0){
    sel=select;
  }

  cout << com << endl;
  
  pwa->Draw(com.Data(),sel.Data(),"goff");
  int np=pwa->GetSelectedRows();

  double* xerr=new double[np];
  for(int i=0;i<np;++i){
    pwa->GetV3()[i]*=0.001;  //  convert to GeV
    pwa->GetV1()[i]*=norm;
    pwa->GetV2()[i]*=norm;
    xerr[i]=0.060;
  }

  TGraphErrors* g=new TGraphErrors(np,
				   pwa->GetV3(), // mass
				   pwa->GetV1(), // intensity
				   xerr,pwa->GetV2());
  g->SetName(tag);
  g->SetTitle(tag);
  if(title!=0){
    g->SetName(title);
    g->SetTitle(title);
  }

  double maxi=0;
  for(int i=0;i<np;++i){
    if(maxi<g->GetY()[i])maxi=g->GetY()[i];
  }
  cout << "Maximum="<<maxi << endl;
  
  
  g->SetMarkerStyle(21);
  g->SetMarkerSize(0.5);
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  if(title!=0)g->SetTitle(title);
  g->SetMaximum(maxi*1.1);
  g->SetMinimum(-maxi*0.1);
  g->GetXaxis()->SetTitle("mass / GeV");
  
  TGraphErrors* clone=(TGraphErrors*)g->DrawClone(opt);
  clone->GetYaxis()->SetRangeUser(-maxi*0.1,maxi*1.1);
  
  //clone->SetLineColor(kGray);

  if(save){
    TString outfile(tag);
    outfile+=".eps";
    gPad->SaveAs(outfile);
  }

  delete g;

  clone->SetName(tag);
  

  return;
}
