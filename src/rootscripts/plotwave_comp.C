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
void plotwave_comp(TString tag, TString qtag, double scale=1, const char* select=0, const char* title=0, bool save=false){

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
    com+="\"):err(\"";
    com+=tag;
    com+="\"):_mass";
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
    pwa->GetV1()[i]*=scale;
    pwa->GetV2()[i]*=scale;
    pwa->GetV3()[i]*=0.001;  //  convert to GeV
    xerr[i]=0.02;
  }

  TGraphErrors* g=new TGraphErrors(np,
				   pwa->GetV3(), // mass
				   pwa->GetV1(), // intensity
				   xerr,pwa->GetV2());
  

  g->SetName(tag);
  g->SetTitle(tag);

  

  double maxi=0;
  for(int i=0;i<np;++i){
    if(maxi<g->GetY()[i])maxi=g->GetY()[i];
  }
  cout << maxi << endl;
  
  g->SetTitle(tit);
  //g->SetMarkerStyle(20);
  g->SetMarkerSize(0.8);
  if(title!=0)g->SetTitle(title);
  g->SetMaximum(maxi*1.1);
  g->SetMinimum(-maxi*0.1);
  TGraphErrors* clone=(TGraphErrors*)g->DrawClone("AP");
  clone->GetYaxis()->SetRangeUser(-maxi*0.1,maxi*1.1);
  clone->SetLineColor(kRed);
  clone->SetMarkerColor(kRed);
  clone->GetXaxis()->SetTitle("Mass / GeV");
  
  clone->GetYaxis()->SetTitle("Intensity");
  clone->GetYaxis()->SetTitleOffset(1.3);

  TString pathname = "/lustre/e18/user/sneubert/Q3PiData/";
  TString filename1="hfit_t_0p1_1p0_zemach_42waves_tdep_30a_2.root";
  
  TFile file1(pathname+filename1);


  TH1F* hist1 = (TH1F*)file1.Get(qtag);//results from mass-independent fit
  if(hist1!=NULL){
    hist1->DrawCopy("same");
  }
  else { cout << "Reference Histo not found!" << endl; }

  double xmax=clone->GetXaxis()->GetXmax();
  TLegend* leg=new TLegend(0.65,0.75,0.85,0.85);
  leg->AddEntry(hist1,"Weitzel/Ryabch","LE");
  leg->AddEntry(clone,"pwaroot (scaled)","LE");

  leg->Draw();

  if(save){
    TString outfile(tag);
    outfile+=".eps";
    gPad->SaveAs(outfile);
  }

  
  delete g;
  return;
}
