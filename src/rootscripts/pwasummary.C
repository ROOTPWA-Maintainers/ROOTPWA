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
void pwasummary(int color=kBlack, double norm=1){

  TString opt("AP");
  
  
  TCanvas* c=(TCanvas*)gROOT->FindObject("PWA");
  if(c==NULL){
    c=new TCanvas("PWA","PWA Summary",10,10,1000,900);
    c->Divide(2,3);
  }
  else opt="P same";

gStyle->SetOptStat(0);
gStyle->SetMarkerSize(0.5);
gStyle->SetFillColor(0);
gStyle->SetPadColor(0);
gStyle->SetLineColor(kGray);


c->cd(1);
pwa->Draw("intens():_mass","","");
TGraph* g2=(TGraph*)gPad->FindObject("Graph");
g2->SetMarkerStyle(23);
g2->SetMarkerSize(0.5);g2->SetMarkerColor(color);
g2->SetTitle("Total Intensity");
g2->DrawClone(opt);

c->cd(2);
pwa->Draw("intens(\"flat\"):_mass","","");
TGraph* g2=(TGraph*)gPad->FindObject("Graph");
g2->SetMarkerStyle(23);
g2->SetMarkerSize(0.5);g2->SetMarkerColor(color);
g2->SetTitle("flat");
g2->DrawClone(opt);


c->cd(3);
pwa->Draw("intens(\"0++0-\"):err(\"0++0-\"):_mass","","goff");
TGraphErrors* g=new TGraphErrors(pwa->GetSelectedRows(),
				 pwa->GetV3(), // mass
				 pwa->GetV1(), // intensity
				 0,pwa->GetV2());
g->SetMarkerStyle(23);
g->SetMarkerSize(0.5);g->SetMarkerColor(color);
g->SetTitle("0++0-");
 g->Draw(opt);


c->cd(4);
 pwa->Draw("intens(\"0-+0+\"):err(\"0-+0+\"):_mass","","goff");
 TGraphErrors* g=new TGraphErrors(pwa->GetSelectedRows(),
 				 pwa->GetV3(), // mass
 				 pwa->GetV1(), // intensity
				 0,pwa->GetV2());
g->SetMarkerStyle(23);
g->SetMarkerSize(0.5);g->SetMarkerColor(color);
g->SetTitle("0-+0+");
g->Draw(opt);




c->cd(5);
pwa->Draw("intens(\"1++0+\"):err(\"1++0+\"):_mass","","goff");
TGraphErrors* g=new TGraphErrors(pwa->GetSelectedRows(),
				 pwa->GetV3(), // mass
				 pwa->GetV1(), // intensity
				 0,pwa->GetV2());
g->SetMarkerStyle(23);
g->SetMarkerSize(0.5);g->SetMarkerColor(color);
g->SetTitle("1++0+");
g->Draw(opt);



c->cd(6);
pwa->Draw("intens(\"2-+0+\"):err(\"2-+0+\"):_mass","","goff");
TGraphErrors* g=new TGraphErrors(pwa->GetSelectedRows(),
				 pwa->GetV3(), // mass
				 pwa->GetV1(), // intensity
				 0,pwa->GetV2());
g->SetMarkerStyle(23);
g->SetMarkerSize(0.5);g->SetMarkerColor(color);
g->SetTitle("2-+0+");
g->Draw(opt);


 TCanvas* c2=(TCanvas*)gROOT->FindObject("PWA2");
 if(c2==NULL){
   c2=new TCanvas("PWA2","PWA Summary",15,15,1000,900);
   c2->Divide(2,3);
   opt="AP";
 }
 else opt="P same";

 c2->cd(1);
 pwa->Draw("intens(\"2++0-\"):err(\"2++0-\"):_mass","","goff");
 TGraphErrors* g=new TGraphErrors(pwa->GetSelectedRows(),
 				 pwa->GetV3(), // mass
 				 pwa->GetV1(), // intensity
 				 0,pwa->GetV2());
 g->SetMarkerStyle(23);
 g->SetMarkerSize(0.5);g->SetMarkerColor(color);
 g->SetTitle("2++0-");
 g->Draw(opt);

c2->cd(2);
 pwa->Draw("intens(\"2++1+\"):err(\"2++1+\"):_mass","","goff");
 TGraphErrors* g=new TGraphErrors(pwa->GetSelectedRows(),
 				 pwa->GetV3(), // mass
 				 pwa->GetV1(), // intensity
 				 0,pwa->GetV2());
 g->SetMarkerStyle(23);
 g->SetMarkerSize(0.5);g->SetMarkerColor(color);
 g->SetTitle("2++1+");
 g->Draw(opt);


c2->cd(3);
pwa->Draw("intens(\"1-+0-\"):err(\"1-+0-\"):_mass","","goff");
TGraphErrors* g=new TGraphErrors(pwa->GetSelectedRows(),
				 pwa->GetV3(), // mass
				 pwa->GetV1(), // intensity
				 0,pwa->GetV2());
g->SetMarkerStyle(23);
g->SetMarkerSize(0.5);g->SetMarkerColor(color);
g->SetTitle("1-+0-");
g->Draw(opt);

c2->cd(4);
pwa->Draw("intens(\"1-+1+\"):err(\"1-+1+\"):_mass","","goff");
TGraphErrors* g=new TGraphErrors(pwa->GetSelectedRows(),
				 pwa->GetV3(), // mass
				 pwa->GetV1(), // intensity
				 0,pwa->GetV2());
g->SetMarkerStyle(23);
g->SetMarkerSize(0.5);g->SetMarkerColor(color);
g->SetTitle("1-+1+");
g->Draw(opt);


c2->cd(5);
pwa->Draw("intens(\"3++0+\"):err(\"3++0+\"):_mass","","goff");
TGraphErrors* g=new TGraphErrors(pwa->GetSelectedRows(),
				 pwa->GetV3(), // mass
				 pwa->GetV1(), // intensity
				 0,pwa->GetV2());
g->SetMarkerStyle(23);
g->SetMarkerSize(0.5);g->SetMarkerColor(color);
g->SetTitle("3++0+");
g->Draw(opt);


c2->cd(6);
pwa->Draw("intens(\"4++1+\"):err(\"4++1+\"):_mass","","goff");
TGraphErrors* g=new TGraphErrors(pwa->GetSelectedRows(),
				 pwa->GetV3(), // mass
				 pwa->GetV1(), // intensity
				 0,pwa->GetV2());
g->SetMarkerStyle(23);
g->SetMarkerSize(0.5);g->SetMarkerColor(color);
g->SetTitle("4++1+");
g->Draw(opt);

c2->Update();


TCanvas* c3=(TCanvas*)gROOT->FindObject("PWA3");
 if(c3==NULL){
   c3=new TCanvas("PWA3","PWA Summary",20,20,1000,900);
   c3->Divide(2,3);
   opt="AP";
 }
 else opt="P same";



c3->cd(1);
pwa->Draw("intens(\"3-+1+\"):err(\"3-+1+\"):_mass","","goff");
TGraphErrors* g=new TGraphErrors(pwa->GetSelectedRows(),
				 pwa->GetV3(), // mass
				 pwa->GetV1(), // intensity
				 0,pwa->GetV2());
g->SetMarkerStyle(23);
g->SetMarkerSize(0.5);g->SetMarkerColor(color);
g->SetTitle("3-+1+");
g->Draw(opt);

c3->cd(2);
pwa->Draw("intens(\"3-+1-\"):err(\"3-+1-\"):_mass","","goff");
TGraphErrors* g=new TGraphErrors(pwa->GetSelectedRows(),
				 pwa->GetV3(), // mass
				 pwa->GetV1(), // intensity
				 0,pwa->GetV2());
g->SetMarkerStyle(23);
g->SetMarkerSize(0.5);g->SetMarkerColor(color);
g->SetTitle("3-+1-");
g->Draw(opt);

c3->cd(3);
pwa->Draw("intens(\"3-+0-\"):err(\"3-+0-\"):_mass","","goff");
TGraphErrors* g=new TGraphErrors(pwa->GetSelectedRows(),
				 pwa->GetV3(), // mass
				 pwa->GetV1(), // intensity
				 0,pwa->GetV2());
g->SetMarkerStyle(23);
g->SetMarkerSize(0.5);g->SetMarkerColor(color);
g->SetTitle("3-+0-");
g->Draw(opt);

c3->Update();

}
