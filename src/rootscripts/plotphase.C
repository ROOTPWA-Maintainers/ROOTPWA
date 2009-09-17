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
void plotphase(int i, int j, TString sel=""){


  TString com="phase(";
  com+=i;
  com+=",";
  com+=j;
  com+="):phaseErr(";
  com+=i;
  com+=",";
  com+=j;
  com+="):_mass";
  
  TFitBin* bin=new TFitBin();
  pwa->SetBranchAddress("fitbin",&bin);
  pwa->GetEntry(0);
  bin->listwaves();

  TString title("#Delta#Phi(");
  title+=bin->waveDesignator(i);
  title+=",";
  title+=bin->waveDesignator(j);
  title+=")";

  cout << title << endl;

  pwa->Draw(com.Data(),sel,"goff");
   TGraphErrors* g=new TGraphErrors(pwa->GetSelectedRows(),
                                   pwa->GetV3(), // mass
                                   pwa->GetV1(), // phase
                                   0,pwa->GetV2());

		 
  g->SetTitle(title);
  g->SetMarkerStyle(23);
  g->SetMarkerSize(0.5);
  g->SetMaximum(200);
  g->SetMinimum(-200);
  if(title!=0)g->SetTitle(title);
  g->DrawClone("AP");

  return;
}
