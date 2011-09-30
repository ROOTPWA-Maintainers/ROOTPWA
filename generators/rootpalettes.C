// workaround for Kpipi script production to be loaded by root

#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"
#include <iostream>

void setDiffColorStyle(const unsigned int NCont) {
  static Int_t *colors =0;
  static Bool_t initialized = kFALSE;

  Double_t Red[3]    = { 0.0, 1.0, 1.0 };
  Double_t Green[3]  = { 0.0, 1.0, 0.0 };
  Double_t Blue[3]   = { 1.0, 1.0, 0.0 };
  Double_t Length[3] = { 0.0, 0.50, 1.0 };

  if(!initialized){
    colors = new Int_t[NCont];
    Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,NCont);
    for (unsigned int i=0; i<NCont; i++) colors[i] = FI+i;
    initialized = kTRUE;
    return;
  }
  gStyle->SetPalette(NCont,colors);
}


void setNormalColorStyle(const unsigned int NCont) {
  static Int_t *colors =0;
  static Bool_t initialized = kFALSE;

  Double_t Red[3]    = { 0.0, 0.0, 1.0 };
  Double_t Green[3]  = { 0.0, 1.0, 0.0 };
  Double_t Blue[3]   = { 1.0, 0.0, 0.0 };
  Double_t Length[3] = { 0.0, 0.50, 1.0 };

  if(!initialized){
    colors = new Int_t[NCont];
    Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,NCont);
    for (unsigned int i=0; i<NCont; i++) colors[i] = FI+i;
    initialized = kTRUE;
    return;
  }
  gStyle->SetPalette(NCont,colors);
}

/*
void setNormalColorStyle() {
  static Int_t  colors2[NCont];
  static Bool_t initialized2 = kFALSE;

  Double_t Red[5]    = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t Green[5]  = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t Blue[5]   = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t Length[5] = { 0.51, 1.00, 0.12, 0.00, 0.00 };


  if(!initialized2){
    Int_t FI = TColor::CreateGradientColorTable(5,Length,Red,Green,Blue,NCont);
    for (int i=0; i<NCont; i++) colors2[i] = FI+i;
    initialized2 = kTRUE;
    return;
  }
  gStyle->SetPalette(NCont,colors2);
}

void setNormalColorStyle2() {
  static Int_t  colors2[NCont];
  static Bool_t initialized2 = kFALSE;

  Double_t Red[5]    = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t Green[5]  = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t Blue[5]   = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t Length[5] = { 0.51, 1.00, 0.12, 0.00, 0.00 };


  if(!initialized2){
    Int_t FI = TColor::CreateGradientColorTable(5,Length,Red,Green,Blue,NCont);
    for (int i=0; i<NCont; i++) colors2[i] = FI+i;
    initialized2 = kTRUE;
    return;
  }
  gStyle->SetPalette(NCont,colors2);
}*/

void rootpalettes() {
}
