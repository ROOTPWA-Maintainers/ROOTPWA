void errors() {

//     gROOT->Reset();
//     gStyle->SetPalette(1);
   
    TCanvas *c1 = new TCanvas("c1","Fehler",1400,700);
    c1->Divide(2,1);
    c1_1->SetTopMargin(0.05);
    c1_1->SetRightMargin(0.05);
    c1_1->SetBottomMargin(0.05);
    c1_1->SetLeftMargin(0.05);
    c1_2->SetTopMargin(0.05);
    c1_2->SetRightMargin(0.05);
    c1_2->SetBottomMargin(0.05);
    c1_2->SetLeftMargin(0.05);


   float w=0.1; // Wichtung da sonst nix erkennbar

   TH1F *h1 = new TH1F("h1","Absolute Fehler",100,-2e-9,2e-9);

   for (int i=0;i<N;i++) {
     
     h1->Fill((h_C[i] - CPUpow[i]),w);

   }
   
   TH1F *h2 = new TH1F("h2","Relative Fehler",100,-0.2e-12,0.2e-12);

   for (int i=0;i<N;i++) {
     
      h2->Fill(1-(h_C[i] / CPUpow[i]));

   }
  
    h1->GetXaxis()->SetTitle("Absolute Fehler");
    h1->GetYaxis()->SetTitle("Ereignisse");
    h1->SetLineWidth(1);
    h1->SetLineColor(1);
    h1->SetLineStyle(1);
    h1->SetFillColor(1);
    h1->SetFillStyle(1005);
   
    h2->GetXaxis()->SetTitle("Relative Fehler");
    h2->GetYaxis()->SetTitle("Ereignisse");
    h2->SetLineWidth(1);
    h2->SetLineColor(1);
    h2->SetLineStyle(1);
    h2->SetFillColor(1);
    h2->SetFillStyle(1005);
   
   
    c1->cd(1);
    h1->Draw();
   
    c1->cd(2);
    h2->Draw();
   

}