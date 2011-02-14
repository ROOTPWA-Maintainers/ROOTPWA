{

double scale=1;

THStack* stack=new THStack("hs","Meson States");
pdg->Draw("_mass/1000>>h1(50,0,5)","_mass<5000 && _status>=3 && _R==0 && (_I==1 || _I==-1) && (_G==-1)","GOFF");
  h1->SetFillColor(kGreen-6);
h1->Scale(scale);
stack->Add(h1);

pdg->Draw("_mass/1000>>h3(50,0,5)","_mass<5000 && _status==2 && _R==0 && (_I==1 || _I==-1) && (_G==-1 || _G==0)","GOFF");
  h3->SetFillColor(kAzure-4);
h3->Scale(scale);
stack->Add(h3);

pdg->Draw("_mass/1000>>h4(50,0,5)","_mass<5000 && _status==1 && _R==0 && (_I==1 || _I==-1) && (_G==-1 || _G==0)","GOFF");
  h4->SetFillColor(kRed-7);
h4->Scale(scale);
stack->Add(h4);



//stack->Draw("");


pdg->Draw("(_mass/1000)^2:_J","_mass<3000 && _R==0 && (_I==1 || _I==-1 || _I==0) ","GOFF");



}
