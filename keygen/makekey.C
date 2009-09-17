// Macro to create gamp key files for 5 (charged) pion events

{
  gROOT->ProcessLine(".L keygen.C+");
  
  partkey p1("pi-");
  partkey p2("pi+");
  partkey p3("pi-");
  partkey p4("pi+");
  partkey p5("pi-");
  
  //            name     iso1 iso2 l
  partkey rh1("rho(770)",&p1,&p2,1);
  
  partkey pi1("pi1(1600)",&p3,&rh1,1);
  partkey rh2("rho(770)",&p4,&p5,1);
  
 
  int l=1; 
  int s=2;

  int j=3;
  int m=1;
  int eps=1; // reflectivity
  int p=-1;
  
  if(eps!=1 && eps!=-1) { 
    std::cout<< "Reflectivity must be +-1 !!!" << std::endl;
    return;
  }
  
  partkey X("X",&rh2,&pi1,l,s);
  wavekey mykey(j,p,m,eps,&X);
  
  std::cout<<mykey.wavename(false)<<std::endl;
  X.write(std::cout,0);
  mykey.write();
  
  // test this wave with gamp:
  std::cout<<std::endl<<" TESTING THIS AMPLITUDE: "<<std::endl<<std::endl;
  TString com="cat testEvents.evt | gamp -P pdgTable.txt \"";
  com+=mykey.wavename(true);
  com+="\" | vamp";
  std::cout<<com<<std::endl;
  gSystem->Exec(com);

  std::cout<<std::endl<<" REFLECTIVITY = "<<eps<<" AMPLITUDE SHOULD BE "<< ((eps>0) ? "SAME" : "OPPOSITE") << " SIGN!"<<std::endl<<std::endl;
  com="cat testEvents.evt | reflect | gamp -P pdgTable.txt \"";
  com+=mykey.wavename(true);
  com+="\" | vamp";
  std::cout<<com<<std::endl;
  gSystem->Exec(com);
  
  std::cout<<"Ok?(y/n)"<<endl;
  char dummy;
  std::cin >> dummy;
  if(dummy=='n'){
    com="rm ";
    com+=mykey.wavename(true);
    gSystem->Exec(com);
    std::cout<<"WAVE Removed!"<<std::endl;
  }
  else {
    com="cp makekey.C scripts/";
    com+=mykey.wavename(true);
    com+=".C";
    gSystem->Exec(com);
    std::cout<<"WAVE Accepted!"<<std::endl;
  }

}

/*
f2(2340):       mass=2.339      width=0.319     0+(4++)
f2(2300):       mass=2.297      width=0.149     0+(4++)
f4(2050):       mass=2.034      width=0.222     0+(8++)
a4(2040):       mass=2.014      width=0.361     2-(8++)
f2(2010):       mass=2.011      width=0.202     0+(4++)
phi3(1850):     mass=1.854      width=0.087     0-(6--)
pi(1800):       mass=1.801      width=0.21      2-(0-+)
f0(1700):       mass=1.715      width=0.125     0+(0++)
rho(1700):      mass=1.7        width=0.24      2+(2--)
rho3(1690):     mass=1.691      width=0.161     2+(6--)
phi(1680):      mass=1.68       width=0.15      0-(2--)
pi2(1670):      mass=1.67       width=0.259     2-(4-+)
omega3(1670):   mass=1.667      width=0.168     0-(6--)
omega(1650):    mass=1.649      width=0.22      0-(2--)
f2'(1525):      mass=1.525      width=0.076     0+(4++)
f1(1510):       mass=1.51       width=0.073     0+(2++)
f0(1500):       mass=1.5        width=0.112     0+(0++)
rho(1450):      mass=1.465      width=0.31      2+(2--)
a0(1450):       mass=1.474      width=0.265     2-(0++)
etaH(1475):     mass=1.475      width=0.081     0+(0-+)
etaL(1405):     mass=1.405      width=0.056     0+(0-+)
eta(1440):      mass=1.42       width=0.06      0+(0-+)
omega(1420):    mass=1.419      width=0.174     0-(2--)
f1(1420):       mass=1.4263     width=0.0555    0+(2++)
f0(1370):       mass=1.35       width=0.35      0+(0++)
a2(1320):       mass=1.318      width=0.107     2-(4++)
pi(1300):       mass=1.3        width=0.4       2-(0-+)
eta(1295):      mass=1.297      width=0.053     0+(0-+)
f1(1285):       mass=1.2819     width=0.024     0+(2++)
f2(1270):       mass=1.2754     width=0.1851    0+(4++)
a1(1269):       mass=1.23       width=0.425     2-(2++)
b1(1235):       mass=1.229      width=0.142     2+(2+-)
h1(1170):       mass=1.17       width=0.36      0-(2+-)
phi(1020):      mass=1.01942    width=0.004468  0-(2--)
a0(980):        mass=0.9931     width=0.071     2-(0++)
f0(980):        mass=0.98       width=0.07      0+(0++)
eta'(958):      mass=0.95778    width=0.000202  0+(0-+)
omega(782):     mass=0.78257    width=0.00844   0-(2--)
rho(770):       mass=0.7693     width=0.1502    2+(2--)
sigma:  mass=0.8        width=0.8       0+(0++)
eta:    mass=0.5473     width=1.18e-06  0+(0-+)
pi0:    mass=0.134977   width=0 2-(0-+)
pi:     mass=0.13957    width=0 2-(0-+)
gamma:  mass=0  width=0 0-(2--)
e:      mass=0.00051    width=0 1-(1+-)
*/
