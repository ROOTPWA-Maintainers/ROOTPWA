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


#include "TDiffractivePhaseSpace.h"
#include <iomanip>
#include <fstream>

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TFile.h"


using namespace std;


// /////////////////////////////////////////////////////////////////////////////////
// // global constants


// // target position
// const double gTargetZPos = -300.0; // [cm]

// // beam parameters:
// const double gBeamMomSigma = 1.2;  // [GeV/c]
// const double gBeamMom      = 189;  // [GeV/c]
// // 2004 beam:
// // const double gBeamDxDz      = 0.00026; // tilt from Quirin was in mrad
// // const double gBeamDxDzSigma = 0.00010;
// // const double gBeamDyDz      = 0.00001; // tilt from Quirin was in mrad
// // const double gBeamDyDzSigma = 0.00018;
// // ideal beam:
// const double gBeamDxDz      = 0.0;
// const double gBeamDxDzSigma = 0.0;
// const double gBeamDyDz      = 0.0;
// const double gBeamDyDzSigma = 0.0;

// // cut on t-distribution
// const double tMin = 0.001;  // [(GeV/c)^2]



TDiffractivePhaseSpace::TDiffractivePhaseSpace() :
  tMin(0.001), daughterMasses(NULL),gProtonMass(0.938272013),gPionMass(0.13957018),gPionMass2(gPionMass * gPionMass)
{}



// /////////////////////////////////////////////////////////////////////////////////
// helper functions

// constructs beam Lorentz vector
TLorentzVector
TDiffractivePhaseSpace::makeBeam()
{
  // throw magnituide of beam momentum
  const double pBeam = gRandom->Gaus(gBeamMom, gBeamMomSigma);
  // throw beam inclination
  const double dxdz = gRandom->Gaus(gBeamDxDz, gBeamDxDzSigma);
  const double dydz = gRandom->Gaus(gBeamDyDz, gBeamDyDzSigma);
  // construct tilted beam momentum Lorentz vector
  const double pz    = pBeam / sqrt(1 + dxdz * dxdz + dydz * dydz);
  const double px    = dxdz * pz;
  const double py    = dydz * pz;
  const double EBeam = sqrt(pBeam * pBeam + gPionMass2);
  return TLorentzVector(px, py, pz, EBeam);
}


// writes event to ascii file read by gamp
bool
TDiffractivePhaseSpace::writePwa2000Ascii(ostream&              out,
					  const TLorentzVector& beam,
					  TGenPhaseSpace&       event)
{
  if (!out) {
    cerr << "Output stream is not writable." << endl;
    return false;
  }
    
   // total number of particles
  unsigned int nfspart=decayProducts.size();
  out << nfspart + 1 << endl;
  // beam particle: geant ID, charge, p_x, p_y, p_z, E
  out << setprecision(numeric_limits<double>::digits10 + 1)
      << "9 -1 " << beam.Px() << " " << beam.Py() << " " << beam.Pz() << " " << beam.E() << endl;
  for (unsigned int i = 0; i < nfspart; ++i) {
    TLorentzVector* hadron = event.GetDecay(i);
    if (!hadron) {
      cerr << "genbod returns NULL pointer to Lorentz vector for daughter " << i << "." << endl;
      continue;
    }
    // hadron: geant ID, charge, p_x, p_y, p_z, E
    out << setprecision(numeric_limits<double>::digits10 + 1)
	<< decayProducts[i].gid << " " << decayProducts[i].charge << " " << hadron->Px() << " " << hadron->Py() << " " << hadron->Pz() << " " << hadron->E() << endl;
    }
  return true;
}


ostream&
TDiffractivePhaseSpace::progressIndicator(const long currentPos,
					  const long nmbTotal,
					  const int  nmbSteps,
					  ostream&   out)
{
  const double step = nmbTotal / (double)nmbSteps;
  if ((nmbTotal >= 0) && ((int)(currentPos / step) - (int)((currentPos - 1) / step) != 0))
    out << "    " << setw(3) << (int)(currentPos / step) * nmbSteps << " %" << endl;
  return out;
}

void 
TDiffractivePhaseSpace::SetSeed(int seed){
  gRandom->SetSeed(seed);
}


void 
TDiffractivePhaseSpace::SetDecayProducts(const std::vector<particleinfo>& info)
{
  decayProducts.clear();
  decayProducts=info;
  BuildDaughterList();
}

void
TDiffractivePhaseSpace::AddDecayProduct(const particleinfo& info){
  decayProducts.push_back(info);
  BuildDaughterList();
}

void
TDiffractivePhaseSpace::BuildDaughterList(){
  if(daughterMasses!=NULL)delete[] daughterMasses;
  daughterMasses=new double[decayProducts.size()];
  for(unsigned int i=0;i<decayProducts.size();++i){
    daughterMasses[i]=decayProducts[i].mass;
  }
}


/////////////////////////////////////////////////////////////////////////////////
// main routine
// void
// TDiffractivePhaseSpace::genPhaseSpaceData(const double   xMassMin          = 2.100,  // lower bound of mass bin [GeV/c^2]
// 		  const double   xMassMax          = 2.140,  // upper bound of mass bin [GeV/c^2]
// 		  const TString& outFileName       = "2100.2140.genbod.evt",
// 		  const TString& thetaHistFileName = "./hTheta.root",  // histogram with experimental distribution of scattering angle
// 		  const int      nmbEvent          = 2000,
// 		  const bool     plot              = false)
// {
//   Double_t daughterMasses[3] = {gPionMass, gPionMass, gPionMass};

//   gRandom->SetSeed(12345);


//   // setup histograms
//   TH1D* ht;
//   TH1D* hm;
//   TH1D* hTheta;
//   TH3D* hVertex3D;
//   TH2D* hVertex;
//   TH1D* hVz;
//   TH1D* hE;
//   if (plot) {
//     ht        = new TH1D("ht", "t", 10000, -0.1, 1);
//     hm        = new TH1D("hm", "3pi mass", 1000, 0.5, 2.5);
//     hTheta    = new TH1D("hThetaGen", "cos theta", 100, 0.99985, 1);
//     hVertex3D = new TH3D("hVertex3D", "Vertex xyz", 100, -2, 2, 100, -2, 2, 200, gTargetZPos - 5, gTargetZPos + 40);
//     hVertex   = new TH2D("hVertex", "Vertex xy", 100, -2, 2, 100, -2, 2);
//     hVz       = new TH1D("hVz","Vertex z", 1000, gTargetZPos - 40, gTargetZPos + 40);
//     hE        = new TH1D("hE", "E", 100, 180, 200);
//   }

//   // open output file
//   ofstream outFile(outFileName);
//   cout << "Writing " << nmbEvent << " events to file '" << outFileName << "'." << endl;

//   // get theta histogram
//   TH1* thetaDist = NULL;
//   {
//     TFile* thetaHistFile = TFile::Open(thetaHistFileName, "READ");
//     if (!thetaHistFile || thetaHistFile->IsZombie()) {
//       cerr << "Cannot open histogram file '" << thetaHistFileName << "'. exiting." << endl;
//       return;
//     }
//     thetaHistFile->GetObject("h1", thetaDist);
//     if (!thetaDist) {
//       cerr << "Cannot find theta histogram in file '" << thetaHistFileName << "'. exiting." << endl;
//       return;
//     }
//   }

//   int countEvent = 0;
//   int attempts   = 0;
//   int tenpercent = (int)(nmbEvent * 0.1);
//   while (countEvent < nmbEvent) { // loop over events
//     ++attempts;


//   }
// }






unsigned int 
TDiffractivePhaseSpace::event(TLorentzVector& beamresult){
  
  
    // construct primary vertex and beam
  const TVector3 vertexPos(0, 0, 
			   gTargetZPos+gRandom->Uniform(-gTargetZLength*0.5,
							gTargetZLength*0.5));
  gbeam = makeBeam();
  
  bool done=false;
  unsigned int attempts=0;
  while(!done){
    // sample theta directly:
    const double theta  = thetaDistribution->GetRandom();
    const double xMass  = gRandom->Uniform(xMassMin, xMassMax);
    const double xMass2 = xMass * xMass;
    
    const double Ea = gbeam.E();
    // account for recoil assume proton recoil
    const double eps  = (xMass2 - gPionMass2) / (2 * gProtonMass * Ea);
    const double epso = 1 - eps;
    const double eat  = Ea * theta;
    //eat *= eat;
    const double mpeps = gProtonMass * eps;
    const double e = Ea - 1 / (2 * gProtonMass * epso) * (eat * eat + mpeps * mpeps);
    const double pw = sqrt(e * e - xMass2); // three momentum magnitude
    const double phi= gRandom->Uniform(0., TMath::TwoPi());
    
    TVector3 p3(1, 0, 0);
    p3.SetMagThetaPhi(pw, theta, phi);
    // rotate to beamdirection:
    p3.RotateUz(gbeam.Vect().Unit());
    
    // build resonance
    TLorentzVector X(p3, e);
    const TLorentzVector q = gbeam - X;
    
    // apply t cut
    const double tGen = -q.M2();
    if (tGen < tMin)
      continue;
    
    // generate phase space distribution
    bool           allowed = phaseSpace.SetDecay(X, decayProducts.size(), daughterMasses,"");
    if (!allowed) {
      cerr << "Decay of M = " << X.M() << " into this final state is kinematically not allowed!" << endl;
      continue;
    }
    
    double weight    = phaseSpace.Generate();
    double maxWeight = phaseSpace.GetWtMax();
    //cerr << "wMax=" << maxWeight << endl;
    //cerr << "w=" << weight << endl;
    double d=gRandom->Uniform();
    //cerr << "d="  << d << endl;
    ++attempts;
    if ( (weight/maxWeight) < d){  // recjection sampling
      continue;
    }
    
    done=true;
    beamresult=gbeam;
  }
  // event is accepted
  //writePwa2000Ascii(outFile, beam, phaseSpace);
  //cerr << attempts << " attempts needed for Phase Space" << endl;
 
  return attempts;
}

unsigned int 
TDiffractivePhaseSpace::event(ostream& stream){
  unsigned int attempts=event(gbeam);
  writePwa2000Ascii(stream, gbeam, phaseSpace);
  return attempts;
}


void 
TDiffractivePhaseSpace::SetBeam(double Mom, double MomSigma,
	double DxDz, double DxDzSigma,
	double DyDz, double DyDzSigma){
  
  gBeamMomSigma=MomSigma;
  gBeamMom=Mom;
  gBeamDxDz=DxDz;
  gBeamDxDzSigma=DxDzSigma;
  gBeamDyDz=DyDz;
  gBeamDyDzSigma=DyDzSigma;
}

