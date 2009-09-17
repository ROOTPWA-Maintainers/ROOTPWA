#include <iostream>
#include <iomanip>
#include <fstream>

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TGenPhaseSpace.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TFile.h"


using namespace std;


/////////////////////////////////////////////////////////////////////////////////
// global constants

// particle masses
const double gProtonMass = 0.938272013;
const double gPionMass   = 0.13957018;
const double gPionMass2  = gPionMass * gPionMass;

// target position
const double gTargetZPos = -300.0; // [cm]

// beam parameters:
const double gBeamMomSigma = 1.2;  // [GeV/c]
const double gBeamMom      = 189;  // [GeV/c]
// 2004 beam:
// const double gBeamDxDz      = 0.00026; // tilt from Quirin was in mrad
// const double gBeamDxDzSigma = 0.00010;
// const double gBeamDyDz      = 0.00001; // tilt from Quirin was in mrad
// const double gBeamDyDzSigma = 0.00018;
// ideal beam:
const double gBeamDxDz      = 0.0;
const double gBeamDxDzSigma = 0.0;
const double gBeamDyDz      = 0.0;
const double gBeamDyDzSigma = 0.0;

// cut on t-distribution
const double tMin = 0.001;  // [(GeV/c)^2]


/////////////////////////////////////////////////////////////////////////////////
// helper functions

// constructs beam Lorentz vector
TLorentzVector
makeBeam()
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
writePwa2000Ascii(ostream&              out,
		  const TLorentzVector& beam,
		  TGenPhaseSpace&       event)
{
  if (!out) {
    cerr << "Output stream is not writable." << endl;
    return false;
  }
    
  const int geantIds[3] = { 8,  9,  9};
  const int charges[3]  = {+1, -1, -1};
  // total number of particles
  out << 4 << endl;
  // beam particle: geant ID, charge, p_x, p_y, p_z, E
  out << setprecision(numeric_limits<double>::digits10 + 1)
      << "9 -1 " << beam.Px() << " " << beam.Py() << " " << beam.Pz() << " " << beam.E() << endl;
  for (unsigned int i = 0; i < 3; ++i) {
    TLorentzVector* hadron = event.GetDecay(i);
    if (!hadron) {
      cerr << "genbod returns NULL pointer to Lorentz vector for daughter " << i << "." << endl;
      continue;
    }
    // hadron: geant ID, charge, p_x, p_y, p_z, E
    out << setprecision(numeric_limits<double>::digits10 + 1)
	<< geantIds[i] << " " << charges[i] << " " << hadron->Px() << " " << hadron->Py() << " " << hadron->Pz() << " " << hadron->E() << endl;
    }
  return true;
}


ostream&
progressIndicator(const long currentPos,
                  const long nmbTotal,
                  const int  nmbSteps = 10,
                  ostream&   out      = cout)
{
  const double step = nmbTotal / (double)nmbSteps;
  if ((nmbTotal >= 0) && ((int)(currentPos / step) - (int)((currentPos - 1) / step) != 0))
    out << "    " << setw(3) << (int)(currentPos / step) * nmbSteps << " %" << endl;
  return out;
}



/////////////////////////////////////////////////////////////////////////////////
// main routine
void
genPhaseSpaceData(const double   xMassMin          = 2.100,  // lower bound of mass bin [GeV/c^2]
		  const double   xMassMax          = 2.140,  // upper bound of mass bin [GeV/c^2]
		  const TString& outFileName       = "2100.2140.genbod.evt",
		  const TString& thetaHistFileName = "./hTheta.root",  // histogram with experimental distribution of scattering angle
		  const int      nmbEvent          = 2000,
		  const bool     plot              = false)
{
  Double_t daughterMasses[3] = {gPionMass, gPionMass, gPionMass};

  gRandom->SetSeed(12345);


  // setup histograms
  TH1D* ht;
  TH1D* hm;
  TH1D* hTheta;
  TH3D* hVertex3D;
  TH2D* hVertex;
  TH1D* hVz;
  TH1D* hE;
  if (plot) {
    ht        = new TH1D("ht", "t", 10000, -0.1, 1);
    hm        = new TH1D("hm", "3pi mass", 1000, 0.5, 2.5);
    hTheta    = new TH1D("hThetaGen", "cos theta", 100, 0.99985, 1);
    hVertex3D = new TH3D("hVertex3D", "Vertex xyz", 100, -2, 2, 100, -2, 2, 200, gTargetZPos - 5, gTargetZPos + 40);
    hVertex   = new TH2D("hVertex", "Vertex xy", 100, -2, 2, 100, -2, 2);
    hVz       = new TH1D("hVz","Vertex z", 1000, gTargetZPos - 40, gTargetZPos + 40);
    hE        = new TH1D("hE", "E", 100, 180, 200);
  }

  // open output file
  ofstream outFile(outFileName);
  cout << "Writing " << nmbEvent << " events to file '" << outFileName << "'." << endl;

  // get theta histogram
  TH1* thetaDist = NULL;
  {
    TFile* thetaHistFile = TFile::Open(thetaHistFileName, "READ");
    if (!thetaHistFile || thetaHistFile->IsZombie()) {
      cerr << "Cannot open histogram file '" << thetaHistFileName << "'. exiting." << endl;
      return;
    }
    thetaHistFile->GetObject("h1", thetaDist);
    if (!thetaDist) {
      cerr << "Cannot find theta histogram in file '" << thetaHistFileName << "'. exiting." << endl;
      return;
    }
  }

  int countEvent = 0;
  int attempts   = 0;
  int tenpercent = (int)(nmbEvent * 0.1);
  while (countEvent < nmbEvent) { // loop over events
    ++attempts;

    // construct primary vertex and beam
    const TVector3       vertexPos(0, 0, gTargetZPos);
    const TLorentzVector beam = makeBeam();

    // sample theta directly:
    const double theta  = thetaDist->GetRandom();
    const double xMass  = gRandom->Uniform(xMassMin, xMassMax);
    const double xMass2 = xMass * xMass;

    const double Ea = beam.E();
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
    p3.RotateUz(beam.Vect().Unit());

    // build resonance
    TLorentzVector X(p3, e);
    const TLorentzVector q = beam - X;

    // apply t cut
    const double tGen = -q.M2();
    if (tGen < tMin)
      continue;

    // generate phase space distribution
    TGenPhaseSpace phaseSpace;
    bool           allowed = phaseSpace.SetDecay(X, 3, daughterMasses);
    if (!allowed) {
      cerr << "Decay of M = " << X.M() << " into 3 pions is not allowed!" << endl;
      continue;
    }
    double maxWeight = phaseSpace.GetWtMax();
    double weight    = phaseSpace.Generate();
    if (weight / maxWeight < gRandom->Uniform())  // recjection sampling
      continue;

    // event is accepted
    ++countEvent;
    progressIndicator(countEvent, nmbEvent);
    if (plot) {
      ht->Fill(tGen);
      hm->Fill(X.M());
      hVertex3D->Fill(vertexPos.X(), vertexPos.Y(), vertexPos.Z());
      hVertex->Fill(vertexPos.X(), vertexPos.Y());
      hVz->Fill(vertexPos.Z());
      hE->Fill(e);
    }
     
    writePwa2000Ascii(outFile, beam, phaseSpace);
  }

  cout << endl << "Needed " << attempts << " attempts to generate " << countEvent << " events." << endl;
  outFile.close();

  if (plot) {
    TCanvas* c = new TCanvas("c", "c", 10, 10, 1000, 1000);
    c->Divide(3, 2);
    c->cd(1);
    hm->Draw();
    c->cd(2);
    ht->Draw();
    c->cd(3);
    hVz->Draw();
    c->cd(4);
    hVertex->Draw();
    c->cd(5);
    hE->Draw();
    c->cd(6);
    hVertex3D->Draw();
    c->Update();
  }
}




