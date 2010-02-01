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


#include <iomanip>
#include <fstream>

#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TFile.h"

#include "TDiffractivePhaseSpace.h"


using namespace std;
using namespace rpwa;


// /////////////////////////////////////////////////////////////////////////////////
// // global constants


// // target position
// const double _targetZPos = -300.0; // [cm]

// // beam parameters:
// const double _beamMomSigma = 1.2;  // [GeV/c]
// const double _beamMom      = 189;  // [GeV/c]
// // 2004 beam:
// // const double _beamDxDz      = 0.00026; // tilt from Quirin was in mrad
// // const double _beamDxDzSigma = 0.00010;
// // const double _beamDyDz      = 0.00001; // tilt from Quirin was in mrad
// // const double _beamDyDzSigma = 0.00018;
// // ideal beam:
// const double _beamDxDz      = 0.0;
// const double _beamDxDzSigma = 0.0;
// const double _beamDyDz      = 0.0;
// const double _beamDyDzSigma = 0.0;

// // cut on t-distribution
// const double _tMin = 0.001;  // [(GeV/c)^2]


TDiffractivePhaseSpace::TDiffractivePhaseSpace()
  : _tMin(0.001),
    _xMassMin(0),
    _xMassMax(0),
    _daughterMasses(NULL),
    _protonMass(0.938272013),
    _pionMass(0.13957018),
    _pionMass2(_pionMass * _pionMass)
{
  _phaseSpace.setWeightType    (nBodyPhaseSpaceGen::S_U_CHUNG);
  _phaseSpace.setKinematicsType(nBodyPhaseSpaceGen::BLOCK);
}


// /////////////////////////////////////////////////////////////////////////////////
// helper functions

// constructs beam Lorentz vector
TLorentzVector
TDiffractivePhaseSpace::makeBeam()
{
  // throw magnituide of beam momentum
  const double pBeam = gRandom->Gaus(_beamMom, _beamMomSigma);
  // throw beam inclination
  const double dxdz = gRandom->Gaus(_beamDxDz, _beamDxDzSigma);
  const double dydz = gRandom->Gaus(_beamDyDz, _beamDyDzSigma);
  // construct tilted beam momentum Lorentz vector
  const double pz    = pBeam / sqrt(1 + dxdz * dxdz + dydz * dydz);
  const double px    = dxdz * pz;
  const double py    = dydz * pz;
  const double EBeam = sqrt(pBeam * pBeam + _pionMass2);
  return TLorentzVector(px, py, pz, EBeam);
}


// writes event to ascii file read by gamp
bool
TDiffractivePhaseSpace::writePwa2000Ascii(ostream&              out,
					  const TLorentzVector& beam,
					  nBodyPhaseSpaceGen&   event)
{
  if (!out) {
    cerr << "Output stream is not writable." << endl;
    return false;
  }
    
   // total number of particles
  unsigned int nfspart=_decayProducts.size();
  out << nfspart + 1 << endl;
  // beam particle: geant ID, charge, p_x, p_y, p_z, E
  out << setprecision(numeric_limits<double>::digits10 + 1)
      << "9 -1 " << beam.Px() << " " << beam.Py() << " " << beam.Pz() << " " << beam.E() << endl;
  for (unsigned int i = 0; i < nfspart; ++i) {
    const TLorentzVector& hadron = event.daughter(i);
    // if (!hadron) {
//       cerr << "genbod returns NULL pointer to Lorentz vector for daughter " << i << "." << endl;
//       continue;
//     }
    // hadron: geant ID, charge, p_x, p_y, p_z, E
    out << setprecision(numeric_limits<double>::digits10 + 1)
	<< _decayProducts[i].gid << " " << _decayProducts[i].charge << " " << hadron.Px() << " " << hadron.Py() << " " << hadron.Pz() << " " << hadron.E() << endl;
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
TDiffractivePhaseSpace::SetSeed(int seed)
{
  gRandom->SetSeed(seed);
  _phaseSpace.setSeed(seed);
}


void 
TDiffractivePhaseSpace::SetDecayProducts(const vector<particleinfo>& info)
{
  _decayProducts.clear();
  _decayProducts=info;
  BuildDaughterList();
}


void
TDiffractivePhaseSpace::AddDecayProduct(const particleinfo& info)
{
  _decayProducts.push_back(info);
  BuildDaughterList();
}


void
TDiffractivePhaseSpace::BuildDaughterList()
{
  if(_daughterMasses!=NULL)delete[] _daughterMasses;
  _daughterMasses=new double[_decayProducts.size()];
  const unsigned int n   = _decayProducts.size();          
  for(unsigned int i=0;i<n;++i){
    _daughterMasses[i]=_decayProducts[i].mass;
  }
  if(n>2){
    _phaseSpace.setDecay((int)n, _daughterMasses);
    if(_xMassMax==0){
      cerr << "TDiffractivePhaseSpace::Please set Mass Range before Decay Products!" << endl;
      throw;
    }
    else {

      cerr << "Calculating max wheight ("<<n<<" fs particles) for m="<<_xMassMax<< endl;
      _phaseSpace.setMaxWeight(1.01 * _phaseSpace.estimateMaxWeight(_xMassMax,1000000));
      cerr << "Max weight:" << _phaseSpace.maxWeight() << endl;

    }
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
//   Double_t daughterMasses[3] = {_pionMass, _pionMass, _pionMass};

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
//     hVertex3D = new TH3D("hVertex3D", "Vertex xyz", 100, -2, 2, 100, -2, 2, 200, _targetZPos - 5, _targetZPos + 40);
//     hVertex   = new TH2D("hVertex", "Vertex xy", 100, -2, 2, 100, -2, 2);
//     hVz       = new TH1D("hVz","Vertex z", 1000, _targetZPos - 40, _targetZPos + 40);
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
TDiffractivePhaseSpace::event(TLorentzVector& beamresult)
{
  // construct primary vertex and beam
  const TVector3 vertexPos(0, 0, 
			   _targetZPos+gRandom->Uniform(-_targetZLength*0.5,
				 			 _targetZLength*0.5));
  _beam = makeBeam();
  
  bool done=false;
  unsigned int attempts=0;
  while(!done){
    // sample theta directly:
    //const double theta  = thetaDistribution->GetRandom();

    const double tprime = gRandom->Exp(_invSlopePar);
    const double xMass  = gRandom->Uniform(_xMassMin, _xMassMax);
    const double xMass2 = xMass * xMass;
    const double xMass4 = xMass2 * xMass2;
    const double m0  = _recoilMass;
    const double ma = _beam.M();
    const double pa= _beam.Vect().Mag();
    const double ma2 = ma * ma;
    const double ma4 = ma2 * ma2;

    const double Ea = _beam.E();
    const double Ea2 = Ea*Ea;
    const double g=2*Ea*m0 + ma2 - xMass2;
    const double EcN=4*Ea2*(g - tprime)
                     - ma4 + 2*ma2*xMass2 - xMass4 ;
    const double EcD=4*Ea*g;
    const double Ec=EcN/EcD;

    // cerr << xMass << endl;
    //cerr << tprime << endl;
    //cerr << "Ea=" << Ea << "    Ec=" << Ec << endl;


    // // account for recoil assume proton recoil
//     const double eps  = (xMass2 - _pionMass2) / (2 * _protonMass * Ea);
//     const double epso = 1 - eps;
//     const double eat  = Ea * theta;
//     //eat *= eat;
//     const double mpeps = _protonMass * eps;
//     const double e = Ea - 1 / (2 * _protonMass * epso) * (eat * eat + mpeps * mpeps);


    const double pw = sqrt(Ec * Ec - xMass2); // three momentum magnitude
    const double t=tprime-(pw*pw+pa*pa-2.*pw*pa);
    // calculate theta following suh-urks paper (formulas 20 and 21)
    const double term1=(xMass2-ma2)/(2*Ec);
    const double term=t-term1*term1;
    if(term<0) {
      //cout << "neg" << endl;
      continue;
     }

    const double theta=1/Ec*sqrt(term);
    const double phi= gRandom->Uniform(0., TMath::TwoPi());
    
    TVector3 p3(1, 0, 0);
    p3.SetMagThetaPhi(pw, theta, phi);
    // rotate to beamdirection:
    p3.RotateUz(_beam.Vect().Unit());
    
    // build resonance
    TLorentzVector X(p3, Ec);
    const TLorentzVector q = _beam - X;
    
    // apply t cut
    const double tGen = -q.M2();
    if (tGen < _tMin){
      //cerr << "tGen < _tMin " << endl;
      continue;
    }
    
    // generate n-body phase space for X system
    ++attempts;
    {
      const double xMass = X.M();
      _phaseSpace.pickMasses(xMass);
      // correct weight for phase space splitting
      const double maxPsWeight = _phaseSpace.maxWeight()  * _xMassMax;
      const double psWeight    = _phaseSpace.calcWeight() * xMass;
      if ((psWeight / maxPsWeight) < _phaseSpace.random())
    	continue;
      _phaseSpace.pickAngles();
      _phaseSpace.calcEventKinematics(X);
    }

    done=true;
    beamresult=_beam;
  } // end while !done
  // event is accepted
   
  return attempts;
}


unsigned int 
TDiffractivePhaseSpace::event(ostream& stream)
{
  unsigned int attempts=event(_beam);
  writePwa2000Ascii(stream, _beam, _phaseSpace);
  return attempts;
}


void 
TDiffractivePhaseSpace::SetBeam(double Mom,  double MomSigma,
				double DxDz, double DxDzSigma,
				double DyDz, double DyDzSigma)
{
  _beamMomSigma=MomSigma;
  _beamMom=Mom;
  _beamDxDz=DxDz;
  _beamDxDzSigma=DxDzSigma;
  _beamDyDz=DyDz;
  _beamDyDzSigma=DyDzSigma;
}

