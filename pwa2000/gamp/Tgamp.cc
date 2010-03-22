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
//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#include "keyfile.h"
#include "particle.h"
#include "lorentz.h"
#include "Vec.h"
#include "event.h"
#include "wave.h"

#include "Tgamp.h"


using namespace std;


//extern int               keydebug;
extern particleDataTable PDGtable;
extern wave              gWave;


Tgamp::Tgamp(const string& pdgTableFileName)
  : _reflect       (false), _mirror (false),
    _suppressOutput(true)
{
  if (pdgTableFileName != "")
    PDGtable.initialize(pdgTableFileName.c_str());
  else
    PDGtable.initialize();
}


Tgamp::~Tgamp()
{
}


complex<double>
Tgamp::Amp(const string& keyFileName,
	   event&        ev) const
{
  keyfile key;
  key.open(keyFileName);

  complex<double> amplitude;
  if (_reflect)
    ev = reflectEvent(ev);
  if(_mirror)
    ev = mirrorEvent(ev);
  key.run(ev, amplitude, _suppressOutput);
  key.rewind();
  key.close();
  return amplitude;
}


complex<double>
Tgamp::Amp(const unsigned int iKey,
	   event&             ev) const
{
  if (iKey >= _keyFileNames.size()) {
    cerr << "invalid index " << iKey << "! returning (0,0)" << endl;
    return complex<double>(0,0);
  }

  return Amp(_keyFileNames[iKey], ev);
}


vector<complex<double> >
Tgamp::Amp(const string& keyFileName,
	   istream&      eventData,
	   const bool    testMode) const  // enables test mode
{
  if (!testMode) {
    vector<complex<double> > amplitudes;
    event                    ev;
    ev.setIOVersion(1);
    while (!(eventData >> ev).eof())
      amplitudes.push_back(Amp(keyFileName, ev));
    return amplitudes;
  } else {
    cout << "!!! test mode!" << endl;
    const int                debug = 0;
    vector<complex<double> > amplitudes;
    event                    ev;
    ev.setIOVersion(1);
    bool first = true;
    while (!(eventData >> ev).eof())
      if (first) {
	amplitudes.push_back(Amp(keyFileName, ev));
	first = false;
      } else {
	if (debug) {
	  cout << "@@Found a wave" << endl;
	  gWave.print();
	  cout << "@@Filling wave" << endl;
	}
	gWave.fill(ev, debug);
	if (debug) {
	  cout << "@@Wave before boosts" << endl;
	  gWave.print();
	}
	gWave.setupFrames(debug);
	if (debug) {
	  cout << "@@Wave after boosts" << endl;
	  gWave.print();
	}
	amplitudes.push_back(gWave.decayAmp(debug));
      }
    return amplitudes;
  }
}


vector<complex<double> >
Tgamp::Amp(const unsigned int iKey,
	   istream&           eventData) const
{
  if (iKey >= _keyFileNames.size()) {
    cerr << "invalid index " << iKey << "!" << endl;
    return vector<complex<double> >();
  }

  return Amp(_keyFileNames[iKey], eventData);
}


// performs reflection through production plane
event
Tgamp::reflectEvent(const event& evIn)
{
  // // get X Lorentz-vector in lab frame
//   const list<particle> daughtersLab = evIn.f_mesons();
//   fourVec              XLab;
//   for (list<particle>::const_iterator i = daughtersLab.begin(); i != daughtersLab.end(); ++i)
//     XLab += i->get4P();

//   // rotation into scattering plane
//   const fourVec          beamLab         = evIn.beam().get4P();
//   const threeVec         prodPlaneNormal = beamLab.V() / XLab.V();
//   static const double    piHalf          = acos(-1.0) / 2.0;
//   rotation               R(prodPlaneNormal.phi(),
//   			   prodPlaneNormal.theta() - piHalf,
//   			   -piHalf);
//   const lorentzTransform rotScatPlane(R);
//   fourVec                X    = XLab;
//   fourVec                beam = beamLab;
//   X    *= rotScatPlane;
//   beam *= rotScatPlane;

//   // boost into X-RF
//   const lorentzTransform boostRf(X);
//   X    *= boostRf;
//   beam *= boostRf;

//   // rotation so that beam is along z-axis
//   R.set(0, signof(beam.x()) * beam.V().theta(), 0);
//   const lorentzTransform rotBeam(R);

//   // total transformation
//   matrix<double>         transMatrix = rotBeam * (boostRf * rotScatPlane);
//   const lorentzTransform L           = lorentzTransform(transMatrix);
//   const lorentzTransform Linv(transMatrix.inv());
     
  // transform event into Gottfried-Jackson frame
  lorentzTransform L=toGottfriedJackson(evIn);
  lorentzTransform Linv(L.inv());
  event evOut = L * evIn;
    
  // reflect final state through production plane
  list<particle> daughtersRefl = evOut.f_mesons();
  for (list<particle>::iterator i = daughtersRefl.begin(); i != daughtersRefl.end(); ++i) {
    fourVec daughter = i->get4P();
    daughter.set(daughter.t(), daughter.x(), -daughter.y(), daughter.z());
    i->set4P(daughter);
  }
  evOut.set_f_mesons(daughtersRefl);
     
  // transform event back to lab frame
  evOut = Linv * evOut;
  
  return evOut;
}

// performs parity reflection
event
Tgamp::mirrorEvent(const event& evIn)
{

  lorentzTransform L=toGottfriedJackson(evIn);
  lorentzTransform Linv(L.inv());
  event evOut = L * evIn;
  
  list<particle> daughtersLab = evOut.f_mesons();
  for (list<particle>::iterator i = daughtersLab.begin(); i != daughtersLab.end(); ++i){
    fourVec daughter = i->get4P();
    daughter.set(daughter.t(), -daughter.x(), -daughter.y(), -daughter.z());
    i->set4P(daughter);
  }
  evOut.set_f_mesons(daughtersLab);
    
  // transform event back to lab frame
  evOut = Linv * evOut;
  return evOut;
}


lorentzTransform
Tgamp::toGottfriedJackson(const event& evIn){
// get X Lorentz-vector in lab frame
  const list<particle> daughtersLab = evIn.f_mesons();
  fourVec              XLab;
  for (list<particle>::const_iterator i = daughtersLab.begin(); i != daughtersLab.end(); ++i)
    XLab += i->get4P();

  // rotation into scattering plane
  const fourVec          beamLab         = evIn.beam().get4P();
  const threeVec         prodPlaneNormal = beamLab.V() / XLab.V();
  static const double    piHalf          = acos(-1.0) / 2.0;
  rotation               R(prodPlaneNormal.phi(),
  			   prodPlaneNormal.theta() - piHalf,
  			   -piHalf);
  const lorentzTransform rotScatPlane(R);
  fourVec                X    = XLab;
  fourVec                beam = beamLab;
  X    *= rotScatPlane;
  beam *= rotScatPlane;

  // boost into X-RF
  const lorentzTransform boostRf(X);
  X    *= boostRf;
  beam *= boostRf;

  // rotation so that beam is along z-axis
  R.set(0, signof(beam.x()) * beam.V().theta(), 0);
  const lorentzTransform rotBeam(R);

  // total transformation
  matrix<double>         transMatrix = rotBeam * (boostRf * rotScatPlane);
  return lorentzTransform(transMatrix);

}
