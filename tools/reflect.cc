
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


#include <complex>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <unistd.h>
#include <assert.h>
#include "event.h"
#include "particle.h"
#include "lorentz.h"
#include "Vec.h"

using namespace std;

#define MASS_PI 3.14159

void printUsage(char* prog) {
        cerr << "usage:" << endl;
        cerr << "  " << prog << " [-h] " << endl;
        cerr << "  -- reflects events" << endl;
}

int main(int argc, char** argv) {
  
  event e;
  threeVec p;
  threeVec prefl;
  fourVec p4refl;
  threeVec beam;
  threeVec N;
  list<particle> f_mesons;
  while(!(cin>>e).eof()) { // begin event loop
    // boost event into X-Rest frame
    f_mesons=e.f_mesons();
    fourVec pX;
    list<particle>::iterator it = f_mesons.begin();
    while (it != f_mesons.end() ) {
      pX+=it->get4P();
      ++it;
    }
    
    // boost into the Gottfried Jackson Frame
    lorentzTransform L, T;
    matrix < double >X (4, 4);
    rotation R;
    fourVec tempBeam, tempX;
    threeVec N;

    tempBeam=e.beam().get4P();
    tempX=pX;

    // Normal vector
    N = tempBeam.V () / tempX.V ();
    
    // rotate system into scattering plane
    T.set (R.set (N.phi (), N.theta () - MASS_PI / 2.0, - MASS_PI / 2.0));
    L = T;

    
    tempX *= T;
    tempBeam *= T;
      
    //tempX.print();
    //tempBeam.print();
    // boost to X rest frame
    T.set (tempX);
    X = T * L;
    L = lorentzTransform(X);
    tempX *= T;
    tempBeam *= T;
           
    // put beam along z
    // T.set (R.set (tempBeam.V ()));
    T.set ( R.set (0.0, signof(tempBeam.x())*tempBeam.V().theta(), 0.0 ) );
    X = T * L;
    L = lorentzTransform(X);
    lorentzTransform Linv(X.inv());
     
    // boost the event
     e=L*e;
    
    // fetch transformed particles and perform parity op
    fourVec pX2;
    f_mesons=e.f_mesons();
    N=e.mesonPlane();
    list<particle>::iterator it2 = f_mesons.begin();
    while (it2 != f_mesons.end() ) {
      //cout << it->Name() << endl;
      //it->print();
      p4refl=it2->get4P();
      //p=p4refl.V();
      //prefl=p-2*(N*p)*N; // reflection on scattering plane
      //p4refl.set(p4refl.t(),prefl);

      //assert(it->get4P().lenSq()-p4refl.lenSq() < 1E-10);
      p=p4refl.V();
      p.set(p.x(),-p.y(),p.z());
      p4refl.set(p4refl.t(),p);
      it2->set4P(p4refl);
      pX2+=p4refl;
      ++it2;
    }
     
    //pX2.print();
    //cerr << pX.len() << " | " << pX2.len() << endl; 
    
    e.set_f_mesons(f_mesons);
     
     // transform back to labframe
     e=Linv*e;

    cout << setprecision(10) << e;
  } // end event loop
  

  return 0;
}
