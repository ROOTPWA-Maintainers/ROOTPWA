
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
#include <string>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include <unistd.h>
#include <assert.h>
#include "event.h"
#include "particle.h"
#include "lorentz.h"
#include "Vec.h"

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

using namespace std;


void printUsage(char* prog) {
  cerr << "Converts pwa2000 evt format to ROOT tree" << endl; 
  cerr << "usage:" << endl;
  cerr << "cat myevents.evt | " << prog << " myevents.root " << endl;
 
}

int main(int argc, char** argv) {
  
  if(argc<2){
    printUsage(argv[0]);
    return 1;
  }

  event e;


  threeVec evtbeam;
  list<particle> f_mesons;
  
  string outfilename(argv[1]);

  TFile* outfile=TFile::Open(outfilename.c_str(),"RECREATE");
  if(outfile==NULL){
    cerr << "Could not open output file '"<<outfilename << "'. Aborting." << endl;
    return 1;
  }
  TTree* outtree=new TTree("events","events");
  TClonesArray* p=new TClonesArray("TLorentzVector");
  TLorentzVector beam;
  double qbeam;
  std::vector<int> q; // array of charges

  outtree->Branch("p",&p);
  outtree->Branch("beam",&beam);
  outtree->Branch("q",&q);
  outtree->Branch("qbeam",&qbeam,"qbeam/i");


  while(!(cin>>e).eof()) { // begin event loop
    p->Delete(); // clear output arrays
    q.clear();
    f_mesons=e.f_mesons();
    fourVec pX;
    list<particle>::iterator it = f_mesons.begin();
    while (it != f_mesons.end() ) {
      pX=it->get4P();
      new ((*p)[p->GetEntries()]) TLorentzVector(pX.x(),pX.y(),pX.z(),pX.t());
      q.push_back(it->Charge());
      ++it;
    }
    fourVec evtbeam=e.beam().get4P();
    beam.SetPxPyPzE(evtbeam.x(),evtbeam.y(),evtbeam.z(),evtbeam.t()); 
    qbeam=e.beam().Charge();
    outtree->Fill();
  }

  outtree->Write();
  outfile->Close();
  return 0;
}
