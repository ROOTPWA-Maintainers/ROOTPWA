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

// This Class' Header ------------------
#include "TPDGDB.h"

// C/C++ Headers ----------------------
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;
using std::ifstream;

// Collaborating Class Headers --------
#include "TString.h"
#include "TPDGEntry.h"
#include "TTree.h"
#include "TFile.h"

// Class Member definitions -----------

ClassImp(TPDGDB)


void 
TPDGDB::Draw(char* com, char* sel, char* opt, int n, int s){
  if(_tree!=NULL)_tree->Draw(com,sel,opt,n,s);
}

unsigned int 
TPDGDB::read(const TString& filename, int num){
  ifstream infile(filename.Data());

  TString outname=filename+".root";
  TFile* outfile=TFile::Open(outname,"RECREATE");
  _tree=new TTree("pdg","pdg");

  TPDGEntry* entry=new TPDGEntry();
  _tree->Branch("TPDGEntry",&entry);

  
  int counter=0;
  char line[500];
  while(infile.good()){
    
    // strip comment lines
    if(infile.peek()=='*'){
      cout << "Stripping comment" << endl;
      
      infile.getline(line,500);
      continue;
    }
    
    
    infile >> (*entry) ;
    entry->Print();
    

    _tree->Fill();

    
    // discard rest of line
    infile.getline(line,500);
    if(num>0 && counter++>=num) break;
  }
  _tree->Write();
  outfile->Close();
  return 0;//_tree->GetEntries();
  
}
