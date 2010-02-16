
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

#include "TSystem.h"

using namespace std;


void printUsage(char* prog) {
  cerr << "Bins evt data in mass bins" << endl; 
  cerr << "options:" << endl;
  cerr << "-n    number of bins" << endl;
  cerr << "-s    mass of start bin" << endl;
  cerr << "-b    bin width in mass" << endl;
  cerr << "-o    output path" << endl;
 
}

int main(int argc, char** argv) {
  
  if(argc<2){
    printUsage(argv[0]);
    return 1;
  }

  TString path;
  double mstart;
  double mbin;
  int nbins;

  int c;
  while ((c = getopt(argc, argv, "n:s:b:o:h")) != -1)
    switch (c) {
    case 'n':
      nbins = atoi(optarg);
      break;
    case 's':
      mstart = atof(optarg);
      break;
    case 'o':
      path = optarg;
      break;
   case 'b':
      mbin = atof(optarg);
      break;
   case 'h':
     printUsage(argv[0]);
     return 1;
    }                                          
  
  
 
  unsigned int nselected=0;
  double m=mstart;

  vector<std::ofstream*> outfiles;// setup output file

  // open outfiles
  for(int ibin=0;ibin<nbins;++ibin){

    double mmin=m;
    double mmax=m+mbin;

    TString binS;
    binS+=mmin*1000;
    binS+=".";
    binS+=mmax*1000;
    // create directory
    TString com;
    com+=path;
    com+=binS;
    com.ReplaceAll(" ","");
    if(!gSystem->mkdir(com.Data()))
      std::cout<<"Directory "
               <<com
               <<" could not be created. Already existent?"<<std::endl;

    TString outfile=path;
    outfile+=binS;
    outfile+="/";
    outfile+=binS;
    outfile+=".evt";
    outfile.ReplaceAll(" ","");
    outfiles.push_back(new std::ofstream(outfile.Data()));
    m=mmax; // step further in m
  }




  list<particle> f_mesons;
  event e;
  while(!(cin>>e).eof()) { // begin event loop
    f_mesons=e.f_mesons();
    fourVec pX;
    fourVec p;
    list<particle>::iterator it = f_mesons.begin();
    while (it != f_mesons.end() ) {
      pX=it->get4P();
      p+=pX;
      ++it;
    }
    // sort by mass
    double m=p.len();
     // select appropriate mass bin
    unsigned int bin=(unsigned int)floor((m-mstart)/mbin);
    if(bin<outfiles.size()){
      std::ofstream* out=outfiles.at(bin);
      ++nselected;
      *out << e ;
    }



  }

 
  return 0;
}
