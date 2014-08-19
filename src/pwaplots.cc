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

// Reads in the results of N fits and creates intensity plots using
// the pwaPlotter class

#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>

int atoi ( const char * str );


#include "pwaPlotter.h"

using namespace std;
using namespace rpwa;


int
main(int argc, char** argv){

  if(argc<2){
    cerr<<"Usage: pwaplot nbins outputfile fit1 fit2 fit3 ..."<<endl;
    return 1;
  }

  unsigned int nbins=atoi(argv[1]);
  string outfilename=argv[2];
  vector<string> inputfiles;
  for(int i=3; i<argc; ++i){
    inputfiles.push_back(argv[i]);
  }

  pwaPlotter plotter;

  for(unsigned int i=0; i<inputfiles.size();++i){
    plotter.addFit(inputfiles[i],inputfiles[i],1,"pwa","fitResult_v2",nbins);
  }


  plotter.produceDensityPlots();
  plotter.printStats();
  plotter.writeAll(outfilename);

  return 0;
}
