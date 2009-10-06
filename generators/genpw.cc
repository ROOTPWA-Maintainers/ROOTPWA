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


/** @brief Simple partial wave event generator (homogeneous in m)
 */


#include <complex>
#include <iostream>
//#include <stringstream>

#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include "TPWWeight.h"
#include "TFile.h"
#include "TH1.h"
#include "TDiffractivePhaseSpace.h"
#include <event.h>

using namespace std;

extern particleDataTable PDGtable;

void printUsage(char* prog) {

}


int main(int argc, char** argv) {

  PDGtable.initialize();
  TPWWeight weighter;
  weighter.addWave("keyfiles/key3pi/SET1/1-1++0+rho770_01_pi-.key",
		   std::complex<double>(15.23,0),0);
  weighter.addWave("keyfiles/key3pi/SET1/1-2++1+rho770_21_pi-.key",
		   std::complex<double>(12.23,3.4),0);
  weighter.loadIntegrals("src/pwafitTest/amplitudes/norm.int");

  TDiffractivePhaseSpace difPS;
  difPS.SetBeam();
  difPS.SetTarget(-300,0.2,2);
  difPS.SetMassRange(1.2,1.4);			
  TFile* infile=TFile::Open(argv[1]);
  difPS.SetThetaDistribution((TH1*)infile->Get("h1"));
  const double mpi=0.13957018;
  difPS.AddDecayProduct(particleinfo(9,-1,mpi));  
  difPS.AddDecayProduct(particleinfo(8,1,mpi)); 
  difPS.AddDecayProduct(particleinfo(9,-1,mpi));


  for(unsigned int i=0;i<100;++i)
    {
      ofstream str("/tmp/event.evt");
      difPS.event(str);
      
      
      event e;
      e.setIOVersion(1);
      //str.close();
      
      ifstream istr("/tmp/event.evt");
      istr >> e;
      
      //e.print();
      
      double val=weighter.weight(e);
      
      
      cerr << val << endl;
      
      val+=1;
    }

  return 0;

}

