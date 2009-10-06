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


  event e;
  e.setIOVersion(1);

  ofstream str("/tmp/event.evt");
  str << "4\n"
      << "9 -1 0. 0. 191.0323749689385 191.0324259546224 \n"  
      << "8 1 0.1776415401059793 -0.9281616207073099 100.6615453374411 100.6660778518155\n"
      << "9 -1 -0.2910619303999746 0.5354572747041667 44.8237555788041 44.82811590468507\n"
      << "9 -1 0.1670942931439715 0.4716546614904428 45.53027444487987 45.53323785416246\n";

  str.close();
  ifstream istr("/tmp/event.evt");
  istr >> e;

  //e.print();

  double val=weighter.weight(e);
  
  
  cerr << val << endl;

  val+=1;

  return 0;

}

