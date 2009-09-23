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
#include "../TFitBin.h"
#include "../TCMatrix.h"
#include "TMatrixD.h"
#include "TString.h"
#include <vector>
#include <iostream>
#include <map>
#include "TComplex.h"


using namespace std;


void testFitBin(){
  TFitBin bin;
  double mass=1000;

  vector<TString> wavenames; // contains rank information 
  //wavenames.push_back("V0_B");
  wavenames.push_back("V0_A");
  wavenames.push_back("V1_A");
  wavenames.push_back("V_flat");

  double logli=10;
  vector<TComplex> amps;

  //amps.push_back(TComplex(1,0));
  amps.push_back(TComplex(2,0));
  amps.push_back(TComplex(1.4,0.6));
  amps.push_back(TComplex(1.2,0));

  int npar=5;
  TMatrixD errm(npar,npar);
  for(int i=0;i<npar;++i){
    for(int j=i;j<npar;++j){
      errm[i][j]=0.1;
      errm[j][i]=0.1;
      //if(i==j)cout << "Sig"<< i << "=" << sqrt(errm[i][i]) << endl;
    }
  }
  errm.Print();
  //cout << " Created error matrix." << endl;
 
 
  vector<pair<int,int> > indices; 
  indices.push_back(pair<int,int>(0,1));
  //indices.push_back(pair<int,int>(1,2));
  indices.push_back(pair<int,int>(2,3));
  indices.push_back(pair<int,int>(4,-1));

  TCMatrix integr(1,1);
  integr.set(0,0,2);
  //integr.set(1,1,2);
 
  cout << "filling TFitBin" << endl;
  cout << "npar=" << npar << endl;
  cout << "Number of amps=" << amps.size()<< endl;
  cout << "Number of indices=" << indices.size() << endl;
  cout << "Number of wavenames=" << wavenames.size()<< endl;

  cout << "Dimension of errormatrix=" << errm.GetNrows() << endl;
  cout << "Dimension of integral matrix=" << integr.nrows() << endl;

  bin.fill(amps,
	   indices,
	   wavenames,
	   1,
	   1,
	   mass,
	   integr,
	   errm,
	   logli,
	   2);
  bin.PrintParameters();

  cout << "Error of A: " << bin.err("A") << endl; 

  return;
}
