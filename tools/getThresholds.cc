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

#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <exception>
#include <string>
#include "TString.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TFile.h"
#include "integral.h"



bool SameJPCEPS(TString s1, TString s2){
  TString ss1=s1(2,5);
  TString ss2=s2(2,5);
  //cout << ss1 << " " << ss2 << endl;
  return ss1==ss2;
}




int
main(int argc, char** argv){

  // read in config file
  ifstream config(argv[1]);
  vector<integral> integrals;
  vector<double> masses;

  while(config.good() && !config.eof()){
    double mass;
    TString normfilename;
    config >> mass >> normfilename;
    ifstream file(normfilename.Data());
    integrals.push_back(integral());
    integrals.back().scan(file);
    masses.push_back(mass);
    cout << mass << endl;
    //getline(0);
    config.ignore(255,'\n');
  }

  list<string> waves=integrals[0].files();
  //unsigned int nwaves=waves.size();
  list<string>::iterator it=waves.begin();
  // wavemap: contains one entry per wave
  // each entry records mass and maximum integral value
  map<TString,pair<double,double> > wavemap; 
  map<TString,TGraph*> graphs;
  TMultiGraph* mg=new TMultiGraph();
  while(it!=waves.end()){
    wavemap[*it]=pair<double,double>(0,0);
    graphs[*it]=new TGraph(masses.size());
    graphs[*it]->SetName(it->c_str());
    mg->Add(graphs[*it]);
    ++it;
  }

  map<TString,TGraph*> jpcgraphs;
  TMultiGraph* mg2=new TMultiGraph();
  //jpcgraphs["0-+0+"]=new TGraph(1000);
  //jpcgraphs["0-+0+"]->SetName("g0-+0+");
  //mg2->Add(jpcgraphs["0-+0+"]);


   map<TString,pair<double,double> >::iterator it2=wavemap.begin(); 
    map<TString,pair<double,double> >::iterator it3=wavemap.begin(); 

  // loop over mass bins
    // unsigned int k=0;
  for(unsigned int i=0; i<masses.size(); ++i){
    cerr << "Processing mass bin "<< masses[i] << endl;
    it2=wavemap.begin(); 
    while(it2!=wavemap.end()){
      it3=it2;
      while(it3!=wavemap.end()){
	try{
	  double val=integrals[i].val(it2->first.Data(),it3->first.Data()).real();
	  if(it2==it3){
	    graphs[it2->first]->SetPoint(i,masses[i],val);
	    if(it2->second.second<val){
	      it2->second=pair<double,double>(masses[i],val);
	    }
	  }
	  else if(SameJPCEPS(it2->first,it3->first)){
	    TString name("g");
	    name.Append(it2->first);
	    name.Append(it3->first);
	    if(jpcgraphs[name]==0){
	      jpcgraphs[name]=new TGraph(masses.size());
	      jpcgraphs[name]->SetName(name);
	      mg2->Add(jpcgraphs[name]);
	    }
	    jpcgraphs[name]->SetPoint(i,masses[i],val);
	    
	  }
	}
	catch(std::exception e){
	  cerr<< "Cought exception" << endl;
	}
	++it3;
      }
      ++it2;
    }
  }
  

  // print results:
  it2=wavemap.begin(); 
  while(it2!=wavemap.end()){
    cout << it2->first << " " << it2->second.first*0.8 << endl;
    ++it2;
  }

  TFile* file=TFile::Open("normgraphs.root","RECREATE");
  mg->Write("graphs");
  mg2->Write("jpcgraphs");
  file->Close();

  return 0;
}
