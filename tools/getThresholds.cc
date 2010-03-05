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
#include "TH1D.h"


using namespace std;


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
    cerr << mass << endl;
    //getline(0);
    config.ignore(255,'\n');
  }
  double weights[masses.size()];
  for(unsigned int i=0;i<masses.size();++i)weights[i]=1;


  list<string> waves=integrals[0].files();
  //unsigned int nwaves=waves.size();
  list<string>::iterator it=waves.begin();
  // wavemap: contains one entry per wave
  // each entry records mass and maximum integral value
  map<TString,pair<double,double> > wavemap; 
  map<TString,TGraph*> graphs;
  map<TString,TGraph*> diaggraphs_re;
  map<TString,TGraph*> diaggraphs_im;
  

  TMultiGraph* mg=new TMultiGraph();
  TMultiGraph* mg_re=new TMultiGraph();
  TMultiGraph* mg_im=new TMultiGraph();
  
  while(it!=waves.end()){
    wavemap[*it]=pair<double,double>(0,0);
    graphs[*it]=new TGraph(masses.size());
    graphs[*it]->SetName(it->c_str());
    mg->Add(graphs[*it]);
    list<string>::iterator itb=it;
    ++itb;
    while(itb!=waves.end()){
      TString s=(*it);
      s+=(*itb);
      
      diaggraphs_re[s]=new TGraph(masses.size());
      diaggraphs_im[s]=new TGraph(masses.size());
      mg_re->Add(diaggraphs_re[s]);
      mg_im->Add(diaggraphs_im[s]);
      
      ++itb;
    }
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
	  std::complex<double> c=integrals[i].val(it2->first.Data(),it3->first.Data());
	  //double n1=sqrt(integrals[i].val(it2->first.Data(),it2->first.Data()).real());
	  //double n2=sqrt(integrals[i].val(it3->first.Data(),it3->first.Data()).real());
	  //double norm=n1*n2;
	  double re=c.real();///norm;
	  double im=c.imag();///norm;
	  double val=integrals[i].val(it2->first.Data(),it3->first.Data()).real();///norm;
	  if(it2!=it3){
	    TString s=(it2->first);
	    s+=it3->first;
	    if(diaggraphs_re[s]!=NULL){
	      diaggraphs_re[s]->SetPoint(i,masses[i],re);
	      diaggraphs_re[s]->SetName(s);
	      diaggraphs_im[s]->SetPoint(i,masses[i],im);
	      diaggraphs_im[s]->SetName(s);
	    }
	  }
	  else if(it2==it3){
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
  // create histos
  map<TString,TGraph*>::iterator git=diaggraphs_re.begin();
  while(git!=diaggraphs_re.end()){
    TGraph* g=git->second;
    double max=-1E6;
    double min=1E6;
    double mean=0;
    double val[masses.size()];
    for(unsigned int i=1; i<masses.size(); ++i){
      double x,y;
      g->GetPoint(i,x,y);
      if(max<y)max=y;
      if(min>y)min=y;
      mean+=y;
      val[i]=y;
    }
    mean/=(double)masses.size();
    double range=(max-min)*0.2;
    //min/
    //double var=0;
    // for(unsigned int i=1; i<masses.size(); ++i){
//       val[i]=(mean-val[i])/mean;
//       var+=val[i]*val[i];
//     }
//     var=sqrt(var);
    
    //cerr << min << "   " << max << endl;
    TString name=g->GetName();
    name.Prepend("h");
    TH1D* histo=new TH1D(name,name,20,min-range,max+range);
    histo->FillN(masses.size(),val,weights);
    histo->Write();
    ++git;
  }



  
  mg->Write("graphs");
  mg2->Write("jpcgraphs");
  mg_re->Write("greal");
  mg_im->Write("gimag");
  file->Close();

  return 0;
}
