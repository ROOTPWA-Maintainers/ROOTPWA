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
#include <iterator>
#include <vector>
#include <string>
#include <map>
#include <set>
#include "TString.h"
#include "TRandom3.h"

using namespace std;

// helper class
class wsetentry {
public:
  TString name;
  double threshold;
  unsigned int index; // in waveset

  void operator =(const wsetentry& w){
    name=w.name; threshold=w.threshold; index=w.index;
  }

  friend bool operator <(const wsetentry& lhs,const wsetentry& rhs);
  friend bool operator ==(const wsetentry& lhs,const wsetentry& rhs);
  friend bool operator !=(const wsetentry& lhs,const wsetentry& rhs);
  bool operator <(const wsetentry& e){return this->index<e.index;}
  bool operator ==(const wsetentry& e){return this->name==e.name;}
  bool operator !=(const wsetentry& e){return !(*this==e);}
};



bool operator <(const wsetentry& lhs,const wsetentry& rhs){
  return lhs.index<rhs.index;
}
bool operator ==(const wsetentry& lhs,const wsetentry& rhs){
  return lhs.name==rhs.name;
}
bool operator !=(const wsetentry& lhs,const wsetentry& rhs){
  return !(lhs==rhs);
}

set<wsetentry>::iterator findByName(set<wsetentry>& aset, wsetentry entry){
   set<wsetentry>::iterator sit=aset.begin();
  while(sit!= aset.end()){
    set<wsetentry>::iterator sit2=sit;
    sit++;
    if(entry.name==sit2->name)return sit2;
  }
  return aset.end();
}

void removeFromSet(set<wsetentry>& aset, wsetentry entry){
  set<wsetentry>::iterator sit=findByName(aset,entry);
  if(sit!=aset.end())aset.erase(sit);
}


///                           Name    Threshold
void readWavelist(set<wsetentry>& result, const TString& input){
  ifstream file(input.Data());
  string line;
  unsigned int index=0;
  while(file.good()){
    getline(file,line);
    unsigned int pos=line.find(" ");
    string name=line.substr(0,pos);
    if(name.length()<2)continue;
    double thres;
    if(pos<line.length()){
      thres=atof(line.substr(pos,line.length()).c_str());
    } else thres=0;
    if(line.length()>1){
      wsetentry entry;
      entry.name=name;entry.threshold=thres;
      entry.index=index++;
      result.insert(entry);
    }
  }
}

set<wsetentry>::iterator
randomEntry(set<wsetentry>& myset, unsigned int start=0){
  if(start>=myset.size()){
    cerr<<" Requesting randomEntry " << start <<" out of range!" << endl;
    throw;
  }
  unsigned int x=start;
  if(start<myset.size()-1){
    unsigned int n=myset.size()-start;
    if(n!=0)x=gRandom->Integer(n)+start;
  }
  //cerr << "x=" <<x << endl;
  set<wsetentry>::iterator it=myset.begin();
  for(unsigned i=0;i<x;++i)++it;
  return it;
}
