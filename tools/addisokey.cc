
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


///////////////////////////////////////////////////////////////////////////
//    This utility compares two keys and returns an isospin-added key if
//    possible. Otherwise it exits with an error
//   
//    The two keys to compare are the only arguments
///////////////////////////////////////////////////////////////////////////

#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include <vector>
#include <iostream>

using namespace std;

int main(int argc, char** argv) {

  if(strlen(argv[1])!=strlen(argv[2])){
    cerr << "Different Key Lengths" << endl;
    return 1;
  }

  TString key1(argv[1]);
  TString key2(argv[2]);

  //cout << key1 << endl;
  //cout << key2 << endl;

  //if(key1.Length()!=key2.Length()){
  //  cerr << "Different Key Lengths" << endl;
  //  return 1;
  //}

  TString key1head=key1(0,7);
  key1.Remove(0,7);
  TString key2head=key2(0,7);
  key2.Remove(0,7);

  if(key1head!=key2head){
    cerr << "Different Headers" << endl;
    return 2;
  }

  TObjArray* key1tokens=key1.Tokenize("_=.");
  TObjArray* key2tokens=key2.Tokenize("_=.");

  unsigned int l=key1tokens->GetEntries();
  if(l!=key2tokens->GetEntries()){
    delete key1tokens;
    delete key2tokens;
    cerr << "Different Number of Tokens" << endl;
    return 3;
  }
  
  TString output(key1head);
  TString rest(key1);

  // loop through tokens, check if everything agrees except pi+ pi- 
  for(unsigned int i=0;i<l;++i){
    TObjString* objtok1=(TObjString*)key1tokens->At(i);
    TObjString* objtok2=(TObjString*)key2tokens->At(i);
    TString tok1=objtok1->GetString();
    TString tok2=objtok2->GetString();

    if(tok1!=tok2){
      if(tok1=="pi+" && tok2=="pi-"){ 
	output.Append("pi+-");
      }
      else if(tok1=="pi-" && tok2=="pi+"){
	output.Append("pi-+");
      }
      else {
	delete key1tokens;
	delete key2tokens;
	cerr << "Different Tokens" << endl;
	return 4;
      }
    }
    else output.Append(tok1);

    rest.Remove(0,tok1.Length());
    // append correct token
    output.Append(rest(0,1));
    rest.Remove(0,1);

  }

  delete key1tokens;
  delete key2tokens;

  cout << output << endl;

  return 0;
}
