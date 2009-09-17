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
// script to calculate binning with equal amount of data
void equalbins(TGraph* h, int n){
  int nbinsOrig=h->GetN();
  double ntot=0;
  double* data=h->GetY();
  double* x=h->GetX();
  for(int j=0; j<nbinsOrig; ++j)ntot+=data[j];

  int perbin=(int)floor((double)ntot/(double)n);
  
 
  
  std::cout<<"Integral: "<<ntot
	   <<"   binsOrig: "
	   <<nbinsOrig<<"   PerBin: "<<perbin<<endl;

  int content=0;
  int laststart=1;
  for(int i=1;i<nbinsOrig; ++i){
    content+=data[i];
    if(content>perbin){
      std::cout<<laststart<<"."<<i<<endl;
      laststart=i+1;
      content=0;
    }
  }
  std::cout<<laststart<<"."<<i<<endl;
}
