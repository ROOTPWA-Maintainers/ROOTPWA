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


// programm to add up two amps with a specific branching and phase
#include <complex>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <unistd.h>
#include <stdlib.h>


using namespace std;


void printUsage(char* prog)
{
  cerr << "Add amplitude files with a phase between and a branching ratio." << endl;
  cerr << "usage single amp mode:" << endl;
  cerr << "    " << prog << " inputfile1 inputfile2 outputfile [phase [ratio]]" << endl;
  cerr << "usage multiple amp mode:" << endl;
  cerr << "    " << prog << " filelist <backup dir>" << endl;
  cerr << "        filelist format: " << endl;
  cerr << "        <file1>" << endl;
  cerr << "        <file2>" << endl;
  cerr << "        <outfile>" << endl;
  cerr << "        <phase>" << endl;
  cerr << "        <amplitude ratio>" <<endl;
}


int main(int argc, char** argv)
{

  vector<string> file1v;
  vector<string> file2v;
  vector<string> outfilev;
  vector<double> phasev;
  vector<double> branchv;
  string backupDir("./");

  if (argc == 3) { // multiple amp mode
    // open filelist
    ifstream flist(argv[1]);
    backupDir=argv[2];
    while(flist.good()){
      char f1[200];
      flist.getline(f1,200);
      char f2[200];
      flist.getline(f2,200);
      char f3[200];
      flist.getline(f3,200);

      file1v.push_back(f1);
      file2v.push_back(f2);
      outfilev.push_back(f3);

      char ph[60];
      double p;
      flist.getline(ph,60);
      if(string(ph)=="pi")p=3.14159592654;
      else p=atof(ph);

      char br[60];
      flist.getline(br,60);

      phasev.push_back(p);
      branchv.push_back(atof(br));
    } // read filelist
  } // end if multiple file mode
  else if (argc > 3 && argc <= 6) {
    file1v.push_back(argv[1]);
    file2v.push_back(argv[2]);
    outfilev.push_back(argv[3]);
    if (argc > 4) {
      double phase;
      if (string(argv[4]) == "pi")
	phase = 3.14159592654;
      else
	phase = atof(argv[4]);
      phasev.push_back(phase);
    }
    else
      phasev.push_back(0);
    if (argc > 5)
      branchv.push_back(atof(argv[5]));
    else
      branchv.push_back(1);
  }
  else {
    printUsage(argv[0]);
    return 1;
  }


  // loop through all amplitudes in list
  unsigned int n=file1v.size();
  cerr << n << " amplitudes in list" << endl;
  for(unsigned int i=0;i<n;++i){

    cout << endl;

    ifstream file1;
    file1.open(file1v[i].c_str());
    cout << file1v[i] <<  endl;
    if(file1.fail()){
      cerr << "*** Cannot open inputfile! Skipping!" << endl;
      continue;
    }
    cout << file2v[i] << endl;
    ifstream file2;
    file2.open(file2v[i].c_str());
    if(file2.fail()){
      cerr << "*** Cannot open inputfile! Skipping!" << endl;
      continue;
    }

    if(outfilev[i]==file1v[i] || outfilev[i]==file2v[i]){
      cerr << "Overwriting of files not allowed! Skipping" << endl;
      continue;
    }

    ofstream out(outfilev[i].c_str());

    cout << "Phase=" <<phasev[i] << "   Br=" << branchv[i] << endl;
    cout << "---> " << outfilev[i] << endl;

    complex<double> a1;
    complex<double> a2;

    complex<double> phase(1,0);
    double phi=phasev[i];
    phase=complex<double>(cos(phi),sin(phi));

    double R=branchv[i]; // branching ratio

    complex<double> amp;

    while (file1.read((char*) &a1,sizeof(complex<double>))){
      //cout << a1 << endl;
      file2.read((char*) &a2,sizeof(a2));
      amp=1./sqrt(1+R)*(a1+R*phase*a2);
      out.write((char*)&amp,sizeof(amp));
    }

    file1.close();
    file2.close();
    out.close();

    // move original amplitudes to backup directory
    if (backupDir != "./") {
      string com("mv ");
      com.append(file1v[i].c_str());
      com.append(" ");
      com.append(backupDir);
      int ret = system(com.c_str());
      if (ret != 0)
	cerr << "command '" << com << "' was not successful." << endl;

      com=("mv ");
      com.append(file2v[i].c_str());
      com.append(" ");
      com.append(backupDir);
      ret = system(com.c_str());
      if (ret != 0)
	cerr << "command '" << com << "' was not successful." << endl;
    }

  } // end loop over files

  return 0;

}
