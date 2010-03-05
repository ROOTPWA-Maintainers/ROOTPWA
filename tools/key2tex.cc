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


// creates latex output from wavefile name
// reades from stdin


#include <iostream>
#include <map>

#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"


using namespace std;


int main(int argc, char** argv)
{
  const unsigned int maxNmbWavesPerPage = 15;
  
  // setup isobar dictionary key->tex
  map<TString, TString> isobars;
  isobars["pi+"     ] = "\\pi^+";
  isobars["pi-"     ] = "\\pi^-";
  isobars["pi+-"    ] = "\\pi^\\pm";
  isobars["pi-+"    ] = "\\pi^\\mp";
  isobars["sigma"   ] = "\\sigma";
  isobars["rho770"  ] = "\\rho(770)";
  isobars["a11269"  ] = "a_1(1269)";
  isobars["a21320"  ] = "a_2(1320)";
  isobars["rho1450" ] = "\\rho(1450)";
  isobars["rho1700" ] = "\\rho(1700)";
  isobars["pi1300"  ] = "\\pi(1300)";
  isobars["pi1800"  ] = "\\pi(1800)";
  isobars["pi21670" ] = "\\pi_2(1670)";
  isobars["f01370"  ] = "f_0(1370)";
  isobars["f01500"  ] = "f_0(1500)";
  isobars["f01700"  ] = "f_0(1700)";
  isobars["f11285"  ] = "f_1(1285)";
  isobars["f11420"  ] = "f_1(1420)";
  isobars["b11235"  ] = "b_1(1235)";
  isobars["b11800"  ] = "b_1(1800)";
  isobars["b11500"  ] = "b_1(1500)";
  isobars["f21270"  ] = "f_2(1270)";
  isobars["f21950"  ] = "f_2(1950)";
  isobars["f21565"  ] = "f_2(1565)";
  isobars["f21270"  ] = "f_2(1270)";
  isobars["f22010"  ] = "f_2(2010)";
  isobars["f11420"  ] = "f_1(1420)";
  isobars["eta1440" ] = "\\eta(1420)";
  isobars["eta21645"] = "\\eta_2(1645)";
  isobars["rho31690"] = "\\rho_3(1690)";
  isobars["a21320"  ] = "a_2(1320)";

  // print latex header
  cout << "\\documentclass[12pt,a4paper]{article}" << endl
       << "\\usepackage{amsmath,amsthm,amssymb}"   << endl
       << "\\begin{document}"                      << endl
       << "\\begin{align*}"                        << endl
       << "\\begin{aligned}"                       << endl;
  string line;
  int    countWave = 0;
  while (!(cin >> line).eof()) { // begin event loop

    TString l(line);
    if (l.IsAlnum()) {
      cout << l;
      continue;
    }

    if (countWave > 0) {
      cout << " \\\\" << endl;
      // pagebreak
      if (countWave % maxNmbWavesPerPage == 0)
	cout << "\\end{aligned}"   << endl
	     << "\\end{align*}"    << endl
	     << "\\pagebreak"      << endl
	     << "\\begin{align*}"  << endl
	     << "\\begin{aligned}" << endl;
    }
    ++countWave;

    cerr << countWave << ": " << line << endl;
    // remove file extension
    l.Remove(l.Length() - 4);
    // extract X quantum numbers
    const TString head = l(0, 7);
    const TString I    = head(0, 1); 
    const TString G    = head(1, 1);
    const TString J    = head(2, 1);
    const TString P    = head(3, 1);
    const TString C    = head(4, 1);
    const TString M    = head(5, 1);
    const TString refl = head(6, 1);
    l.Remove(0, 7);
    // print X quantum numbers
    cout << I << "^" << G << J << "^{" << P << C << "}" << M << "^" << refl << "\\quad & ";

    // tokenize input
    TObjArray* tokens = l.Tokenize("_=");
    int        mode   = 0;
    for (int i = 0; i < tokens->GetEntries(); ++i) {
      const TString token = ((TObjString*)tokens->At(i))->GetString();
      cerr << "    " << mode << ": '" << token << "'" << endl;
      if (mode == 0) {  // isobar mode
	if (isobars.find(token) != isobars.end())
	  cout << isobars[token];
	else
	  cout << token;
	cout << " ";
	// check which mode to switch to depending whether we get _ or =
	l.Remove(0, token.Length());
	if (l(0, 1) == "_")
	  mode = 1;
	else
	  mode = 2;
      } else if (mode == 1) {  // ls mode
	if (token.Length() == 1)  // only l
	  cout << "[" << token << "] ";
	else
	  cout << "\\left[\\begin{array}{c}" << token(0, 1) << "\\\\"
	       << token(1, 1) << "\\end{array}\\right]";
	l.Remove(0, token.Length());
	mode = 0;
      } else if (mode == 2) {
	cout << "\\rightarrow ";
	if (isobars.find(token) != isobars.end())
	  cout << isobars[token];
	else
	  cout << token;
	cout << " ";
	l.Remove(0, token.Length());
	if (l(0, 1) == "_" )
	  mode = 1;
	else
	  mode = 2;
      }
      l.Remove(0, 1); // remove delimiter
    }
    cout << " & ";    
  }
  cout << "\\end{aligned}"  << endl
       << "\\end{align*}"   << endl
       << "\\end{document}" << endl;
};
