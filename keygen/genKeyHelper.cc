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
//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <complex>
#include <limits>

#include "TSystem.h"

#include "Tgamp.h"
#include "particleKey.h"
#include "waveKey.h"

#include "utilities.h"
#include "genKeyHelper.h"


using namespace std;
using namespace rpwa;


bool
rpwa::testKeyFile(const string& keyFileName,       // file name of key file under test
		  const string& dataFileName,      // file with test data in .evt format
		  const string& pdgTableFileName,  // path to PDG table file
		  const double  precision,         // warn threshold for comparison |1 - amp. / refl. amp.|
		  const bool    printAmp)          // if true amplitude values are printed in addition to ratios amp. / refl. amp.
{
  ifstream dataFile(dataFileName.c_str());
  if (!dataFile.good()) {
    printErr << "problems reading data file." << endl;
    return false;
  }
  Tgamp gamp(pdgTableFileName);
  gamp.reflect(false);
  vector<complex<double> > amps = gamp.Amp(keyFileName, dataFile);
  //vector<complex<double> > amps = calcAmplitudes(keyFileName, dataFile, pdgTableFileName, false);
  dataFile.clear();
  dataFile.seekg(0, ios::beg);
  if (!dataFile.good()) {
    printErr << "problems rereading data file." << endl;
    return false;
  }
  gamp.reflect(true);
  vector<complex<double> > ampsRefl = gamp.Amp(keyFileName, dataFile);
  if (amps.size() != ampsRefl.size()) {
    cerr << "different number of events for reflected and unmodified data." << endl;
    return false;
  }
  printInfo << "comparing amplitudes of " << amps.size() << " events "
	    << "with the amplitudes of the respective reflected events:" << endl;
  unsigned int countAmpRatioNotOk = 0;
  for (unsigned int i = 0; i < amps.size(); ++i) {
    const unsigned int nmbDigits = numeric_limits<double>::digits10 + 1;
    ostringstream s;
    s.precision(nmbDigits);
    s.setf(ios_base::scientific, ios_base::floatfield);
    s << "event " << i << ": "
      << "real part / real part refl. = " << amps[i].real() / ampsRefl[i].real() << ", "
      << "imag. part / imag. part refl. = " << amps[i].imag() / ampsRefl[i].imag();
    bool ampRatioOk = true;
    if (   (fabs(amps[i].real() / ampsRefl[i].real()) - 1 > precision)
	|| (fabs(amps[i].imag() / ampsRefl[i].imag()) - 1 > precision)) {
      ampRatioOk = false;
      ++countAmpRatioNotOk;
    }
    cout << "    " << s.str() << ((!ampRatioOk) ? " <!!!" : "") << endl;
    if (printAmp) {
      s.str("");
      s << "amp. = " << amps[i] << ", amp. refl. = " << ampsRefl[i];
      cout << "        " << s.str() << endl;
    }
  }
  if (countAmpRatioNotOk == 0) {
    printInfo << "relative differences of absolute values of all real and imaginary parts "
	      << "are smaller than " << precision << endl;
    return true;
  } else {
    printWarn << "for " << countAmpRatioNotOk << " amplitude(s) the relative differences "
	      << "of the abolute values of real or imaginary part "
	      << "differ by more than " << precision << endl;
    return false;
  }
}


void
rpwa::generateKeyFile(const waveKey& wave,              // complete isobar decay spec
		      const string&  srcMacroFileName,  // path to macro that will be copied to <wave name>.C
		      const bool     testKey,           // if true amplitude behavior under reflectivity is tested
		      const string&  dataFileName,      // file with test data in .evt format
		      const string&  pdgTableFileName)  // path to PDG table file
{
  const string       waveName    = wave.waveName(false).Data();
  const string       keyFileName = wave.waveName(true ).Data();
  const particleKey& X           = wave.mother();
  printInfo << "generating .key file for wave " << waveName << endl
	    << "    isobar decomposition:" << endl;
  X.write(cout, 4);
  cout << "    number of negative FS particles = " << X.countFsCharge(-1) << endl
       << "    number of neutral  FS particles = " << X.countFsCharge( 0) << endl
       << "    number of positive FS particles = " << X.countFsCharge(+1) << endl;

  // test amplitude
  bool keyFileOk = true;
  if (testKey) {
    wave.write("none");
    keyFileOk = testKeyFile(keyFileName, dataFileName, pdgTableFileName, 1e-9, true);
    printInfo << "reflectivity = " << wave.refl() << ": amplitude of reflected event should have the "
	      << ((wave.refl() > 0) ? "SAME" : "OPPOSITE") << " sign." << endl;
  }
  
  cout <<"keep .key file? (y/n)" << endl;
  char keepFile;
  cin >> keepFile;
  if ((keepFile == 'n') || (keepFile == 'N')) {
    const string command = "rm '" + keyFileName + "'";
    printInfo << "executing " << command << endl;
    gSystem->Exec(command.c_str());
    cout << ".key file removed!" << endl;
  } else {
    // regenerate keyfile with output mode set to binary
    wave.write("binary");
    // copy macro file
    const string command = "cp '" + srcMacroFileName + "' '" + keyFileName + ".C'";
    printInfo << "executing " << command << endl;
    gSystem->Exec(command.c_str());
    cout << "copied macro to file '" << keyFileName << ".C'" << endl;
  }
}
