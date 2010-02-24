///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert, Boris Grube
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
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cassert>

#include "TMath.h"

#include "utilities.h"
#include "waveKey.h"


using namespace std;
using namespace rpwa;


waveKey::waveKey(particleKey* mother,  // decay chain of mother particle
		 const int    J,       // total angular momentum
		 const int    P,       // parity
		 const int    M,       // z-projection of J
		 const int    refl,    // reflectivity
		 const int    I,       // isospin
		 const int    G,       // G-parity
		 const int    C)       // charge conjugation
  : _I           (I),
    _G           (G),
    _J           (J),
    _P           (P),
    _C           (C),
    _M           (M),
    _refl        (refl),
    _mother      (mother),
    _outputFormat("binary"),
    _waveDebug   (0)
{
  assert(mother);
  if(_M < 0) {
    printWarn << "M = " << _M << " < 0 not allowed in reflectivity basis. "
	      << "setting M = -M > 0" << endl;
    _M *= -1;
  }
  if ((_refl != +1) && (_refl != -1)) { 
    printErr << "reflectivity " << _refl << " does not make sense; must be +-1. aborting." << endl;
    throw;
  }
}


waveKey::~waveKey()
{
}


// constructs wave name according to IGJPCME<isobar 1>_LS_<isobar 2>
TString
waveKey::waveName(const bool fileSuffix) const
{
  stringstream name;
  name << _I << (_G > 0 ? "+" : "-")
       << _J << (_P > 0 ? "+" : "-") << (_C > 0 ? "+" : "-")
       << _M << (_refl >= 0 ? "+" : "-")
       << _mother->isobar(0)->waveName()
       << "_" << _mother->L() << _mother->S() << "_"
       << _mother->isobar(1)->waveName();
  if (fileSuffix)
    name << ".key";
  return name.str().c_str();
}


void
waveKey::write(ostream&      out,
	       const string& outputFormat) const
{
  // get final state
  vector<const particleKey*> fsParticles;
  _mother->fsParticles(fsParticles);
  string fsName;
  for (unsigned int i = 0; i < fsParticles.size(); ++i)
    fsName += fsParticles[i]->name();
  // insert header
  out << "########################################################" << endl
      << "# gamp key file for                                    #" << endl
      << "# " << fsName << " final state" << endl
      << "########################################################" << endl
      << "# IGJPCME<Isobar1>_LS_<Isobar2>                        #" << endl
      << "#                                                      #" << endl
      << "# "<< waveName() << endl
      << "#                                                      #" << endl
      << "########################################################" << endl
      << "# !!! all angular momenta are in units of 1/2 hbar !!! #" << endl
      << "########################################################" << endl
      << endl;
  // instert basic parameters
  out << "debug = " << _waveDebug << ";" << endl
      << "channel = t;" << endl
      << "mode = " << ((outputFormat != "") ? outputFormat : _outputFormat) << ";" << endl
      << endl;
  // insert amplitude using reflectivity basis
  writeBoseSymmetrizedAmp(out);
  out << ";" << endl;
}


void
waveKey::write(const string& outputFormat) const
{
  const string keyFileName = waveName(true).Data();
  printInfo << "creating keyfile '" << keyFileName << "'" << endl;
  ofstream out(keyFileName.c_str());
  write(out, outputFormat);
  out.close();
}


// sets name of _mother based on given quantum numbers
void
waveKey::setStateName(const int J,
		      const int P,
		      const int M) const
{
  stringstream name;
  name << "J = " << J * 2 << " P = " << P << " M = " << M * 2 << " ";
  _mother->setName(name.str().c_str());
}


void 
waveKey::assignFsIds(vector<int>&       fsIds,
		     const vector<int>& fsNegIds,
		     const vector<int>& fsPosIds,
		     const vector<int>& fsNeutIds,
		     const vector<int>& fsCharges) const
{
  assert(fsCharges.size() == fsNegIds.size() + fsPosIds.size() + fsNeutIds.size());
  unsigned int countFsNeg  = 0;
  unsigned int countFsPos  = 0;
  unsigned int countFsNeut = 0;
  for (unsigned int i = 0; i < fsCharges.size(); ++i) {
    if (fsCharges[i] == -1)
      fsIds[i] = fsNegIds[countFsNeg++];
    else if (fsCharges[i] == +1)
      fsIds[i] = fsPosIds[countFsPos++];
    else if (fsCharges[i] == 0)
      fsIds[i] = fsNeutIds[countFsNeut++];
    assert(abs(fsCharges[i]) <= 1);
  }
}


// builds amplitude in reflectivity basis following the notation in
// Suh-Urk's note "Formulas for Partial-Wave Analysis - Version V"
// |A M refl> = theta(M) [ |A M> - refl P (-)^(J - M) |A -M> ] (eq. 38)
void 
waveKey::writeReflectivityBasisAmp(ostream&           out,
				   const unsigned int offset) const
{
  const double theta = (_M == 0 ? 0.5 : 1 / sqrt(2));
  indent(out, offset);
  out << showpoint << setprecision(18)
      << theta << " * (" << endl;
  setStateName(_J, _P, _M);
  _mother->write(out, offset + 4);
  const double factor = -(double)_refl * (double)_P * pow(-1, _J - _M);
  indent(out, offset);
  out << (factor >= 0 ? " +" : " -") << endl;
  setStateName(_J, _P, -_M);
  _mother->write(out, offset + 4);
  indent(out, offset);
  out << ")" << endl;
}


void 
waveKey::writeBoseSymmetrizedAmp(ostream&           out,
				 const unsigned int offset) const
{
  // get final state particle charges
  vector<int> fsCharges;
  _mother->fsCharges(fsCharges);
  const unsigned int nmbFsNeg  = _mother->countFsCharge(-1);
  const unsigned int nmbFsNeut = _mother->countFsCharge( 0);
  const unsigned int nmbFsPos  = _mother->countFsCharge(+1);

  // initialize IDs for final state particles
  vector<int> fsNegIds (nmbFsNeg);
  vector<int> fsNeutIds(nmbFsNeut);
  vector<int> fsPosIds (nmbFsPos);
  for (unsigned int i = 0; i < nmbFsNeg; ++i)
    fsNegIds[i]  = i + 1;
  for (unsigned int i = 0; i < nmbFsNeut; ++i)
    fsNeutIds[i] = i + 1;
  for (unsigned int i = 0; i < nmbFsPos; ++i)
    fsPosIds[i]  = i + 1;

  // calculate normalization factor for symmetrization of FS particles
  const double nmbCombinations =   TMath::Factorial(nmbFsNeg)
                                 * TMath::Factorial(nmbFsNeut)
                                 * TMath::Factorial(nmbFsPos);
  const double normFactor      = 1 / sqrt(nmbCombinations);

  // Bose symmetrize amplitudes in reflectivity basis
  indent(out, offset);
  out << showpoint << setprecision(18)
      << normFactor << " * (" << endl << endl;
  bool         first     = true;
  unsigned int countTerm = 0;
  do {  // permute the negative FS particles
    do {  // permute the neutral FS particles
      do {  // permute the positive FS particles
	vector<int> fsIds(nmbFsPos + nmbFsNeg + nmbFsNeut);
	assignFsIds(fsIds, fsNegIds, fsPosIds, fsNeutIds, fsCharges);
	_mother->setFsIds(fsIds);
	indent(out, offset);
	if (first)
	  first = false;
	else {
	  indent(out, offset);
	  out << " +" << endl;
	}
	indent(out, offset + 4);
	out << "# Bose symmetrization term " << ++countTerm
	    << " of " << (unsigned int)nmbCombinations << ": A(";
	for (unsigned int i = 0; i < fsIds.size(); ++i) {
	  out << fsIds[i];
	  if (fsCharges[i] == -1)
	    out << "-";
	  else if (fsCharges[i] == +1)
	    out << "+";
	  if (i == fsIds.size() - 1)
	    out << ")" << endl;
	  else
	    out << ", ";
	}
	writeReflectivityBasisAmp(out, offset + 4);
	out << endl;
      } while (next_permutation(fsPosIds.begin(), fsPosIds.end()));
    } while (next_permutation(fsNeutIds.begin(), fsNeutIds.end()));
  } while (next_permutation(fsNegIds.begin(), fsNegIds.end()));
  indent(out, offset);
  out << ")";
}
