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
//#include <iomanip>
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
		 const int    C,       // charge conjugation
		 const bool  doRefl)
  : _I           (I),
    _G           (G),
    _J           (J),
    _P           (P),
    _C           (C),
    _M           (M),
    _refl        (refl),
    _mother      (mother),
    _outputFormat("binary"),
    _waveDebug   (0),
    _doReflectivity (doRefl)
{
  assert(mother);
  if(_refl==0)_doReflectivity=false;

  if(_doReflectivity && _M < 0 ) {
    printWarn << "M = " << _M << " < 0 not allowed in reflectivity basis. "
	      << "setting M = -M > 0" << endl;
    _M *= -1;
  }
  if (_doReflectivity && (_refl != +1) && (_refl != -1)) { 
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
  name << _I << sign(_G) << _J << sign(_P) << sign(_C) << _M << sign(_refl)
       << _mother->isobar(0)->waveName()
       << "_" << _mother->L() << _mother->S() << "_"
       << _mother->isobar(1)->waveName();
  if (fileSuffix)
    name << ".key";
  return name.str();
}


void
waveKey::write(ostream&      out,
	       const string& outputFormat) const
{
  // get final state
  vector<const particleKey*> fsParts;
  _mother->fsParticles(fsParts);
  string fsName;
  for (unsigned int i = 0; i < fsParts.size(); ++i)
    fsName += fsParts[i]->name();
  // insert header
  out << "########################################################" << endl
      << "# gamp key file for                                    #" << endl
      << "# " << fsName << " final state" << endl
      << "########################################################" << endl
      << "# IGJPCMR<Isobar1>_LS_<Isobar2>                        #" << endl
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
waveKey::assignFsIds(vector<int>&                      fsIds,
		     const map<TString, vector<int> >& fsPartIds) const
{
  unsigned int nmbFsPart = 0;
  for (map<TString, vector<int> >::const_iterator i = fsPartIds.begin(); i != fsPartIds.end(); ++i)
    nmbFsPart += i->second.size();
  vector<const particleKey*> fsParts;
  _mother->fsParticles(fsParts);
  assert(fsParts.size() == nmbFsPart);
  fsIds.clear();
  fsIds.resize(nmbFsPart, 0);
  map<TString, unsigned int> countFs;
  for (unsigned int i = 0; i < fsParts.size(); ++i) {
    const TString partName = fsParts[i]->name();
    map<TString, vector<int> >::const_iterator ids = fsPartIds.find(partName);
    assert(ids != fsPartIds.end());
    fsIds[i] = ids->second[countFs[partName]];
    ++countFs[partName];
  }
}


// builds amplitude in reflectivity basis following the notation in
// Suh-Urk's note "Formulas for Partial-Wave Analysis - Version V"
// |A M refl> = theta(M) [ |A M> - refl P (-)^(J - M) |A -M> ] (eq. 38)
void 
waveKey::writeReflectivityBasisAmp(ostream&           out,
				   const unsigned int offset) const
{
  if(_doReflectivity){
    const double theta = (_M == 0 ? 0.5 : 1 / sqrt(2));
    indent(out, offset);
    out << "# reflectivity eigenstate for epsilon = " << sign(_refl) << endl;
    indent(out, offset);
    out << maxPrecision(theta) << " * (" << endl;
    setStateName(_J, _P, _M);
    _mother->write(out, offset + 4);
    const double factor = -(double)_refl * (double)_P * pow((double)-1, _J - _M);
    indent(out, offset);
    out << sign(factor) << endl;
    setStateName(_J, _P, -_M);
    _mother->write(out, offset + 4);
    indent(out, offset);
    out << ")" << endl;
  }
  else {
    indent(out, offset);
    setStateName(_J, _P, _M);
    _mother->write(out, offset + 4);
  }
}


// permutes indistinguishable FS particles
void 
waveKey::permuteFsParts(map<TString, vector<int> >&           fsPartIds,
			map<TString, vector<int> >::iterator& thisPart,
			const vector<const particleKey*>      fsParts,
			bool&                                 firstCall,
			unsigned int&                         countTerm,
			const unsigned int                    nmbCombinations,
			ostream&                              out,
			const unsigned int                    offset) const
{
  do {
    map<TString, vector<int> >::iterator nextPart = thisPart;
    if (++nextPart != fsPartIds.end())
      permuteFsParts(fsPartIds, nextPart, fsParts, firstCall,
		     countTerm, nmbCombinations, out, offset);
    else {
      vector<int> fsIds;
      assignFsIds(fsIds, fsPartIds);
      _mother->setFsIds(fsIds);
      indent(out, offset);
      if (firstCall)
	firstCall = false;
      else {
	indent(out, offset);
	out << " +" << endl;
      }
      indent(out, offset + 4);
      out << "# Bose symmetrization term " << ++countTerm
	  << " of " << nmbCombinations << ": A(";
      for (unsigned int j = 0; j < fsIds.size(); ++j) {
	out << fsParts[j]->name() << "_" << fsIds[j];
	if (j == fsIds.size() - 1)
	  out << ")" << endl;
	else
	  out << ", ";
      }
      writeReflectivityBasisAmp(out, offset + 4);
      out << endl;
    }
  } while (next_permutation(thisPart->second.begin(), thisPart->second.end()));
}


void 
waveKey::writeBoseSymmetrizedAmp(ostream&           out,
				 const unsigned int offset) const
{
  // get final state particle mutliplicities
  const map<TString, unsigned int> fsPartMult = _mother->fsPartMult();

  // initialize IDs for final state particles
  map<TString, vector<int> > fsPartIds;
  for (map<TString, unsigned int>::const_iterator i = fsPartMult.begin(); i != fsPartMult.end(); ++i) {
    const TString      partName = i->first;
    const unsigned int partMult = i->second;
    fsPartIds[partName].resize(partMult, 0);
    for (unsigned int j = 0; j < partMult; ++j)
      fsPartIds[partName][j] = j + 1;
  }

  // calculate normalization factor for symmetrization of FS particles
  double nmbCombinations = 1;
  for (map<TString, unsigned int>::const_iterator i = fsPartMult.begin(); i != fsPartMult.end(); ++i)
    nmbCombinations *= TMath::Factorial(i->second);
  const double normFactor = 1 / sqrt(nmbCombinations);

  // Bose symmetrize amplitudes in reflectivity basis
  indent(out, offset);
  out << maxPrecision(normFactor) << " * (" << endl << endl;
  bool                       first     = true;
  unsigned int               countTerm = 0;
  vector<const particleKey*> fsParts;
  _mother->fsParticles(fsParts);
  map<TString, vector<int> >::iterator firstPart = fsPartIds.begin();
  permuteFsParts(fsPartIds, firstPart, fsParts, first,
		  countTerm, (unsigned int)nmbCombinations, out, offset);
  indent(out, offset);
  out << ")";
}
