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

#include "utilities.h"
#include "particleKey.h"


using namespace std;
using namespace rpwa;


particleKey::particleKey(const TString& name,     // particle name
			 particleKey*   isobar1,  // first isobar
			 particleKey*   isobar2,  // second isobar (or bachelor)
			 int            L,        // relative angular momentum between isobars
			 int            S,        // total spin of isobars
			 const string&  massDep)  // mass dependence of amplitude (i.e. for sigma)
  : _name        (name),
    _L           (L),
    _S           (S),
    _isFsParticle(false),
    _id          (0),
    _massDep     (massDep)
{
  _isobars[0] = isobar1;
  _isobars[1] = isobar2;
  if (!_isobars[0] || !_isobars[1])
    _isFsParticle = true;
  for (unsigned int i = 0; i < 2; ++i)
    if (_isobars[i])
      _nmbFsParticles[i] = _isobars[i]->nmbFsParticles();
    else
      _nmbFsParticles[i] = 0;
}


particleKey::~particleKey()
{
}


TString 
particleKey::waveName() const
{
  TString name(_name);
  if (nmbFsParticles() > 2) {  // append decay chain according to <mother>=<isobar1>_LS_<isobar2>
    name += "=";
    name += _isobars[0]->waveName();
    name += "_";
    name += _L;
    if(_S >= 0)  // put S only when explicitly specified
      name += _S;
    name += "_";
    name += _isobars[1]->waveName();
  }
  name.ReplaceAll("(", "");
  name.ReplaceAll(")", "");
  return name;
}


unsigned int 
particleKey::nmbFsParticles(const int index) const
{
  if (index > 1) {
    printWarn << _name << ": index = " << index << " out of range. returning 0." << endl;
    return 0;
  }
  if (_isFsParticle)
    return 1;
  else
    if (index < 0)
      return _nmbFsParticles[0] + _nmbFsParticles[1];
    else
      return _nmbFsParticles[index];
}


void 
particleKey::fsParticles(vector<const particleKey*>& particles) const
{
  if (_isFsParticle)
    particles.push_back(this);
  else {
    _isobars[0]->fsParticles(particles);
    _isobars[1]->fsParticles(particles);
  }
}


void 
particleKey::fsCharges(vector<int>& charges) const
{
  if (_isFsParticle) {
    if (_name.Contains("pi-"))
      charges.push_back(-1);
    else if (_name.Contains("pi+"))
      charges.push_back(+1);
    else if (_name.Contains("pi0"))
      charges.push_back(0);
  } else {
    _isobars[0]->fsCharges(charges);
    _isobars[1]->fsCharges(charges);
  }
}


unsigned int
particleKey::countFsCharge(const int charge) const
{
  string fsName;
  switch (charge) {
  case -1:
    fsName = "pi-";
    break;
  case 0:
    fsName = "pi0";
    break;
  case +1:
    fsName = "pi+";
    break;
  default:
    return 0;
  }
  unsigned count = 0;
  if (_isFsParticle)
    if (_name.Contains(fsName.c_str()))
      count = 1;
    else
      count = 0;
  else {
    count += _isobars[0]->countFsCharge(charge);
    count += _isobars[1]->countFsCharge(charge);
  }
  return count;
}


void 
particleKey::setFsIds(const vector<int>& fsIds)
{
  if (_isFsParticle) {
    if (fsIds.size() != 1)
      printErr << _name << ": size of ID vector does not match "
	       << "number of final state particles." << endl;
    else
      _id = fsIds[0];
  } else {
    // split ID vector and pass subvectors to isobars
    vector<int> isoFsIds[2];
    isoFsIds[0].assign(fsIds.begin(),                      fsIds.begin() + _nmbFsParticles[0]);
    isoFsIds[1].assign(fsIds.begin() + _nmbFsParticles[0], fsIds.begin() + nmbFsParticles());
    for (unsigned int i = 0; i < 2; ++i)
      _isobars[i]->setFsIds(isoFsIds[i]);
  }
}


void
particleKey::write(ostream&           out,
		   const unsigned int offset) const
{
  indent(out, offset);
  int indentation = _name.Length() + 1;
  if (!_isFsParticle) {  // write decay chain
    // mother particle name and opening brace
    out << _name << "{" << endl;
    // isobars
    _isobars[0]->write(out, offset + indentation);
    _isobars[1]->write(out, offset + indentation);
    // LS
    indent(out, offset + indentation);
    out << "l=" << 2 * _L << endl;
    if (_S >= 0) {  // put S only when explicitly specified
      indent(out, offset + indentation);
      out << "s=" << 2 * _S << endl;
    }
    // closing brace
    indent(out, offset);
    out << "}" ;
    if (_massDep != "")
      out << " massdep=" << _massDep;
    out << endl;
  } else  // write final state particle with index
    out << _name << "[" << _id << "]" << endl;
}
