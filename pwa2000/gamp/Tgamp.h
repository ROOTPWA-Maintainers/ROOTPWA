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


#ifndef TGAMP_HH
#define TGAMP_HH


#include <complex>
#include <vector>
#include <iostream>
#include <string>
#include <lorentz.h>

class event;


class Tgamp {

public:

  Tgamp(const std::string& pdgTableFileName = "");
  virtual ~Tgamp();
  
  std::complex<double> Amp(const std::string& keyFileName,
			   event&             ev) const;
  std::complex<double> Amp(const unsigned int iKey,
			   event&             ev) const;
  std::vector<std::complex<double> > Amp(const std::string& keyFileName,
					 std::istream&      eventData,
					 const bool         testMode = false) const;
  std::vector<std::complex<double> > Amp(const unsigned int iKey,
					 std::istream&      eventData) const;

  void addWave(const std::string& keyFileName) { _keyFileNames.push_back(keyFileName); }

  // performs reflection through production plane
  static event reflectEvent  (const event& evIn);
  static event mirrorEvent   (const event& evIn);
  void         reflect       (const bool   flag = true) { _reflect        = flag; }
  void         mirror        (const bool   flag = true) { _mirror        = true; }
  void         suppressOutput(const bool   flag = true) { _suppressOutput = flag; }

private:

  std::vector<std::string> _keyFileNames;
  bool                     _reflect;         // if true events are reflected through production plane
  bool                     _mirror;         // if true events are parity-mirrored
  bool                     _suppressOutput;  // if true output from keyparse() is suppressed


  static lorentzTransform toGottfriedJackson(const event& ev);

};


#endif  // TGAMP_HH


//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
